#pragma once

// DamageProfile and the LSD helpers live in libtaph (taph/damage_profile.hpp) with the
// library profiler that fits them. fqdup's derep/extend/split use them through these
// aliases; the estimators below are fqdup's own.

#include "fqdup/logger.hpp"
#include "taph/damage_profile.hpp"

using taph::kOxChannelZDetect;
using taph::LsdLlrBinAccum;
using taph::LsdClassifyParams;
using taph::LengthBinOptions;
using taph::lsd_llr_score;
using taph::lsd_classify_read;
using taph::lsd_accumulate;
using taph::lsd_accumulate_soft;
using taph::compute_lsd_edges;
using taph::LSD_L_MIN;
using taph::LSD_L_MAX;
using taph::LSD_HIST_BINS;
using taph::LSD_MIN_READS_PER_BIN;
using taph::lsd_hist_bin;
using taph::Damage3PrimeMode;
using taph::make_lsd_classify_params;
using taph::LsdPrebuilt;
using DamageProfile                 = taph::BulkDamageProfile;
using LengthBinDamageProfile        = taph::LengthBinDamageProfile;
using LengthStratifiedDamageProfile = taph::LengthStratifiedDamageProfile;

// Run deamination damage estimation on a FASTQ file and return a populated DamageProfile.
// max_reads: stop after this many reads (0 = use all reads).
DamageProfile estimate_damage(
    const std::string& path,
    double mask_threshold,
    taph::SampleDamageProfile::LibraryType forced_lib =
        taph::SampleDamageProfile::LibraryType::UNKNOWN,
    int64_t max_reads = 0);

// ── Unified estimation + adapter/hexamer QC -----------------------------------
// T6.2: shared damage path used by `damage`, `derep`, and `extend` subcommands.

enum class DamageClipPolicy { Off, ReportOnly, Refit };

struct DamageEstimateOptions {
    double  mask_threshold      = 0.05;
    taph::SampleDamageProfile::LibraryType forced_lib =
        taph::SampleDamageProfile::LibraryType::UNKNOWN;
    int64_t max_reads           = 0;
    size_t  threads             = 1;
    bool    qc_enabled          = true;
    int64_t adapter_scan_reads  = 1'000'000;
    DamageClipPolicy clip_policy = DamageClipPolicy::ReportOnly;
    // T6.2 fix (Bug 2): when true, run the input scan and populate QC tables
    // but skip the deamination damage-fit emission — DamageProfile returned has
    // enabled=false and all zeros. Callers that want QC without committing
    // to damage-aware hashing (the common derep default) should set this.
    bool    skip_deam_fit       = false;
};

struct DamageQcReport {
    bool                       enabled              = false;
    bool                       profile_clipped      = false;
    int64_t                    adapter_reads_scanned = 0;
    taph::AdapterStubs         adapter;
    taph::HexStats             hex_stats;
    taph::LibraryQcFlags       flags;
    taph::PreservationSummary  preservation;
    double                     short_read_frac      = 0.0;
    std::vector<std::string>   flag_names;
};

struct DamageEstimate {
    DamageProfile  profile;
    DamageQcReport qc;
    // Log-length histogram accumulated during the primary scan.
    // Passed as prebuilt_hist to estimate_damage_by_length() to skip its
    // own histogram sub-pass (saves one full file read).
    std::vector<uint64_t> lsd_hist;
    // Total records READ in the scan (incl. length-filtered short reads). Equals max_reads
    // when the cap bound the scan -> more reads exist; < max_reads means EOF was reached.
    // Used to decide full-depth escalation without being fooled by short-read filtering.
    int64_t reads_seen = 0;
};

DamageEstimate estimate_damage_with_qc(const std::string& path,
                                       const DamageEstimateOptions& opts);
   // when set, skip FASTQ reader entirely

// Stripped-down per-bin T/TC accumulator for DamageSplitModel.
// Avoids all CpG/oxoG/base-composition tracking and eliminates log() calls from
// the read loop by precomputing additive LLR coefficients from the bulk model.
// Use instead of estimate_damage_by_length() when only the split model is needed.
LengthStratifiedDamageProfile estimate_damage_split_model(
    const std::string& path,
    const DamageProfile& bulk,
    const std::vector<uint64_t>* prebuilt_hist = nullptr,
    int n_workers = 4);

// Precomputed per-bin LOD tables for per-read ancient/modern classification.
// Built once from a LengthStratifiedDamageProfile; scoring is O(n_pos) with
// no log/division per read. Falls back to bulk exponential when bins is empty.
struct DamageSplitModel {
    struct Bin {
        int   lo = 0, hi = 0;
        bool  ss_mode = false;
        // lod_T[i] = log(p_dam[i] / p_und[i])   (observe T at 5' pos i → ancient)
        // lod_C[i] = log((1-p_dam[i])/(1-p_und[i])) (observe C → modern)
        std::array<float, LengthBinDamageProfile::N_POS> lod5_T{};
        std::array<float, LengthBinDamageProfile::N_POS> lod5_C{};
        std::array<float, LengthBinDamageProfile::N_POS> lod3_T{};
        std::array<float, LengthBinDamageProfile::N_POS> lod3_C{};
        // CpG-context 5' C->T lod (used when the next base is G): methyl-C deaminates faster, so the
        // damage excess is scaled by the bin's CpG/bulk amplitude ratio. Falls back to the bulk lod5
        // when the CpG amplitude is ill-estimated. Recovers the CpG channel single-rate scoring misses.
        std::array<float, LengthBinDamageProfile::N_POS> lod5_T_cpg{};
        std::array<float, LengthBinDamageProfile::N_POS> lod5_C_cpg{};
    };

    std::vector<Bin> bins;
    DamageProfile    fallback;

    // UNIFIED SCORER (integration): score() delegates to the shared closed-form
    // taph::lsd_llr_score — the SAME per-read math `fqdup profile` uses via
    // lsd_llr_from_sig. The per-bin empirical `bins` LOD tables are retained only
    // for BIN model (de)serialization compat; they no longer drive scoring. This
    // makes split≡profile "one mechanism" (user hard constraint) and, on synthetic
    // validation, is neutral for ss (ΔAUC≈−0.001) and +0.04 AUC for ds vs the old
    // empirical-bootstrap tables (which underweight the 3′ G→A because their 3′ LOD
    // is trained on a 5′-only-labeled subset). cls_ is derived lazily from fallback
    // so every construction path (build + read_split_model_bin) is covered.
    mutable taph::LsdClassifyParams cls_{};
    mutable bool cls_ready_ = false;

    bool valid() const { return !bins.empty(); }

    // Build from a finalized LengthStratifiedDamageProfile (must have n_damaged
    // and per_pos_5prime_ct_damaged populated via bulk_for_llr in estimate_damage_by_length).
    static DamageSplitModel build(const LengthStratifiedDamageProfile& lsd,
                                  const DamageProfile& bulk);

    // Score one read. Returns LLR > 0 for damaged, < 0 for undamaged.
    // UNIFIED: delegates to the shared taph::lsd_llr_score (== profile's math).
    float score(const std::string& seq, int n_pos = 14) const;

    // Legacy per-bin empirical-LOD scorer (pre-unification). Retained for A/B
    // comparison and BIN-model diagnostics; not used by split/profile routing.
    float score_empirical_legacy(const std::string& seq, int n_pos = 14) const;
};
