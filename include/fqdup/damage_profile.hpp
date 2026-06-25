#pragma once

// DamageProfile — ancient DNA damage model + per-position empirical mask.
//
// Global scope (not anonymous namespace) so it can be shared across TUs

// Oxidative channel detection threshold (binomial z-score)
static constexpr float kOxChannelZDetect = 3.0f;
// (derep.cpp, extend.cpp).  All methods are inline; no ODR risk.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "fqdup/logger.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/length_bin_damage_profile.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/sample_damage_profile.hpp"

inline constexpr int LSD_L_MIN             = 30;
inline constexpr int LSD_L_MAX             = 500;
inline constexpr int LSD_HIST_BINS         = 128;
inline constexpr int LSD_MIN_READS_PER_BIN = 100;

inline int lsd_hist_bin(int L) {
    static const double kLogMin = std::log(static_cast<double>(LSD_L_MIN));
    static const double kStep   = (std::log(static_cast<double>(LSD_L_MAX + 1)) - kLogMin)
                                  / LSD_HIST_BINS;
    int Lcap = L < LSD_L_MAX ? L : LSD_L_MAX;
    int hb   = static_cast<int>((std::log(static_cast<double>(Lcap)) - kLogMin) / kStep);
    if (hb < 0)              hb = 0;
    if (hb >= LSD_HIST_BINS) hb = LSD_HIST_BINS - 1;
    return hb;
}

// 3' damage prior mode. ds_ga = use 3' G→A signal (DS libraries),
// ss_ct = use 3' C→T signal (SS libraries), neutral = no 3' damage prior
// (UNKNOWN library type — drop the 3' channel from masking entirely).
enum class Damage3PrimeMode : uint8_t { ds_ga, ss_ct, neutral };

struct DamageProfile {
    // 3' damage prior mode — set from taph::SampleDamageProfile::library_type. UNKNOWN maps to
    // neutral so we don't silently apply DS or SS priors. ss_mode (below) is
    // kept for compatibility and derived from this field.
    Damage3PrimeMode mode_3prime = Damage3PrimeMode::ds_ga;

    // Fitted exponential model parameters (for logging and expected_mismatches only).
    double d_max_5prime   = 0.0;
    double d_max_3prime   = 0.0;
    double lambda_5prime  = 0.5;
    double lambda_3prime  = 0.5;
    double background     = 0.02;
    // Interior baseline fractions used by per-read LLR (null-distribution models).
    // bg_5_tc = 5' interior T/(T+C); bg_3_channel = 3' interior T/(T+C) for SS
    // libraries, A/(A+G) for DS libraries.
    double bg_5_tc        = 0.5;
    double bg_3_channel   = 0.5;
    // Bulk CpG-split 5' C→T d_max (from taph::dmax_ct5_cpg_like / _noncpg_like).
    // Negative (NaN) = not available. Used to weight the per-read 5' LLR: CpG
    // C/T sites are stronger evidence of ancient than non-CpG sites.
    double d_cpg_5prime    = -1.0;
    double d_noncpg_5prime = -1.0;
    double mask_threshold = 0.05;
    double pcr_error_rate = 0.0;
    double mixture_d_damaged = 0.0;
    double mixture_pi_damaged = 0.0;
    double mixture_d_population_highgc = 0.0;
    int    mixture_K = 0;
    int    mixture_n_components = 0;
    bool   mixture_converged = false;
    bool   mixture_identifiable = false;
    // Validated reference-free ancient fraction (libtaph SampleDamageProfile::pi, SOLUTION §6.3).
    // Shadow-mode step 2: carried for the contract-vs-mixture divergence shadow; not yet a verdict source.
    double pi_point    = -1.0;
    double pi_lo       = -1.0;
    double pi_hi       = -1.0;
    taph::DamageConfidence pi_state = taph::DamageConfidence::UNDETERMINED;
    bool   pi_detected = false;
    bool   enabled        = false;
    bool   ss_mode        = false;  // true = single-stranded library (C→T at both ends)

    // Empirical per-position mask (primary masking authority).
    //
    // mask_pos[p] = true means position p (measured from either terminus) sits
    // in the damage zone and will be replaced with a neutral byte before hashing.
    //
    // Symmetry invariant: canonical_hash(seq) == canonical_hash(revcomp(seq)).
    // apply_damage_mask() uses max(mask_pos[i], mask_pos[j]) semantics so that
    // position p from the 5' end and position p from the 3' end are treated
    // identically regardless of which strand is processed.
    // Upper bound on terminal positions tracked for damage masking.
    // 32 covers extreme aDNA scenarios (very heavy deamination > 20 bp deep)
    // while staying tiny in memory (~32 B per profile). All loops use this
    // constant directly — bump and recompile to extend further.
    static constexpr int MASK_POSITIONS = 32;
    bool mask_pos[MASK_POSITIONS] = {};

    // Populate mask_pos from the fitted exponential model.
    // Used when damage parameters are given manually rather than estimated.
    void populate_mask_from_model() {
        for (int p = 0; p < MASK_POSITIONS; ++p) {
            double exc5 = d_max_5prime * std::exp(-lambda_5prime * p);
            double exc3 = d_max_3prime * std::exp(-lambda_3prime * p);
            mask_pos[p] = std::max(exc5, exc3) > mask_threshold;
        }
    }

    // Helpers retained for expected_mismatches() and informational logging only.
    double dmg_5(int p) const { return d_max_5prime * std::exp(-lambda_5prime * p); }
    double dmg_3(int p) const { return d_max_3prime * std::exp(-lambda_3prime * p); }

    // P(damage-driven mismatch at pos) for a read of length L. Used by Phase 3
    // LR scoring (T5.2) to weigh damage-channel mismatches. ss_mode applies
    // C→T at both ends; ds_mode applies C→T at 5' and G→A at 3'. Caller is
    // responsible for checking that the mismatch alt is in the damage channel.
    double p_damage_at(int pos, int L) const {
        if (!enabled || L <= 0) return 0.0;
        double d5 = dmg_5(pos);
        double d3 = dmg_3(L - 1 - pos);
        return std::max(d5, d3);  // either end can drive damage at this position
    }

    double expected_mismatches(int L) const {
        // d_max and the exponential model are defined on the conditional T/(T+C)
        // scale: dmg_5(p) is the excess deamination rate among C/T bases only.
        // The per-position mismatch probability between two reads is therefore
        // dmg × (1-dmg) × f_C, where f_C = 1-background is the C fraction within
        // the T+C pool (background = T/(T+C) at undamaged positions). Without this
        // factor the expectation is inflated by ~1/(1-background) ≈ 2× for typical
        // aDNA base composition. The same factor applies to the 3' G→A channel by
        // complement symmetry (G fraction in A+G pool ≈ 1-background).
        const double f_suscep = 1.0 - background;
        double e = 0.0;
        for (int p = 0; p < L; ++p) {
            double d5 = dmg_5(p);
            double d3 = dmg_3(L - 1 - p);
            e += 2.0 * d5 * (1.0 - d5) * f_suscep;
            e += 2.0 * d3 * (1.0 - d3) * f_suscep;
        }
        e += 2.0 * L * pcr_error_rate;
        return e;
    }

    int mismatch_tolerance(int L) const {
        double lam = expected_mismatches(L);
        if (lam < 1e-9) return 0;
        double cumP = 0.0, pois = std::exp(-lam);
        int k = 0;
        while (cumP + pois < 0.99 && k < 100) {
            cumP += pois;
            ++k;
            pois *= lam / k;
        }
        return k;
    }

    void print_info(int typical_read_length) const {
        log_info("--- Damage-Aware Deduplication ---");
        log_info("  5'-end d_max:  " + std::to_string(d_max_5prime));
        log_info("  3'-end d_max:  " + std::to_string(d_max_3prime));
        log_info("  5'-end lambda: " + std::to_string(lambda_5prime));
        log_info("  3'-end lambda: " + std::to_string(lambda_3prime));
        log_info("  Background:    " + std::to_string(background));
        log_info("  Mask threshold:" + std::to_string(mask_threshold));
        if (pcr_error_rate > 0.0)
            log_info("  PCR error rate:" + std::to_string(pcr_error_rate) + " per bp");
        int n_masked = 0;
        std::string pos_str;
        for (int p = 0; p < MASK_POSITIONS; ++p) {
            if (mask_pos[p]) {
                if (n_masked > 0) pos_str += ',';
                pos_str += std::to_string(p);
                ++n_masked;
            }
        }
        if (n_masked > 0)
            log_info("  Masked positions: " + pos_str + " (" +
                     std::to_string(n_masked) + " bp each end)");
        else
            log_info("  Masked positions: none");
        if (typical_read_length > 0) {
            double e   = expected_mismatches(typical_read_length);
            int    tol = mismatch_tolerance(typical_read_length);
            log_info("  Expected mismatches (L=" + std::to_string(typical_read_length) +
                     "): " + std::to_string(e) +
                     ", 99th-pct tolerance: " + std::to_string(tol));
        }
    }
};

struct LengthBinOptions {
    enum class Mode { DISABLED, AUTO, QUANTILE, EXPLICIT };

    Mode mode = Mode::DISABLED;
    int quantile_bins = 1;
    std::vector<int> explicit_edges;

    bool enabled() const { return mode != Mode::DISABLED; }
};

// Use canonical definitions from libtaph (taph/length_bin_damage_profile.hpp).
using LengthBinDamageProfile       = taph::LengthBinDamageProfile;
using LengthStratifiedDamageProfile = taph::LengthStratifiedDamageProfile;

// Per-bin ancient-subset accumulator used when fusing the LSD pass into the
// oxoG second pass (Fix E).  Field layout is identical to the local
// LlrBinAccum struct inside estimate_damage_by_length so that prebuilt data
// can be injected without copying.
struct LsdLlrBinAccum {
    int64_t n_damaged   = 0;
    int64_t n_undamaged = 0;
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> h_3_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> n_3_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_cpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc_cpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_noncpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc_noncpg{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_g{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tg_5_anc{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_mod_g{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tg_5_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> h_3_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> n_3_mod{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> a_5_anc_all{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> c_5_anc_all{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> g_5_anc_all{};
    std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_all{};
    // Soft-EM posterior-weighted accumulators for ancient-fraction d_max.
    // Soft-EM posterior-weighted C→T counts at N_SOFT_POS terminal positions.
    // Using multiple positions lets d_anc estimation avoid adapter-artifact
    // contamination at pos 0 by picking the peak over non-masked positions.
    static constexpr int N_SOFT_POS = 4;
    double sw_t5_anc[N_SOFT_POS]  = {};
    double sw_tc5_anc[N_SOFT_POS] = {};
    double sw_h3_anc[N_SOFT_POS]  = {};
    double sw_n3_anc[N_SOFT_POS]  = {};
    double sw_sum = 0.0;
};

// Parameters for per-read ancient/modern classification (LLR scorer).
struct LsdClassifyParams {
    double d_anc        = 0.0;
    double lam_5        = 0.3;
    double lam_3        = 0.3;
    double bg_5         = 0.5;
    double bg_3         = 0.5;
    bool   is_ss           = false;
    bool   skip_pos0_5prime = false;  // skip pos 0 in 5' LLR (hexamer artifact)
    double cpg_scale    = 1.0;
    double noncpg_scale = 1.0;
    // Step-2 shadow (SOLUTION §6.7): contract per-read amplitude (D_MAX_CONSERVED) + UNDETERMINED gate.
    // When d_anc_contract >= 0 the oxoG pass computes the contract class alongside the current one and
    // logs divergence; it never drives a verdict. contract_gated_off = library is not pi-DETECTED.
    double d_anc_contract     = -1.0;
    bool   contract_gated_off = false;
};

// Build classify params from a finalized bulk DamageProfile.
LsdClassifyParams make_lsd_classify_params(const DamageProfile& bulk);

// Return raw LLR score for one read (positive = ancient-consistent).
double lsd_llr_score(const std::string& seq, const LsdClassifyParams& p);

// Classify one read as ancient (true) or modern (false) using the LLR.
bool lsd_classify_read(const std::string& seq, const LsdClassifyParams& p);

// Accumulate terminal stats for one read into an LsdLlrBinAccum.
// ancient=true: accumulate damage-channel + composition + oxoG counts.
// ancient=false: accumulate only oxoG marker (T/G) for the modern subset.
void lsd_accumulate(const std::string& seq, LsdLlrBinAccum& acc,
                    bool ancient, bool is_ss);

// Soft-EM accumulation: add posterior-weighted pos-0 counts to sw_* fields.
// w = P(ancient | read, π), computed as sigmoid(llr + log_prior_odds).
void lsd_accumulate_soft(const std::string& seq, LsdLlrBinAccum& acc,
                         double w, bool is_ss);

// Compute bin edges from a log-length histogram + options (same logic as the
// edge-picking block inside estimate_damage_by_length, extracted so damage.cpp
// can call it before the oxoG pass to pre-configure worker LSD state).
std::vector<int> compute_lsd_edges(const std::vector<uint64_t>& hist,
                                   const LengthBinOptions& options);

// Prebuilt LSD data produced by a fused oxoG+LSD pass.  Pass via the
// prebuilt parameter of estimate_damage_by_length to skip the FASTQ reader.
struct LsdPrebuilt {
    std::vector<int>            edges;         // from compute_lsd_edges
    taph::LengthBinStats        merged_stats;  // accumulated but NOT finalized
    std::vector<LsdLlrBinAccum> llr_bins;      // per-bin, merged across threads
};

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

LengthStratifiedDamageProfile estimate_damage_by_length(
    const std::string& path,
    taph::SampleDamageProfile::LibraryType forced_lib,
    const LengthBinOptions& options,
    const std::vector<uint64_t>* prebuilt_hist = nullptr,
    size_t reader_threads = 0,
    int64_t max_reads = 0,
    const DamageProfile* bulk_for_llr = nullptr,
    const LsdPrebuilt* prebuilt = nullptr);   // when set, skip FASTQ reader entirely

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

    bool valid() const { return !bins.empty(); }

    // Build from a finalized LengthStratifiedDamageProfile (must have n_damaged
    // and per_pos_5prime_ct_damaged populated via bulk_for_llr in estimate_damage_by_length).
    static DamageSplitModel build(const LengthStratifiedDamageProfile& lsd,
                                  const DamageProfile& bulk);

    // Score one read. Returns LLR > 0 for damaged, < 0 for undamaged.
    // Uses per-bin empirical tables when valid(), else bulk exponential fallback.
    float score(const std::string& seq, int n_pos = 14) const;
};
