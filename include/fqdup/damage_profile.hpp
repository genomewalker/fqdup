#pragma once

// DamageProfile — ancient DNA damage model + per-position empirical mask.
//
// Global scope (not anonymous namespace) so it can be shared across TUs
// (derep.cpp, extend.cpp).  All methods are inline; no ODR risk.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "fqdup/logger.hpp"
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

struct DamageProfile {
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
    double mixture_d_reference = 0.0;
    int    mixture_K = 0;
    int    mixture_n_components = 0;
    bool   mixture_converged = false;
    bool   mixture_identifiable = false;
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
    static constexpr int MASK_POSITIONS = 15;
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

struct LengthBinDamageProfile {
    static constexpr int N_POS = 15;
    int length_lo = 0;
    int length_hi = 0;
    int64_t n_reads = 0;
    double d_max_5prime = 0.0;
    double d_max_3prime = 0.0;
    double lambda_5prime = 0.0;
    double lambda_3prime = 0.0;
    double bg_5prime = 0.0;
    double bg_3prime = 0.0;
    double cpg_contrast = std::numeric_limits<double>::quiet_NaN();
    bool validated = false;
    bool ss_mode = false;  // true if per_pos_3prime is C→T (SS), false if G→A (DS)
    std::string source = "none";
    // Per-position damage frequencies (positions 0..N_POS-1 from the 5'/3' end).
    // -1.0 marks positions with insufficient coverage (<100 reads).
    std::array<double, N_POS> per_pos_5prime_ct{};
    std::array<double, N_POS> per_pos_3prime{};

    // Per-length-bin K-component GC mixture (populated by libdart finalize).
    static constexpr int N_GC_BINS = 10;
    double mixture_d_damaged   = 0.0;
    double mixture_d_reference = 0.0;
    double mixture_d_population = 0.0;
    double mixture_pi_damaged  = 0.0;
    int    mixture_n_components = 0;
    bool   mixture_converged   = false;
    bool   mixture_identifiable = false;
    // Per-GC-bin d_max (10 bins, 0-10% ... 90-100%) within this length bin.
    // -1.0 if bin not valid (insufficient reads / C-sites).
    std::array<double,  N_GC_BINS> gc_d_max{};
    std::array<int64_t, N_GC_BINS> gc_n_reads{};
    std::array<double,  N_GC_BINS> gc_p_damaged{};

    // Per-read LLR unmixing results (populated when bulk params are supplied).
    // ancient = LLR(read | bulk ancient class) > LLR(read | bulk modern class).
    int64_t n_damaged = 0;
    int64_t n_undamaged  = 0;
    double  d_max_5_damaged = -1.0;
    double  d_max_3_damaged = -1.0;
    // Per-position T/(T+C) at 5' end (always) and at 3' end (T/(T+C) for SS,
    // A/(A+G) for DS) within the ancient-classified read subset.
    std::array<double, N_POS> per_pos_5prime_ct_damaged{};
    std::array<double, N_POS> per_pos_3prime_damaged{};
    // CpG-split 5' C→T in the ancient subset (metaDMG-style split). CpG = C
    // followed by G in the 5'→3' direction. -1 marks insufficient coverage.
    std::array<double, N_POS> per_pos_5prime_ct_cpg_damaged{};
    std::array<double, N_POS> per_pos_5prime_ct_noncpg_damaged{};
    double d_max_5_cpg_damaged    = -1.0;
    double d_max_5_noncpg_damaged = -1.0;
    // 8-oxoG marker: T/(T+G) at 5' terminal positions, split by ancient/modern.
    // s_gt_5_damaged_vs_undamaged = f_damaged(pos0) - f_modern(pos0); NaN if coverage low.
    // CAUTION: at pos 0 this signal is contaminated by the C→T elevation that
    // defines the ancient class. The clean 8-oxoG metric is g_to_t_5_damaged
    // below (derived from G-base depletion, independent of C→T).
    std::array<double, N_POS> per_pos_5prime_gt_damaged{};
    std::array<double, N_POS> per_pos_5prime_gt_undamaged{};
    double s_gt_5_damaged_vs_undamaged = std::numeric_limits<double>::quiet_NaN();

    // Clean 8-oxoG rate: G-fraction depletion at terminal vs interior baseline
    // within the ancient subset. Under pure C→T deamination, G-fraction at
    // terminal equals G-fraction at interior. Under 8-oxoG, terminal G-fraction
    // drops because G→T. g_to_t_5_damaged = p_G(interior) - p_G(terminal).
    double g_to_t_5_damaged = std::numeric_limits<double>::quiet_NaN();
    double pG_terminal_5_damaged = std::numeric_limits<double>::quiet_NaN();
    double pG_interior_5_damaged = std::numeric_limits<double>::quiet_NaN();

    // Reference-free trinucleotide spectrum (bulk, from libdart).
    // 64 contexts = prev*16 + mid*4 + next; A=0,C=1,G=2,T=3.
    // Terminal = read pos 1..4 from that end; interior = pos 10..14 (null baseline).
    // Downstream post-processing can contrast the terminal counters against
    // their interior counterparts for reference-free context analysis.
    std::array<int64_t, 64> tri_5prime_terminal{};
    std::array<int64_t, 64> tri_5prime_interior{};
    std::array<int64_t, 64> tri_3prime_terminal{};
    std::array<int64_t, 64> tri_3prime_interior{};
};

struct LengthStratifiedDamageProfile {
    std::string method = "single";
    std::vector<int> edges;
    std::vector<LengthBinDamageProfile> bins;
    int64_t reads_scanned = 0;
    int min_length = 0;
    int max_length = 0;

    // Joint length × GC 2-component mixture:
    // shared d_ancient, shared pi_ancient, per-cell w_ancient(b,g).
    // Treats each (length_bin, gc_bin) cell's d_max as an observation from
    // a 2-component Gaussian (undamaged at μ=0 vs damaged at μ=d_joint_ancient).
    double  d_joint_ancient   = 0.0;
    double  pi_joint_ancient  = 0.0;
    double  d_joint_population = 0.0;  // c_sites-weighted mean over all cells
    bool    joint_converged   = false;
    bool    joint_separated   = false;
    // w_ancient(b,g) = posterior P(damaged | cell_d_max_{b,g});
    // rows = length bins, cols = 10 GC bins. Empty if fit didn't run.
    std::vector<std::array<double, LengthBinDamageProfile::N_GC_BINS>> cell_w_ancient;
};

// Run DART damage estimation on a FASTQ file and return a populated DamageProfile.
// max_reads: stop after this many reads (0 = use all reads).
DamageProfile estimate_damage(
    const std::string& path,
    double mask_threshold,
    taph::SampleDamageProfile::LibraryType forced_lib =
        taph::SampleDamageProfile::LibraryType::UNKNOWN,
    int64_t max_reads = 0);

LengthStratifiedDamageProfile estimate_damage_by_length(
    const std::string& path,
    taph::SampleDamageProfile::LibraryType forced_lib,
    const LengthBinOptions& options,
    const std::vector<uint64_t>* prebuilt_hist = nullptr,
    size_t reader_threads = 0,
    int64_t max_reads = 0,
    const DamageProfile* bulk_for_llr = nullptr);
