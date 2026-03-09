#pragma once

// DamageProfile — ancient DNA damage model + per-position empirical mask.
//
// Global scope (not anonymous namespace) so it can be shared across TUs
// (derep.cpp, extend.cpp).  All methods are inline; no ODR risk.

#include <algorithm>
#include <cmath>
#include <string>

#include "fqdup/logger.hpp"
#include "dart/sample_damage_profile.hpp"

struct DamageProfile {
    // Fitted exponential model parameters (for logging and expected_mismatches only).
    double d_max_5prime   = 0.0;
    double d_max_3prime   = 0.0;
    double lambda_5prime  = 0.5;
    double lambda_3prime  = 0.5;
    double background     = 0.02;
    double mask_threshold = 0.05;
    double pcr_error_rate = 0.0;
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
        double e = 0.0;
        for (int p = 0; p < L; ++p) {
            double d5 = dmg_5(p);
            double d3 = dmg_3(L - 1 - p);
            e += 2.0 * d5 * (1.0 - d5);
            e += 2.0 * d3 * (1.0 - d3);
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

// Run DART damage estimation on a FASTQ file and return a populated DamageProfile.
// max_reads: stop after this many reads (0 = use all reads).
DamageProfile estimate_damage(
    const std::string& path,
    double mask_threshold,
    dart::SampleDamageProfile::LibraryType forced_lib =
        dart::SampleDamageProfile::LibraryType::UNKNOWN,
    int64_t max_reads = 0);
