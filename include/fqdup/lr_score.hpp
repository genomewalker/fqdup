#pragma once

// Per-edge log-likelihood ratio for Phase 3 absorption decisions.
//
// Hypotheses (likelihoods over the observed mismatches at positions m_i):
//   H_pcr  : child is a PCR/sequencing error of parent
//            P(obs | H_pcr) = ∏ p_err(Q_i)
//   H_real : child is a real distinct molecule (variant or co-located capture)
//            P(obs | H_real) ≈ ∏ p_real(m_i, alt, base, bundle_occ)
//
// LR = log[ P(obs | H_pcr) * π_pcr / (P(obs | H_real) * π_real) ]
//
// Positive LR ⇒ favor PCR ⇒ absorb. Threshold gates the call (T5.4).
//
// `compute_lr` is intentionally independent of derep_detail internals so it
// can be unit-tested in isolation. Callers pass packed mismatch positions +
// quality bytes; the scorer does the math.
//
// Notes
// -----
// * Quality is integer Phred (1..60). FASTQ inputs are guaranteed to carry
//   per-base qualities; `phred_to_perr` is strict and clamps q to [1, 60] to
//   bound p_err away from 0/1 without inventing a neutral fallback.
// * Damage-aware path is delegated to T5.2 (damage_profile.hpp::p_damage)
//   and stitched in as a separate term. Here we keep the math purely
//   quality- and bundle-driven.
// * Bundle occupancy enters as a *prior* on H_real: rare bundles ⇒
//   independent capture events at this locus are unlikely ⇒ small π_real.
//   We pass it through `bundle_occ` (n clusters at this bundle_key).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstddef>

namespace fqdup::lr {

constexpr double  kP_real_uniform  = 0.75;   // baseline P(specific alt | real variant at pos)
                                             // (1/3 for the alt × 3-fold inflation for typical bias)
constexpr double  kBundlePrior_k   = 0.5;    // smoothing for bundle prior π_pcr/π_real

// Phred → p_err. Strict: FASTQ inputs are guaranteed to carry per-base
// qualities, so q==0 is treated as an upstream bug and clamped to Q1 rather
// than silently substituted with a "neutral" value. Upper clamp at Q60 keeps
// p_err away from underflow.
inline double phred_to_perr(uint8_t q) {
    if (q < 1)  q = 1;
    if (q > 60) q = 60;
    return std::pow(10.0, -static_cast<double>(q) / 10.0);
}

// Per-mismatch log-likelihood ratio. Positive ⇒ favor PCR error.
//
// p_damage in [0,1] is the prior probability that a mismatch at this position
// is driven by aDNA deamination on a real molecule (T5.2). It is added to the
// H_real likelihood for damage-channel mismatches, weakening the PCR favouring.
// For non-damage-channel mismatches the caller passes 0.0.
inline double per_mismatch_log_lr(uint8_t qual_at_pos,
                                  double  p_damage = 0.0) {
    double p_err  = phred_to_perr(qual_at_pos);
    // P(obs | H_pcr)  = p_err / 3
    // P(obs | H_real) = kP_real_uniform + p_damage * (1 - kP_real_uniform)
    //                 = baseline_real + damage_lift
    double num = p_err / 3.0;
    double den = kP_real_uniform + p_damage * (1.0 - kP_real_uniform);
    return std::log(num) - std::log(den);
}

// Bundle prior: log(π_pcr / π_real). For capture-enriched aDNA, HIGH bundle
// occupancy = locus was heavily captured = many independent molecules pulled
// in = mismatches more likely to be real co-located variants, not PCR siblings.
// So the prior penalises absorption (negative) as occupancy grows. Singleton
// bundle returns 0 (neutral; no local evidence either way).
//
// Magnitude calibrated so that a 100-cluster bundle subtracts ~kBundlePrior_alpha
// nats from the LR — typically 0.5 nats by default, comparable to one quality
// step at Q30. T5.7 calibration may tune this.
inline double bundle_log_prior(uint64_t bundle_occupancy) {
    if (bundle_occupancy <= 1) return 0.0;
    constexpr double kBundlePrior_alpha = 0.10;
    return -kBundlePrior_alpha *
           std::log1p(static_cast<double>(bundle_occupancy - 1));
}

// Compute the full LR for an edge with up to two mismatches (H≤2).
// Returns log( P_pcr * π_pcr / (P_real * π_real) ).
//
//   quals     : pointer into QualArena bytes for the *child* sequence
//   qual_len  : child quality length (bounds for safety)
//   m_pos     : pointer to mismatch positions (1 or 2 entries)
//   n_mm      : 1 (H=1) or 2 (H=2)
//   bundle_occ: occupancy of the parent's bundle key
// `quals` is indexed parallel to mismatches: quals[i] is the Phred at the i-th
// mismatch position (caller pre-extracts via QualArena::q_at). `qual_len` is
// the number of mismatch entries, NOT a sequence length. m_pos is retained for
// future per-position priors but not used to index quals.
inline double compute_lr(const uint8_t*  quals,
                         int             /*qual_len*/,
                         const uint16_t* /*m_pos*/,
                         int             n_mm,
                         uint64_t        bundle_occ,
                         const double*   p_damage = nullptr) {
    double lr = bundle_log_prior(bundle_occ);
    for (int i = 0; i < n_mm; ++i) {
        uint8_t q  = quals[i];
        double  pd = p_damage ? p_damage[i] : 0.0;
        lr += per_mismatch_log_lr(q, pd);
    }
    return lr;
}

} // namespace fqdup::lr
