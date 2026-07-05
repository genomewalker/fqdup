#pragma once

// Empirical posterior-odds model for Phase 3 PCR-error absorption (T5.8).
//
// Replaces the three hard-coded constants from lr_score.hpp + the dataset-
// specific lr_threshold knob:
//   - kP_real_uniform=0.75   → per-bin P(real mismatch | bin) fit from B1.
//   - kBundlePrior_alpha=0.10 → per-occ-bin log(π_pcr/π_real) fit from B1.
//   - lr_threshold=-9.0      → gone. Decision is "S > 0" by definition of
//                              posterior odds.
//
// The model is fit ONCE per run from candidate edges collected in B1, then
// queried in B2. Bins are coarse (subst_class × terminal-distance × occ ×
// damage-channel) so each cell sees enough mass to fit; sparse cells fall
// back to a Beta(α,β) shrinkage prior toward the global rate.
//
// Cross-bundle recurrence is the key evidence: if substitution X at position
// Y appears in N distinct bundles, that is N independent sightings of "this
// looks like a real variant" rather than N copies of one PCR family. P_real
// for a bin is fit from the bundle-distinct count, not raw mismatch count.

#include "fqdup/lr_score.hpp"  // phred_to_perr

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_set>
#include <vector>

namespace fqdup::errcor_emp {

// ── Bin axes ────────────────────────────────────────────────────────────────

// Substitution class collapsed by strand complement (6 classes):
//   0: C>T / G>A   (deamination channel for DS aDNA)
//   1: A>G / T>C   (transition, non-damage)
//   2: C>A / G>T   (oxidation-related transversion)
//   3: A>C / T>G   (transversion)
//   4: C>G / G>C   (transversion)
//   5: A>T / T>A   (transversion)
constexpr int kNumSubstClasses = 6;
constexpr int kNumTermDistBins = 4;   // d in [0..2], [3..5], [6..10], [11..]
constexpr int kNumOccBins      = 6;   // [1], [2..3], [4..7], [8..15], [16..63], [64..]
constexpr int kNumDamageBins   = 2;   // 0=not in damage zone / not damage channel, 1=both
constexpr int kNumBins =
    kNumSubstClasses * kNumTermDistBins * kNumOccBins * kNumDamageBins;

// 2-bit nucleotide encoding (A=0,C=1,G=2,T=3) → substitution class.
// Returns -1 if parent==alt (no mismatch) or invalid input.
inline int subst_class(uint8_t parent_2b, uint8_t alt_2b) {
    if (parent_2b == alt_2b) return -1;
    if (parent_2b > 3 || alt_2b > 3) return -1;
    // Canonicalize by complement: A<->T (0<->3), C<->G (1<->2).
    auto comp = [](uint8_t x) -> uint8_t { return static_cast<uint8_t>(3u - x); };
    uint8_t p = parent_2b, a = alt_2b;
    // Map each unordered complement-equivalent pair to one canonical pair.
    // Pick the pair where p is the "lower-by-canonical-key" of {p, comp(p)}.
    // Encoding table (canonical p,a → class):
    //   (C,T)=>0  (1,3) ; (G,A) maps via comp to (C,T)
    //   (A,G)=>1  (0,2) ; (T,C) maps via comp to (A,G)
    //   (C,A)=>2  (1,0) ; (G,T) maps via comp to (C,A)
    //   (A,C)=>3  (0,1) ; (T,G) maps via comp to (A,C)
    //   (C,G)=>4  (1,2) ; (G,C) maps via comp to (C,G)
    //   (A,T)=>5  (0,3) ; (T,A) maps via comp to (A,T)
    // To canonicalize: if p in {G,T} (i.e. p>=2) flip both via comp.
    if (p >= 2) { p = comp(p); a = comp(a); }
    // Now p in {A=0, C=1}. 8 valid (p,a) cells.
    static constexpr int kTable[2][4] = {
        // p=A: a=A,C,G,T
        { -1,  3,  1,  5 },
        // p=C: a=A,C,G,T
        {  2, -1,  4,  0 },
    };
    return kTable[p][a];
}

inline int term_dist_bin(int pos, int len) {
    if (len <= 0) return 0;
    int d = std::min(pos, len - 1 - pos);
    if (d < 0)  d = 0;
    if (d <= 2)  return 0;
    if (d <= 5)  return 1;
    if (d <= 10) return 2;
    return 3;
}

inline int occ_bin(uint64_t occ) {
    if (occ <= 1)  return 0;
    if (occ <= 3)  return 1;
    if (occ <= 7)  return 2;
    if (occ <= 15) return 3;
    if (occ <= 63) return 4;
    return 5;
}

inline int bin_index(int subst, int term, int occ, int dmg) {
    return ((subst * kNumTermDistBins + term) * kNumOccBins + occ) * kNumDamageBins + dmg;
}

// ── Feature records ─────────────────────────────────────────────────────────

struct EdgeMmFeatures {
    uint16_t pos;          // full-length position
    uint8_t  parent_2b;
    uint8_t  alt_2b;
    uint8_t  qual;
    uint8_t  damage_chan;  // 1 if (parent^alt)==2 AND in damage zone, else 0
    double   p_damage;     // p_damage_at(pos, len) when damage_chan, else 0
};

struct EdgeCandidate {
    uint32_t child_id;
    uint32_t parent_id;
    uint64_t bundle_key;   // raw bundle hash; used for distinct-bundle counting
    uint16_t child_len;
    uint64_t bundle_occ;   // occupancy at bundle_key
    uint8_t  n_mm;         // 1 or 2
    EdgeMmFeatures mm[2];
    // T8.6: indel channel for the rescue path. Defaults are zero so existing
    // (substitution-only) paths are unaffected. `indel_in_mask` is set when
    // any indel event lies within the damage-mask zone.
    uint8_t  n_ins        = 0;
    uint8_t  n_del        = 0;
    uint8_t  indel_in_mask = 0;
};

// T8.6: tunable α weights for the indel channel. Defaults match the GPT-5.5
// review: interior indels are heavily penalised (likely real); damage-mask
// indels (slippage near deaminated termini) get a lift.
struct IndelScoreParams {
    double alpha_ins = 2.0;
    double alpha_del = 2.0;
    double mask_bonus = 1.0;  // added to S when indel_in_mask is set
};

// ── Shrinkage prior ─────────────────────────────────────────────────────────
// Beta(α=1, β=9999): weak prior, mean 1e-4 — matches the per-base sequencing
// error rate at Phred-40 (10^-4) and the typical baseline rate at which a
// (subst x term x occ x dmg) bin signature recurs across DISTINCT bundles
// in a short-read library. With β=99 the prior was 100× too generous to
// "real" — every bin defaulted to p_real≈0.01 even with zero recurrence
// evidence, which gave per-mismatch log-odds of ≈-6 nats and stalled all
// absorption on small/synthetic inputs (errcor and pipeline-P8 tests).
// High-evidence bins still move toward their empirical rate; rare bins
// now sit near 1e-4, matching the noise floor.
constexpr double kBetaA = 1.0;
constexpr double kBetaB = 9999.0;

// Floor on p_real to keep log() well-defined and bound the score from blowing
// up on a near-empty bin. 1e-6 ⇒ at most +13.8 nats per mismatch from H_real.
constexpr double kPRealFloor = 1e-6;

// ── Cached log lookups (avoid per-mismatch pow/log in score()) ──────────────
// log(P(specific err base) / 1) = log(p_err / 3) for q in [0..60]. Singleton
// table built at first call, never freed. phred_to_perr clamps q to [1,60];
// q=0 maps to neutral (treated as q=1 for safety).
inline const std::array<double, 61>& log_p_err_over3_table() {
    static const std::array<double, 61> t = []() {
        std::array<double, 61> a{};
        for (int q = 0; q < 61; ++q) {
            double p_err = lr::phred_to_perr(static_cast<uint8_t>(q == 0 ? 1 : q));
            a[q] = std::log(p_err / 3.0);
        }
        return a;
    }();
    return t;
}
inline double log_p_err_over3(uint8_t q) {
    if (q > 60) q = 60;
    return log_p_err_over3_table()[q];
}

// ── Model ───────────────────────────────────────────────────────────────────

struct ErrCorEmpiricalModel {
    // Per-bin estimate of P(mismatch is "real-variant-like" | bin), defined
    // operationally as P(this bin's substitution recurs across ≥1 other
    // distinct bundle). Higher ⇒ stronger evidence the mismatch is real ⇒
    // less willing to absorb.
    std::array<double, kNumBins> p_real_bin{};
    // Cached log(p_real_bin[b]); populated at end of fit(). score() reads
    // this directly, eliminating per-mismatch std::log().
    std::array<double, kNumBins> log_p_real_bin{};
    double log_p_real_global = std::log(0.0001);

    // Per-occ-bin log(π_pcr / π_real), fit empirically from the proportion
    // of cross-bundle recurring vs singleton mismatches at each occupancy.
    std::array<double, kNumOccBins> log_pi_ratio_occ{};

    // Global fallbacks for empty/sparse bins.
    double p_real_global       = 0.0001;
    double log_pi_ratio_global = 0.0;

    bool fitted = false;

    // Cross-bundle accumulator (WIN 1): bundles_per_bin_[b] = the set of
    // distinct bundle_keys observed in bin b. fit_accumulate() streams edges
    // into it so the full edge vector never has to be materialized;
    // fit_finalize() consumes it. Sized kNumBins on first use.
    std::vector<std::unordered_set<uint64_t>> bundles_per_bin_;

    // Build from collected edges. `total_distinct_bundles` is the number of
    // distinct bundle_keys observed across all edges (used as the recurrence
    // denominator).
    void fit(const std::vector<EdgeCandidate>& edges,
             uint64_t total_distinct_bundles);

    // WIN 1 streaming fit. accumulate_edge() folds one edge's mismatches into
    // `acc` (kNumBins sets); it is const and takes a caller-owned accumulator
    // so parallel workers can each fill a private accumulator and then merge —
    // set-union is order-independent, so the finalized tables are bit-identical
    // to the single-pass fit().
    void accumulate_edge(std::vector<std::unordered_set<uint64_t>>& acc,
                         const EdgeCandidate& e) const;
    void fit_accumulate(const EdgeCandidate& e);
    void fit_merge(const std::vector<std::unordered_set<uint64_t>>& other);
    void fit_finalize(uint64_t total_distinct_bundles);

    // Posterior log-odds score for an edge.
    //   S = log P(obs|PCR) − log P(obs|real) + log π_pcr − log π_real
    // Absorb iff S > 0.
    double score(const EdgeCandidate& e) const;

    // T8.6 (provisional): same as score(), then adjusts for the indel channel.
    //   S' = S - α_ins·n_ins - α_del·n_del + (indel_in_mask ? mask_bonus : 0)
    // The substitution channel uses fitted bins; the indel channel is a
    // hand-set additive penalty until enough rescue edges accumulate to fit
    // p_ins / p_del per bin.
    double score_with_indel(const EdgeCandidate& e,
                            const IndelScoreParams& ip) const {
        double S = score(e);
        S -= ip.alpha_ins * static_cast<double>(e.n_ins);
        S -= ip.alpha_del * static_cast<double>(e.n_del);
        if (e.indel_in_mask) S += ip.mask_bonus;
        return S;
    }

    // JSON snapshot of the fitted model for run logs.
    std::string to_json() const;

private:
    double p_real_for(int subst, int term, int occ, int dmg) const {
        if (subst < 0) return p_real_global;
        double v = p_real_bin[bin_index(subst, term, occ, dmg)];
        if (!(v > 0.0)) v = p_real_global;     // unset / NaN → fallback
        if (v < kPRealFloor) v = kPRealFloor;
        if (v > 1.0 - kPRealFloor) v = 1.0 - kPRealFloor;
        return v;
    }
    double log_pi_for(int occ) const {
        if (!fitted) return log_pi_ratio_global;
        return log_pi_ratio_occ[occ];
    }
};

inline double ErrCorEmpiricalModel::score(const EdgeCandidate& e) const {
    const int o_bin = occ_bin(e.bundle_occ);
    double S = log_pi_for(o_bin);
    for (int i = 0; i < e.n_mm; ++i) {
        const auto& m = e.mm[i];
        int s = subst_class(m.parent_2b, m.alt_2b);
        int t = term_dist_bin(m.pos, e.child_len);
        int d = m.damage_chan ? 1 : 0;
        // log(p_err / 3) — table lookup, no pow()/log().
        const double log_p_pcr = log_p_err_over3(m.qual);
        // log(p_real) — cached unless damage lift applies (rare path).
        double log_p_real;
        if (m.damage_chan && m.p_damage > 0.0) {
            double p_real = p_real_for(s, t, o_bin, d);
            double pd = m.p_damage > 1.0 ? 1.0 : m.p_damage;
            p_real = p_real + pd * (1.0 - p_real);
            log_p_real = std::log(p_real);
        } else if (s < 0) {
            log_p_real = log_p_real_global;
        } else {
            log_p_real = fitted ? log_p_real_bin[bin_index(s, t, o_bin, d)]
                                : log_p_real_global;
        }
        S += log_p_pcr - log_p_real;
    }
    return S;
}

// ── Fit ─────────────────────────────────────────────────────────────────────
//
// For each bin b (subst × term × occ × dmg):
//   n_obs[b]      = total mismatch observations falling into b
//   n_bundles[b]  = number of DISTINCT bundle_keys in which b was observed
// Then:
//   p_real_bin[b] = (n_bundles[b] + α) / (total_distinct_bundles + α + β)
//
// For the per-occupancy prior:
//   For each edge, count whether ANY of its mismatches' (subst×term) cell is
//   "recurring" (seen in ≥2 distinct bundles globally) vs "singleton".
//   log_pi_ratio_occ[o] = log( (n_singleton[o]+α) / (n_recurring[o]+α) )
// Singleton-dominated occupancies favor PCR; recurring-dominated favor real.
// ── Per-bin P(real | bin) from cross-bundle recurrence ─────────────────────
// For each bin b, count the number of DISTINCT bundles in which a mismatch
// falling in that bin was observed. A signature that recurs across many
// independent bundles is evidence of a systematic real pattern (damage,
// sequencing bias, biological variant); a signature confined to few bundles is
// more consistent with random PCR error.
//   p_real_bin[b] = (n_bundles_b + α) / (total_distinct_bundles + α + β)
// Always in (0, 1) because n_bundles_b ≤ total_distinct_bundles.
inline void ErrCorEmpiricalModel::accumulate_edge(
        std::vector<std::unordered_set<uint64_t>>& acc,
        const EdgeCandidate& e) const {
    for (int i = 0; i < e.n_mm; ++i) {
        const auto& m = e.mm[i];
        int s = subst_class(m.parent_2b, m.alt_2b);
        if (s < 0) continue;
        int t = term_dist_bin(m.pos, e.child_len);
        int o = occ_bin(e.bundle_occ);
        int d = m.damage_chan ? 1 : 0;
        acc[bin_index(s, t, o, d)].insert(e.bundle_key);
    }
}

inline void ErrCorEmpiricalModel::fit_accumulate(const EdgeCandidate& e) {
    if (bundles_per_bin_.size() != static_cast<size_t>(kNumBins))
        bundles_per_bin_.assign(kNumBins, {});
    accumulate_edge(bundles_per_bin_, e);
}

inline void ErrCorEmpiricalModel::fit_merge(
        const std::vector<std::unordered_set<uint64_t>>& other) {
    if (bundles_per_bin_.size() != static_cast<size_t>(kNumBins))
        bundles_per_bin_.assign(kNumBins, {});
    const int n = std::min<int>(kNumBins, static_cast<int>(other.size()));
    for (int b = 0; b < n; ++b)
        bundles_per_bin_[b].insert(other[b].begin(), other[b].end());
}

inline void ErrCorEmpiricalModel::fit(const std::vector<EdgeCandidate>& edges,
                                      uint64_t total_distinct_bundles) {
    bundles_per_bin_.assign(kNumBins, {});
    for (const auto& e : edges) accumulate_edge(bundles_per_bin_, e);
    fit_finalize(total_distinct_bundles);
}

inline void ErrCorEmpiricalModel::fit_finalize(uint64_t total_distinct_bundles) {
    auto& bundles_per_bin = bundles_per_bin_;
    // Adaptive shrinkage: scale β with total_distinct_bundles so the Beta
    // prior dominates on small inputs. With β fixed at 9999, a single
    // observation in a 1-bundle dataset would push p_real_bin toward
    // 2/10001 ≈ 2e-4 — an order of magnitude above the true PCR rate.
    // Scaling β as max(9999, 100·N) keeps the prior dominant unless we have
    // hundreds of independent bundles supporting the bin.
    const double beta_eff =
        std::max(kBetaB, 100.0 * static_cast<double>(total_distinct_bundles));
    double denom = static_cast<double>(total_distinct_bundles) + kBetaA + beta_eff;
    if (denom <= 0) denom = 1.0;
    for (int b = 0; b < kNumBins; ++b) {
        double nb = static_cast<double>(bundles_per_bin[b].size());
        double v = (nb + kBetaA) / denom;
        if (v < kPRealFloor)        v = kPRealFloor;
        if (v > 1.0 - kPRealFloor)  v = 1.0 - kPRealFloor;
        p_real_bin[b] = v;
    }
    // Global fallback = pure Beta prior mean (used only for unset bins).
    p_real_global = kBetaA / (kBetaA + beta_eff);

    // ── Per-occ-bin log(π_pcr / π_real) from within-bundle recurrence ──────
    // PCR amplification produces FAMILIES of children that share the same
    // mismatch signature relative to their common parent (the error
    // happened once during amplification, then was copied). A parent whose
    // children's edges have many duplicate (mm_pos, alt_base) signatures
    // is exhibiting a PCR family; one whose edges have all-unique
    // signatures is more likely seeing independent real variation.
    // log_pi_ratio_occ[o] = log( (recurring + α) / (unique + α) )
    // Positive ⇒ favor PCR ⇒ higher absorption at this occupancy.
    //
    // The within-parent recurrence test is a poor PCR signal: PCR errors
    // randomise position across copies, so a 30x parent with 30 H=1 children
    // typically has all-singleton signatures (sig_count=1 each) — which the
    // recurrence test misreads as "favor real" and crushes absorption.
    //
    // Use a fixed prior log(π_pcr/π_real) reflecting the operational
    // baseline: in any near-clone candidate edge from amplified short-read
    // data, the PCR-error hypothesis is at least an order of magnitude more
    // likely a-priori than the real-variant hypothesis (real biological
    // variation between siblings of the same molecule is essentially zero;
    // sequencing errors and PCR errors dominate at ~10⁻³–10⁻⁴ per base).
    // log(10) ≈ 2.30 nats — modest enough that strongly-recurring real
    // signatures (large p_real_bin) still flip the decision to "protect".
    constexpr double kLogPiPriorPCR = 3.00;
    log_pi_ratio_global = kLogPiPriorPCR;
    for (int o = 0; o < kNumOccBins; ++o) log_pi_ratio_occ[o] = kLogPiPriorPCR;

    fitted = true;

    // Cache log(p_real_bin[b]) so score() doesn't call std::log() per mismatch.
    for (int b = 0; b < kNumBins; ++b) {
        double v = p_real_bin[b];
        if (!(v > 0.0)) v = p_real_global;
        if (v < kPRealFloor)        v = kPRealFloor;
        if (v > 1.0 - kPRealFloor)  v = 1.0 - kPRealFloor;
        log_p_real_bin[b] = std::log(v);
    }
    log_p_real_global = std::log(p_real_global > 0.0 ? p_real_global : kPRealFloor);
}

inline std::string ErrCorEmpiricalModel::to_json() const {
    char buf[64];
    std::string s = "{\"fitted\":";
    s += fitted ? "true" : "false";
    auto fmt = [&](double v) {
        std::snprintf(buf, sizeof(buf), "%.6g", v);
        return std::string(buf);
    };
    s += ",\"p_real_global\":" + fmt(p_real_global);
    s += ",\"log_pi_ratio_global\":" + fmt(log_pi_ratio_global);
    s += ",\"log_pi_ratio_occ\":[";
    for (int o = 0; o < kNumOccBins; ++o) {
        if (o) s += ",";
        s += fmt(log_pi_ratio_occ[o]);
    }
    s += "],\"p_real_bin\":[";
    for (int b = 0; b < kNumBins; ++b) {
        if (b) s += ",";
        s += fmt(p_real_bin[b]);
    }
    s += "]}";
    return s;
}

} // namespace fqdup::errcor_emp
