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
};

// ── Shrinkage prior ─────────────────────────────────────────────────────────
// Beta(α=1, β=99): weak prior, mean 0.01 — i.e. "by default we expect ~1% of
// mismatches in a bin to recur across distinct bundles". High-evidence bins
// move toward their empirical rate; rare bins stay near the prior.
constexpr double kBetaA = 1.0;
constexpr double kBetaB = 99.0;

// Floor on p_real to keep log() well-defined and bound the score from blowing
// up on a near-empty bin. 1e-6 ⇒ at most +13.8 nats per mismatch from H_real.
constexpr double kPRealFloor = 1e-6;

// ── Model ───────────────────────────────────────────────────────────────────

struct ErrCorEmpiricalModel {
    // Per-bin estimate of P(mismatch is "real-variant-like" | bin), defined
    // operationally as P(this bin's substitution recurs across ≥1 other
    // distinct bundle). Higher ⇒ stronger evidence the mismatch is real ⇒
    // less willing to absorb.
    std::array<double, kNumBins> p_real_bin{};

    // Per-occ-bin log(π_pcr / π_real), fit empirically from the proportion
    // of cross-bundle recurring vs singleton mismatches at each occupancy.
    std::array<double, kNumOccBins> log_pi_ratio_occ{};

    // Global fallbacks for empty/sparse bins.
    double p_real_global       = 0.01;
    double log_pi_ratio_global = 0.0;

    bool fitted = false;

    // Build from collected edges. `total_distinct_bundles` is the number of
    // distinct bundle_keys observed across all edges (used as the recurrence
    // denominator).
    void fit(const std::vector<EdgeCandidate>& edges,
             uint64_t total_distinct_bundles);

    // Posterior log-odds score for an edge.
    //   S = log P(obs|PCR) − log P(obs|real) + log π_pcr − log π_real
    // Absorb iff S > 0.
    double score(const EdgeCandidate& e) const;

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
    double S = log_pi_for(occ_bin(e.bundle_occ));
    for (int i = 0; i < e.n_mm; ++i) {
        const auto& m = e.mm[i];
        int s = subst_class(m.parent_2b, m.alt_2b);
        int t = term_dist_bin(m.pos, e.child_len);
        int o = occ_bin(e.bundle_occ);
        int d = m.damage_chan ? 1 : 0;
        double p_err  = lr::phred_to_perr(m.qual);
        double p_pcr  = p_err / 3.0;
        double p_real = p_real_for(s, t, o, d);
        // Damage-channel mismatches near terminal positions get a bounded
        // lift: a fraction p_damage of their mass moves from H_pcr to H_real.
        if (m.damage_chan && m.p_damage > 0.0) {
            double pd = m.p_damage;
            if (pd > 1.0) pd = 1.0;
            p_real = p_real + pd * (1.0 - p_real);
        }
        S += std::log(p_pcr) - std::log(p_real);
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
inline void ErrCorEmpiricalModel::fit(const std::vector<EdgeCandidate>& edges,
                                      uint64_t total_distinct_bundles) {
    std::array<uint64_t, kNumBins> n_obs{};
    std::vector<std::unordered_set<uint64_t>> bundles_per_bin(kNumBins);

    // Reduced-cell recurrence (subst × term) — used for the per-occ prior.
    constexpr int kRedCells = kNumSubstClasses * kNumTermDistBins;
    auto red_idx = [](int s, int t) { return s * kNumTermDistBins + t; };
    std::vector<std::unordered_set<uint64_t>> bundles_per_red(kRedCells);

    // Pass 1: accumulate.
    for (const auto& e : edges) {
        for (int i = 0; i < e.n_mm; ++i) {
            const auto& m = e.mm[i];
            int s = subst_class(m.parent_2b, m.alt_2b);
            if (s < 0) continue;
            int t = term_dist_bin(m.pos, e.child_len);
            int o = occ_bin(e.bundle_occ);
            int d = m.damage_chan ? 1 : 0;
            int b = bin_index(s, t, o, d);
            n_obs[b]++;
            bundles_per_bin[b].insert(e.bundle_key);
            bundles_per_red[red_idx(s, t)].insert(e.bundle_key);
        }
    }

    // Global rate from total recurrences across all bins.
    uint64_t total_n_obs = 0, total_n_bundles_sum = 0;
    for (int b = 0; b < kNumBins; ++b) {
        total_n_obs += n_obs[b];
        total_n_bundles_sum += bundles_per_bin[b].size();
    }
    double denom = static_cast<double>(total_distinct_bundles) + kBetaA + kBetaB;
    if (denom <= 0) denom = 1.0;
    p_real_global = (static_cast<double>(total_n_bundles_sum) + kBetaA) / denom;
    if (p_real_global < kPRealFloor) p_real_global = kPRealFloor;

    // Per-bin shrunk rate.
    for (int b = 0; b < kNumBins; ++b) {
        double nb = static_cast<double>(bundles_per_bin[b].size());
        p_real_bin[b] = (nb + kBetaA) / denom;
    }

    // Per-occ-bin π_pcr/π_real from singleton-vs-recurring split.
    std::array<uint64_t, kNumOccBins> n_singleton{}, n_recurring{};
    for (const auto& e : edges) {
        bool any_recurring = false;
        for (int i = 0; i < e.n_mm; ++i) {
            const auto& m = e.mm[i];
            int s = subst_class(m.parent_2b, m.alt_2b);
            if (s < 0) continue;
            int t = term_dist_bin(m.pos, e.child_len);
            if (bundles_per_red[red_idx(s, t)].size() >= 2) {
                any_recurring = true;
                break;
            }
        }
        int o = occ_bin(e.bundle_occ);
        if (any_recurring) n_recurring[o]++;
        else               n_singleton[o]++;
    }
    uint64_t tot_sing = 0, tot_rec = 0;
    for (int o = 0; o < kNumOccBins; ++o) { tot_sing += n_singleton[o]; tot_rec += n_recurring[o]; }
    log_pi_ratio_global = std::log((static_cast<double>(tot_sing) + kBetaA) /
                                    (static_cast<double>(tot_rec)  + kBetaA));
    for (int o = 0; o < kNumOccBins; ++o) {
        double s_o = static_cast<double>(n_singleton[o]) + kBetaA;
        double r_o = static_cast<double>(n_recurring[o]) + kBetaA;
        log_pi_ratio_occ[o] = std::log(s_o / r_o);
    }

    fitted = true;
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
