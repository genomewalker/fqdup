#pragma once
// Oxidation score worker for fqdup damage second pass.
#include "damage_worker.hpp"
#include "taph/sample_damage_profile.hpp"
#include "taph/llr_table.hpp"
#include "taph/oxog_score.hpp"
#include "fqdup/damage_profile.hpp"
#include <array>
#include <cmath>
#include <string>
#include <vector>
// ---- second pass: oxidation-compatible composition score (S_DS / S_SS) -----
//
// For each read, compute q = P(ancient|read) via 5' C→T LLR from pass-1 priors.
// For DS, dual-orient: score both seq and rc(seq), use the orientation with higher q.
// After orientation, accumulate weighted T/(T+G) in the middle third by GC bin:
//   p_anc = Σ(q·T) / Σ(q·M)      (ancient-enriched)
//   p_bg  = Σ((1-q)·T) / Σ((1-q)·M)  (background)
//   S = Σ_b α_b (p_anc,b − p_bg,b)    (GC-weighted across bins)
//
// Unlike D = T/(T+G) − A/(A+C), S does not cancel in DS because both strands
// contribute to the same T/(T+G) numerator after orientation.

inline static float per_read_damaged_weight(
    const std::string& seq,
    const taph::SampleDamageProfile& dp)
{
    int L = static_cast<int>(seq.size());
    int bin = taph::SampleDamageProfile::get_gc_bin(seq);
    const auto& gc = dp.gc_bins[bin];
    float bg   = (gc.valid && gc.baseline_tc > 0.0f) ? gc.baseline_tc : dp.fit_baseline_5prime;
    bg = std::clamp(bg, 1e-6f, 1.0f - 1e-6f);   // guard log(pd/bg) from -inf on sparse GC bins
    float dmax = (gc.valid && gc.d_max > 0.005f)     ? gc.d_max       : dp.d_max_5prime;
    float lam  = dp.lambda_5prime;
    float pi   = (gc.valid && gc.p_damaged > 0.0f)   ? gc.p_damaged   : dp.pi_damaged;
    pi = std::clamp(pi, 0.001f, 0.999f);
    float llr = 0.0f;
    for (int p = 0; p < L; ++p) {
        float expected_excess = dmax * std::exp(-lam * static_cast<float>(p));
        if (expected_excess < 0.01f) break;
        char b = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[p])));
        float pd = std::clamp(bg + (1.0f - bg) * expected_excess, 1e-6f, 1.0f - 1e-6f);
        if      (b == 'T') llr += std::log(pd / bg);
        else if (b == 'C') llr += std::log((1.0f - pd) / (1.0f - bg));
    }
    float logit_pi = std::log(pi / (1.0f - pi));
    return 1.0f / (1.0f + std::exp(-(logit_pi + llr)));
}

static constexpr int N_GC = taph::SampleDamageProfile::N_GC_BINS;
using LLRTable = taph::LLRTable;

// OxBinAcc, update_ox_bins, merge_ox_bins, compute_ox_scores live in libtaph/oxog_score.hpp
#include "taph/oxog_score.hpp"
using taph::OxBinAcc;

struct OxogWorkerState {
    OxBinAcc bins[N_GC];
    // Fix E: LSD data accumulated when fuse_lsd=true in oxog_worker.
    taph::LengthBinStats        lbs;
    std::vector<LsdLlrBinAccum> llr_acc;  // per length-bin, size = n_bins
    // Step-2 shadow counters (SOLUTION §6.7): contract per-read class vs current d_anc class. Merged and
    // logged in damage.cpp; never drives a verdict.
    uint64_t shadow_n = 0, shadow_anc_old = 0, shadow_anc_new = 0, shadow_flip = 0;
};

inline static void oxog_worker(WorkQueue& queue,
                        OxogWorkerState& state,
                        const taph::SampleDamageProfile& dp,
                        const LLRTable* llr,
                        bool is_ss,
                        bool fuse_lsd,
                        const LsdClassifyParams* cls_params,
                        const std::vector<int>* lsd_edges,
                        double lsd_log_prior_odds)
{
    const int n_lsd_bins = fuse_lsd && lsd_edges
        ? 1 + static_cast<int>(lsd_edges->size()) : 0;

    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (const std::string& seq : batch) {
            int L = static_cast<int>(seq.size());
            if (L < 30) continue;
            int bin = taph::SampleDamageProfile::get_gc_bin(seq);
            const LLRTable& t = llr[bin];
            double q;
            bool rev = false;
            if (is_ss) {
                q = taph::llr_table_weight_fwd(seq.data(), static_cast<int>(seq.size()), t);
            } else {
                float w_fwd = taph::llr_table_weight_fwd(seq.data(), static_cast<int>(seq.size()), t);
                float w_rev = taph::llr_table_weight_rev(seq.data(), static_cast<int>(seq.size()), t);
                if (w_rev > w_fwd) { rev = true; q = w_rev; }
                else               {              q = w_fwd; }
            }
            int beg = L / 3, end = L - (L / 3);
            if (end <= beg) continue;
            int T = 0, G = 0, Cv = 0, Av = 0;
            if (!rev) {
                for (int j = beg; j < end; ++j) {
                    char b = static_cast<char>(seq[j] & ~0x20u);
                    if      (b == 'T') ++T;
                    else if (b == 'G') ++G;
                    else if (b == 'C') ++Cv;
                    else if (b == 'A') ++Av;
                }
            } else {
                // DS read oriented in reverse: complement-map the middle third
                // without allocating an rc string (A→T, C→G, G→C, T→A).
                for (int j = beg; j < end; ++j) {
                    char b = static_cast<char>(seq[j] & ~0x20u);
                    if      (b == 'A') ++T;
                    else if (b == 'C') ++G;
                    else if (b == 'G') ++Cv;
                    else if (b == 'T') ++Av;
                }
            }
            int M = T + G;
            if (M == 0) continue;
            taph::update_ox_bins(state.bins[bin], q, T, G, Cv, Av);

            // Fix E: fused LSD accumulation.
            if (fuse_lsd && L >= LSD_L_MIN && cls_params && n_lsd_bins > 0) {
                int L_binned = std::min(L, LSD_L_MAX);
                state.lbs.update(seq, L_binned);
                // bin_index: count how many edges L_binned is >= 
                int lb = 0;
                for (int e : *lsd_edges) { if (L_binned >= e) ++lb; else break; }
                if (lb >= 0 && lb < n_lsd_bins) {
                    double llr_val = lsd_llr_score(seq, *cls_params);
                    bool   anc     = llr_val > 0.0;
                    // Step-2 shadow (SOLUTION §6.7): contract class = LLR under the cohort amplitude
                    // (D_MAX_CONSERVED), gated off when the library is not pi-DETECTED. Counted only.
                    if (cls_params->d_anc_contract >= 0.0) {
                        bool anc_new = false;
                        if (!cls_params->contract_gated_off) {
                            LsdClassifyParams qc = *cls_params;
                            qc.d_anc = cls_params->d_anc_contract;
                            anc_new = lsd_llr_score(seq, qc) > 0.0;
                        }
                        ++state.shadow_n;
                        if (anc)               ++state.shadow_anc_old;
                        if (anc_new)           ++state.shadow_anc_new;
                        if (anc != anc_new)    ++state.shadow_flip;
                    }
                    if (anc) ++state.llr_acc[lb].n_damaged;
                    else     ++state.llr_acc[lb].n_undamaged;
                    lsd_accumulate(seq, state.llr_acc[lb], anc, cls_params->is_ss);
                    double post = 1.0 / (1.0 + std::exp(-(llr_val + lsd_log_prior_odds)));
                    lsd_accumulate_soft(seq, state.llr_acc[lb], post, cls_params->is_ss);
                }
            }
        }
    }
}

