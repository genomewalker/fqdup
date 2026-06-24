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

// ---- Post-pass1: reweight OxoG sig histograms with fitted dp -----------
// Replaces the reader2 FASTQ loop. Numerically identical to the per-read
// pass for first moments; second moments are exact via stored Σ T_i^2 etc.
inline static void oxog_from_sig_hist(
    taph::OxBinAcc* ox_bins,
    int n_gc,
    const std::vector<OxSigCell>& merged_fwd,
    const std::vector<OxSigCell>& merged_rev,
    const LLRTable* llr_tables,
    bool is_ss)
{
    static constexpr int N_SIG = N_SIG10;
    char buf[11] = {};
    auto apply = [&](const std::vector<OxSigCell>& hist) {
        for (int gc = 0; gc < n_gc; ++gc) {
            const LLRTable& tbl = llr_tables[gc];
            for (int sig = 0; sig < N_SIG; ++sig) {
                const auto& cell = hist[gc * N_SIG + sig];
                if (cell.n == 0) continue;
                uint32_t s = static_cast<uint32_t>(sig);
                for (int i = 0; i < 10; ++i) {
                    const int v = s % 3; s /= 3;
                    buf[i] = (v == 1) ? 'T' : (v == 2) ? 'C' : 'A';
                }
                // ceiling: LLRTable has up to 15 positions; sig encodes 10.
                // Positions 10-14 contribute ~0.5% to q for typical lambda values.
                // 3^15=18 GB per thread — extending the sig is rejected.
                const double q  = taph::llr_table_weight_fwd(buf, 10, tbl);
                const double qm = 1.0 - q;
                auto& x = ox_bins[gc];
                const double a  = q  * cell.t,  b_ = q  * (cell.t + cell.g);
                const double c  = qm * cell.t,  d  = qm * (cell.t + cell.g);
                x.A += a; x.B += b_; x.C += c; x.D += d;
                x.AA += q*q   * cell.t2; x.AB += q*q   * cell.tm; x.BB += q*q   * cell.m2;
                x.CC += qm*qm * cell.t2; x.CD += qm*qm * cell.tm; x.DD += qm*qm * cell.m2;
                x.AC += q*qm  * cell.t2; x.AD += q*qm  * cell.tm;
                x.BC += q*qm  * cell.tm; x.BD += q*qm  * cell.m2;
                x.M_total += cell.t + cell.g; x.reads += cell.n;
                x.T_all   += cell.t; x.TG_all += cell.t + cell.g;
                x.C_all   += cell.cv; x.CA_all += cell.cv + cell.av;
                x.A_ca += q  * cell.av; x.B_ca += q  * (cell.cv + cell.av);
                x.C_ca += qm * cell.av; x.D_ca += qm * (cell.cv + cell.av);
            }
        }
    };
    apply(merged_fwd);
    if (!is_ss) apply(merged_rev);
}

// ---- LSD single-pass: reconstruct LsdLlrBinAccum from sig5 histograms ---
// lsd_llr_from_sig: LLR for a read encoded as 5' ternary (T=1,C=2,other=0)
// and 3' quinary (A=1,G=2,T=3,C=4,other=0). Uses the same decay model as
// lsd_llr_score but without per-read CpG context (scale fixed at 1.0).
// ceiling: cpg_scale/noncpg_scale not applied — needs seq[p+1] which the sig
//          does not encode; upgrade: extend sig with next-base ternary context.
inline static float lsd_llr_from_sig(
    uint32_t sig5f,
    uint32_t sig5r_4,
    const LsdClassifyParams& p,
    bool is_ss)
{
    float llr = 0.0f;
    // 5' end: advance past skipped position 0, then decode remaining trits
    uint32_t s5 = sig5f;
    if (p.skip_pos0_5prime) s5 /= 3;
    for (int i = p.skip_pos0_5prime ? 1 : 0; i < 5; ++i) {
        const int tv = s5 % 3; s5 /= 3;
        const float d = static_cast<float>(p.d_anc * std::exp(-p.lam_5 * i));
        if (d < 1e-4f) break;
        const float bg = static_cast<float>(p.bg_5);
        const float pd = std::clamp(bg + (1.0f - bg) * d, 1e-6f, 1.0f - 1e-6f);
        if      (tv == 1) llr += std::log(pd / bg);
        else if (tv == 2) llr += std::log((1.0f - pd) / (1.0f - bg));
    }
    // 3' end
    uint32_t s3 = sig5r_4;
    for (int i = 0; i < 5; ++i) {
        const int rv = s3 % 5; s3 /= 5;
        const float d = static_cast<float>(p.d_anc * std::exp(-p.lam_3 * i));
        if (d < 1e-4f) break;
        const float bg = static_cast<float>(p.bg_3);
        const float pd = std::clamp(bg + (1.0f - bg) * d, 1e-6f, 1.0f - 1e-6f);
        if (!is_ss) {  // DS: A=hit (rv==1), G=eligible (rv==2)
            if      (rv == 1) llr += std::log(pd / bg);
            else if (rv == 2) llr += std::log((1.0f - pd) / (1.0f - bg));
        } else {        // SS: T=hit (rv==3), C=eligible (rv==4)
            if      (rv == 3) llr += std::log(pd / bg);
            else if (rv == 4) llr += std::log((1.0f - pd) / (1.0f - bg));
        }
    }
    return llr;
}

// Reconstruct LsdLlrBinAccum for n_bins bins from per-(bin,sig) counts.
// Positions 0..4 are exact; positions 5..14 are approximated from merged_lbs
// aggregate counts scaled by the per-bin damaged fraction.
// Shadow counter output: sh_n/sh_old/sh_new/sh_flip filled when cls.d_anc_contract >= 0.
inline static std::vector<LsdLlrBinAccum> reconstruct_lsd_llr_accum(
    const std::vector<int32_t>& merged_cnt,  // [n_bins × N_LSD_SIG5]
    int n_bins,
    const taph::LengthBinStats& merged_lbs,
    const LsdClassifyParams& cls,
    bool is_ss,
    double lsd_log_prior_odds,
    uint64_t* sh_n, uint64_t* sh_old, uint64_t* sh_new, uint64_t* sh_flip)
{
    std::vector<LsdLlrBinAccum> out(n_bins);
    const bool do_shadow = (sh_n != nullptr) && (cls.d_anc_contract >= 0.0);
    LsdClassifyParams contract_cls = cls;
    contract_cls.d_anc = cls.d_anc_contract;

    for (int b = 0; b < n_bins; ++b) {
        auto& acc = out[b];
        const int32_t* cnt = merged_cnt.data() + b * N_LSD_SIG5;
        for (int sig = 0; sig < N_LSD_SIG5; ++sig) {
            const int32_t n = cnt[sig];
            if (n == 0) continue;
            const uint32_t sig5f   = static_cast<uint32_t>(sig) / N_LSD_SIG5_3R;
            const uint32_t sig5r_4 = static_cast<uint32_t>(sig) % N_LSD_SIG5_3R;
            const float llr_val = lsd_llr_from_sig(sig5f, sig5r_4, cls, is_ss);
            const bool  anc     = llr_val > 0.0f;
            const double post   = 1.0 / (1.0 + std::exp(-(static_cast<double>(llr_val) + lsd_log_prior_odds)));
            if (anc) acc.n_damaged   += n; else acc.n_undamaged += n;
            acc.sw_sum += post * n;
            if (do_shadow && !cls.contract_gated_off) {
                const bool anc_c = lsd_llr_from_sig(sig5f, sig5r_4, contract_cls, is_ss) > 0.0f;
                *sh_n += n;
                if (anc)          *sh_old += n;
                if (anc_c)        *sh_new += n;
                if (anc != anc_c) *sh_flip += n;
            } else if (do_shadow) {
                *sh_n += n;
                if (anc) *sh_old += n;
            }
            // 5' positions 0..4
            uint32_t s5 = sig5f;
            for (int p = 0; p < 5; ++p) {
                const int fv = s5 % 3; s5 /= 3;
                const bool is_t = (fv == 1), is_tc = (fv != 0);
                if (anc) {
                    acc.t_5_anc[p]  += (int64_t)n * is_t;
                    acc.tc_5_anc[p] += (int64_t)n * is_tc;
                } else {
                    acc.t_5_mod[p]  += (int64_t)n * is_t;
                    acc.tc_5_mod[p] += (int64_t)n * is_tc;
                }
                if (p < LsdLlrBinAccum::N_SOFT_POS) {
                    acc.sw_t5_anc[p]  += post * n * is_t;
                    acc.sw_tc5_anc[p] += post * n * is_tc;
                }
            }
            // 3' positions 0..4
            uint32_t s3 = sig5r_4;
            for (int p = 0; p < 5; ++p) {
                const int rv = s3 % 5; s3 /= 5;
                const bool hit  = is_ss ? (rv == 3) : (rv == 1);
                const bool elig = is_ss ? (rv == 3 || rv == 4) : (rv == 1 || rv == 2);
                if (anc) {
                    acc.h_3_anc[p] += (int64_t)n * hit;
                    acc.n_3_anc[p] += (int64_t)n * elig;
                } else {
                    acc.h_3_mod[p] += (int64_t)n * hit;
                    acc.n_3_mod[p] += (int64_t)n * elig;
                }
                if (p < LsdLlrBinAccum::N_SOFT_POS) {
                    acc.sw_h3_anc[p] += post * n * hit;
                    acc.sw_n3_anc[p] += post * n * elig;
                }
            }
        }
        // Positions 5..14: approximate from aggregate lbs × bulk damaged fraction.
        const int64_t n_total = acc.n_damaged + acc.n_undamaged;
        const int bsz = static_cast<int>(merged_lbs.profiles.size());
        if (n_total > 0 && b < bsz) {
            const auto& lp = merged_lbs.profiles[static_cast<size_t>(b)];
            const double fa = static_cast<double>(acc.n_damaged) / n_total;
            const double fm = 1.0 - fa;
            for (int p = 5; p < LengthBinDamageProfile::N_POS; ++p) {
                acc.t_5_anc[p]  = static_cast<int64_t>(lp.t_freq_5prime[p]   * fa);
                acc.tc_5_anc[p] = static_cast<int64_t>(lp.tc_total_5prime[p] * fa);
                acc.t_5_mod[p]  = static_cast<int64_t>(lp.t_freq_5prime[p]   * fm);
                acc.tc_5_mod[p] = static_cast<int64_t>(lp.tc_total_5prime[p] * fm);
                if (!is_ss) {
                    const double ag = lp.a_freq_3prime[p] + lp.g_freq_3prime[p];
                    acc.h_3_anc[p] = static_cast<int64_t>(lp.a_freq_3prime[p] * fa);
                    acc.n_3_anc[p] = static_cast<int64_t>(ag * fa);
                    acc.h_3_mod[p] = static_cast<int64_t>(lp.a_freq_3prime[p] * fm);
                    acc.n_3_mod[p] = static_cast<int64_t>(ag * fm);
                } else {
                    acc.h_3_anc[p] = static_cast<int64_t>(lp.t_freq_3prime[p]    * fa);
                    acc.n_3_anc[p] = static_cast<int64_t>(lp.tc_total_3prime[p]  * fa);
                    acc.h_3_mod[p] = static_cast<int64_t>(lp.t_freq_3prime[p]    * fm);
                    acc.n_3_mod[p] = static_cast<int64_t>(lp.tc_total_3prime[p]  * fm);
                }
            }
        }
    }
    return out;
}

