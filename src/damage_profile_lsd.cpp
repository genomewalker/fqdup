// LSD helpers: compute_lsd_edges, make_lsd_classify_params, lsd_accumulate*.
// damage_profile.cpp
// estimate_damage() — deamination + QC full pipeline damage estimation.
// Called by both 'derep' (Pass 0) and 'extend' (Pass 0).

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"

#include "taph/frame_selector_decl.hpp"
#include "taph/length_gc_joint_mixture.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/log_length_gmm.hpp"
// taph/sample_damage_profile.hpp is included transitively via damage_profile.hpp

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

// ---- helpers for Fix-E LSD fusion (also used by damage.cpp) -----------

std::vector<int> compute_lsd_edges(const std::vector<uint64_t>& hist,
                                   const LengthBinOptions& options)
{
    const double log_min = std::log(static_cast<double>(LSD_L_MIN));
    const double log_max = std::log(static_cast<double>(LSD_L_MAX + 1));
    std::vector<int> edges;
    if (options.mode == LengthBinOptions::Mode::EXPLICIT) {
        for (int e : options.explicit_edges)
            if (e > LSD_L_MIN && e < LSD_L_MAX) edges.push_back(e);
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    } else if (options.mode == LengthBinOptions::Mode::QUANTILE) {
        int nb = std::clamp(options.quantile_bins, 1,
                            static_cast<int>(taph::LengthBinStats::MAX_BINS));
        if (nb > 1)
            edges = taph::detect_quantile_length_edges(hist, log_min, log_max,
                                                       LSD_L_MIN, LSD_L_MAX, nb);
    } else if (options.mode == LengthBinOptions::Mode::AUTO) {
        taph::LogLengthGmmResult gr = taph::detect_log_length_gmm_edges(
            hist, log_min, log_max, LSD_L_MIN, LSD_L_MAX,
            static_cast<int>(taph::LengthBinStats::MAX_BINS));
        edges = gr.edges;
    }
    const int max_edges = static_cast<int>(taph::LengthBinStats::MAX_BINS) - 1;
    if (static_cast<int>(edges.size()) > max_edges)
        edges.resize(max_edges);
    return edges;
}

LsdClassifyParams make_lsd_classify_params(const DamageProfile& bulk)
{
    LsdClassifyParams p;
    p.d_anc = (bulk.mixture_converged && bulk.mixture_d_damaged > 0.01)
        ? bulk.mixture_d_damaged
        : std::max(bulk.d_max_5prime, bulk.d_max_3prime);
    p.lam_5 = bulk.lambda_5prime;
    p.lam_3 = bulk.lambda_3prime;
    p.bg_5  = bulk.bg_5_tc;
    p.bg_3  = bulk.bg_3_channel;
    p.is_ss = bulk.ss_mode;
    const double d_cpg = bulk.d_cpg_5prime;
    const double d_ncp = bulk.d_noncpg_5prime;
    if (d_cpg > 0.01 && d_ncp > 0.01 && !std::isnan(d_cpg) && !std::isnan(d_ncp)) {
        const double d_avg = 0.5 * (d_cpg + d_ncp);
        if (d_avg > 1e-6) {
            p.cpg_scale    = d_cpg / d_avg;
            p.noncpg_scale = d_ncp / d_avg;
        }
    }
    return p;
}

double lsd_llr_score(const std::string& seq, const LsdClassifyParams& p)
{
    constexpr int    P_MAX = 5;
    constexpr double EPS   = 1e-10;
    const int L = static_cast<int>(seq.size());
    double llr = 0.0;
    const int np5 = std::min(P_MAX, L);
    for (int i = (p.skip_pos0_5prime ? 1 : 0); i < np5; ++i) {
        char c = seq[i];
        int hit = (c == 'T' || c == 't') ? 1
                : (c == 'C' || c == 'c') ? 0 : -1;
        if (hit < 0) continue;
        double site_scale = 1.0;
        if (i + 1 < L) {
            char n = seq[i + 1];
            bool cpg = (n == 'G' || n == 'g');
            site_scale = cpg ? p.cpg_scale : p.noncpg_scale;
        }
        double d_site = p.d_anc * site_scale;
        double pa = std::clamp(p.bg_5 + d_site * std::exp(-p.lam_5 * i) * (1.0 - p.bg_5),
                               EPS, 1.0 - EPS);
        double pm = std::clamp(p.bg_5, EPS, 1.0 - EPS);
        llr += hit ? (std::log(pa) - std::log(pm))
                   : (std::log(1.0 - pa) - std::log(1.0 - pm));
    }
    const int np3 = std::min(P_MAX, L);
    for (int i = 0; i < np3; ++i) {
        char c = seq[L - 1 - i];
        int hit;
        if (p.is_ss)
            hit = (c == 'T' || c == 't') ? 1 : (c == 'C' || c == 'c') ? 0 : -1;
        else
            hit = (c == 'A' || c == 'a') ? 1 : (c == 'G' || c == 'g') ? 0 : -1;
        if (hit < 0) continue;
        double pa = std::clamp(p.bg_3 + p.d_anc * std::exp(-p.lam_3 * i) * (1.0 - p.bg_3),
                               EPS, 1.0 - EPS);
        double pm = std::clamp(p.bg_3, EPS, 1.0 - EPS);
        llr += hit ? (std::log(pa) - std::log(pm))
                   : (std::log(1.0 - pa) - std::log(1.0 - pm));
    }
    return llr;
}

bool lsd_classify_read(const std::string& seq, const LsdClassifyParams& p)
{
    return lsd_llr_score(seq, p) > 0.0;
}

void lsd_accumulate(const std::string& seq, LsdLlrBinAccum& acc,
                    bool ancient, bool is_ss)
{
    const int L  = static_cast<int>(seq.size());
    const int np = std::min<int>(LengthBinDamageProfile::N_POS, L);
    for (int p = 0; p < np; ++p) {
        char c  = seq[p];
        char cU = (c == 'A' || c == 'C' || c == 'G' || c == 'T') ? c
                : (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : 'N';
        bool is_t = (cU == 'T');
        bool is_c = (cU == 'C');
        bool is_g = (cU == 'G');
        if (ancient) {
            if (is_t) { ++acc.t_5_anc[p]; ++acc.tc_5_anc[p]; }
            else if (is_c) { ++acc.tc_5_anc[p]; }
            if ((is_t || is_c) && (p + 1) < L) {
                char n  = seq[p + 1];
                char nU = (n == 'A' || n == 'C' || n == 'G' || n == 'T') ? n
                        : (n >= 'a' && n <= 'z') ? static_cast<char>(n - 32) : 'N';
                bool cpg = (nU == 'G');
                if (cpg) { if (is_t) ++acc.t_5_anc_cpg[p]; ++acc.tc_5_anc_cpg[p]; }
                else     { if (is_t) ++acc.t_5_anc_noncpg[p]; ++acc.tc_5_anc_noncpg[p]; }
            }
        }
        if (is_t || is_g) {
            if (ancient) { if (is_t) ++acc.t_5_anc_g[p]; ++acc.tg_5_anc[p]; }
            else         { if (is_t) ++acc.t_5_mod_g[p]; ++acc.tg_5_mod[p]; }
        }
        if (!ancient) {
            if (is_t) { ++acc.t_5_mod[p]; ++acc.tc_5_mod[p]; }
            else if (is_c) { ++acc.tc_5_mod[p]; }
        }
        if (ancient) {
            if      (cU == 'A') ++acc.a_5_anc_all[p];
            else if (cU == 'C') ++acc.c_5_anc_all[p];
            else if (cU == 'G') ++acc.g_5_anc_all[p];
            else if (cU == 'T') ++acc.t_5_anc_all[p];
        }
    }
    for (int p = 0; p < np; ++p) {
        char c = seq[L - 1 - p];
        if (is_ss) {
            bool is_t3 = (c == 'T' || c == 't');
            bool is_c3 = (c == 'C' || c == 'c');
            if (ancient) {
                if (is_t3) { ++acc.h_3_anc[p]; ++acc.n_3_anc[p]; }
                else if (is_c3) { ++acc.n_3_anc[p]; }
            } else {
                if (is_t3) { ++acc.h_3_mod[p]; ++acc.n_3_mod[p]; }
                else if (is_c3) { ++acc.n_3_mod[p]; }
            }
        } else {
            bool is_a3 = (c == 'A' || c == 'a');
            bool is_g3 = (c == 'G' || c == 'g');
            if (ancient) {
                if (is_a3) { ++acc.h_3_anc[p]; ++acc.n_3_anc[p]; }
                else if (is_g3) { ++acc.n_3_anc[p]; }
            } else {
                if (is_a3) { ++acc.h_3_mod[p]; ++acc.n_3_mod[p]; }
                else if (is_g3) { ++acc.n_3_mod[p]; }
            }
        }
    }
}

void lsd_accumulate_soft(const std::string& seq, LsdLlrBinAccum& acc,
                         double w, bool is_ss)
{
    if (seq.empty()) return;
    const int L = static_cast<int>(seq.size());
    // 5' end: positions 0..N_SOFT_POS-1
    const int np5 = std::min(LsdLlrBinAccum::N_SOFT_POS, L);
    for (int p = 0; p < np5; ++p) {
        char c = seq[p];
        bool is_t = (c == 'T' || c == 't');
        bool is_c = (c == 'C' || c == 'c');
        if (is_t) { acc.sw_t5_anc[p] += w; acc.sw_tc5_anc[p] += w; }
        else if (is_c) { acc.sw_tc5_anc[p] += w; }
    }
    // 3' end: positions 0..N_SOFT_POS-1 (pos 0 = last base, pos 1 = second-to-last)
    const int np3 = std::min(LsdLlrBinAccum::N_SOFT_POS, L);
    for (int p = 0; p < np3; ++p) {
        char c = seq[L - 1 - p];
        bool hit;
        if (is_ss)
            hit = (c == 'T' || c == 't');
        else
            hit = (c == 'A' || c == 'a');
        bool eligible = is_ss ? (c == 'T' || c == 't' || c == 'C' || c == 'c')
                               : (c == 'A' || c == 'a' || c == 'G' || c == 'g');
        if (eligible) {
            if (hit) acc.sw_h3_anc[p] += w;
            acc.sw_n3_anc[p] += w;
        }
    }
    acc.sw_sum += w;
}

