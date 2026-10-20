// Split model: estimate_damage_split_model and DamageSplitModel (estimate_damage_by_length is in libtaph).
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

static constexpr double MIN_COV_DP = 100.0;

// ---- estimate_damage_split_model -------------------------------------------
// Stripped accumulation pass: only T/TC counts at N_POS positions per bin,
// classified by precomputed additive LLR (no log() in the read loop).

LengthStratifiedDamageProfile estimate_damage_split_model(
    const std::string& path,
    const DamageProfile& bulk,
    const std::vector<uint64_t>* prebuilt_hist,
    int n_workers)
{
    LengthStratifiedDamageProfile out;
    out.min_length = LSD_L_MIN;
    out.max_length = LSD_L_MAX;

    const double log_min = std::log(static_cast<double>(LSD_L_MIN));
    const double log_max = std::log(static_cast<double>(LSD_L_MAX + 1));
    const double step    = (log_max - log_min) / static_cast<double>(LSD_HIST_BINS);

    // ---- Pass 1: histogram (skipped when prebuilt_hist supplied) ----
    std::vector<uint64_t> hist;
    int64_t reads_scanned = 0;
    if (prebuilt_hist) {
        hist = *prebuilt_hist;
        for (auto c : hist) reads_scanned += static_cast<int64_t>(c);
    } else {
        hist.assign(LSD_HIST_BINS, 0);
        auto reader = make_fastq_reader(path, static_cast<size_t>(std::max(1, n_workers)));  // parallel inflate
        FastqRecord rec;
        while (reader->read(rec)) {
            int L = static_cast<int>(rec.seq.size());
            if (L < LSD_L_MIN || L > LSD_L_MAX) continue;
            double lx = std::log(static_cast<double>(L));
            int bi = std::clamp(static_cast<int>((lx - log_min) / step), 0, LSD_HIST_BINS - 1);
            ++hist[bi];
            ++reads_scanned;
        }
    }
    out.reads_scanned = reads_scanned;

    // ---- pick edges via GMM (same logic as estimate_damage_by_length AUTO) ----
    taph::LogLengthGmmResult gr = taph::detect_log_length_gmm_edges(
        hist, log_min, log_max, LSD_L_MIN, LSD_L_MAX,
        static_cast<int>(taph::LengthBinStats::MAX_BINS));
    out.edges  = gr.edges;
    out.method = gr.edges.empty() ? "auto_gmm_single" : "auto_gmm";
    const std::vector<int>& edges = out.edges;
    const int n_bins = 1 + static_cast<int>(edges.size());

    // ---- Precompute additive LLR classifier coefficients ----
    // log-ratio for observing T(hit=1) or C(hit=0) at each terminal position.
    // Classification: sum over positions; >0 → ancient.
    constexpr int CLASS_POS = 5;
    constexpr double EPS    = 1e-10;
    const bool   ss       = bulk.enabled && bulk.ss_mode;
    const double d_anc    = bulk.enabled
        ? (bulk.mixture_converged && bulk.mixture_d_damaged > 0.01
           ? bulk.mixture_d_damaged
           : std::max(bulk.d_max_5prime, bulk.d_max_3prime))
        : 0.0;
    const double lam5     = bulk.enabled ? bulk.lambda_5prime  : 0.3;
    const double lam3     = bulk.enabled ? bulk.lambda_3prime  : 0.3;
    const double bg5      = bulk.enabled ? bulk.bg_5_tc        : 0.5;
    const double bg3      = bulk.enabled ? bulk.bg_3_channel   : 0.5;

    double coeff5_T[CLASS_POS]{}, coeff5_C[CLASS_POS]{};
    double coeff3_T[CLASS_POS]{}, coeff3_C[CLASS_POS]{};
    for (int p = 0; p < CLASS_POS; ++p) {
        double p_anc5 = std::clamp(bg5 + d_anc * std::exp(-lam5 * p) * (1.0 - bg5), EPS, 1.0 - EPS);
        double p_mod5 = std::clamp(bg5, EPS, 1.0 - EPS);
        coeff5_T[p] = std::log(p_anc5) - std::log(p_mod5);
        coeff5_C[p] = std::log(1.0 - p_anc5) - std::log(1.0 - p_mod5);
        double p_anc3 = std::clamp(bg3 + d_anc * std::exp(-lam3 * p) * (1.0 - bg3), EPS, 1.0 - EPS);
        double p_mod3 = std::clamp(bg3, EPS, 1.0 - EPS);
        coeff3_T[p] = std::log(p_anc3) - std::log(p_mod3);
        coeff3_C[p] = std::log(1.0 - p_anc3) - std::log(1.0 - p_mod3);
    }

    // ---- Stripped per-bin accumulator ----
    // No CpG fields: the CpG/bulk amplitude ratio (r_cpg) is taken from the first LSD pass
    // (per_pos_5prime_ct_cpg_damaged). This pass only refines the bulk damaged-class profile.
    struct SplitBinAccum {
        int64_t n_dam = 0, n_und = 0;
        // all reads (bulk T/(T+C)) and damaged-class only
        std::array<int64_t, LengthBinDamageProfile::N_POS> all5_t{}, all5_tc{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> dam5_t{}, dam5_tc{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> all3_t{}, all3_tc{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> dam3_t{}, dam3_tc{};
    };

    auto bin_of = [&](int L) {
        int k = 0;
        while (k < static_cast<int>(edges.size()) && L >= edges[k]) ++k;
        return k;
    };

    if (n_workers < 1) n_workers = 1;
    if (n_workers > 8) n_workers = 8;
    constexpr int BATCH = 4096;

    struct Queue {
        std::mutex              mtx;
        std::condition_variable cv_fill, cv_drain;
        std::vector<std::vector<std::string>> q;
        bool done = false;
        int  cap;
        explicit Queue(int c) : cap(c) {}
    } queue(n_workers * 2);

    std::vector<std::vector<SplitBinAccum>> wacc(n_workers);
    for (auto& v : wacc) v.assign(n_bins, {});

    std::vector<std::thread> workers;
    workers.reserve(n_workers);
    for (int t = 0; t < n_workers; ++t) {
        workers.emplace_back([&, t] {
            std::vector<std::string> batch;
            for (;;) {
                {
                    std::unique_lock<std::mutex> lk(queue.mtx);
                    queue.cv_fill.wait(lk, [&]{ return !queue.q.empty() || queue.done; });
                    if (queue.q.empty()) break;
                    batch = std::move(queue.q.back());
                    queue.q.pop_back();
                    queue.cv_drain.notify_one();
                }
                for (const auto& seq : batch) {
                    int L = static_cast<int>(seq.size());
                    if (L < LSD_L_MIN) continue;
                    int Lb = std::min(L, LSD_L_MAX);
                    int b  = bin_of(Lb);
                    if (b < 0 || b >= n_bins) continue;
                    auto& acc = wacc[t][b];

                    // classify: additive LLR, no transcendentals.
                    // ds classifies on 5' only; the 3' G->A enrichment is then measured against
                    // this class. Classifying on 3' as well would make the 3' lod self-referential.
                    double llr = 0.0;
                    int np5 = std::min(CLASS_POS, L);
                    for (int p = 0; p < np5; ++p) {
                        char c = seq[p];
                        if      (c == 'T' || c == 't') llr += coeff5_T[p];
                        else if (c == 'C' || c == 'c') llr += coeff5_C[p];
                    }
                    if (ss) {
                        int np3 = std::min(CLASS_POS, L);
                        for (int p = 0; p < np3; ++p) {
                            char c = seq[L - 1 - p];
                            if      (c == 'T' || c == 't') llr += coeff3_T[p];
                            else if (c == 'C' || c == 'c') llr += coeff3_C[p];
                        }
                    }
                    const bool anc = (llr > 0.0);
                    if (anc) ++acc.n_dam; else ++acc.n_und;

                    // accumulate T/(T+C) at 5' positions
                    int np = std::min<int>(LengthBinDamageProfile::N_POS, L);
                    for (int p = 0; p < np; ++p) {
                        char c = seq[p];
                        bool is_t = (c == 'T' || c == 't');
                        bool is_c = (c == 'C' || c == 'c');
                        if (is_t || is_c) {
                            if (is_t) ++acc.all5_t[p];
                            ++acc.all5_tc[p];
                            if (anc) {
                                if (is_t) ++acc.dam5_t[p];
                                ++acc.dam5_tc[p];
                            }
                        }
                    }
                    // accumulate 3' positions (SS: T/(T+C), DS: A/(A+G))
                    if (ss) {
                        for (int p = 0; p < np; ++p) {
                            char c = seq[L - 1 - p];
                            bool is_t = (c == 'T' || c == 't');
                            bool is_c = (c == 'C' || c == 'c');
                            if (is_t || is_c) {
                                if (is_t) ++acc.all3_t[p];
                                ++acc.all3_tc[p];
                                if (anc) {
                                    if (is_t) ++acc.dam3_t[p];
                                    ++acc.dam3_tc[p];
                                }
                            }
                        }
                    } else {
                        for (int p = 0; p < np; ++p) {
                            char c = seq[L - 1 - p];
                            bool is_a = (c == 'A' || c == 'a');
                            bool is_g = (c == 'G' || c == 'g');
                            if (is_a || is_g) {
                                if (is_a) ++acc.all3_t[p];
                                ++acc.all3_tc[p];
                                if (anc) {
                                    if (is_a) ++acc.dam3_t[p];
                                    ++acc.dam3_tc[p];
                                }
                            }
                        }
                    }
                }
            }
        });
    }

    {
        auto reader = make_fastq_reader(path, static_cast<size_t>(std::max(1, n_workers)));  // parallel inflate
        FastqRecord rec;
        std::vector<std::string> batch;
        batch.reserve(BATCH);
        while (reader->read(rec)) {
            batch.push_back(std::move(rec.seq));
            if (static_cast<int>(batch.size()) == BATCH) {
                std::unique_lock<std::mutex> lk(queue.mtx);
                queue.cv_drain.wait(lk, [&]{ return static_cast<int>(queue.q.size()) < queue.cap; });
                queue.q.push_back(std::move(batch));
                queue.cv_fill.notify_one();
                batch.clear();
                batch.reserve(BATCH);
            }
        }
        if (!batch.empty()) {
            std::unique_lock<std::mutex> lk(queue.mtx);
            queue.cv_drain.wait(lk, [&]{ return static_cast<int>(queue.q.size()) < queue.cap; });
            queue.q.push_back(std::move(batch));
            queue.cv_fill.notify_one();
        }
        {
            std::unique_lock<std::mutex> lk(queue.mtx);
            queue.done = true;
            queue.cv_fill.notify_all();
        }
    }
    for (auto& w : workers) w.join();

    // ---- merge and populate output ----
    std::vector<SplitBinAccum> merged(n_bins);
    for (int t = 0; t < n_workers; ++t)
        for (int b = 0; b < n_bins; ++b) {
            merged[b].n_dam += wacc[t][b].n_dam;
            merged[b].n_und += wacc[t][b].n_und;
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                merged[b].all5_t[p]  += wacc[t][b].all5_t[p];
                merged[b].all5_tc[p] += wacc[t][b].all5_tc[p];
                merged[b].dam5_t[p]  += wacc[t][b].dam5_t[p];
                merged[b].dam5_tc[p] += wacc[t][b].dam5_tc[p];
                merged[b].all3_t[p]  += wacc[t][b].all3_t[p];
                merged[b].all3_tc[p] += wacc[t][b].all3_tc[p];
                merged[b].dam3_t[p]  += wacc[t][b].dam3_t[p];
                merged[b].dam3_tc[p] += wacc[t][b].dam3_tc[p];
            }
        }

    out.bins.reserve(n_bins);
    for (int b = 0; b < n_bins; ++b) {
        const auto& m = merged[b];
        LengthBinDamageProfile lb;
        lb.length_lo   = (b == 0) ? LSD_L_MIN : edges[b - 1];
        lb.length_hi   = (b == n_bins - 1) ? LSD_L_MAX : (edges[b] - 1);
        lb.n_damaged   = m.n_dam;
        lb.n_undamaged = m.n_und;
        lb.n_reads     = m.n_dam + m.n_und;
        lb.ss_mode     = ss;
        constexpr double MIN_COV = 100.0;
        for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
            lb.per_pos_5prime_ct[p] =
                m.all5_tc[p] >= MIN_COV ? static_cast<double>(m.all5_t[p]) / m.all5_tc[p] : -1.0;
            lb.per_pos_5prime_ct_damaged[p] =
                m.dam5_tc[p] >= MIN_COV ? static_cast<double>(m.dam5_t[p]) / m.dam5_tc[p] : -1.0;
            lb.per_pos_3prime[p] =
                m.all3_tc[p] >= MIN_COV ? static_cast<double>(m.all3_t[p]) / m.all3_tc[p] : -1.0;
            lb.per_pos_3prime_damaged[p] =
                m.dam3_tc[p] >= MIN_COV ? static_cast<double>(m.dam3_t[p]) / m.dam3_tc[p] : -1.0;
        }
        out.bins.push_back(lb);
    }
    log_info("Pass 0: damage split model (" + std::to_string(reads_scanned) +
             " reads; method=" + out.method +
             "; n_bins=" + std::to_string(n_bins) +
             "; workers=" + std::to_string(n_workers) + ") in " + path);
    return out;
}

// ---- DamageSplitModel -------------------------------------------------------

DamageSplitModel DamageSplitModel::build(const LengthStratifiedDamageProfile& lsd,
                                          const DamageProfile& bulk)
{
    DamageSplitModel m;
    m.fallback = bulk;

    for (const auto& lb : lsd.bins) {
        const int64_t n_tot = lb.n_damaged + lb.n_undamaged;
        if (n_tot == 0 || lb.n_undamaged == 0) continue;

        const double pi_d = static_cast<double>(lb.n_damaged) / n_tot;
        const double pi_u = static_cast<double>(lb.n_undamaged) / n_tot;

        Bin bin;
        bin.lo      = lb.length_lo;
        bin.hi      = lb.length_hi;
        bin.ss_mode = lb.ss_mode;

        // CpG amplitude ratio: scale the per-position 5' C->T damage EXCESS (not the baseline) by
        // d_max_5_cpg / d_max_5_bulk when both are well-estimated; else r_cpg=1 (CpG lod == bulk lod).
        double r_cpg = 1.0;
        if (std::isfinite(lb.d_max_5_cpg_damaged) && std::isfinite(lb.d_max_5_damaged)
                && lb.d_max_5_damaged > 1e-4f)
            r_cpg = std::clamp(static_cast<double>(lb.d_max_5_cpg_damaged) / lb.d_max_5_damaged, 0.0, 5.0);
        for (int i = 0; i < LengthBinDamageProfile::N_POS; ++i) {
            double pb = lb.per_pos_5prime_ct[i];
            double pd = lb.per_pos_5prime_ct_damaged[i];
            if (pb < 0.0 || pd < 0.0) {
                bin.lod5_T[i] = bin.lod5_C[i] = bin.lod5_T_cpg[i] = bin.lod5_C_cpg[i] = 0.0f; continue;
            }
            double pu = std::max(0.0001, std::min(0.9999, (pb - pi_d * pd) / pi_u));
            pd = std::max(0.0001, std::min(0.9999, pd));
            bin.lod5_T[i] = static_cast<float>(std::log(pd / pu));
            bin.lod5_C[i] = static_cast<float>(std::log((1.0 - pd) / (1.0 - pu)));
            // CpG-elevated: scale the excess (pd-pu) by r_cpg, keep the baseline pu.
            double pd_cpg = std::max(0.0001, std::min(0.9999, pu + (pd - pu) * r_cpg));
            bin.lod5_T_cpg[i] = static_cast<float>(std::log(pd_cpg / pu));
            bin.lod5_C_cpg[i] = static_cast<float>(std::log((1.0 - pd_cpg) / (1.0 - pu)));
        }
        // Build the 3' G->A lod for BOTH ds and ss (was ss-only). per_pos_3prime holds the 3' G->A
        // damage (A-fraction) for ds and the 3' deamination for ss; the lod stays 0 where there is no
        // 3' signal, so the score()'s now-ungated 3' term is harmless when absent. Recovers the joint end.
        for (int i = 0; i < LengthBinDamageProfile::N_POS; ++i) {
            double pb = lb.per_pos_3prime[i];
            double pd = lb.per_pos_3prime_damaged[i];
            if (pb < 0.0 || pd < 0.0) { bin.lod3_T[i] = bin.lod3_C[i] = 0.0f; continue; }
            double pu = std::max(0.0001, std::min(0.9999, (pb - pi_d * pd) / pi_u));
            pd = std::max(0.0001, std::min(0.9999, pd));
            bin.lod3_T[i] = static_cast<float>(std::log(pd / pu));
            bin.lod3_C[i] = static_cast<float>(std::log((1.0 - pd) / (1.0 - pu)));
        }
        m.bins.push_back(std::move(bin));
    }
    return m;
}

float DamageSplitModel::score(const std::string& seq, int n_pos) const
{
    // UNIFIED per-read scorer: delegate to the shared closed-form LLR that
    // `fqdup profile` uses (taph::lsd_llr_score ≡ lsd_llr_from_sig math), so
    // split and profile select damaged/undamaged reads with ONE mechanism.
    // cls_ is derived once from the (always-populated) fallback bulk fit; the
    // per-library d_anc/λ/bg + ss|ds channel + CpG scale all come from data.
    if (!cls_ready_) { cls_ = make_lsd_classify_params(fallback); cls_ready_ = true; }
    return static_cast<float>(taph::lsd_llr_score(seq, cls_));
}

float DamageSplitModel::score_empirical_legacy(const std::string& seq, int n_pos) const
{
    const int L = static_cast<int>(seq.size());
    float s = 0.0f;

    if (!bins.empty()) {
        const Bin* bin = &bins.back();
        for (const auto& b : bins)
            if (L >= b.lo && L <= b.hi) { bin = &b; break; }

        const int lim5 = std::min(n_pos, LengthBinDamageProfile::N_POS - 1);
        // i=0 is the 5' terminus — the MOST discriminative position (max C->T amplitude). The classifier
        // (coeff5, built from p=0 with exp(-lam*p) peaking at p=0) already weights it; score() masked it
        // out by starting at i=1, discarding the single strongest feature. Score from the terminus.
        for (int i = 0; i <= lim5 && i < L; ++i) {
            char c = seq[i];
            // CpG context (next base G): methyl-C deaminates faster -> use the CpG-elevated lod.
            const bool cpg = (i + 1 < L && seq[i + 1] == 'G');
            if      (c == 'T') s += cpg ? bin->lod5_T_cpg[i] : bin->lod5_T[i];
            else if (c == 'C') s += cpg ? bin->lod5_C_cpg[i] : bin->lod5_C[i];
        }
        // 3' G->A joint-ends LLR for BOTH ds and ss (was ss-only). The lod3 table encodes the 3' G->A
        // log-odds (built from per_pos_3prime) and is 0 wherever the bin has no 3' damage, so adding it
        // for ds is harmless when absent and recovers the joint-end signal single-end scoring drops.
        // Validated: composition-matched synthetic AUC +0.015-0.03 (5'-only -> 5'+3') at FLB-KapK damage.
        // Library-aware 3' transition (must match the per_pos_3prime build): ds 3' = G->A (damaged=A,
        // ref=G); ss 3' = C->T (damaged=T, ref=C). lod3_T is the DAMAGED-base lod, lod3_C the REF-base
        // lod, regardless of transition. The previous code scored A/G unconditionally, which mismatched
        // the C->T-trained ss lod3 (scoring the wrong bases). ss synthetic AUC gain doubles when correct.
        for (int i = 0; i <= lim5 && i < L; ++i) {   // i=0 = 3' terminus (max G->A), see 5' note above
            char c = seq[L - 1 - i];
            if (bin->ss_mode) {              // ss: 3' C->T
                if      (c == 'T') s += bin->lod3_T[i];
                else if (c == 'C') s += bin->lod3_C[i];
            } else {                         // ds: 3' G->A
                if      (c == 'A') s += bin->lod3_T[i];
                else if (c == 'G') s += bin->lod3_C[i];
            }
        }
        return s;
    }

    // Bulk exponential fallback
    // No damage signal → all reads score below threshold → undamaged
    if (!fallback.enabled) return -1.0f;
    // Bail only when NEITHER end carries signal. The 3' G->A live end is scored for BOTH ds and ss:
    // when the 5' is ->G-overcall artifact-dead (d_max_5prime~0, e.g. FLB57md) the damage lives entirely
    // on the 3' end, so 5'-only ds scoring is blind to it. A d_max_5prime=0 5' loop contributes log(1)=0
    // per position, so keeping it is harmless; the 3' term supplies the discrimination.
    if (fallback.d_max_5prime <= 0.0 && fallback.d_max_3prime <= 0.0)
        return -1.0f;
    for (int i = 1; i <= n_pos && i < L; ++i) {
        double pd = fallback.bg_5_tc + fallback.d_max_5prime *
                    std::exp(-fallback.lambda_5prime * i);
        double pu = fallback.bg_5_tc;
        pd = std::max(0.0001, std::min(0.9999, pd));
        pu = std::max(0.0001, std::min(0.9999, pu));
        char c = seq[i];
        if      (c == 'T') s += static_cast<float>(std::log(pd / pu));
        else if (c == 'C') s += static_cast<float>(std::log((1 - pd) / (1 - pu)));
    }
    if (fallback.d_max_3prime > 0.0) {
        for (int i = 1; i <= n_pos && i < L; ++i) {
            double pd = fallback.bg_3_channel + fallback.d_max_3prime *
                        std::exp(-fallback.lambda_3prime * i);
            double pu = fallback.bg_3_channel;
            pd = std::max(0.0001, std::min(0.9999, pd));
            pu = std::max(0.0001, std::min(0.9999, pu));
            char c = seq[L - 1 - i];
            // Library-aware 3' transition: ss 3' = C->T (damaged=T, ref=C), ds 3' = G->A (damaged=A, ref=G).
            const char dam = fallback.ss_mode ? 'T' : 'A';
            const char ref = fallback.ss_mode ? 'C' : 'G';
            if      (c == dam) s += static_cast<float>(std::log(pd / pu));
            else if (c == ref) s += static_cast<float>(std::log((1 - pd) / (1 - pu)));
        }
    }
    return s;
}

