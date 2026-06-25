// Length-stratified damage: estimate_damage_by_length, estimate_damage_split_model, DamageSplitModel.
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

// ---- length-stratified damage estimation -------------------------------
// Two-pass: pass 1 builds a log-length histogram and picks bin edges
// (explicit | quantile | auto GMM+BIC); pass 2 routes each read to its
// per-bin SampleDamageProfile accumulator and finalizes independently.

LengthStratifiedDamageProfile estimate_damage_by_length(
    const std::string& path,
    taph::SampleDamageProfile::LibraryType forced_lib,
    const LengthBinOptions& options,
    const std::vector<uint64_t>* prebuilt_hist,
    size_t reader_threads,
    int64_t max_reads,
    const DamageProfile* bulk_for_llr,
    const LsdPrebuilt* prebuilt)
{
    LengthStratifiedDamageProfile out;
    out.min_length = LSD_L_MIN;
    out.max_length = LSD_L_MAX;

    const double log_min = std::log(static_cast<double>(LSD_L_MIN));
    const double log_max = std::log(static_cast<double>(LSD_L_MAX + 1));
    const double step    = (log_max - log_min) / static_cast<double>(LSD_HIST_BINS);

    // ---- Pass 1: log-length histogram (skipped when caller supplies prebuilt_hist) -----
    std::vector<uint64_t> hist;
    int64_t reads_scanned = 0;
    if (prebuilt_hist) {
        hist = *prebuilt_hist;
        for (auto c : hist) reads_scanned += static_cast<int64_t>(c);
    } else {
        hist.assign(LSD_HIST_BINS, 0);
        auto reader = make_fastq_reader(path, reader_threads);
        FastqRecord rec;
        int64_t record_pos = 0;
        while (reader->read(rec)) {
            ++record_pos;
            if (max_reads > 0 && record_pos > max_reads) break;
            int L = static_cast<int>(rec.seq.size());
            if (L < LSD_L_MIN) continue;
            if (L > LSD_L_MAX) L = LSD_L_MAX;
            double lx = std::log(static_cast<double>(L));
            int bi = static_cast<int>((lx - log_min) / step);
            if (bi < 0) bi = 0;
            if (bi >= LSD_HIST_BINS) bi = LSD_HIST_BINS - 1;
            ++hist[bi];
            ++reads_scanned;
        }
    }
    out.reads_scanned = reads_scanned;

    // ---- pick edges ----------------------------------------------------
    std::vector<int> edges;
    if (options.mode == LengthBinOptions::Mode::EXPLICIT) {
        out.method = "explicit";
        int dropped_oor = 0;
        for (int e : options.explicit_edges) {
            if (e > LSD_L_MIN && e < LSD_L_MAX) edges.push_back(e);
            else ++dropped_oor;
        }
        if (dropped_oor > 0)
            log_warn("damage_profile: " + std::to_string(dropped_oor) +
                     " explicit length edges out of range (" +
                     std::to_string(LSD_L_MIN) + "," + std::to_string(LSD_L_MAX) +
                     ") were dropped");
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    } else if (options.mode == LengthBinOptions::Mode::QUANTILE) {
        int nb = std::clamp(options.quantile_bins, 1,
                            static_cast<int>(taph::LengthBinStats::MAX_BINS));
        out.method = "quantile";
        if (nb > 1)
            edges = taph::detect_quantile_length_edges(hist, log_min, log_max,
                                                      LSD_L_MIN, LSD_L_MAX, nb);
    } else if (options.mode == LengthBinOptions::Mode::AUTO) {
        taph::LogLengthGmmResult gr = taph::detect_log_length_gmm_edges(
            hist, log_min, log_max, LSD_L_MIN, LSD_L_MAX,
            static_cast<int>(taph::LengthBinStats::MAX_BINS));
        edges = gr.edges;
        out.method = edges.empty() ? "auto_gmm_single" : "auto_gmm";
    } else {
        out.method = "single";
    }
    const int max_edges = static_cast<int>(taph::LengthBinStats::MAX_BINS) - 1;
    if (static_cast<int>(edges.size()) > max_edges) {
        log_warn("damage_profile: " +
                 std::to_string(edges.size() - max_edges) +
                 " length edges truncated (taph::LengthBinStats::MAX_BINS=" +
                 std::to_string(taph::LengthBinStats::MAX_BINS) + ")");
        edges.resize(max_edges);
    }
    out.edges = edges;

    int n_bins = 1 + static_cast<int>(edges.size());

    // ---- Pass 2: multi-threaded per-bin accumulation via LengthBinStats -----
    int n_threads = reader_threads > 0
        ? static_cast<int>(reader_threads)
        : static_cast<int>(std::thread::hardware_concurrency());
    if (n_threads < 1) n_threads = 1;
    constexpr int BATCH_SZ = 8192;

    struct SeqQueue {
        std::mutex              mtx;
        std::condition_variable cv_not_empty, cv_not_full;
        std::vector<std::vector<std::string>> batches;
        bool done = false;
        int  max_depth;
        explicit SeqQueue(int d) : max_depth(d) {}
    };
    SeqQueue queue(2 * n_threads);

    auto push_batch = [&](std::vector<std::string>&& b) {
        std::unique_lock<std::mutex> lk(queue.mtx);
        queue.cv_not_full.wait(lk, [&]{ return (int)queue.batches.size() < queue.max_depth; });
        queue.batches.push_back(std::move(b));
        queue.cv_not_empty.notify_one();
    };
    auto pop_batch = [&](std::vector<std::string>& b) {
        std::unique_lock<std::mutex> lk(queue.mtx);
        queue.cv_not_empty.wait(lk, [&]{ return !queue.batches.empty() || queue.done; });
        if (queue.batches.empty()) return false;
        b = std::move(queue.batches.back());
        queue.batches.pop_back();
        queue.cv_not_full.notify_one();
        return true;
    };
    auto set_done = [&]() {
        std::unique_lock<std::mutex> lk(queue.mtx);
        queue.done = true;
        queue.cv_not_empty.notify_all();
    };

    std::vector<taph::LengthBinStats> worker_stats(n_threads);
    for (auto& w : worker_stats) {
        w.forced_library_type = forced_lib;
        w.configure(edges);
    }

    // ---- Per-read LLR accumulator (optional, when bulk params supplied) -----
    // Use the public type so prebuilt data from damage.cpp can be injected
    // directly without a copy (Fix E).
    using LlrBinAccum = LsdLlrBinAccum;
    const bool do_llr = (bulk_for_llr != nullptr) && bulk_for_llr->enabled;
    const double llr_d_anc = do_llr
        ? (bulk_for_llr->mixture_converged && bulk_for_llr->mixture_d_damaged > 0.01
           ? bulk_for_llr->mixture_d_damaged
           : std::max(bulk_for_llr->d_max_5prime, bulk_for_llr->d_max_3prime))
        : 0.0;
    const double llr_lam_5 = do_llr ? bulk_for_llr->lambda_5prime : 0.3;
    const double llr_lam_3 = do_llr ? bulk_for_llr->lambda_3prime : 0.3;
    const double llr_bg_5  = do_llr ? bulk_for_llr->bg_5_tc       : 0.5;
    const double llr_bg_3  = do_llr ? bulk_for_llr->bg_3_channel  : 0.5;
    const bool   llr_ss    = do_llr ? bulk_for_llr->ss_mode       : false;

    // CpG-weighted classifier scaling. If bulk CpG / non-CpG d_max estimates
    // are available, scale d_anc per-site by the ratio of that context's
    // damage rate to the average. Keeps the site-averaged rate at d_anc while
    // up-weighting CpG hits (stronger ancient evidence per site).
    double cpg_scale    = 1.0;
    double noncpg_scale = 1.0;
    if (do_llr) {
        const double d_cpg = bulk_for_llr->d_cpg_5prime;
        const double d_ncp = bulk_for_llr->d_noncpg_5prime;
        if (d_cpg > 0.01 && d_ncp > 0.01 && !std::isnan(d_cpg) && !std::isnan(d_ncp)) {
            const double d_avg = 0.5 * (d_cpg + d_ncp);
            if (d_avg > 1e-6) {
                cpg_scale    = d_cpg / d_avg;
                noncpg_scale = d_ncp / d_avg;
            }
        }
    }

    // Per-thread × per-bin ancient-subset accumulators.
    std::vector<std::vector<LlrBinAccum>> worker_llr(n_threads);
    for (auto& v : worker_llr) v.assign(n_bins, {});

    auto bin_index = [&edges](int L) -> int {
        int k = 0;
        while (k < static_cast<int>(edges.size()) && L >= edges[k]) ++k;
        return k;
    };

    // Per-read LLR using shared bulk mixture parameters.
    // P(obs=hit | ancient, p) = bg + d_anc * exp(-λ*p) * (1 - bg)
    // P(obs=hit | modern)     = bg   (position-independent)
    // "hit" = T at C-sites (5' any lib; 3' SS) or A at G-sites (3' DS).
    auto classify_read = [&](const std::string& seq) -> bool {
        constexpr int P_MAX  = 5;
        constexpr double EPS = 1e-10;
        int L = static_cast<int>(seq.size());
        double llr = 0.0;
        int np5 = std::min(P_MAX, L);
        for (int p = 0; p < np5; ++p) {
            char c = seq[p];
            int hit = (c == 'T' || c == 't') ? 1
                    : (c == 'C' || c == 'c') ? 0 : -1;
            if (hit < 0) continue;
            // CpG-context weighting: if next base is G, this C/T site is a CpG site;
            // use the CpG-scaled d_anc. Otherwise non-CpG scale.
            double site_scale = 1.0;
            if (p + 1 < L) {
                char n = seq[p + 1];
                bool cpg = (n == 'G' || n == 'g');
                site_scale = cpg ? cpg_scale : noncpg_scale;
            }
            double d_anc_site = llr_d_anc * site_scale;
            double p_anc = llr_bg_5 + d_anc_site * std::exp(-llr_lam_5 * p) * (1.0 - llr_bg_5);
            double p_mod = llr_bg_5;
            p_anc = std::clamp(p_anc, EPS, 1.0 - EPS);
            p_mod = std::clamp(p_mod, EPS, 1.0 - EPS);
            llr += hit ? (std::log(p_anc) - std::log(p_mod))
                       : (std::log(1.0 - p_anc) - std::log(1.0 - p_mod));
        }
        int np3 = std::min(P_MAX, L);
        for (int p = 0; p < np3; ++p) {
            char c = seq[L - 1 - p];
            int hit;
            if (llr_ss) {
                hit = (c == 'T' || c == 't') ? 1 : (c == 'C' || c == 'c') ? 0 : -1;
            } else {
                hit = (c == 'A' || c == 'a') ? 1 : (c == 'G' || c == 'g') ? 0 : -1;
            }
            if (hit < 0) continue;
            double p_anc = llr_bg_3 + llr_d_anc * std::exp(-llr_lam_3 * p) * (1.0 - llr_bg_3);
            double p_mod = llr_bg_3;
            p_anc = std::clamp(p_anc, EPS, 1.0 - EPS);
            p_mod = std::clamp(p_mod, EPS, 1.0 - EPS);
            llr += hit ? (std::log(p_anc) - std::log(p_mod))
                       : (std::log(1.0 - p_anc) - std::log(1.0 - p_mod));
        }
        return llr > 0.0;
    };

    auto accumulate_5prime_terminals = [&](const std::string& seq, LlrBinAccum& acc, bool ancient) {
        int L = static_cast<int>(seq.size());
        int np = std::min<int>(LengthBinDamageProfile::N_POS, L);
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
            }
            // CpG split: need next base to decide; CpG = C followed by G in 5'→3' direction.
            // For C/T sites, apply the same CpG marker (interpreting T as T-at-CpG = damaged CpG-C).
            if ((is_t || is_c) && (p + 1) < L) {
                char n  = seq[p + 1];
                char nU = (n == 'A' || n == 'C' || n == 'G' || n == 'T') ? n
                        : (n >= 'a' && n <= 'z') ? static_cast<char>(n - 32) : 'N';
                bool cpg = (nU == 'G');
                if (ancient) {
                    if (cpg) {
                        if (is_t) ++acc.t_5_anc_cpg[p];
                        ++acc.tc_5_anc_cpg[p];
                    } else {
                        if (is_t) ++acc.t_5_anc_noncpg[p];
                        ++acc.tc_5_anc_noncpg[p];
                    }
                }
            }
            // 8-oxoG marker: T/(T+G) at 5' terminal positions.
            if (is_t || is_g) {
                if (ancient) {
                    if (is_t) ++acc.t_5_anc_g[p];
                    ++acc.tg_5_anc[p];
                } else {
                    if (is_t) ++acc.t_5_mod_g[p];
                    ++acc.tg_5_mod[p];
                }
            }
            // Per-position base composition in the ancient subset — feeds the
            // clean G-depletion 8-oxoG metric (p_G(interior) - p_G(terminal)).
            if (ancient) {
                if      (cU == 'A') ++acc.a_5_anc_all[p];
                else if (cU == 'C') ++acc.c_5_anc_all[p];
                else if (cU == 'G') ++acc.g_5_anc_all[p];
                else if (cU == 'T') ++acc.t_5_anc_all[p];
            }
        }
    };

    auto accumulate_damaged_terminals = [&](const std::string& seq, LlrBinAccum& acc) {
        accumulate_5prime_terminals(seq, acc, /*damaged=*/true);
        int L = static_cast<int>(seq.size());
        int np = std::min<int>(LengthBinDamageProfile::N_POS, L);
        for (int p = 0; p < np; ++p) {
            char c = seq[L - 1 - p];
            if (llr_ss) {
                if (c == 'T' || c == 't') { ++acc.h_3_anc[p]; ++acc.n_3_anc[p]; }
                else if (c == 'C' || c == 'c') { ++acc.n_3_anc[p]; }
            } else {
                if (c == 'A' || c == 'a') { ++acc.h_3_anc[p]; ++acc.n_3_anc[p]; }
                else if (c == 'G' || c == 'g') { ++acc.n_3_anc[p]; }
            }
        }
    };

    // ---- Fix E: prebuilt path — skip FASTQ reader entirely ---------------
    // When the caller has already accumulated LSD data during the oxoG second
    // pass, inject it here and skip the worker launch + file read entirely.
    if (prebuilt != nullptr) {
        // Override edges and n_bins with whatever was computed before the oxoG pass.
        edges  = prebuilt->edges;
        n_bins = 1 + static_cast<int>(edges.size());
        out.edges = edges;
        // Inject prebuilt LLR bins as a single synthetic "thread" so the
        // merge loop below works unchanged (iterates worker_llr.size()).
        if (do_llr && !prebuilt->llr_bins.empty()) {
            worker_llr.assign(1, prebuilt->llr_bins);
        }
        // No workers were launched; just mark the queue done for cleanliness.
        set_done();
        // worker_stats is empty (no threads started); master is configured below.
    } else {
        std::vector<std::thread> workers;
        workers.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t) {
            workers.emplace_back([&, t]{
                std::vector<std::string> batch;
                while (pop_batch(batch)) {
                    for (const auto& seq : batch) {
                        int L = static_cast<int>(seq.size());
                        if (L < LSD_L_MIN) continue;
                        int L_binned = std::min(L, LSD_L_MAX);
                        worker_stats[t].update(seq, L_binned);
                        if (do_llr) {
                            int b = bin_index(L_binned);
                            if (b >= 0 && b < n_bins) {
                                bool anc = classify_read(seq);
                                if (anc) {
                                    ++worker_llr[t][b].n_damaged;
                                    accumulate_damaged_terminals(seq, worker_llr[t][b]);
                                } else {
                                    ++worker_llr[t][b].n_undamaged;
                                    accumulate_5prime_terminals(seq, worker_llr[t][b], /*damaged=*/false);
                                }
                            }
                        }
                    }
                }
            });
        }

        try {
            auto reader = make_fastq_reader(path, reader_threads);
            FastqRecord rec;
            std::vector<std::string> batch;
            batch.reserve(BATCH_SZ);
            int64_t record_pos = 0;
            while (reader->read(rec)) {
                ++record_pos;
                if (max_reads > 0 && record_pos > max_reads) break;
                batch.push_back(std::move(rec.seq));
                if (static_cast<int>(batch.size()) == BATCH_SZ) {
                    push_batch(std::move(batch));
                    batch.clear();
                    batch.reserve(BATCH_SZ);
                }
            }
            if (!batch.empty()) push_batch(std::move(batch));
        } catch (...) {
            set_done();
            for (auto& w : workers) w.join();
            throw;
        }
        set_done();
        for (auto& w : workers) w.join();
    }

    // ---- merge + finalize ------------------------------------------------
    taph::LengthBinStats master;
    master.forced_library_type = forced_lib;
    if (prebuilt != nullptr) {
        // Prebuilt stats were accumulated in damage.cpp workers; take them as-is.
        master = prebuilt->merged_stats;
        master.forced_library_type = forced_lib;
    } else {
        master.configure(edges);
        for (auto& w : worker_stats) master.merge(w);
    }
    master.finalize_all();

    // ---- extract LengthBinDamageProfile ----------------------------------
    out.bins.reserve(n_bins);
    for (int b = 0; b < n_bins; ++b) {
        LengthBinDamageProfile lb;
        lb.length_lo = (b == 0) ? LSD_L_MIN : edges[b - 1];
        lb.length_hi = (b == n_bins - 1) ? LSD_L_MAX : (edges[b] - 1);
        const auto& pr = master.profiles[b];
        lb.n_reads   = static_cast<int64_t>(pr.n_reads);
        if (pr.n_reads < static_cast<size_t>(LSD_MIN_READS_PER_BIN)) {
            lb.source = "insufficient_reads";
            out.bins.push_back(lb);
            continue;
        }
        lb.d_max_5prime  = pr.d_max_5prime;
        lb.d_max_3prime  = pr.d_max_3prime;
        lb.lambda_5prime = pr.lambda_5prime;
        lb.lambda_3prime = pr.lambda_3prime;
        lb.bg_5prime     = pr.fit_baseline_5prime;
        lb.bg_3prime     = pr.fit_baseline_3prime;
        lb.validated     = pr.damage_validated;
        lb.source        = pr.d_max_source_str();
        bool bin_is_ss = (pr.library_type ==
                          taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
        if (pr.library_type ==
                taph::SampleDamageProfile::LibraryType::UNKNOWN &&
            pr.library_type_evaluable) {
            bin_is_ss = (pr.library_p_ss > pr.library_p_ds);
        }
        lb.ss_mode = bin_is_ss;
        for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
            double cov5 = pr.tc_total_5prime[p];
            lb.per_pos_5prime_ct[p] =
                (cov5 >= MIN_COV_DP) ? pr.t_freq_5prime[p] : -1.0;
            double v3 = -1.0;
            if (bin_is_ss) {
                double cov3 = pr.tc_total_3prime[p];
                if (cov3 >= MIN_COV_DP)
                    v3 = pr.t_freq_3prime[p] / cov3;
            } else {
                double cov3 = pr.ag_total_3prime[p];
                if (cov3 >= MIN_COV_DP)
                    v3 = pr.a_freq_3prime[p];
            }
            lb.per_pos_3prime[p] = v3;
        }
        lb.mixture_d_damaged    = pr.mixture_d_damaged;
        lb.mixture_d_population_highgc = pr.mixture_d_population_highgc;
        lb.mixture_d_population = pr.mixture_d_population;
        lb.mixture_pi_damaged   = pr.mixture_pi_damaged;
        lb.mixture_n_components = pr.mixture_n_components;
        lb.mixture_converged    = pr.mixture_converged;
        lb.mixture_identifiable = pr.mixture_identifiable;
        for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
            const auto& gb = pr.gc_bins[g];
            lb.gc_d_max[g]     = gb.valid ? gb.d_max : -1.0;
            lb.gc_n_reads[g]   = static_cast<int64_t>(gb.n_reads);
            // C1: gc_p_damaged is a probability in [0,1]; -1.0 is the
            // "no-coverage / not-computed" sentinel for an invalid bin and is
            // emitted as JSON null by profile_json.cpp (guarded `< 0.0`). A
            // genuine P(damaged)=0 stays numeric and is preserved.
            lb.gc_p_damaged[g] = gb.valid ? gb.p_damaged : -1.0;
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                uint64_t t = gb.t_counts[p], c = gb.c_counts[p];
                lb.gc_per_pos_ct[g][p] = (t + c >= 100) ? static_cast<double>(t) / (t + c) : -1.0;
            }
        }
        for (int i = 0; i < 64; ++i) {
            lb.tri_5prime_terminal[i] = static_cast<int64_t>(pr.tri_5prime_terminal[i]);
            lb.tri_5prime_interior[i] = static_cast<int64_t>(pr.tri_5prime_interior[i]);
            lb.tri_3prime_terminal[i] = static_cast<int64_t>(pr.tri_3prime_terminal[i]);
            lb.tri_3prime_interior[i] = static_cast<int64_t>(pr.tri_3prime_interior[i]);
        }
        lb.ox_stop_rate_baseline = pr.ox_stop_conversion_rate_baseline;
        lb.ox_stop_rate_terminal = pr.ox_stop_rate_terminal;
        lb.ox_uniformity_ratio   = pr.ox_uniformity_ratio;
        lb.channel_c_valid       = pr.channel_c_valid;
        lb.ca_stop_rate_baseline = pr.ca_stop_rate_baseline;
        lb.ca_stop_rate_terminal = pr.ca_stop_rate_terminal;
        lb.ca_uniformity_ratio   = pr.ca_uniformity_ratio;
        lb.channel_f_z           = pr.channel_f_z;
        lb.channel_f_mh_z        = pr.channel_f_mh_z;
        lb.channel_f_common_or   = pr.channel_f_common_or;
        lb.channel_f_valid       = pr.channel_f_valid;
        lb.cg_stop_rate_baseline = pr.cg_stop_rate_baseline;
        lb.cg_stop_rate_terminal = pr.cg_stop_rate_terminal;
        lb.cg_uniformity_ratio   = pr.cg_uniformity_ratio;
        lb.channel_g_z           = pr.channel_g_z;
        lb.channel_g_valid       = pr.channel_g_valid;
        lb.at_stop_rate_baseline = pr.at_stop_rate_baseline;
        lb.at_stop_rate_terminal = pr.at_stop_rate_terminal;
        lb.at_uniformity_ratio   = pr.at_uniformity_ratio;
        lb.channel_h_z           = pr.channel_h_z;
        lb.channel_h_z_p2plus    = pr.channel_h_z_p2plus;
        lb.channel_h_valid       = pr.channel_h_valid;
        out.bins.push_back(lb);
    }

    // ---- Merge per-thread LLR accumulators → per-bin ancient-subset stats ---
    if (do_llr) {
        for (int b = 0; b < n_bins; ++b) {
            LlrBinAccum merged;
            for (size_t t = 0; t < worker_llr.size(); ++t) {
                const auto& w = worker_llr[t][b];
                merged.n_damaged += w.n_damaged;
                merged.n_undamaged  += w.n_undamaged;
                for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                    merged.t_5_anc[p]  += w.t_5_anc[p];
                    merged.tc_5_anc[p] += w.tc_5_anc[p];
                    merged.h_3_anc[p]  += w.h_3_anc[p];
                    merged.n_3_anc[p]  += w.n_3_anc[p];
                    merged.t_5_anc_cpg[p]     += w.t_5_anc_cpg[p];
                    merged.tc_5_anc_cpg[p]    += w.tc_5_anc_cpg[p];
                    merged.t_5_anc_noncpg[p]  += w.t_5_anc_noncpg[p];
                    merged.tc_5_anc_noncpg[p] += w.tc_5_anc_noncpg[p];
                    merged.t_5_anc_g[p]       += w.t_5_anc_g[p];
                    merged.tg_5_anc[p]        += w.tg_5_anc[p];
                    merged.t_5_mod_g[p]       += w.t_5_mod_g[p];
                    merged.tg_5_mod[p]        += w.tg_5_mod[p];
                    merged.a_5_anc_all[p]     += w.a_5_anc_all[p];
                    merged.c_5_anc_all[p]     += w.c_5_anc_all[p];
                    merged.g_5_anc_all[p]     += w.g_5_anc_all[p];
                    merged.t_5_anc_all[p]     += w.t_5_anc_all[p];
                }
            }
            auto& lb = out.bins[b];
            lb.n_damaged = merged.n_damaged;
            lb.n_undamaged  = merged.n_undamaged;
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                lb.per_pos_5prime_ct_damaged[p] = (merged.tc_5_anc[p] >= MIN_COV_DP)
                    ? static_cast<double>(merged.t_5_anc[p]) / merged.tc_5_anc[p]
                    : -1.0;
                lb.per_pos_3prime_damaged[p] = (merged.n_3_anc[p] >= MIN_COV_DP)
                    ? static_cast<double>(merged.h_3_anc[p]) / merged.n_3_anc[p]
                    : -1.0;
                lb.per_pos_5prime_ct_cpg_damaged[p] = (merged.tc_5_anc_cpg[p] >= MIN_COV_DP)
                    ? static_cast<double>(merged.t_5_anc_cpg[p]) / merged.tc_5_anc_cpg[p]
                    : -1.0;
                lb.per_pos_5prime_ct_noncpg_damaged[p] = (merged.tc_5_anc_noncpg[p] >= MIN_COV_DP)
                    ? static_cast<double>(merged.t_5_anc_noncpg[p]) / merged.tc_5_anc_noncpg[p]
                    : -1.0;
                lb.per_pos_5prime_gt_damaged[p] = (merged.tg_5_anc[p] >= MIN_COV_DP)
                    ? static_cast<double>(merged.t_5_anc_g[p]) / merged.tg_5_anc[p]
                    : -1.0;
                lb.per_pos_5prime_gt_undamaged[p] = (merged.tg_5_mod[p] >= MIN_COV_DP)
                    ? static_cast<double>(merged.t_5_mod_g[p]) / merged.tg_5_mod[p]
                    : -1.0;
            }
            auto pick_dmax = [](const std::array<double, LengthBinDamageProfile::N_POS>& pp,
                                double bg) -> double {
                if (pp[0] < 0 && pp[1] < 0) return -1.0;
                double v0 = pp[0];
                double v1 = pp[1];
                double use = (v0 >= 0 && v1 >= 0 && v0 < v1) ? v1 : (v0 >= 0 ? v0 : v1);
                return std::max(0.0, use - bg);
            };
            lb.d_max_5_damaged = pick_dmax(lb.per_pos_5prime_ct_damaged, llr_bg_5);
            lb.d_max_3_damaged = pick_dmax(lb.per_pos_3prime_damaged,    llr_bg_3);
            // CpG / non-CpG d_max (5' only — 3' CpG requires context past read end).
            // Use a neutral 0.5 baseline for the T/(T+C) scale in the damaged subset;
            // if the ancient-subset non-CpG baseline is known from later interior
            // positions we could refine, but p≥5 values in the ancient subset are
            // close to 0.5 in practice so this is fine for reporting.
            lb.d_max_5_cpg_damaged    = pick_dmax(lb.per_pos_5prime_ct_cpg_damaged,    llr_bg_5);
            lb.d_max_5_noncpg_damaged = pick_dmax(lb.per_pos_5prime_ct_noncpg_damaged, llr_bg_5);
            // 8-oxoG excess: Δ(ancient - modern) of T/(T+G) at pos 0 (negative = ancient is cleaner).
            if (merged.tg_5_anc[0] >= MIN_COV_DP && merged.tg_5_mod[0] >= MIN_COV_DP) {
                double f_anc = static_cast<double>(merged.t_5_anc_g[0]) / merged.tg_5_anc[0];
                double f_mod = static_cast<double>(merged.t_5_mod_g[0]) / merged.tg_5_mod[0];
                lb.s_gt_5_damaged_vs_undamaged = f_anc - f_mod;
            } else {
                lb.s_gt_5_damaged_vs_undamaged = std::numeric_limits<double>::quiet_NaN();
            }
            // Clean 8-oxoG metric: G-fraction depletion at 5' terminal vs interior
            // baseline. Under pure C→T deamination the G fraction is unchanged, so
            // any drop at terminal positions isolates G→T oxidative damage.
            // Interior baseline = mean p_G over positions 5..N_POS-1.
            auto pG_at = [&](int p) -> double {
                int64_t tot = merged.a_5_anc_all[p] + merged.c_5_anc_all[p]
                            + merged.g_5_anc_all[p] + merged.t_5_anc_all[p];
                if (tot < static_cast<int64_t>(MIN_COV_DP)) return -1.0;
                return static_cast<double>(merged.g_5_anc_all[p]) / tot;
            };
            {
                double pG0 = pG_at(0);
                double sum = 0.0;
                int    n   = 0;
                for (int p = 5; p < LengthBinDamageProfile::N_POS; ++p) {
                    double v = pG_at(p);
                    if (v >= 0.0) { sum += v; ++n; }
                }
                double pG_int = (n > 0) ? sum / n : -1.0;
                if (pG0 >= 0.0 && pG_int >= 0.0) {
                    lb.pG_terminal_5_damaged = pG0;
                    lb.pG_interior_5_damaged = pG_int;
                    lb.g_to_t_5_damaged      = pG_int - pG0;
                } else {
                    lb.pG_terminal_5_damaged = std::numeric_limits<double>::quiet_NaN();
                    lb.pG_interior_5_damaged = std::numeric_limits<double>::quiet_NaN();
                    lb.g_to_t_5_damaged      = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }

    // ---- Joint length × GC 2-component mixture (shared d_damaged) ----------
    {
        auto jr = taph::fit_length_gc_joint_mixture(master);
        out.d_joint_damaged    = jr.d_damaged;
        out.pi_joint_damaged   = jr.pi_damaged;
        out.d_joint_population = jr.d_population;
        out.joint_converged    = jr.converged;
        out.joint_separated    = jr.separated;
        out.cell_w_damaged     = std::move(jr.cell_w_damaged);
    }

    log_info("Pass 0: deamination damage estimation — length stratified (" +
             std::to_string(reads_scanned) + " reads; method=" + out.method +
             "; n_bins=" + std::to_string(n_bins) +
             "; threads=" + std::to_string(n_threads) + ") in " + path);
    for (int b = 0; b < n_bins; ++b) {
        const auto& lb = out.bins[b];
        log_info("  bin[" + std::to_string(b) + "] L=[" +
                 std::to_string(lb.length_lo) + "," + std::to_string(lb.length_hi) +
                 "] n=" + std::to_string(lb.n_reads) +
                 " d5=" + std::to_string(lb.d_max_5prime) +
                 " d3=" + std::to_string(lb.d_max_3prime) +
                 " lam5=" + std::to_string(lb.lambda_5prime) +
                 " lam3=" + std::to_string(lb.lambda_3prime) +
                 " src=" + lb.source);
    }
    return out;
}

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
        auto reader = make_fastq_reader(path);
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

                    // classify: additive LLR, no transcendentals
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
        auto reader = make_fastq_reader(path);
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
    const int L = static_cast<int>(seq.size());
    float s = 0.0f;

    if (!bins.empty()) {
        const Bin* bin = &bins.back();
        for (const auto& b : bins)
            if (L >= b.lo && L <= b.hi) { bin = &b; break; }

        const int lim5 = std::min(n_pos, LengthBinDamageProfile::N_POS - 1);
        for (int i = 1; i <= lim5 && i < L; ++i) {
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
        for (int i = 1; i <= lim5 && i < L; ++i) {
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

