// damage_profile.cpp
// estimate_damage() — DART full pipeline damage estimation.
// Called by both 'derep' (Pass 0) and 'extend' (Pass 0).

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"

#include "taph/frame_selector_decl.hpp"
#include "taph/length_gc_joint_mixture.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/log_length_gmm.hpp"
// dart/sample_damage_profile.hpp is included transitively via damage_profile.hpp

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

static DamageProfile estimate_damage_impl(
    FastqReaderBase& reader,
    const std::string& path,
    double mask_threshold,
    taph::SampleDamageProfile::LibraryType forced_lib,
    int64_t max_reads)
{
    FastqRecord rec;
    taph::SampleDamageProfile dart_profile;
    dart_profile.forced_library_type = forced_lib;

    int      reads_scanned = 0;
    int      typical_len   = 0;
    int64_t  record_pos    = 0;
    int64_t  len_sum       = 0;

    while (reader.read(rec)) {
        record_pos++;
        if (max_reads > 0 && record_pos > max_reads) break;
        int L = static_cast<int>(rec.seq.size());
        if (L < 30) continue;
        len_sum += L;
        taph::FrameSelector::update_sample_profile(dart_profile, rec.seq);
        reads_scanned++;
    }

    if (reads_scanned > 0)
        typical_len = static_cast<int>(len_sum / reads_scanned);

    taph::FrameSelector::finalize_sample_profile(dart_profile);

    double d_max_5   = dart_profile.d_max_5prime;
    double d_max_3   = dart_profile.d_max_3prime;
    double lambda_5  = dart_profile.lambda_5prime;
    double lambda_3  = dart_profile.lambda_3prime;
    double bg_5      = dart_profile.fit_baseline_5prime;
    double bg_3      = dart_profile.fit_baseline_3prime;
    double background = (bg_5 + bg_3) / 2.0;
    double d_max_combined = dart_profile.d_max_combined;

    const bool is_ss = (dart_profile.library_type ==
                        taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);

    DamageProfile profile;
    profile.d_max_5prime   = d_max_5;
    profile.d_max_3prime   = d_max_3;
    profile.lambda_5prime  = lambda_5;
    profile.lambda_3prime  = lambda_3;
    profile.background     = background;
    profile.mask_threshold = mask_threshold;
    profile.bg_5_tc        = dart_profile.fit_baseline_5prime;
    profile.bg_3_channel   = dart_profile.fit_baseline_3prime;
    profile.mixture_d_damaged = dart_profile.mixture_d_ancient;
    profile.mixture_pi_damaged = dart_profile.mixture_pi_ancient;
    profile.mixture_d_reference = dart_profile.mixture_d_reference;
    profile.mixture_K = dart_profile.mixture_K;
    profile.mixture_n_components = dart_profile.mixture_n_components;
    profile.mixture_converged = dart_profile.mixture_converged;
    profile.mixture_identifiable = dart_profile.mixture_identifiable;
    profile.ss_mode        = is_ss;
    profile.enabled        = (d_max_combined > 0.02 || d_max_5 > 0.02 || d_max_3 > 0.02
                              || dart_profile.damage_validated);

    // Empirical per-position mask from DART's normalized frequencies.
    // DS: max(5' C→T excess, 3' G→A excess) > threshold.
    // SS: max(5' C→T excess, 3' T/(T+C) excess) > threshold — both ends have C→T.
    double bg_tc = dart_profile.baseline_t_freq /
                   (dart_profile.baseline_t_freq + dart_profile.baseline_c_freq + 1e-9);
    for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dart_profile.tc_total_5prime[p] >= MIN_COV_DP)
            excess_5 = dart_profile.t_freq_5prime[p] - bg_5;
        if (is_ss) {
            if (dart_profile.tc_total_3prime[p] >= MIN_COV_DP) {
                double tc3 = dart_profile.tc_total_3prime[p];
                double ct3_freq = (tc3 > 0)
                    ? dart_profile.t_freq_3prime[p] / tc3
                    : bg_tc;
                excess_3 = ct3_freq - bg_tc;
            }
        } else {
            if (dart_profile.ag_total_3prime[p] >= MIN_COV_DP)
                excess_3 = dart_profile.a_freq_3prime[p] - bg_3;
        }
        profile.mask_pos[p] = (excess_5 > mask_threshold) || (excess_3 > mask_threshold);
    }

    // Fill short gaps in the empirical damage mask (≤1 position) to handle
    // sampling noise creating non-contiguous masks.
    {
        constexpr int MAX_GAP = 1;
        int last_p = -1;
        for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
            if (profile.mask_pos[p]) {
                if (last_p >= 0 && (p - last_p - 1) <= MAX_GAP) {
                    for (int g = last_p + 1; g < p; ++g)
                        profile.mask_pos[g] = true;
                }
                last_p = p;
            }
        }
    }

    std::string sample_note = (max_reads > 0 && record_pos > max_reads)
        ? " [sampled " + std::to_string(max_reads) + " reads]" : "";
    log_info("Pass 0: DART damage estimation — " +
             std::to_string(reads_scanned) + " reads in " + path + sample_note);
    log_info("  Library type: " + std::string(dart_profile.library_type_str()) +
             (dart_profile.library_type_auto_detected ? " (auto-detected)" : " (forced)"));
    log_info("  5'-end: d_max=" + std::to_string(d_max_5) +
             " lambda=" + std::to_string(lambda_5) +
             " bg=" + std::to_string(bg_5));
    log_info("  3'-end: d_max=" + std::to_string(d_max_3) +
             " lambda=" + std::to_string(lambda_3) +
             " bg=" + std::to_string(bg_3));
    log_info("  d_max_combined=" + std::to_string(d_max_combined) +
             " (source=" + dart_profile.d_max_source_str() + ")" +
             (dart_profile.mixture_converged ? " [mixture]" : "") +
             (dart_profile.damage_validated  ? " [validated]" : "") +
             (dart_profile.damage_artifact   ? " [ARTIFACT]" : ""));
    if (profile.enabled && typical_len > 0) {
        profile.print_info(typical_len);
    } else {
        log_info("  Damage below threshold — standard exact hashing will be used");
    }
    return profile;
}

DamageProfile estimate_damage(
    const std::string& path,
    double mask_threshold,
    taph::SampleDamageProfile::LibraryType forced_lib,
    int64_t max_reads)
{
    auto reader = make_fastq_reader(path);
    return estimate_damage_impl(*reader, path, mask_threshold, forced_lib, max_reads);
}

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
    const DamageProfile* bulk_for_llr)
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
                 " length edges truncated (libtaph LengthBinStats::MAX_BINS=" +
                 std::to_string(taph::LengthBinStats::MAX_BINS) + ")");
        edges.resize(max_edges);
    }
    out.edges = edges;

    const int n_bins = 1 + static_cast<int>(edges.size());

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
    struct LlrBinAccum {
        int64_t n_damaged = 0;
        int64_t n_undamaged  = 0;
        std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> h_3_anc{};  // T (SS) or A (DS)
        std::array<int64_t, LengthBinDamageProfile::N_POS> n_3_anc{};  // T+C or A+G
        // CpG-split C→T at 5' in the ancient subset (metaDMG-style split).
        std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_cpg{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc_cpg{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_noncpg{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> tc_5_anc_noncpg{};
        // 8-oxoG marker (G→T) at 5' in both subsets (contrast = oxidation signal).
        std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_g{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> tg_5_anc{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_mod_g{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> tg_5_mod{};
        // Per-position 4-base histograms at 5' end in the ancient subset.
        // Used to compute the composition-based 8-oxoG rate via G-depletion
        // (pG_terminal - pG_interior). Independent of the C→T signal.
        std::array<int64_t, LengthBinDamageProfile::N_POS> a_5_anc_all{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> c_5_anc_all{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> g_5_anc_all{};
        std::array<int64_t, LengthBinDamageProfile::N_POS> t_5_anc_all{};
    };
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

    // ---- merge + finalize ------------------------------------------------
    taph::LengthBinStats master;
    master.forced_library_type = forced_lib;
    master.configure(edges);
    for (auto& w : worker_stats) master.merge(w);
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
        lb.cpg_contrast  = pr.cpg_contrast;
        lb.validated     = pr.damage_validated;
        lb.source        = pr.d_max_source_str();
        const bool bin_is_ss = (pr.library_type ==
                                taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
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
        lb.mixture_d_damaged    = pr.mixture_d_ancient;
        lb.mixture_d_reference  = pr.mixture_d_reference;
        lb.mixture_d_population = pr.mixture_d_population;
        lb.mixture_pi_damaged   = pr.mixture_pi_ancient;
        lb.mixture_n_components = pr.mixture_n_components;
        lb.mixture_converged    = pr.mixture_converged;
        lb.mixture_identifiable = pr.mixture_identifiable;
        for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
            const auto& gb = pr.gc_bins[g];
            lb.gc_d_max[g]     = gb.valid ? gb.d_max : -1.0;
            lb.gc_n_reads[g]   = static_cast<int64_t>(gb.n_reads);
            lb.gc_p_damaged[g] = gb.valid ? gb.p_damaged : -1.0;
        }
        for (int i = 0; i < 64; ++i) {
            lb.tri_5prime_terminal[i] = static_cast<int64_t>(pr.tri_5prime_terminal[i]);
            lb.tri_5prime_interior[i] = static_cast<int64_t>(pr.tri_5prime_interior[i]);
            lb.tri_3prime_terminal[i] = static_cast<int64_t>(pr.tri_3prime_terminal[i]);
            lb.tri_3prime_interior[i] = static_cast<int64_t>(pr.tri_3prime_interior[i]);
        }
        out.bins.push_back(lb);
    }

    // ---- Merge per-thread LLR accumulators → per-bin ancient-subset stats ---
    if (do_llr) {
        for (int b = 0; b < n_bins; ++b) {
            LlrBinAccum merged;
            for (int t = 0; t < n_threads; ++t) {
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

    // ---- Joint length × GC 2-component mixture (shared d_ancient) ----------
    {
        auto jr = taph::fit_length_gc_joint_mixture(master);
        out.d_joint_ancient    = jr.d_ancient;
        out.pi_joint_ancient   = jr.pi_ancient;
        out.d_joint_population = jr.d_population;
        out.joint_converged    = jr.converged;
        out.joint_separated    = jr.separated;
        out.cell_w_ancient     = std::move(jr.cell_w_ancient);
    }

    log_info("Pass 0: DART damage estimation — length stratified (" +
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

