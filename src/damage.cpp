// fqdup damage — multi-threaded single-pass FASTQ damage profiler
//
// Each worker thread owns its own SampleDamageProfile (no sharing, no locks).
// The producer fills a bounded work queue; workers drain it.
// All per-thread profiles are merged with merge_sample_profiles() at the end,
// then finalize_sample_profile() is called once on the combined result.
//
// Usage: fqdup damage -i FILE [-p THREADS] [--library-type auto|ds|ss]
//                             [--mask-threshold FLOAT] [--tsv FILE]

#include "fqdup/fastq_common.hpp"

#include "dart/frame_selector_decl.hpp"
#include "dart/sample_damage_profile.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

static constexpr int    N_POS     = 15;
static constexpr double MIN_COV   = 100.0;
static constexpr int    BATCH_SZ  = 8192;   // sequences per work unit

// ---- bounded work queue ------------------------------------------------
struct WorkQueue {
    std::mutex              mtx;
    std::condition_variable cv_not_empty;
    std::condition_variable cv_not_full;
    std::vector<std::vector<std::string>> batches;
    bool   done      = false;
    int    max_depth;

    explicit WorkQueue(int max_depth) : max_depth(max_depth) {}

    void push(std::vector<std::string>&& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_full.wait(lk, [&]{ return (int)batches.size() < max_depth; });
        batches.push_back(std::move(batch));
        cv_not_empty.notify_one();
    }

    bool pop(std::vector<std::string>& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_empty.wait(lk, [&]{ return !batches.empty() || done; });
        if (batches.empty()) return false;
        batch = std::move(batches.back());
        batches.pop_back();
        cv_not_full.notify_one();
        return true;
    }

    void set_done() {
        std::unique_lock<std::mutex> lk(mtx);
        done = true;
        cv_not_empty.notify_all();
    }
};

// ---- per-thread state --------------------------------------------------
struct WorkerState {
    dart::SampleDamageProfile profile;
    int64_t reads_scanned = 0;
    int64_t reads_skipped = 0;
    int     len_min       = INT_MAX;
    int     len_max       = 0;
    int64_t len_sum       = 0;
};

static void worker_fn(WorkQueue& queue, WorkerState& state) {
    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (const std::string& seq : batch) {
            int L = static_cast<int>(seq.size());
            if (L < 30) { ++state.reads_skipped; continue; }
            dart::FrameSelector::update_sample_profile(state.profile, seq);
            if (L < state.len_min) state.len_min = L;
            if (L > state.len_max) state.len_max = L;
            state.len_sum += L;
            ++state.reads_scanned;
        }
    }
}

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

static float per_read_ancient_weight(
    const std::string& seq,
    const dart::SampleDamageProfile& dp)
{
    int L = static_cast<int>(seq.size());
    int bin = dart::SampleDamageProfile::get_gc_bin(seq);
    const auto& gc = dp.gc_bins[bin];
    float bg   = (gc.valid && gc.baseline_tc > 0.0f) ? gc.baseline_tc : dp.fit_baseline_5prime;
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

static constexpr int N_GC = dart::SampleDamageProfile::N_GC_BINS;

struct OxBinAcc {
    double A = 0, B = 0, C = 0, D = 0;           // first moments
    double AA = 0, AB = 0, BB = 0;                 // second moments (anc)
    double CC = 0, CD = 0, DD = 0;                 // second moments (bg)
    double AC = 0, AD = 0, BC = 0, BD = 0;         // cross moments
    double M_total = 0;
    uint64_t reads = 0;
};

struct OxogWorkerState {
    OxBinAcc bins[N_GC];
};

static void oxog_worker(WorkQueue& queue,
                        OxogWorkerState& state,
                        const dart::SampleDamageProfile& dp,
                        bool is_ss)
{
    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (const std::string& seq : batch) {
            int L = static_cast<int>(seq.size());
            if (L < 30) continue;
            int bin = dart::SampleDamageProfile::get_gc_bin(seq);
            std::string oriented;
            double q;
            if (is_ss) {
                oriented = seq;
                q = per_read_ancient_weight(seq, dp);
            } else {
                float w_fwd = per_read_ancient_weight(seq, dp);
                std::string rc = revcomp(seq);
                float w_rev = per_read_ancient_weight(rc, dp);
                if (w_rev > w_fwd) { oriented = std::move(rc); q = w_rev; }
                else               { oriented = seq;            q = w_fwd; }
            }
            int beg = L / 3, end = L - (L / 3);
            if (end <= beg) continue;
            int T = 0, G = 0;
            for (int j = beg; j < end; ++j) {
                char b = static_cast<char>(std::toupper(
                             static_cast<unsigned char>(oriented[j])));
                if      (b == 'T') ++T;
                else if (b == 'G') ++G;
            }
            int M = T + G;
            if (M == 0) continue;
            double a = q * T, b_ = q * M;
            double c = (1.0 - q) * T, d = (1.0 - q) * M;
            auto& x = state.bins[bin];
            x.A += a;  x.B += b_;
            x.C += c;  x.D += d;
            x.AA += a*a; x.AB += a*b_; x.BB += b_*b_;
            x.CC += c*c; x.CD += c*d;  x.DD += d*d;
            x.AC += a*c; x.AD += a*d;  x.BC += b_*c; x.BD += b_*d;
            x.M_total += M;
            ++x.reads;
        }
    }
}

// ---- helpers -----------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " damage -i FILE [options]\n"
        << "\nOptions:\n"
        << "  -i FILE                    Input FASTQ (.gz or plain)\n"
        << "  -p N                       Worker threads (default: all cores)\n"
        << "  --library-type auto|ds|ss  Library type for 3'-end interpretation (default: auto)\n"
        << "  --mask-threshold FLOAT     Mask when excess P(deam) > T (default: 0.05)\n"
        << "  --tsv FILE                 Write per-position table as TSV\n"
        << "  --json FILE                Write full damage profile as JSON\n";
}

// ---- main --------------------------------------------------------------
int damage_main(int argc, char** argv) {
    std::string in_path;
    int         n_threads     = static_cast<int>(std::thread::hardware_concurrency());
    if (n_threads < 1) n_threads = 1;
    double      mask_threshold = 0.05;
    std::string tsv_path;
    std::string json_path;
    bool        run_oxog      = true;
    dart::SampleDamageProfile::LibraryType forced_lib =
        dart::SampleDamageProfile::LibraryType::UNKNOWN;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            in_path = argv[++i];
        } else if ((arg == "-p" || arg == "--threads") && i + 1 < argc) {
            n_threads = std::stoi(argv[++i]);
            if (n_threads < 1) n_threads = 1;
        } else if (arg == "--library-type" && i + 1 < argc) {
            std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                forced_lib = dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                forced_lib = dart::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt != "auto") {
                std::cerr << "Error: unknown --library-type: " << lt << "\n";
                return 1;
            }
        } else if (arg == "--mask-threshold" && i + 1 < argc) {
            mask_threshold = std::stod(argv[++i]);
        } else if (arg == "--tsv" && i + 1 < argc) {
            tsv_path = argv[++i];
        } else if (arg == "--json" && i + 1 < argc) {
            json_path = argv[++i];
        } else if (arg == "--no-oxog") {
            run_oxog = false;
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Error: unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (in_path.empty()) {
        std::cerr << "Error: -i FILE required\n";
        print_usage(argv[0]);
        return 1;
    }
    if (mask_threshold <= 0.0 || mask_threshold >= 1.0) {
        std::cerr << "Error: --mask-threshold must be in (0, 1), got " << mask_threshold << "\n";
        return 1;
    }

    // ---- set up workers ------------------------------------------------
    // Queue depth = 2 * n_threads so the producer stays ahead without OOM.
    WorkQueue queue(2 * n_threads);
    std::vector<WorkerState> states(n_threads);
    std::vector<std::thread> workers;
    workers.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t)
        workers.emplace_back(worker_fn, std::ref(queue), std::ref(states[t]));

    // ---- producer: stream FASTQ into batches ---------------------------
    try {
        auto reader = make_fastq_reader(in_path, static_cast<size_t>(n_threads));
        FastqRecord rec;
        std::vector<std::string> batch;
        batch.reserve(BATCH_SZ);
        while (reader->read(rec)) {
            batch.push_back(std::move(rec.seq));
            if ((int)batch.size() == BATCH_SZ) {
                queue.push(std::move(batch));
                batch.clear();
                batch.reserve(BATCH_SZ);
            }
        }
        if (!batch.empty())
            queue.push(std::move(batch));
    } catch (...) {
        queue.set_done();
        for (auto& w : workers) w.join();
        throw;
    }
    queue.set_done();

    for (auto& w : workers) w.join();

    // ---- merge all per-thread results ----------------------------------
    dart::SampleDamageProfile dp;
    dp.forced_library_type = forced_lib;
    int64_t reads_scanned = 0, reads_skipped = 0;
    int     len_min = INT_MAX, len_max = 0;
    int64_t len_sum = 0;

    for (auto& s : states) {
        dart::FrameSelector::merge_sample_profiles(dp, s.profile);
        reads_scanned += s.reads_scanned;
        reads_skipped += s.reads_skipped;
        if (s.len_min < len_min) len_min = s.len_min;
        if (s.len_max > len_max) len_max = s.len_max;
        len_sum += s.len_sum;
    }
    dp.forced_library_type = forced_lib;  // ensure forced_lib survives merge
    dart::FrameSelector::finalize_sample_profile(dp);

    const bool is_ss = (dp.library_type == dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED);

    // ---- second pass: oxidation-compatible composition score -----------
    // Run when pass-1 found non-trivial 5' damage (gives meaningful q weights).
    double S_oxog = 0.0, SE_oxog = 0.0;
    bool has_oxog_score = false;

    if (run_oxog && dp.d_max_5prime > 0.01f) {
        WorkQueue queue2(2 * n_threads);
        std::vector<OxogWorkerState> ox_states(n_threads);
        std::vector<std::thread> workers2;
        workers2.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t)
            workers2.emplace_back(oxog_worker, std::ref(queue2),
                                   std::ref(ox_states[t]), std::cref(dp), is_ss);
        try {
            auto reader2 = make_fastq_reader(in_path, static_cast<size_t>(n_threads));
            FastqRecord rec;
            std::vector<std::string> batch;
            batch.reserve(BATCH_SZ);
            while (reader2->read(rec)) {
                batch.push_back(std::move(rec.seq));
                if ((int)batch.size() == BATCH_SZ) {
                    queue2.push(std::move(batch));
                    batch.clear();
                    batch.reserve(BATCH_SZ);
                }
            }
            if (!batch.empty()) queue2.push(std::move(batch));
        } catch (...) {
            queue2.set_done();
            for (auto& w : workers2) w.join();
            throw;
        }
        queue2.set_done();
        for (auto& w : workers2) w.join();

        // merge per-thread accumulators
        OxBinAcc merged[N_GC] = {};
        for (auto& os : ox_states)
            for (int b = 0; b < N_GC; ++b) {
                merged[b].A += os.bins[b].A;  merged[b].B += os.bins[b].B;
                merged[b].C += os.bins[b].C;  merged[b].D += os.bins[b].D;
                merged[b].AA += os.bins[b].AA; merged[b].AB += os.bins[b].AB; merged[b].BB += os.bins[b].BB;
                merged[b].CC += os.bins[b].CC; merged[b].CD += os.bins[b].CD; merged[b].DD += os.bins[b].DD;
                merged[b].AC += os.bins[b].AC; merged[b].AD += os.bins[b].AD;
                merged[b].BC += os.bins[b].BC; merged[b].BD += os.bins[b].BD;
                merged[b].M_total += os.bins[b].M_total;
                merged[b].reads   += os.bins[b].reads;
            }

        // aggregate across GC bins
        double M_all = 0.0;
        for (int b = 0; b < N_GC; ++b) {
            if (merged[b].B > 0 && merged[b].D > 0)
                M_all += merged[b].M_total;
        }

        if (M_all > 0) {
            double var_S = 0.0;
            for (int b = 0; b < N_GC; ++b) {
                const auto& x = merged[b];
                if (x.B <= 0 || x.D <= 0) continue;
                double alpha = x.M_total / M_all;
                double p_anc = x.A / x.B;
                double p_bg  = x.C / x.D;
                S_oxog += alpha * (p_anc - p_bg);

                double v_anc = (x.AA - 2*p_anc*x.AB + p_anc*p_anc*x.BB) / (x.B*x.B);
                double v_bg  = (x.CC - 2*p_bg *x.CD + p_bg *p_bg *x.DD) / (x.D*x.D);
                double cov   = (x.AC - p_bg*x.AD - p_anc*x.BC + p_anc*p_bg*x.BD) / (x.B*x.D);
                double var_b = v_anc + v_bg - 2*cov;
                if (var_b < 0) var_b = 0;
                var_S += alpha * alpha * var_b;
            }
            SE_oxog = std::sqrt(var_S);
            has_oxog_score = true;
        }
    }

    // ---- compute per-position mask -------------------------------------
    // After finalize_sample_profile():
    //   t_freq_5prime[p]  = T/(T+C) frequency (already normalised)
    //   a_freq_3prime[p]  = A/(A+G) frequency (already normalised)
    //   t_freq_3prime[p]  = raw T count at 3'  (NOT normalised)
    const double bg_5  = dp.fit_baseline_5prime;
    const double bg_3  = dp.fit_baseline_3prime;
    const double bg_tc = dp.baseline_t_freq /
                         (dp.baseline_t_freq + dp.baseline_c_freq + 1e-9);

    bool mask_pos[N_POS] = {};
    for (int p = 0; p < N_POS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dp.tc_total_5prime[p] >= MIN_COV)
            excess_5 = dp.t_freq_5prime[p] - bg_5;
        if (is_ss) {
            if (dp.tc_total_3prime[p] >= MIN_COV)
                excess_3 = dp.t_freq_3prime[p] / dp.tc_total_3prime[p] - bg_tc;
        } else {
            if (dp.ag_total_3prime[p] >= MIN_COV)
                excess_3 = dp.a_freq_3prime[p] - bg_3;
        }
        mask_pos[p] = (excess_5 > mask_threshold) || (excess_3 > mask_threshold);
    }

    int n_masked = 0;
    std::string masked_str;
    for (int p = 0; p < N_POS; ++p) {
        if (mask_pos[p]) {
            if (n_masked) masked_str += ',';
            masked_str += std::to_string(p);
            ++n_masked;
        }
    }

    // ---- comma-formatted integer helper --------------------------------
    auto fmt_count = [](int64_t n) {
        std::string s = std::to_string(n);
        int pos = static_cast<int>(s.size()) - 3;
        while (pos > 0) { s.insert(pos, ","); pos -= 3; }
        return s;
    };

    // ---- human-readable report ----------------------------------------
    std::cout << "=== fqdup damage ===\n";
    std::cout << "Input:   " << in_path << "\n";
    std::cout << "Threads: " << n_threads << "\n";
    std::cout << "Reads:   " << fmt_count(reads_scanned) << " scanned";
    if (reads_skipped)
        std::cout << "  (" << fmt_count(reads_skipped) << " skipped, <30 bp)";
    std::cout << "\n";
    if (reads_scanned > 0) {
        double mean_len = static_cast<double>(len_sum) / reads_scanned;
        std::cout << "Length:  min=" << len_min
                  << "  mean=" << std::fixed << std::setprecision(1) << mean_len
                  << "  max=" << len_max << "\n";
    }
    std::cout << "\n";

    std::cout << "Library: " << dp.library_type_str()
              << (dp.library_type_rescued ? " [rescued]" : "")
              << (dp.library_type_auto_detected ? " (auto-detected)" : " (forced)")
              << "\n";
    if (dp.composition_bias_5prime || dp.composition_bias_3prime) {
        std::cout << "  WARN:  composition bias at ";
        if (dp.composition_bias_5prime) std::cout << "5'";
        if (dp.composition_bias_5prime && dp.composition_bias_3prime) std::cout << " + ";
        if (dp.composition_bias_3prime) std::cout << "3'";
        std::cout << " — library type call may be unreliable\n";
    }
    std::cout << "  BIC  bias=" << std::fixed << std::setprecision(1) << dp.library_bic_bias
              << "  DS=" << dp.library_bic_ds
              << "  SS=" << dp.library_bic_ss
              << "  SS_full=" << dp.library_bic_mix << "\n";
    std::cout << "  fit  CT5_amp=" << std::setprecision(4) << dp.libtype_amp_ct5
              << "  ΔBIC=" << std::setprecision(1) << dp.libtype_dbic_ct5
              << "  GA3_amp=" << dp.libtype_amp_ga3
              << "  ΔBIC=" << dp.libtype_dbic_ga3
              << "  GA0_amp=" << dp.libtype_amp_ga0
              << "  ΔBIC=" << dp.libtype_dbic_ga0
              << "  CT3_amp=" << dp.libtype_amp_ct3
              << "  ΔBIC=" << dp.libtype_dbic_ct3 << "\n";
    bool artifact_5 = dp.position_0_artifact_5prime || dp.inverted_pattern_5prime;
    bool artifact_3 = dp.position_0_artifact_3prime || dp.inverted_pattern_3prime;
    if ((artifact_5 || artifact_3) && (dp.fit_offset_5prime > 1 || dp.fit_offset_3prime > 1)) {
        std::cout << "  NOTE: adapter artifact — 5' fit start=pos" << dp.fit_offset_5prime
                  << "  3' fit start=pos" << dp.fit_offset_3prime
                  << " (d_max corrected to peak of pos1-5)\n";
    }
    std::cout << "  5' terminal shift: " << std::showpos << std::fixed << std::setprecision(4)
              << dp.terminal_shift_5prime << std::noshowpos
              << "  (z=" << std::setprecision(1) << dp.terminal_z_5prime << ")\n";
    if (is_ss) {
        std::cout << "  3' ctrl shift:     " << std::showpos << std::fixed << std::setprecision(4)
                  << dp.ctrl_shift_3prime << std::noshowpos
                  << "  (z=" << std::setprecision(1) << dp.ctrl_z_3prime << ")  [SS signal]\n";
    } else {
        std::cout << "  3' terminal shift: " << std::showpos << std::fixed << std::setprecision(4)
                  << dp.terminal_shift_3prime << std::noshowpos
                  << "  (z=" << std::setprecision(1) << dp.terminal_z_3prime << ")\n";
    }
    std::cout << "\n";

    std::cout << "5'-end   d_max=" << std::fixed << std::setprecision(4) << dp.d_max_5prime
              << "  lambda=" << std::setprecision(3) << dp.lambda_5prime
              << "  bg=" << std::setprecision(4) << bg_5 << "\n";
    std::cout << "3'-end   d_max=" << std::fixed << std::setprecision(4) << dp.d_max_3prime
              << "  lambda=" << std::setprecision(3) << dp.lambda_3prime
              << "  bg=" << std::setprecision(4) << bg_3 << "\n";
    std::cout << "combined d_max=" << std::fixed << std::setprecision(4) << dp.d_max_combined
              << "  (source=" << dp.d_max_source_str() << ")";
    if (dp.mixture_converged) std::cout << " [mixture]";
    if (dp.damage_validated)  std::cout << " [validated]";
    if (dp.damage_artifact)   std::cout << " [ARTIFACT]";
    std::cout << "\n\n";

    std::cout << "Mask threshold: " << mask_threshold << " → "
              << n_masked << " position" << (n_masked != 1 ? "s" : "") << " masked";
    if (n_masked) std::cout << " (pos " << masked_str << ")";
    std::cout << "\n\n";

    // ---- complement asymmetry (8-oxoG QC) ---------------------------------
    // D = T/(T+G) - A/(A+C) at interior positions.
    // Under a balanced DS library with no G→T oxidation, D ≈ 0 (Chargaff).
    // Positive D is compatible with 8-oxoG but not uniquely diagnostic.
    {
        double gt_total = dp.baseline_g_to_t_count + dp.baseline_g_total;
        double ca_total = dp.baseline_c_to_a_count + dp.baseline_c_ox_total;
        if (gt_total >= 500.0 && ca_total >= 500.0) {
            std::cout << "8-oxoG (G\u2192T):\n";
            // Model-based detection: GT(p) = A*exp(-mu*p) + B; B = uniform background
            std::cout << "  GT fit:  B=" << std::fixed << std::setprecision(4) << dp.g_bg_fitted
                      << "  A=" << dp.g_term_fitted
                      << "  mu=" << std::setprecision(2) << dp.g_decay_fitted
                      << "  s_gt=" << std::showpos << std::setprecision(4) << dp.s_gt << std::noshowpos
                      << "  AC_interior=" << std::setprecision(4) << dp.ox_ca_baseline;
            std::cout << (dp.ox_damage_detected ? "  [detected]" : "  [not detected]") << "\n";
            if (has_oxog_score) {
                std::cout << "  S_oxog=" << std::showpos << std::fixed << std::setprecision(4)
                          << S_oxog << std::noshowpos
                          << "  SE=" << std::setprecision(4) << SE_oxog
                          << (is_ss ? "  [SS: ancient-weighted T/(T+G) enrichment]"
                                    : "  [DS: dual-orient ancient-weighted T/(T+G) enrichment]")
                          << "\n";
            }
            std::cout << "\n";
        }
    }

    // ---- depurination (Channel E) --------------------------------------
    std::cout << "Depurination:";
    if (dp.depurination_detected) {
        std::cout << "  [DETECTED]"
                  << "  enrichment_5=" << std::showpos << std::setprecision(4)
                  << dp.purine_enrichment_5prime << std::noshowpos
                  << "  enrichment_3=" << std::showpos
                  << dp.purine_enrichment_3prime << std::noshowpos;
    } else {
        std::cout << "  not detected"
                  << "  enrichment_5=" << std::showpos << std::setprecision(4)
                  << dp.purine_enrichment_5prime << std::noshowpos
                  << "  enrichment_3=" << std::showpos
                  << dp.purine_enrichment_3prime << std::noshowpos;
    }
    std::cout << "\n\n";

    // ---- per-position table --------------------------------------------
    const char* col3_hdr = is_ss ? "3'_CT  " : "3'_GA  ";
    std::cout << "pos  5'_CT   5'_GA   " << col3_hdr << " 5'_GT  \n";
    std::cout << "---  ------  ------  ------  ------\n";
    for (int p = 0; p < N_POS; ++p) {
        double freq5_ct = (dp.tc_total_5prime[p] >= MIN_COV) ? dp.t_freq_5prime[p] : 0.0;
        double ag5 = dp.a_freq_5prime[p] + dp.g_freq_5prime[p];
        double freq5_ga = (ag5 >= MIN_COV) ? dp.a_freq_5prime[p] / ag5 : 0.0;
        double freq3;
        if (is_ss) {
            freq3 = (dp.tc_total_3prime[p] >= MIN_COV)
                    ? dp.t_freq_3prime[p] / dp.tc_total_3prime[p] : 0.0;
        } else {
            freq3 = (dp.ag_total_3prime[p] >= MIN_COV) ? dp.a_freq_3prime[p] : 0.0;
        }
        double tg_denom = dp.t_from_g_5prime[p] + dp.g_count_5prime[p];
        double freq5_gt = (tg_denom >= MIN_COV) ? dp.t_from_g_5prime[p] / tg_denom : 0.0;
        std::cout << std::setw(3) << p
                  << "  " << std::fixed << std::setprecision(4) << freq5_ct
                  << "  " << freq5_ga
                  << "  " << freq3
                  << "  " << freq5_gt;
        if (mask_pos[p]) std::cout << "  *";
        std::cout << "\n";
    }

    // ---- optional TSV --------------------------------------------------
    if (!tsv_path.empty()) {
        std::ofstream tsv(tsv_path);
        if (!tsv) {
            std::cerr << "Error: cannot write TSV: " << tsv_path << "\n";
            return 1;
        }
        tsv << "pos\tend5\tfreq5\tend3\tfreq3\tfreq5_gt\tmask\tcov5\tcov3\tcov5_gt\n";
        for (int p = 0; p < N_POS; ++p) {
            double cov5  = dp.tc_total_5prime[p];
            double freq5 = (cov5 >= MIN_COV) ? dp.t_freq_5prime[p] : -1.0;
            double cov3, freq3;
            if (is_ss) {
                cov3  = dp.tc_total_3prime[p];
                freq3 = (cov3 >= MIN_COV) ? dp.t_freq_3prime[p] / cov3 : -1.0;
            } else {
                cov3  = dp.ag_total_3prime[p];
                freq3 = (cov3 >= MIN_COV) ? dp.a_freq_3prime[p] : -1.0;
            }
            double cov5_gt  = dp.t_from_g_5prime[p] + dp.g_count_5prime[p];
            double freq5_gt = (cov5_gt >= MIN_COV) ? dp.t_from_g_5prime[p] / cov5_gt : -1.0;
            tsv << p
                << "\t5prime\t" << std::fixed << std::setprecision(6) << freq5
                << "\t3prime\t" << freq3
                << "\t" << freq5_gt
                << "\t" << (mask_pos[p] ? "1" : "0")
                << "\t" << static_cast<int64_t>(cov5)
                << "\t" << static_cast<int64_t>(cov3)
                << "\t" << static_cast<int64_t>(cov5_gt)
                << "\n";
        }
        std::cout << "TSV written: " << tsv_path << "\n";
    }

    // ---- optional JSON -------------------------------------------------
    if (!json_path.empty()) {
        std::ofstream j(json_path);
        if (!j) {
            std::cerr << "Error: cannot write JSON: " << json_path << "\n";
            return 1;
        }
        auto fp4 = [](double v) { return std::to_string(static_cast<int>(v * 10000) / 10000.0); };
        j << std::fixed;
        j << "{\n";
        j << "  \"input\": \"" << in_path << "\",\n";
        j << "  \"n_reads\": " << reads_scanned << ",\n";
        j << "  \"library_type\": \"" << dp.library_type_str() << "\",\n";
        j << "  \"library_type_auto\": " << (dp.library_type_auto_detected ? "true" : "false") << ",\n";
        j << "  \"library_type_rescued\": " << (dp.library_type_rescued ? "true" : "false") << ",\n";
        const char* ds_str = (dp.damage_status == dart::SampleDamageProfile::DamageStatus::PRESENT) ? "present"
                           : (dp.damage_status == dart::SampleDamageProfile::DamageStatus::WEAK)    ? "weak"
                           : "absent";
        j << "  \"damage_status\": \"" << ds_str << "\",\n";
        j << "  \"deamination\": {\n";
        j << "    \"d_max_5prime\": " << std::setprecision(6) << dp.d_max_5prime << ",\n";
        j << "    \"d_max_3prime\": " << dp.d_max_3prime << ",\n";
        j << "    \"d_max_combined\": " << dp.d_max_combined << ",\n";
        j << "    \"d_metamatch\": " << dp.d_metamatch << ",\n";
        j << "    \"source\": \"" << dp.d_max_source_str() << "\",\n";
        j << "    \"lambda_5prime\": " << dp.lambda_5prime << ",\n";
        j << "    \"lambda_3prime\": " << dp.lambda_3prime << ",\n";
        j << "    \"bg_5prime\": " << dp.fit_baseline_5prime << ",\n";
        j << "    \"bg_3prime\": " << dp.fit_baseline_3prime << ",\n";
        j << "    \"validated\": " << (dp.damage_validated ? "true" : "false") << ",\n";
        j << "    \"artifact\": " << (dp.damage_artifact ? "true" : "false") << ",\n";
        // Context-split CpG-like estimates
        auto nan_or = [](float v) -> std::string {
            return std::isnan(v) ? "null" : std::to_string(v);
        };
        j << "    \"cpg_like\": {\n";
        j << "      \"dmax_ct5_cpg\": "         << nan_or(dp.dmax_ct5_cpg_like)    << ",\n";
        j << "      \"dmax_ct5_noncpg\": "      << nan_or(dp.dmax_ct5_noncpg_like) << ",\n";
        j << "      \"cpg_ratio\": "            << nan_or(dp.cpg_ratio)            << ",\n";
        j << "      \"log2_cpg_ratio\": "       << nan_or(dp.log2_cpg_ratio)       << ",\n";
        j << "      \"baseline_cpg\": "         << nan_or(dp.fit_baseline_ct5_cpg_like)    << ",\n";
        j << "      \"baseline_noncpg\": "      << nan_or(dp.fit_baseline_ct5_noncpg_like) << ",\n";
        j << "      \"cov_terminal_cpg\": "     << std::setprecision(0) << dp.cov_ct5_cpg_like_terminal    << ",\n";
        j << "      \"cov_terminal_noncpg\": "  << dp.cov_ct5_noncpg_like_terminal << ",\n";
        j << "      \"effcov_terminal_cpg\": "  << dp.effcov_ct5_cpg_like_terminal    << ",\n";
        j << "      \"effcov_terminal_noncpg\": " << dp.effcov_ct5_noncpg_like_terminal << "\n";
        j << "    },\n";
        // Experimental: upstream-context-aware deamination (AC, CC, GC, TC)
        j << "    \"context_deamination\": {\n";
        j << "      \"dmax_AC\": " << nan_or(dp.dmax_ct5_by_upstream[dart::SampleDamageProfile::CTX_AC]) << ",\n";
        j << "      \"dmax_CC\": " << nan_or(dp.dmax_ct5_by_upstream[dart::SampleDamageProfile::CTX_CC]) << ",\n";
        j << "      \"dmax_GC\": " << nan_or(dp.dmax_ct5_by_upstream[dart::SampleDamageProfile::CTX_GC]) << ",\n";
        j << "      \"dmax_TC\": " << nan_or(dp.dmax_ct5_by_upstream[dart::SampleDamageProfile::CTX_TC]) << ",\n";
        j << "      \"dipyr_contrast\": " << nan_or(dp.dipyr_contrast) << ",\n";
        j << "      \"cpg_contrast\": " << nan_or(dp.cpg_contrast) << ",\n";
        j << "      \"heterogeneity_chi2\": " << std::setprecision(2) << dp.context_heterogeneity_chi2 << ",\n";
        j << "      \"heterogeneity_p\": " << std::setprecision(4) << dp.context_heterogeneity_p << ",\n";
        j << "      \"heterogeneity_detected\": " << (dp.context_heterogeneity_detected ? "true" : "false") << "\n";
        j << "    },\n";
        j << "    \"per_pos_5prime_ct\": [";
        for (int p = 0; p < N_POS; ++p) {
            double v = (dp.tc_total_5prime[p] >= MIN_COV) ? dp.t_freq_5prime[p] : -1.0;
            j << std::setprecision(6) << v; if (p < N_POS-1) j << ",";
        }
        j << "],\n";
        j << "    \"per_pos_3prime\": [";
        for (int p = 0; p < N_POS; ++p) {
            double v;
            if (is_ss) v = (dp.tc_total_3prime[p] >= MIN_COV) ? dp.t_freq_3prime[p] / dp.tc_total_3prime[p] : -1.0;
            else        v = (dp.ag_total_3prime[p]  >= MIN_COV) ? dp.a_freq_3prime[p] : -1.0;
            j << std::setprecision(6) << v; if (p < N_POS-1) j << ",";
        }
        j << "]\n";
        j << "  },\n";
        // complement_asymmetry: model-based G→T analysis.
        // GT(p) = A*exp(-mu*p) + B; B = uniform background (8-oxoG); A = terminal artifact.
        // s_gt = B - A/(A+C)_interior (Chargaff contrast; meaningful for SS, ~0 for DS).
        j << "  \"complement_asymmetry\": {\n";
        j << "    \"D\": " << std::setprecision(6) << dp.ox_gt_asymmetry << ",\n";
        j << "    \"tg_interior\": " << dp.ox_gt_baseline << ",\n";
        j << "    \"ac_interior\": " << dp.ox_ca_baseline << ",\n";
        j << "    \"tg_terminal\": " << dp.ox_gt_rate_terminal << ",\n";
        j << "    \"ac_terminal\": " << dp.ox_ca_rate_terminal << ",\n";
        j << "    \"gt_bg_fitted\": " << dp.g_bg_fitted << ",\n";
        j << "    \"gt_term_fitted\": " << dp.g_term_fitted << ",\n";
        j << "    \"gt_decay_fitted\": " << dp.g_decay_fitted << ",\n";
        j << "    \"s_gt\": " << dp.s_gt << ",\n";
        j << "    \"per_pos_5prime_gt\": [";
        for (int p = 0; p < N_POS; ++p) {
            double denom = dp.t_from_g_5prime[p] + dp.g_count_5prime[p];
            double v = (denom >= MIN_COV) ? dp.t_from_g_5prime[p] / denom : -1.0;
            j << std::setprecision(6) << v; if (p < N_POS-1) j << ",";
        }
        j << "],\n";
        j << "    \"channel_c_detected\": " << (dp.ox_damage_detected ? "true" : "false") << ",\n";
        j << "    \"s_oxog\": " << std::setprecision(6) << (has_oxog_score ? S_oxog : 0.0) << ",\n";
        j << "    \"se_s_oxog\": " << (has_oxog_score ? SE_oxog : 0.0) << ",\n";
        j << "    \"s_oxog_16ctx\": [";
        for (int i = 0; i < 16; ++i) {
            float v = dp.s_oxog_16ctx[i];
            if (std::isnan(v)) j << "null";
            else j << std::setprecision(6) << v;
            if (i < 15) j << ",";
        }
        j << "],\n";
        j << "    \"cov_oxog_16ctx\": [";
        for (int i = 0; i < 16; ++i) {
            j << std::setprecision(0) << dp.cov_oxog_16ctx[i];
            if (i < 15) j << ",";
        }
        j << "]\n";
        j << "  },\n";
        j << "  \"interior_ct_cluster\": {\n";
        j << "    \"short_asym_log2oe\": " << std::setprecision(6) << dp.interior_ct_cluster_short_asym_log2oe << ",\n";
        j << "    \"short_log2oe\": " << dp.interior_ct_cluster_short_log2oe << ",\n";
        j << "    \"short_z\": " << dp.interior_ct_cluster_short_z << ",\n";
        j << "    \"short_obs\": " << dp.interior_ct_cluster_short_obs << ",\n";
        j << "    \"short_exp\": " << dp.interior_ct_cluster_short_exp << ",\n";
        j << "    \"reads_used\": " << dp.interior_ct_cluster_reads_used << ",\n";
        j << "    \"sep_log2oe\": [";
        for (int i = 0; i < 10; ++i) {
            j << std::setprecision(6) << dp.interior_ct_cluster_sep_log2oe[i];
            if (i < 9) j << ",";
        }
        j << "]\n";
        j << "  },\n";
        j << "  \"depurination\": {\n";
        j << "    \"detected\": " << (dp.depurination_detected ? "true" : "false") << ",\n";
        j << "    \"enrichment_5prime\": " << std::setprecision(6) << dp.purine_enrichment_5prime << ",\n";
        j << "    \"enrichment_3prime\": " << dp.purine_enrichment_3prime << ",\n";
        j << "    \"rate_interior\": " << dp.purine_rate_interior << "\n";
        j << "  },\n";
        j << "  \"bic\": {\n";
        j << "    \"bias\": " << std::setprecision(2) << dp.library_bic_bias << ",\n";
        j << "    \"ds\": " << dp.library_bic_ds << ",\n";
        j << "    \"ss\": " << dp.library_bic_ss << ",\n";
        j << "    \"ct5_amp\": " << std::setprecision(6) << dp.libtype_amp_ct5 << ",\n";
        j << "    \"ga3_amp\": " << dp.libtype_amp_ga3 << ",\n";
        j << "    \"ga0_amp\": " << dp.libtype_amp_ga0 << ",\n";
        j << "    \"ct3_amp\": " << dp.libtype_amp_ct3 << "\n";
        j << "  }\n";
        j << "}\n";
        std::cout << "JSON written: " << json_path << "\n";
    }

    return 0;
}
