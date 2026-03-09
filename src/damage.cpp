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

#include <climits>
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

// ---- helpers -----------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " damage -i FILE [options]\n"
        << "\nOptions:\n"
        << "  -i FILE                    Input FASTQ (.gz or plain)\n"
        << "  -p N                       Worker threads (default: all cores)\n"
        << "  --library-type auto|ds|ss  Library type for 3'-end interpretation (default: auto)\n"
        << "  --mask-threshold FLOAT     Mask when excess P(deam) > T (default: 0.05)\n"
        << "  --tsv FILE                 Write per-position table as TSV\n";
}

// ---- main --------------------------------------------------------------
int damage_main(int argc, char** argv) {
    std::string in_path;
    int         n_threads     = static_cast<int>(std::thread::hardware_concurrency());
    if (n_threads < 1) n_threads = 1;
    double      mask_threshold = 0.05;
    std::string tsv_path;
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
        auto reader = make_fastq_reader(in_path);
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

    // ---- compute per-position mask -------------------------------------
    // After finalize_sample_profile():
    //   t_freq_5prime[p]  = T/(T+C) frequency (already normalised)
    //   a_freq_3prime[p]  = A/(A+G) frequency (already normalised)
    //   t_freq_3prime[p]  = raw T count at 3'  (NOT normalised)
    const bool   is_ss = (dp.library_type == dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
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

    // ---- per-position table --------------------------------------------
    const char* col3_hdr = is_ss ? "3'_CT  " : "3'_GA  ";
    std::cout << "pos  5'_CT   5'_GA   " << col3_hdr << "\n";
    std::cout << "---  ------  ------  ------\n";
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
        std::cout << std::setw(3) << p
                  << "  " << std::fixed << std::setprecision(4) << freq5_ct
                  << "  " << freq5_ga
                  << "  " << freq3;
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
        tsv << "pos\tend5\tfreq5\tend3\tfreq3\tmask\tcov5\tcov3\n";
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
            tsv << p
                << "\t5prime\t" << std::fixed << std::setprecision(6) << freq5
                << "\t3prime\t" << freq3
                << "\t" << (mask_pos[p] ? "1" : "0")
                << "\t" << static_cast<int64_t>(cov5)
                << "\t" << static_cast<int64_t>(cov3)
                << "\n";
        }
        std::cout << "TSV written: " << tsv_path << "\n";
    }

    return 0;
}
