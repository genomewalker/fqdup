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
#include "fqdup/damage_profile.hpp"

#include "taph/frame_selector_decl.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/sample_damage_profile.hpp"

#include <sstream>

#include <algorithm>
#include <array>
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
    taph::SampleDamageProfile           profile;
    std::array<uint64_t, LSD_HIST_BINS> lsd_hist{};
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
            if (L < LSD_L_MIN) { ++state.reads_skipped; continue; }
            taph::FrameSelector::update_sample_profile(state.profile, seq);
            ++state.lsd_hist[lsd_hist_bin(L)];
            if (L < state.len_min) state.len_min = L;
            if (L > state.len_max) state.len_max = L;
            state.len_sum += L;
            ++state.reads_scanned;
        }
    }
}

// ---- clip pass: re-profile after stripping 5' and/or 3' adapter stubs ------
static void clip_worker_fn(WorkQueue& queue, WorkerState& state,
                           const std::vector<std::string>& stubs5,
                           const std::vector<std::string>& stubs3) {
    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (std::string& seq : batch) {
            if (!stubs5.empty() && (int)seq.size() >= 6) {
                for (const auto& stub : stubs5) {
                    if (seq.compare(0, 6, stub) == 0) { seq.erase(0, 6); break; }
                }
            }
            if (!stubs3.empty()) {
                bool trimmed;
                do {
                    trimmed = false;
                    int L = static_cast<int>(seq.size());
                    if (L < 12) break;
                    for (const auto& stub : stubs3) {
                        if (seq.compare(L - 6, 6, stub) == 0) {
                            seq.erase(L - 6, 6);
                            trimmed = true;
                            break;
                        }
                    }
                } while (trimmed);
            }
            int L = static_cast<int>(seq.size());
            if (L < LSD_L_MIN) { ++state.reads_skipped; continue; }
            taph::FrameSelector::update_sample_profile(state.profile, seq);
            ++state.lsd_hist[lsd_hist_bin(L)];
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

static float per_read_damaged_weight(
    const std::string& seq,
    const taph::SampleDamageProfile& dp)
{
    int L = static_cast<int>(seq.size());
    int bin = taph::SampleDamageProfile::get_gc_bin(seq);
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

static constexpr int N_GC = taph::SampleDamageProfile::N_GC_BINS;

struct OxBinAcc {
    double A = 0, B = 0, C = 0, D = 0;           // first moments
    double AA = 0, AB = 0, BB = 0;                 // second moments (anc)
    double CC = 0, CD = 0, DD = 0;                 // second moments (bg)
    double AC = 0, AD = 0, BC = 0, BD = 0;         // cross moments
    double M_total = 0;
    uint64_t reads = 0;
    // Orientation-corrected Chargaff contrast (Option 2: D_oriented).
    // Accumulated unweighted from all oriented reads; not split by q.
    double T_all = 0, TG_all = 0;   // T count and T+G total in middle third
    double C_all = 0, CA_all = 0;   // C count and C+A total in middle third
};

struct OxogWorkerState {
    OxBinAcc bins[N_GC];
};

static void oxog_worker(WorkQueue& queue,
                        OxogWorkerState& state,
                        const taph::SampleDamageProfile& dp,
                        bool is_ss)
{
    std::vector<std::string> batch;
    while (queue.pop(batch)) {
        for (const std::string& seq : batch) {
            int L = static_cast<int>(seq.size());
            if (L < 30) continue;
            int bin = taph::SampleDamageProfile::get_gc_bin(seq);
            std::string oriented;
            double q;
            if (is_ss) {
                oriented = seq;
                q = per_read_damaged_weight(seq, dp);
            } else {
                float w_fwd = per_read_damaged_weight(seq, dp);
                std::string rc = revcomp(seq);
                float w_rev = per_read_damaged_weight(rc, dp);
                if (w_rev > w_fwd) { oriented = std::move(rc); q = w_rev; }
                else               { oriented = seq;            q = w_fwd; }
            }
            int beg = L / 3, end = L - (L / 3);
            if (end <= beg) continue;
            int T = 0, G = 0, Cv = 0, Av = 0;
            for (int j = beg; j < end; ++j) {
                char b = static_cast<char>(std::toupper(
                             static_cast<unsigned char>(oriented[j])));
                if      (b == 'T') ++T;
                else if (b == 'G') ++G;
                else if (b == 'C') ++Cv;
                else if (b == 'A') ++Av;
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
            x.T_all  += T;  x.TG_all += M;
            x.C_all  += Cv; x.CA_all += Cv + Av;
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
        << "  --json FILE                Write full damage profile as JSON\n"
        << "  --length-bins SPEC         Length-stratified damage: auto | N | e1,e2,... (default: off)\n";
}

static LengthBinOptions parse_length_bins(const std::string& spec, bool& ok) {
    LengthBinOptions opt;
    ok = true;
    if (spec.empty() || spec == "off" || spec == "none" || spec == "disabled") {
        opt.mode = LengthBinOptions::Mode::DISABLED;
        return opt;
    }
    if (spec == "auto") {
        opt.mode = LengthBinOptions::Mode::AUTO;
        return opt;
    }
    if (spec.find(',') != std::string::npos) {
        opt.mode = LengthBinOptions::Mode::EXPLICIT;
        std::stringstream ss(spec);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            if (tok.empty()) continue;
            try { opt.explicit_edges.push_back(std::stoi(tok)); }
            catch (...) { ok = false; return opt; }
        }
        if (opt.explicit_edges.empty()) ok = false;
        return opt;
    }
    try {
        int n = std::stoi(spec);
        if (n < 1) { ok = false; return opt; }
        opt.mode = LengthBinOptions::Mode::QUANTILE;
        opt.quantile_bins = n;
        return opt;
    } catch (...) { ok = false; return opt; }
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
    LengthBinOptions lb_opts = []{ bool ok; return parse_length_bins("auto", ok); }();
    taph::SampleDamageProfile::LibraryType forced_lib =
        taph::SampleDamageProfile::LibraryType::UNKNOWN;

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
                forced_lib = taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                forced_lib = taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
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
        } else if (arg == "--length-bins" && i + 1 < argc) {
            bool ok = true;
            lb_opts = parse_length_bins(argv[++i], ok);
            if (!ok) {
                std::cerr << "Error: invalid --length-bins: " << argv[i] << "\n";
                return 1;
            }
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

    // ---- pre-scan: first PRE_SCAN_READS reads to detect adapter stubs ----
    // Single-threaded; only needs hexamer_count_5prime / n_hexamers_*.
    // Cost: <5s on any file. Result: clip_stubs populated before full pass.
    static constexpr int64_t PRE_SCAN_READS = 500000;

    taph::AdapterStubs stubs;

    {
        WorkerState scan_state;
        std::array<uint32_t, 4096> hex3_terminal{};
        uint64_t n_hex3 = 0;
        auto reader_pre = make_fastq_reader(in_path, 1);
        FastqRecord rec;
        int64_t n = 0;
        while (reader_pre->read(rec) && n < PRE_SCAN_READS) {
            int L = static_cast<int>(rec.seq.size());
            if (L < LSD_L_MIN) continue;
            taph::FrameSelector::update_sample_profile(scan_state.profile, rec.seq);
            if (L >= 12) {
                int code = taph::encode_hex_at(rec.seq, L - 6);
                if (code >= 0) { ++hex3_terminal[code]; ++n_hex3; }
            }
            ++n;
        }
        stubs = taph::detect_adapter_stubs(scan_state.profile, hex3_terminal.data(), n_hex3);
    }

    // ---- full pass: profile (with optional 5' clip) --------------------
    WorkQueue queue(2 * n_threads);
    std::vector<WorkerState> states(n_threads);
    std::vector<std::thread> workers;
    workers.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        if (stubs.adapter_clipped || stubs.adapter3_clipped)
            workers.emplace_back(clip_worker_fn, std::ref(queue),
                                 std::ref(states[t]),
                                 std::cref(stubs.stubs5), std::cref(stubs.stubs3));
        else
            workers.emplace_back(worker_fn, std::ref(queue), std::ref(states[t]));
    }

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
    taph::SampleDamageProfile dp;
    dp.forced_library_type = forced_lib;
    int64_t reads_scanned = 0, reads_skipped = 0;
    int     len_min = INT_MAX, len_max = 0;
    int64_t len_sum = 0;

    std::vector<uint64_t> merged_lsd_hist(LSD_HIST_BINS, 0);

    for (auto& s : states) {
        taph::FrameSelector::merge_sample_profiles(dp, s.profile);
        for (int hb = 0; hb < LSD_HIST_BINS; ++hb)
            merged_lsd_hist[hb] += s.lsd_hist[hb];
        reads_scanned += s.reads_scanned;
        reads_skipped += s.reads_skipped;
        if (s.len_min < len_min) len_min = s.len_min;
        if (s.len_max > len_max) len_max = s.len_max;
        len_sum += s.len_sum;
    }
    dp.forced_library_type = forced_lib;  // ensure forced_lib survives merge
    taph::FrameSelector::finalize_sample_profile(dp);

    // Recompute top enriched hexamers from final (possibly clipped) dp for JSON output.
    stubs.top_enriched = taph::compute_hex_enriched_5prime(dp);
    stubs.flag_hex_artifact = !stubs.top_enriched.empty()
        && stubs.top_enriched[0].log2fc > 1.5
        && !stubs.top_enriched[0].damage_consistent;

    const bool is_ss = (dp.library_type == taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);

    // ---- second pass: oxidation-compatible composition score -----------
    // Run when pass-1 found non-trivial 5' damage (gives meaningful q weights).
    double S_oxog = 0.0, SE_oxog = 0.0;
    double D_oriented = 0.0;
    double oxog_score_z = 0.0, oxog_score_p = 1.0;
    bool has_oxog_score = false;
    double depur_score_z_5    = 0.0;
    double depur_score_z_3    = std::numeric_limits<double>::quiet_NaN();
    double depur_score_z      = 0.0;
    double depur_score_p      = 1.0;
    double depur_ctrl_shift_5 = 0.0;
    double depur_ctrl_shift_3 = 0.0;
    double oxog_trinuc_cosine = std::numeric_limits<double>::quiet_NaN();

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
                merged[b].T_all   += os.bins[b].T_all;
                merged[b].TG_all  += os.bins[b].TG_all;
                merged[b].C_all   += os.bins[b].C_all;
                merged[b].CA_all  += os.bins[b].CA_all;
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

        // Orientation-corrected Chargaff contrast D_oriented = T/(T+G) - C/(C+A)
        // on oriented reads (DS rotated by max-q orientation, SS as-is).
        // Unlike the dart D field, this breaks the DS Chargaff cancellation because
        // both strands are measured from the same oriented direction after rotation.
        double w_dor = 0.0;
        for (int b = 0; b < N_GC; ++b) {
            const auto& x = merged[b];
            if (x.TG_all > 0 && x.CA_all > 0) {
                double w = x.TG_all;
                D_oriented += w * (x.T_all / x.TG_all - x.C_all / x.CA_all);
                w_dor += w;
            }
        }
        if (w_dor > 0) D_oriented /= w_dor;
    }

    // ---- Score tests: oxog interior, depurination, damage mask (→ libtaph) ----
    {
        auto oxog_r = taph::compute_oxog_interior_score(dp);
        oxog_score_z = oxog_r.z;
        oxog_score_p = oxog_r.p;
    }
    auto depur = taph::compute_depur_score(dp, is_ss);
    depur_score_z_5    = depur.z5;
    depur_score_z_3    = depur.z3;
    depur_score_z      = depur.z;
    depur_score_p      = depur.p;
    depur_ctrl_shift_5 = depur.shift5;
    depur_ctrl_shift_3 = depur.shift3;

    auto dmask = taph::compute_damage_mask(dp, is_ss, mask_threshold,
                                           static_cast<int>(MIN_COV));

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
              << "  bg=" << std::setprecision(4) << dp.fit_baseline_5prime << "\n";
    std::cout << "3'-end   d_max=" << std::fixed << std::setprecision(4) << dp.d_max_3prime
              << "  lambda=" << std::setprecision(3) << dp.lambda_3prime
              << "  bg=" << std::setprecision(4) << dp.fit_baseline_3prime << "\n";
    std::cout << "combined d_max=" << std::fixed << std::setprecision(4) << dp.d_max_combined
              << "  (source=" << dp.d_max_source_str() << ")";
    if (dp.mixture_converged) std::cout << " [mixture]";
    if (dp.damage_validated)  std::cout << " [validated]";
    if (dp.damage_artifact)   std::cout << " [ARTIFACT]";
    std::cout << "\n\n";

    std::cout << "Mask threshold: " << mask_threshold << " → "
              << dmask.n_masked << " position" << (dmask.n_masked != 1 ? "s" : "") << " masked";
    if (dmask.n_masked) std::cout << " (pos " << dmask.masked_str << ")";
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
    std::cout << "  score_z=" << std::showpos << std::setprecision(2)
              << depur_score_z << std::noshowpos
              << "  (z5=" << std::showpos << depur_score_z_5 << std::noshowpos;
    if (!std::isnan(depur_score_z_3))
        std::cout << "  z3=" << std::showpos << std::setprecision(2)
                  << depur_score_z_3 << std::noshowpos;
    std::cout << ")\n\n";

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
        if (dmask.pos[p]) std::cout << "  *";
        std::cout << "\n";
    }

    // ---- optional length-stratified damage ----------------------------------------
    DamageProfile bulk_dp;
    bulk_dp.enabled              = true;
    bulk_dp.d_max_5prime         = dp.d_max_5prime;
    bulk_dp.d_max_3prime         = dp.d_max_3prime;
    bulk_dp.lambda_5prime        = dp.lambda_5prime;
    bulk_dp.lambda_3prime        = dp.lambda_3prime;
    bulk_dp.background           = dp.fit_baseline_5prime;
    bulk_dp.bg_5_tc              = dp.fit_baseline_5prime;
    bulk_dp.bg_3_channel         = dp.fit_baseline_3prime;
    bulk_dp.ss_mode              = is_ss;
    bulk_dp.mixture_d_damaged    = dp.mixture_d_ancient;
    bulk_dp.mixture_pi_damaged   = dp.mixture_pi_ancient;
    bulk_dp.mixture_d_reference  = dp.mixture_d_reference;
    bulk_dp.mixture_n_components = dp.mixture_n_components;
    bulk_dp.mixture_converged    = dp.mixture_converged;
    bulk_dp.mixture_identifiable = dp.mixture_identifiable;
    bulk_dp.d_cpg_5prime         = dp.dmax_ct5_cpg_like;
    bulk_dp.d_noncpg_5prime      = dp.dmax_ct5_noncpg_like;

    LengthStratifiedDamageProfile lsd;
    if (lb_opts.enabled()) {
        lsd = estimate_damage_by_length(in_path, forced_lib, lb_opts,
                                        &merged_lsd_hist,
                                        static_cast<size_t>(n_threads), 0, &bulk_dp);
        std::cout << "\nLength-stratified damage (method=" << lsd.method
                  << ", n_bins=" << lsd.bins.size() << "):\n";
        std::cout << "bin  L_lo  L_hi  n_reads  d5      d3      lam5    lam3    src\n";
        for (size_t b = 0; b < lsd.bins.size(); ++b) {
            const auto& lb = lsd.bins[b];
            std::cout << std::setw(3) << b
                      << "  " << std::setw(4) << lb.length_lo
                      << "  " << std::setw(4) << lb.length_hi
                      << "  " << std::setw(7) << lb.n_reads
                      << "  " << std::fixed << std::setprecision(4) << lb.d_max_5prime
                      << "  " << lb.d_max_3prime
                      << "  " << std::setprecision(3) << lb.lambda_5prime
                      << "  " << lb.lambda_3prime
                      << "  " << lb.source << "\n";
        }
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
                << "\t" << (dmask.pos[p] ? "1" : "0")
                << "\t" << static_cast<int64_t>(cov5)
                << "\t" << static_cast<int64_t>(cov3)
                << "\t" << static_cast<int64_t>(cov5_gt)
                << "\n";
        }
        std::cout << "TSV written: " << tsv_path << "\n";
    }

    // ---- Pre-compute calibrated score statistics used across JSON sections ----
    auto cpg_score  = taph::compute_cpg_score(dp);
    double cpg_score_z = cpg_score.z, cpg_score_p = cpg_score.p;

    auto hex_stats  = taph::compute_hex_stats(dp);
    double hex_shift_g = hex_stats.shift_g;
    double hex_shift_z = hex_stats.shift_z;
    double hex_shift_p = hex_stats.shift_p;

    // Short-read spike fraction (used by QC flags + JSON).
    double short_read_frac = -1.0;
    if (!lsd.bins.empty()) {
        uint64_t short_reads = 0, total_len_reads = 0;
        for (const auto& lb : lsd.bins) {
            total_len_reads += lb.n_reads;
            if (lb.length_hi <= 50) short_reads += lb.n_reads;
        }
        if (total_len_reads > 0)
            short_read_frac = static_cast<double>(short_reads) / total_len_reads;
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
        const char* ds_str = (dp.damage_status == taph::SampleDamageProfile::DamageStatus::PRESENT) ? "present"
                           : (dp.damage_status == taph::SampleDamageProfile::DamageStatus::WEAK)    ? "weak"
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
        j << "    \"mixture_gc\": {\n";
        j << "      \"d_ancient\": " << dp.mixture_d_ancient << ",\n";
        j << "      \"pi_ancient\": " << dp.mixture_pi_ancient << ",\n";
        j << "      \"d_reference\": " << dp.mixture_d_reference << ",\n";
        j << "      \"identifiable\": " << (dp.mixture_identifiable ? "true" : "false") << ",\n";
        j << "      \"converged\": " << (dp.mixture_converged ? "true" : "false") << ",\n";
        j << "      \"n_components\": " << dp.mixture_n_components << "\n";
        j << "    },\n";
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
        j << "      \"effcov_terminal_noncpg\": " << dp.effcov_ct5_noncpg_like_terminal << ",\n";
        j << "      \"cpg_score_z\": " << std::setprecision(6) << cpg_score_z << ",\n";
        j << "      \"cpg_score_p\": " << cpg_score_p << "\n";
        j << "    },\n";
        // Experimental: upstream-context-aware deamination (AC, CC, GC, TC)
        j << "    \"context_deamination\": {\n";
        j << "      \"dmax_AC\": " << nan_or(dp.dmax_ct5_by_upstream[taph::SampleDamageProfile::CTX_AC]) << ",\n";
        j << "      \"dmax_CC\": " << nan_or(dp.dmax_ct5_by_upstream[taph::SampleDamageProfile::CTX_CC]) << ",\n";
        j << "      \"dmax_GC\": " << nan_or(dp.dmax_ct5_by_upstream[taph::SampleDamageProfile::CTX_GC]) << ",\n";
        j << "      \"dmax_TC\": " << nan_or(dp.dmax_ct5_by_upstream[taph::SampleDamageProfile::CTX_TC]) << ",\n";
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
        j << "],\n";
        j << "    \"by_length_method\": \"" << lsd.method << "\",\n";
        j << "    \"by_length\": [";
        for (size_t b = 0; b < lsd.bins.size(); ++b) {
            const auto& lb = lsd.bins[b];
            if (b > 0) j << ",";
            j << "\n      {";
            j << "\"length_lo\":" << lb.length_lo
              << ",\"length_hi\":" << lb.length_hi
              << ",\"n_reads\":" << lb.n_reads
              << ",\"d_max_5prime\":" << std::setprecision(6) << lb.d_max_5prime
              << ",\"d_max_3prime\":" << lb.d_max_3prime
              << ",\"lambda_5prime\":" << lb.lambda_5prime
              << ",\"lambda_3prime\":" << lb.lambda_3prime
              << ",\"bg_5prime\":" << lb.bg_5prime
              << ",\"bg_3prime\":" << lb.bg_3prime
              << ",\"cpg_contrast\":"
              << (std::isnan(lb.cpg_contrast) ? std::string("null") : std::to_string(lb.cpg_contrast))
              << ",\"validated\":" << (lb.validated ? "true" : "false")
              << ",\"ss_mode\":" << (lb.ss_mode ? "true" : "false")
              << ",\"source\":\"" << lb.source << "\"";
            j << ",\"per_pos_5prime_ct\":[";
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                j << std::setprecision(6) << lb.per_pos_5prime_ct[p];
                if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
            }
            j << "],\"per_pos_3prime\":[";
            for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                j << std::setprecision(6) << lb.per_pos_3prime[p];
                if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
            }
            j << "]";
            j << ",\"mixture\":{\"d_ancient\":" << std::setprecision(6) << lb.mixture_d_damaged
              << ",\"d_reference\":" << lb.mixture_d_reference
              << ",\"d_population\":" << lb.mixture_d_population
              << ",\"pi_ancient\":" << lb.mixture_pi_damaged
              << ",\"n_components\":" << lb.mixture_n_components
              << ",\"converged\":" << (lb.mixture_converged ? "true" : "false")
              << ",\"identifiable\":" << (lb.mixture_identifiable ? "true" : "false")
              << "}";
            j << ",\"gc_d_max\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << std::setprecision(6) << lb.gc_d_max[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_n_reads\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << lb.gc_n_reads[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "],\"gc_p_damaged\":[";
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << std::setprecision(6) << lb.gc_p_damaged[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "]";
            auto dnan = [](double v) -> std::string {
                return std::isnan(v) ? "null" : std::to_string(v);
            };
            auto write_dnan_arr = [&](const auto& arr) {
                for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                    j << dnan(arr[p]);
                    if (p + 1 < LengthBinDamageProfile::N_POS) j << ",";
                }
            };
            j << ",\"llr\":{\"n_damaged\":" << lb.n_damaged
              << ",\"n_undamaged\":" << lb.n_undamaged
              << ",\"d_max_5_damaged\":"          << dnan(lb.d_max_5_damaged)
              << ",\"d_max_3_damaged\":"          << dnan(lb.d_max_3_damaged)
              << ",\"d_max_5_cpg_damaged\":"      << dnan(lb.d_max_5_cpg_damaged)
              << ",\"d_max_5_noncpg_damaged\":"   << dnan(lb.d_max_5_noncpg_damaged)
              << ",\"s_gt_5_damaged_vs_undamaged\":" << dnan(lb.s_gt_5_damaged_vs_undamaged)
              << ",\"g_to_t_5_damaged\":"         << dnan(lb.g_to_t_5_damaged)
              << ",\"pG_terminal_5_damaged\":"    << dnan(lb.pG_terminal_5_damaged)
              << ",\"pG_interior_5_damaged\":"    << dnan(lb.pG_interior_5_damaged)
              << ",\"per_pos_5prime_ct_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_ct_damaged);
            j << "],\"per_pos_3prime_damaged\":[";
            write_dnan_arr(lb.per_pos_3prime_damaged);
            j << "],\"per_pos_5prime_ct_cpg_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_ct_cpg_damaged);
            j << "],\"per_pos_5prime_ct_noncpg_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_ct_noncpg_damaged);
            j << "],\"per_pos_5prime_gt_damaged\":[";
            write_dnan_arr(lb.per_pos_5prime_gt_damaged);
            j << "],\"per_pos_5prime_gt_undamaged\":[";
            write_dnan_arr(lb.per_pos_5prime_gt_undamaged);
            j << "]}";
            j << ",\"trinuc\":{\"tri_5prime_terminal\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_5prime_terminal[i]; if (i < 63) j << ","; }
            j << "],\"tri_5prime_interior\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_5prime_interior[i]; if (i < 63) j << ","; }
            j << "],\"tri_3prime_terminal\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_3prime_terminal[i]; if (i < 63) j << ","; }
            j << "],\"tri_3prime_interior\":[";
            for (int i = 0; i < 64; ++i) { j << lb.tri_3prime_interior[i]; if (i < 63) j << ","; }
            j << "]}";
            j << "}";
        }
        if (!lsd.bins.empty()) j << "\n    ";
        j << "],\n";
        j << "    \"by_length_joint\": {\n";
        j << "      \"d_ancient\": "   << std::setprecision(6) << lsd.d_joint_ancient << ",\n";
        j << "      \"pi_ancient\": "  << lsd.pi_joint_ancient << ",\n";
        j << "      \"d_population\": " << lsd.d_joint_population << ",\n";
        j << "      \"converged\": "   << (lsd.joint_converged ? "true" : "false") << ",\n";
        j << "      \"separated\": "   << (lsd.joint_separated ? "true" : "false") << ",\n";
        j << "      \"cell_w_ancient\": [";
        for (size_t b = 0; b < lsd.cell_w_ancient.size(); ++b) {
            if (b > 0) j << ",";
            j << "[";
            const auto& row = lsd.cell_w_ancient[b];
            for (int g = 0; g < LengthBinDamageProfile::N_GC_BINS; ++g) {
                j << std::setprecision(6) << row[g];
                if (g + 1 < LengthBinDamageProfile::N_GC_BINS) j << ",";
            }
            j << "]";
        }
        j << "]\n    }\n";
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
        j << "    \"channel_c_valid\": " << (dp.channel_c_valid ? "true" : "false") << ",\n";
        j << "    \"channel_c_detected\": " << (dp.ox_damage_detected ? "true" : "false") << ",\n";
        j << "    \"ox_is_artifact\": " << (dp.ox_is_artifact ? "true" : "false") << ",\n";
        j << "    \"ox_d_max\": " << std::setprecision(6) << dp.ox_d_max << ",\n";
        j << "    \"ox_stop_rate_terminal\": " << dp.ox_stop_rate_terminal << ",\n";
        j << "    \"ox_stop_rate_interior\": " << dp.ox_stop_rate_interior << ",\n";
        j << "    \"ox_stop_rate_baseline\": " << dp.ox_stop_conversion_rate_baseline << ",\n";
        j << "    \"ox_uniformity_ratio\": " << dp.ox_uniformity_ratio << ",\n";
        j << "    \"ox_gt_rate_interior\": " << dp.ox_gt_rate_interior << ",\n";
        j << "    \"ox_gt_uniformity\": " << dp.ox_gt_uniformity << ",\n";
        j << "    \"ox_ca_rate_interior\": " << dp.ox_ca_rate_interior << ",\n";
        j << "    \"ox_ca_uniformity\": " << dp.ox_ca_uniformity << ",\n";
        j << "    \"s_oxog\": " << std::setprecision(6) << (has_oxog_score ? S_oxog : 0.0) << ",\n";
        j << "    \"se_s_oxog\": " << (has_oxog_score ? SE_oxog : 0.0) << ",\n";
        j << "    \"ox_d_oriented\": " << std::setprecision(6) << D_oriented << ",\n";
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
        j << "],\n";
        // Multi-context composite (Option 1): authentic 8-oxoG produces broad
        // positive signal across all 16 NxGxN contexts; a prep artifact spikes
        // in 1-2 contexts.  oxog_ctx_z is coverage-weighted z across contexts.
        {
            double ctx_wsum = 0.0, ctx_w2 = 0.0;
            int ctx_pos = 0;
            for (int i = 0; i < 16; ++i) {
                float v = dp.s_oxog_16ctx[i];
                double cov = static_cast<double>(dp.cov_oxog_16ctx[i]);
                if (std::isnan(v) || cov <= 0) continue;
                double w = std::sqrt(cov);
                ctx_wsum += v * w;
                ctx_w2   += cov;
                if (v > 0) ++ctx_pos;
            }
            double ctx_z = (ctx_w2 > 0) ? ctx_wsum / std::sqrt(ctx_w2) : 0.0;
            bool multi_ctx = (ctx_pos >= 10) && (ctx_z > 3.0);
            j << "    \"oxog_ctx_n_positive\": " << ctx_pos << ",\n";
            j << "    \"oxog_ctx_z\": " << std::setprecision(4) << ctx_z << ",\n";
            j << "    \"oxog_score_z\": " << std::setprecision(6) << oxog_score_z << ",\n";
            j << "    \"oxog_score_p\": " << std::setprecision(6) << oxog_score_p << ",\n";
        }
        {
            auto otr = taph::compute_oxog_trinuc(dp);
            oxog_trinuc_cosine = otr.cosine;
            if (std::isnan(otr.cosine))
                j << "    \"oxog_trinuc_cosine\": null,\n";
            else
                j << "    \"oxog_trinuc_cosine\": " << std::setprecision(6) << otr.cosine << ",\n";
            j << "    \"oxog_trinuc_n_context\": " << otr.n_ctx << "\n";
        }
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
        j << "],\n";
        {
            double O = dp.interior_ct_cluster_short_obs;
            double E = dp.interior_ct_cluster_short_exp;
            double sz = 0.0, sp = 1.0;
            if (O > 0 && E > 1e-9) {
                double sign_val = (O >= E) ? 1.0 : -1.0;
                sz = sign_val * std::sqrt(2.0 * (O * std::log(O / E) - (O - E)));
                sp = 0.5 * std::erfc(sz / std::sqrt(2.0));
            }
            j << "    \"short_score_z\": " << std::setprecision(6) << sz << ",\n";
            j << "    \"short_score_p\": " << sp << "\n";
        }
        j << "  },\n";
        j << "  \"depurination\": {\n";
        j << "    \"detected\": " << (dp.depurination_detected ? "true" : "false") << ",\n";
        j << "    \"enrichment_5prime\": " << std::setprecision(6) << dp.purine_enrichment_5prime << ",\n";
        j << "    \"enrichment_3prime\": " << dp.purine_enrichment_3prime << ",\n";
        j << "    \"rate_interior\": " << std::setprecision(6) << dp.purine_rate_interior << ",\n";
        j << "    \"depur_ctrl_shift_5prime\": " << depur_ctrl_shift_5 << ",\n";
        j << "    \"depur_ctrl_shift_3prime\": " << depur_ctrl_shift_3 << ",\n";
        j << "    \"depur_score_z_5prime\": " << depur_score_z_5 << ",\n";
        j << "    \"depur_score_z_3prime\": "
          << (std::isnan(depur_score_z_3) ? "null" : std::to_string(depur_score_z_3)) << ",\n";
        j << "    \"depur_score_z\": " << depur_score_z << ",\n";
        j << "    \"depur_score_p\": " << depur_score_p << "\n";
        j << "  },\n";
        // Reference-free trinucleotide spectrum (bulk, 64 contexts).
        // Index = prev*16 + mid*4 + next, A=0,C=1,G=2,T=3.
        // Terminal = read pos 1..4 from that end; interior = pos 10..14.
        // Post-processing (e.g. in Python) can contrast the 5' and 3' terminal
        // counters against their interior counterparts to produce a 32-entry
        // terminal-minus-interior context contrast per end without any external
        // reference panel or trained model.
        {
            auto emit_arr = [&](const char* name,
                                const std::array<uint64_t, 64>& v,
                                bool trailing_comma) {
                j << "    \"" << name << "\": [";
                for (int i = 0; i < 64; ++i) {
                    j << v[i];
                    if (i < 63) j << ",";
                }
                j << "]" << (trailing_comma ? "," : "") << "\n";
            };
            j << "  \"trinuc_spectrum\": {\n";
            emit_arr("tri_5prime_terminal", dp.tri_5prime_terminal, true);
            emit_arr("tri_5prime_interior", dp.tri_5prime_interior, true);
            emit_arr("tri_3prime_terminal", dp.tri_3prime_terminal, true);
            emit_arr("tri_3prime_interior", dp.tri_3prime_interior, false);
            j << "  },\n";
        }
        auto pres = taph::compute_preservation_summary(dp, is_ss,
            stubs.adapter_clipped, stubs.flag_hex_artifact,
            cpg_score_z, oxog_score_z, oxog_trinuc_cosine, hex_shift_p);

        j << "  \"preservation\": {\n";
        j << "    \"score\": " << std::setprecision(6) << dp.preservation_score << ",\n";
        j << "    \"evidence\": " << dp.preservation_evidence << ",\n";
        j << "    \"reliability\": " << dp.preservation_reliability << ",\n";
        j << "    \"f5\": " << dp.preservation_f5 << ",\n";
        j << "    \"f3\": " << dp.preservation_f3 << ",\n";
        j << "    \"f_coh\": " << dp.preservation_f_coh << ",\n";
        j << "    \"f_cpg\": " << dp.preservation_f_cpg << ",\n";
        {
            j << "    \"authenticity_eff\": " << std::setprecision(6) << pres.authenticity_eff << ",\n";
            j << "    \"authenticity_evidence\": " << pres.authenticity_evidence << ",\n";
            if (pres.d5_was_corrected)
                j << "    \"d5_hexamer_corrected\": " << pres.d5_hexamer_corrected << ",\n";
            j << "    \"oxidation_eff\": " << pres.oxidation_eff << ",\n";
            j << "    \"oxidation_evidence\": " << pres.oxidation_evidence << ",\n";
            j << "    \"qc_risk_eff\": " << pres.qc_risk_eff << ",\n";
            j << "    \"qc_evidence\": " << pres.qc_evidence << ",\n";
            j << "    \"label\": \"" << pres.label << "\"\n";
        }
        j << "  },\n";
        {
            auto flags = taph::compute_library_qc_flags(dp, is_ss,
                stubs.flag_hex_artifact,
                hex_stats.jsd, hex_stats.entropy_terminal,
                short_read_frac);

            j << "  \"library_qc\": {\n";
            if (stubs.adapter_clipped) {
                j << "    \"adapter_stubs_clipped\": [";
                for (int i = 0; i < (int)stubs.stubs5.size(); ++i) {
                    if (i > 0) j << ",";
                    j << "\"" << stubs.stubs5[i] << "\"";
                }
                j << "],\n";
            }
            if (stubs.adapter3_clipped) {
                j << "    \"adapter_stubs_clipped_3prime\": [";
                for (int i = 0; i < (int)stubs.stubs3.size(); ++i) {
                    if (i > 0) j << ",";
                    j << "\"" << stubs.stubs3[i] << "\"";
                }
                j << "],\n";
            }
            j << "    \"adapter_offset_5prime\": " << dp.fit_offset_5prime << ",\n";
            j << "    \"adapter_offset_3prime\": " << dp.fit_offset_3prime << ",\n";
            j << "    \"position0_artifact_5prime\": " << (dp.position_0_artifact_5prime ? "true" : "false") << ",\n";
            j << "    \"position0_artifact_3prime\": " << (dp.position_0_artifact_3prime ? "true" : "false") << ",\n";
            j << "    \"hexamer_entropy_5prime\": " << std::setprecision(4) << hex_stats.entropy_terminal << ",\n";
            j << "    \"hexamer_entropy_interior\": " << hex_stats.entropy_interior << ",\n";
            j << "    \"hexamer_terminal_interior_jsd\": " << hex_stats.jsd << ",\n";
            j << "    \"hex_shift_g\": " << std::setprecision(4) << hex_shift_g << ",\n";
            j << "    \"hex_shift_z\": " << hex_shift_z << ",\n";
            j << "    \"hex_shift_p\": " << std::setprecision(6) << hex_shift_p << ",\n";
            j << "    \"hexamer_excess_tc\": " << std::setprecision(6) << dp.hexamer_excess_tc << ",\n";
            // Top overrepresented 5' hexamers vs interior.
            // damage_consistent=true: starts with T (C→T deamination product) — expected in aDNA.
            // damage_consistent=false with high log2fc: adapter remnant or ligation bias artifact.
            j << "    \"top_hexamers_5prime\": [";
            {
                int n_out = 0;
                for (const auto& hr : stubs.top_enriched) {
                    if (n_out >= 5) break;
                    auto seq = taph::decode_hex(hr.idx);
                    if (n_out > 0) j << ",";
                    j << "{\"seq\":\"" << seq.data() << "\","
                      << "\"log2fc\":" << std::setprecision(3) << hr.log2fc << ","
                      << "\"damage_consistent\":" << (hr.damage_consistent ? "true" : "false") << "}";
                    ++n_out;
                }
            }
            j << "],\n";
            {
                static const std::pair<const char*, const char*> kAdapters[] = {
                    {"ACACTC", "TruSeq/P5"},
                    {"AATGAT", "TruSeq/Universal"},
                    {"GATCGG", "TruSeq/i7"},
                    {"CTGTCT", "Nextera/Tn5"},
                    {"AGATCG", "TruSeq/R1"},
                    {"TGGAAT", "TruSeq/R2"},
                    {"GCGAAT", "TruSeq/R2alt"},
                };
                std::string top_seq = stubs.stubs5.empty() ? "" : stubs.stubs5[0];
                if (!top_seq.empty()) {
                    const char* name = "unknown";
                    for (const auto& kv : kAdapters)
                        if (top_seq == kv.first) { name = kv.second; break; }
                    j << "    \"adapter_prefix_identified\": {\"seq\":\"" << top_seq
                      << "\",\"name\":\"" << name << "\"},\n";
                }
            }
            j << "    \"depurination_detected\": " << (dp.depurination_detected ? "true" : "false") << ",\n";
            j << "    \"short_read_frac\": " << std::setprecision(4)
              << (short_read_frac < 0 ? 0.0 : short_read_frac) << ",\n";

            j << "    \"flags\": [";
            bool first_flag = true;
            auto emit_flag = [&](const char* name) {
                if (!first_flag) j << ",";
                j << "\"" << name << "\"";
                first_flag = false;
            };
            if (flags.adapter_remnant_5prime)   emit_flag("adapter_remnant_5prime");
            if (flags.adapter_remnant_3prime)   emit_flag("adapter_remnant_3prime");
            if (flags.hexamer_composition_bias) emit_flag("hexamer_composition_bias");
            if (flags.hexamer_terminal_shift)   emit_flag("hexamer_terminal_shift");
            if (flags.short_read_spike)         emit_flag("short_read_spike");
            if (flags.depurination)             emit_flag("depurination");
            if (flags.ds_3prime_signal_absent)  emit_flag("ds_3prime_signal_absent");
            if (flags.ga3_inward_displaced)     emit_flag("ga3_inward_displaced");
            if (flags.hexamer_artifact_bias)    emit_flag("hexamer_artifact_bias");
            j << "]\n";
            j << "  },\n";
        }
        {
            auto dcp = taph::compute_damage_context_profile(
                dp, cpg_score_z, hex_shift_z,
                stubs.adapter_clipped, stubs.adapter3_clipped,
                stubs.flag_hex_artifact);
            auto emit_score = [&](const char* name, float v, bool trailing) {
                j << "    \"" << name << "\": ";
                if (!std::isfinite(v)) j << "null";
                else                    j << std::setprecision(4) << v;
                j << (trailing ? ",\n" : "\n");
            };
            auto json_escape = [](const std::string& s) {
                std::string out; out.reserve(s.size() + 8);
                static const char kHex[] = "0123456789abcdef";
                for (unsigned char c : s) {
                    switch (c) {
                        case '"':  out += "\\\""; break;
                        case '\\': out += "\\\\"; break;
                        case '\b': out += "\\b";  break;
                        case '\f': out += "\\f";  break;
                        case '\n': out += "\\n";  break;
                        case '\r': out += "\\r";  break;
                        case '\t': out += "\\t";  break;
                        default:
                            if (c < 0x20) {
                                out += "\\u00";
                                out += kHex[(c >> 4) & 0xF];
                                out += kHex[c & 0xF];
                            } else {
                                out += static_cast<char>(c);
                            }
                    }
                }
                return out;
            };
            j << "  \"damage_context_profile\": {\n";
            j << "    \"method\": \"" << dcp.method << "\",\n";
            j << "    \"reference_required\": "
              << (dcp.reference_required ? "true" : "false") << ",\n";
            j << "    \"alignment_required\": "
              << (dcp.alignment_required ? "true" : "false") << ",\n";
            j << "    \"dominant_process\": \""
              << dcp.dominant_process_str << "\",\n";
            j << "    \"interpretation\": \""
              << json_escape(dcp.interpretation) << "\",\n";
            emit_score("terminal_deamination_score",  dcp.terminal_deamination_score,  true);
            emit_score("cpg_context_score",           dcp.cpg_context_score,           true);
            emit_score("dipyrimidine_context_score",  dcp.dipyrimidine_context_score,  true);
            emit_score("oxidative_context_score",     dcp.oxidative_context_score,     true);
            emit_score("fragmentation_context_score", dcp.fragmentation_context_score, true);
            emit_score("library_artifact_score",      dcp.library_artifact_score,      true);
            auto emit_f = [&](const char* name, float v, bool trailing) {
                j << "      \"" << name << "\": ";
                if (std::isnan(v) || std::isinf(v)) j << "null";
                else                                 j << std::setprecision(6) << v;
                j << (trailing ? ",\n" : "\n");
            };
            j << "    \"evidence\": {\n";
            emit_f("d_max_5",                  dcp.evidence.d_max_5,                  true);
            emit_f("d_max_3",                  dcp.evidence.d_max_3,                  true);
            emit_f("lambda_5",                 dcp.evidence.lambda_5,                 true);
            emit_f("lambda_3",                 dcp.evidence.lambda_3,                 true);
            emit_f("log2_cpg_ratio",           dcp.evidence.log2_cpg_ratio,           true);
            emit_f("cpg_z",                    dcp.evidence.cpg_z,                    true);
            emit_f("dipyr_contrast",           dcp.evidence.dipyr_contrast,           true);
            emit_f("ox_gt_asymmetry",          dcp.evidence.ox_gt_asymmetry,          true);
            emit_f("s_oxog_mean",              dcp.evidence.s_oxog_mean,              true);
            emit_f("s_oxog_max",               dcp.evidence.s_oxog_max,               true);
            emit_f("purine_enrichment_5prime", dcp.evidence.purine_enrichment_5prime, true);
            emit_f("hex_shift_z",              dcp.evidence.hex_shift_z,              true);
            j << "      \"adapter_clipped\": "
              << (dcp.evidence.adapter_clipped ? "true" : "false") << ",\n";
            j << "      \"adapter3_clipped\": "
              << (dcp.evidence.adapter3_clipped ? "true" : "false") << ",\n";
            j << "      \"flag_hex_artifact\": "
              << (dcp.evidence.flag_hex_artifact ? "true" : "false") << ",\n";
            j << "      \"position_0_artifact_5prime\": "
              << (dcp.evidence.position_0_artifact_5prime ? "true" : "false") << ",\n";
            j << "      \"position_0_artifact_3prime\": "
              << (dcp.evidence.position_0_artifact_3prime ? "true" : "false") << ",\n";
            j << "      \"fit_offset_5prime\": " << dcp.evidence.fit_offset_5prime << ",\n";
            j << "      \"fit_offset_3prime\": " << dcp.evidence.fit_offset_3prime << ",\n";
            j << "      \"n_reads\": " << dcp.evidence.n_reads << "\n";
            j << "    }\n";
            j << "  },\n";
        }
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
