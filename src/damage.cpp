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
#include "damage_html_assets.hpp"

#include "taph/frame_selector_decl.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/sample_damage_profile.hpp"
#include "taph/profile_json.hpp"
#include "fqdup/version.hpp"

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

// ---- paired-end work queue + worker -------------------------------------
struct PairedBatch {
    std::vector<std::string> r1;
    std::vector<std::string> r2;
};

struct PairedWorkQueue {
    std::mutex              mtx;
    std::condition_variable cv_not_empty;
    std::condition_variable cv_not_full;
    std::vector<PairedBatch> batches;
    bool   done      = false;
    int    max_depth;

    explicit PairedWorkQueue(int max_depth) : max_depth(max_depth) {}

    void push(PairedBatch&& batch) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_not_full.wait(lk, [&]{ return (int)batches.size() < max_depth; });
        batches.push_back(std::move(batch));
        cv_not_empty.notify_one();
    }
    bool pop(PairedBatch& batch) {
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

static void paired_worker_fn(PairedWorkQueue& queue, WorkerState& state) {
    PairedBatch batch;
    while (queue.pop(batch)) {
        const size_t n = std::min(batch.r1.size(), batch.r2.size());
        for (size_t k = 0; k < n; ++k) {
            const std::string& s1 = batch.r1[k];
            const std::string& s2 = batch.r2[k];
            int L1 = static_cast<int>(s1.size());
            int L2 = static_cast<int>(s2.size());
            if (L1 < LSD_L_MIN || L2 < LSD_L_MIN) { ++state.reads_skipped; continue; }
            const bool ok = taph::FrameSelector::update_sample_profile_paired(
                state.profile, s1, s2);
            if (!ok) { ++state.reads_skipped; continue; }
            // Length tracking: use R1 as the canonical insert proxy.
            ++state.lsd_hist[lsd_hist_bin(L1)];
            if (L1 < state.len_min) state.len_min = L1;
            if (L1 > state.len_max) state.len_max = L1;
            state.len_sum += L1;
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

// Pre-computed per-GC-bin LLR table: replaces exp+2×log+2×clamp inside the
// oxoG worker inner loop. Populated once after pass-1 dp is finalized.
struct LLRTable {
    float addT[15];   // log(pd/bg)          at each position
    float addC[15];   // log((1-pd)/(1-bg))  at each position
    int   max_pos;    // first p where excess < 0.01 (loop break)
    float logit_pi;
};

static void build_llr_tables(LLRTable* tables, const taph::SampleDamageProfile& dp) {
    for (int b = 0; b < N_GC; ++b) {
        auto& t = tables[b];
        const auto& gc = dp.gc_bins[b];
        float bg   = (gc.valid && gc.baseline_tc > 0.0f) ? gc.baseline_tc : dp.fit_baseline_5prime;
        float dmax = (gc.valid && gc.d_max > 0.005f)     ? gc.d_max       : dp.d_max_5prime;
        float lam  = dp.lambda_5prime;
        float pi   = (gc.valid && gc.p_damaged > 0.0f)   ? gc.p_damaged   : dp.pi_damaged;
        pi = std::clamp(pi, 0.001f, 0.999f);
        t.logit_pi = std::log(pi / (1.0f - pi));
        t.max_pos  = 0;
        for (int p = 0; p < 15; ++p) {
            float excess = dmax * std::exp(-lam * static_cast<float>(p));
            if (excess < 0.01f) { t.max_pos = p; break; }
            float pd = std::clamp(bg + (1.0f - bg) * excess, 1e-6f, 1.0f - 1e-6f);
            t.addT[p] = std::log(pd / bg);
            t.addC[p] = std::log((1.0f - pd) / (1.0f - bg));
            t.max_pos  = p + 1;
        }
    }
}

// Forward scoring: T/C at positions 0..max_pos-1 via table lookup.
static float weight_from_table(const std::string& seq, const LLRTable& t) {
    int mp = std::min(t.max_pos, static_cast<int>(seq.size()));
    float llr = 0.0f;
    for (int p = 0; p < mp; ++p) {
        char b = static_cast<char>(seq[p] & ~0x20u);
        if      (b == 'T') llr += t.addT[p];
        else if (b == 'C') llr += t.addC[p];
    }
    return 1.0f / (1.0f + std::exp(-(t.logit_pi + llr)));
}

// Reverse-complement scoring without allocating an rc string.
// RC position p = complement(seq[L-1-p]): A→T, G→C.
static float weight_rev_from_table(const std::string& seq, const LLRTable& t) {
    int L  = static_cast<int>(seq.size());
    int mp = std::min(t.max_pos, L);
    float llr = 0.0f;
    for (int p = 0; p < mp; ++p) {
        char b = static_cast<char>(seq[L - 1 - p] & ~0x20u);
        if      (b == 'A') llr += t.addT[p];   // A in original → T in RC
        else if (b == 'G') llr += t.addC[p];   // G in original → C in RC
    }
    return 1.0f / (1.0f + std::exp(-(t.logit_pi + llr)));
}

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
    // Fix E: LSD data accumulated when fuse_lsd=true in oxog_worker.
    taph::LengthBinStats        lbs;
    std::vector<LsdLlrBinAccum> llr_acc;  // per length-bin, size = n_bins
};

static void oxog_worker(WorkQueue& queue,
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
                q = weight_from_table(seq, t);
            } else {
                float w_fwd = weight_from_table(seq, t);
                float w_rev = weight_rev_from_table(seq, t);
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

// ---- helpers -----------------------------------------------------------
static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " damage (-i FILE | -1 R1.fq -2 R2.fq) [options]\n"
        << "\nInput modes (mutually exclusive):\n"
        << "  -i FILE                    Single-end / merged input FASTQ (.gz or plain)\n"
        << "  -1 FILE  -2 FILE           Paired-end raw reads (un-merged); R1 contributes\n"
        << "                             5' counters, R2 (complement-mapped) contributes\n"
        << "                             3' counters. Read 3' ends ignored (adapter zone).\n"
        << "\nOptions:\n"
        << "  -p N                       Worker threads (default: all cores)\n"
        << "  --library-type auto|ds|ss  Library type for 3'-end interpretation (default: auto)\n"
        << "  --mask-threshold FLOAT     Mask when excess P(deam) > T (default: 0.05)\n"
        << "  --tsv FILE                 Write per-position table as TSV\n"
        << "  --json FILE                Write full damage profile as JSON\n"
        << "  --html FILE                Write interactive damage report as self-contained HTML\n"
        << "  --length-bins SPEC         Length-stratified damage: auto | N | e1,e2,... (default: off)\n"
        << "  --adapter-scan-reads N     Reads sampled (single-thread) for adapter-stub detection\n"
        << "                             (default: 1000000; 0 = scan all reads)\n";
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
    std::string r1_path;
    std::string r2_path;
    int         n_threads     = static_cast<int>(std::thread::hardware_concurrency());
    if (n_threads < 1) n_threads = 1;
    double      mask_threshold = 0.05;
    std::string tsv_path;
    std::string json_path;
    std::string html_path;
    bool        run_oxog      = true;
    LengthBinOptions lb_opts = []{ bool ok; return parse_length_bins("auto", ok); }();
    int64_t     adapter_scan_reads = 1'000'000;  // 0 = scan entire file
    taph::SampleDamageProfile::LibraryType forced_lib =
        taph::SampleDamageProfile::LibraryType::UNKNOWN;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            in_path = argv[++i];
        } else if ((arg == "-1" || arg == "--read1") && i + 1 < argc) {
            r1_path = argv[++i];
        } else if ((arg == "-2" || arg == "--read2") && i + 1 < argc) {
            r2_path = argv[++i];
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
        } else if (arg == "--html" && i + 1 < argc) {
            html_path = argv[++i];
        } else if (arg == "--no-oxog") {
            run_oxog = false;
        } else if (arg == "--adapter-scan-reads" && i + 1 < argc) {
            long long v = std::stoll(argv[++i]);
            if (v < 0) {
                std::cerr << "Error: --adapter-scan-reads must be >= 0 (0 = scan all), got " << v << "\n";
                return 1;
            }
            adapter_scan_reads = static_cast<int64_t>(v);
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

    const bool paired_mode = !r1_path.empty() && !r2_path.empty();
    if (paired_mode && !in_path.empty()) {
        std::cerr << "Error: -i and -1/-2 are mutually exclusive\n";
        print_usage(argv[0]);
        return 1;
    }
    if (!paired_mode && (!r1_path.empty() || !r2_path.empty())) {
        std::cerr << "Error: -1 and -2 must both be provided for paired mode\n";
        print_usage(argv[0]);
        return 1;
    }
    if (!paired_mode && in_path.empty()) {
        std::cerr << "Error: -i FILE or -1/-2 required\n";
        print_usage(argv[0]);
        return 1;
    }
    if (mask_threshold <= 0.0 || mask_threshold >= 1.0) {
        std::cerr << "Error: --mask-threshold must be in (0, 1), got " << mask_threshold << "\n";
        return 1;
    }

    // ---- pre-scan: adapter-stub detection on first --adapter-scan-reads reads ----
    // Single-threaded; only needs hexamer_count_5prime / n_hexamers_*.
    // 0 = unlimited (scan entire file). Default 1M ≈ a few seconds on aDNA inputs.
    const int64_t pre_scan_reads = adapter_scan_reads;

    taph::AdapterStubs stubs;

    std::array<uint32_t, 4096> hex3_terminal{};
    uint64_t n_hex3 = 0;

    // SE only: reader opened here, kept alive through full pass (single-pass).
    // se_scan_buf holds records buffered during hexamer analysis.
    std::unique_ptr<FastqReaderBase> reader_se;
    std::vector<FastqRecord> se_scan_buf;

    if (!paired_mode) {
        WorkerState scan_state;
        if (pre_scan_reads > 0)
            se_scan_buf.reserve(static_cast<size_t>(pre_scan_reads));
        reader_se = make_fastq_reader(in_path, static_cast<size_t>(n_threads));
        FastqRecord rec;
        int64_t n = 0;
        while (reader_se->read(rec) && (pre_scan_reads == 0 || n < pre_scan_reads)) {
            int L = static_cast<int>(rec.seq.size());
            if (L >= LSD_L_MIN) {
                taph::FrameSelector::update_sample_profile(scan_state.profile, rec.seq);
                if (L >= 12) {
                    int code = taph::encode_hex_at(rec.seq, L - 6);
                    if (code >= 0) { ++hex3_terminal[code]; ++n_hex3; }
                }
                ++n;
            }
            se_scan_buf.push_back(std::move(rec));
        }
        // reader_se stays open — remaining records consumed in full pass below.
        stubs = taph::detect_adapter_stubs(scan_state.profile, hex3_terminal.data(), n_hex3);
    } else {
        // Paired pre-scan: 5' hexamer from R1[0..5], 3' terminal hexamer from
        // RC(R2[0..5]) (R2's 5' = molecule 3' end, complement-mapped).
        WorkerState scan_state;
        auto reader_r1 = make_fastq_reader(r1_path, 1);
        auto reader_r2 = make_fastq_reader(r2_path, 1);
        FastqRecord rec1, rec2;
        int64_t n = 0;
        auto rc_base = [](char c) -> char {
            switch (c) { case 'A': return 'T'; case 'T': return 'A';
                         case 'C': return 'G'; case 'G': return 'C';
                         case 'a': return 't'; case 't': return 'a';
                         case 'c': return 'g'; case 'g': return 'c'; }
            return 'N';
        };
        while (reader_r1->read(rec1) && reader_r2->read(rec2)
               && (pre_scan_reads == 0 || n < pre_scan_reads)) {
            int L1 = static_cast<int>(rec1.seq.size());
            int L2 = static_cast<int>(rec2.seq.size());
            if (L1 < LSD_L_MIN || L2 < LSD_L_MIN) continue;
            taph::FrameSelector::update_sample_profile_paired(
                scan_state.profile, rec1.seq, rec2.seq);
            // Molecule 3' terminal hexamer: revcomp(R2[0..5]).
            if (L2 >= 6) {
                std::string h3(6, 'N');
                for (int i = 0; i < 6; ++i) h3[i] = rc_base(rec2.seq[5 - i]);
                int code = taph::encode_hex_at(h3, 0);
                if (code >= 0) { ++hex3_terminal[code]; ++n_hex3; }
            }
            ++n;
        }
        stubs = taph::detect_adapter_stubs(scan_state.profile, hex3_terminal.data(), n_hex3);
    }

    // Per-read adapter stub fraction — fast second pass over pre-scan reads.
    // Only runs when stubs were detected; adds ~same wall time as the pre-scan.
    int64_t n_stub5_hits = 0, n_stub3_hits = 0, n_stub_reads_checked = 0;

    // Build adapter prefix code sets for prefix-conditioned F/G/H recomputation.
    // Populated after detect_adapter_stubs; consumed after finalize_sample_profile.
    std::vector<uint32_t> adapter_pfx_codes_5p, adapter_pfx_codes_3p;
    for (const auto& s : stubs.stubs5) {
        int code = taph::encode_hex_at(s, 0);
        if (code >= 0) adapter_pfx_codes_5p.push_back(static_cast<uint32_t>(code));
    }
    for (const auto& s : stubs.stubs3) {
        int code = taph::encode_hex_at(s, 0);
        if (code >= 0) adapter_pfx_codes_3p.push_back(static_cast<uint32_t>(code));
    }

    // ---- full pass: profile --------------------------------------------
    std::vector<WorkerState> states(n_threads);
    std::vector<std::thread> workers;
    workers.reserve(n_threads);

    if (!paired_mode) {
        WorkQueue queue(2 * n_threads);
        for (int t = 0; t < n_threads; ++t) {
            if (stubs.adapter_clipped || stubs.adapter3_clipped)
                workers.emplace_back(clip_worker_fn, std::ref(queue),
                                     std::ref(states[t]),
                                     std::cref(stubs.stubs5), std::cref(stubs.stubs3));
            else
                workers.emplace_back(worker_fn, std::ref(queue), std::ref(states[t]));
        }
        try {
            FastqRecord rec;
            std::vector<std::string> batch;
            batch.reserve(BATCH_SZ);
            auto enqueue = [&](std::string& seq) {
                if (!stubs.stubs5.empty() || !stubs.stubs3.empty()) {
                    int L = static_cast<int>(seq.size());
                    if (L >= 6) {
                        for (const auto& s : stubs.stubs5)
                            if (seq.compare(0, 6, s) == 0) { ++n_stub5_hits; break; }
                        for (const auto& s : stubs.stubs3)
                            if (seq.compare(L - 6, 6, s) == 0) { ++n_stub3_hits; break; }
                    }
                    ++n_stub_reads_checked;
                }
                batch.push_back(std::move(seq));
                if ((int)batch.size() == BATCH_SZ) {
                    queue.push(std::move(batch));
                    batch.clear();
                    batch.reserve(BATCH_SZ);
                }
            };
            for (auto& r : se_scan_buf) enqueue(r.seq);
            se_scan_buf.clear();
            se_scan_buf.shrink_to_fit();
            while (reader_se->read(rec)) enqueue(rec.seq);
            if (!batch.empty()) queue.push(std::move(batch));
        } catch (...) {
            queue.set_done();
            for (auto& w : workers) w.join();
            throw;
        }
        queue.set_done();
        for (auto& w : workers) w.join();
    } else {
        // Paired-end: clip workers not used (PE-mode skips read 3' ends, so
        // 3' adapter contamination doesn't enter per_pos_3prime; and PE-mode
        // pre-scan uses update_sample_profile_paired which doesn't drive the
        // 5' adapter-stub detector well enough to reliably clip on this path).
        PairedWorkQueue queue(2 * n_threads);
        for (int t = 0; t < n_threads; ++t)
            workers.emplace_back(paired_worker_fn, std::ref(queue), std::ref(states[t]));
        try {
            auto reader_r1 = make_fastq_reader(r1_path, static_cast<size_t>(n_threads));
            auto reader_r2 = make_fastq_reader(r2_path, static_cast<size_t>(n_threads));
            FastqRecord rec1, rec2;
            PairedBatch batch;
            batch.r1.reserve(BATCH_SZ);
            batch.r2.reserve(BATCH_SZ);
            while (reader_r1->read(rec1) && reader_r2->read(rec2)) {
                if (!stubs.stubs5.empty() || !stubs.stubs3.empty()) {
                    if (static_cast<int>(rec1.seq.size()) >= 6)
                        for (const auto& s : stubs.stubs5)
                            if (rec1.seq.compare(0, 6, s) == 0) { ++n_stub5_hits; break; }
                    if (static_cast<int>(rec2.seq.size()) >= 6)
                        for (const auto& s : stubs.stubs3)
                            if (rec2.seq.compare(0, 6, s) == 0) { ++n_stub3_hits; break; }
                    ++n_stub_reads_checked;
                }
                batch.r1.push_back(std::move(rec1.seq));
                batch.r2.push_back(std::move(rec2.seq));
                if ((int)batch.r1.size() == BATCH_SZ) {
                    queue.push(std::move(batch));
                    batch = PairedBatch{};
                    batch.r1.reserve(BATCH_SZ);
                    batch.r2.reserve(BATCH_SZ);
                }
            }
            if (!batch.r1.empty()) queue.push(std::move(batch));
        } catch (...) {
            queue.set_done();
            for (auto& w : workers) w.join();
            throw;
        }
        queue.set_done();
        for (auto& w : workers) w.join();
    }

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

    // Recompute F/G/H z-scores excluding reads whose first (5') or last (3')
    // hexamer matches a detected adapter stub — single-pass prefix-conditioned
    // correction without a second read pass.
    if (!adapter_pfx_codes_5p.empty() || !adapter_pfx_codes_3p.empty())
        taph::FrameSelector::recompute_fgh_excluding_adapter_prefixes(
            dp, adapter_pfx_codes_5p, adapter_pfx_codes_3p);

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
    double oxog_context_cosine = std::numeric_limits<double>::quiet_NaN();
    double oxog_gt_asymmetry_val   = std::numeric_limits<double>::quiet_NaN();
    double oxog_gt_rate_val     = std::numeric_limits<double>::quiet_NaN();

    // Oxog second pass uses single-end reads; skip in paired mode.
    if (paired_mode) run_oxog = false;

    // Fix E: build bulk_dp here (needed by both the oxoG pass and LSD fusion).
    // Fields come directly from the finalized dp; nothing between finalize and
    // this point modifies the values used below.
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

    // Fuse per-read ancient/modern classification into the oxoG pass whenever
    // damage is detectable. When --length-bins is also enabled, compute multi-bin
    // LSD edges; otherwise a single global bin accumulates the ancient fraction.
    const bool fuse_lsd = run_oxog && !paired_mode && dp.d_max_5prime > 0.01f;
    std::vector<int>   lsd_fuse_edges;
    LsdClassifyParams  lsd_cls_params{};
    double             lsd_log_prior_odds = 0.0;  // log(π/(1-π)), default π=0.5
    if (fuse_lsd) {
        lsd_cls_params = make_lsd_classify_params(bulk_dp);
        if (lb_opts.enabled())
            lsd_fuse_edges = compute_lsd_edges(merged_lsd_hist, lb_opts);
        // else lsd_fuse_edges stays empty → one global bin
        // Estimate π (endogenous fraction) from bulk d_max / d_anc.
        // Avoids using the hard-call threshold as the prior — crucial for
        // samples with tiny endogenous fractions (PPV collapses at LLR>0).
        if (lsd_cls_params.d_anc > 0.01) {
            double pi_prior = std::clamp(
                static_cast<double>(dp.d_max_5prime) / lsd_cls_params.d_anc,
                1e-6, 1.0 - 1e-6);
            lsd_log_prior_odds = std::log(pi_prior / (1.0 - pi_prior));
        }
    }

    // Populated by the oxoG pass when fuse_lsd=true; passed to
    // estimate_damage_by_length to skip its FASTQ reader entirely.
    LsdPrebuilt lsd_prebuilt;
    bool        lsd_prebuilt_ready = false;

    if (run_oxog && dp.d_max_5prime > 0.01f) {
        LLRTable llr_tables[N_GC];
        build_llr_tables(llr_tables, dp);
        WorkQueue queue2(2 * n_threads);
        const int n_lsd_bins = fuse_lsd
            ? 1 + static_cast<int>(lsd_fuse_edges.size()) : 0;
        std::vector<OxogWorkerState> ox_states(n_threads);
        if (fuse_lsd) {
            for (auto& s : ox_states) {
                s.lbs.forced_library_type = forced_lib;
                s.lbs.configure(lsd_fuse_edges);
                s.llr_acc.assign(n_lsd_bins, {});
            }
        }
        std::vector<std::thread> workers2;
        workers2.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t)
            workers2.emplace_back(oxog_worker, std::ref(queue2),
                                   std::ref(ox_states[t]), std::cref(dp),
                                   llr_tables, is_ss,
                                   fuse_lsd,
                                   fuse_lsd ? &lsd_cls_params : nullptr,
                                   fuse_lsd ? &lsd_fuse_edges : nullptr,
                                   lsd_log_prior_odds);
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
        // Unlike the deamination D field, this breaks the DS Chargaff cancellation because
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

        // Fix E: merge per-thread LSD state into LsdPrebuilt for injection below.
        if (fuse_lsd) {
            taph::LengthBinStats lsd_master;
            lsd_master.forced_library_type = forced_lib;
            lsd_master.configure(lsd_fuse_edges);
            for (auto& os : ox_states) lsd_master.merge(os.lbs);

            const int n_lsd_bins = 1 + static_cast<int>(lsd_fuse_edges.size());
            std::vector<LsdLlrBinAccum> merged_llr(n_lsd_bins);
            for (auto& os : ox_states) {
                for (int b = 0; b < n_lsd_bins; ++b) {
                    const auto& w = os.llr_acc[b];
                    auto& m = merged_llr[b];
                    m.n_damaged   += w.n_damaged;
                    m.n_undamaged += w.n_undamaged;
                    for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                        m.t_5_anc[p]          += w.t_5_anc[p];
                        m.tc_5_anc[p]         += w.tc_5_anc[p];
                        m.h_3_anc[p]          += w.h_3_anc[p];
                        m.n_3_anc[p]          += w.n_3_anc[p];
                        m.t_5_anc_cpg[p]      += w.t_5_anc_cpg[p];
                        m.tc_5_anc_cpg[p]     += w.tc_5_anc_cpg[p];
                        m.t_5_anc_noncpg[p]   += w.t_5_anc_noncpg[p];
                        m.tc_5_anc_noncpg[p]  += w.tc_5_anc_noncpg[p];
                        m.t_5_anc_g[p]        += w.t_5_anc_g[p];
                        m.tg_5_anc[p]         += w.tg_5_anc[p];
                        m.t_5_mod_g[p]        += w.t_5_mod_g[p];
                        m.tg_5_mod[p]         += w.tg_5_mod[p];
                        m.t_5_mod[p]          += w.t_5_mod[p];
                        m.tc_5_mod[p]         += w.tc_5_mod[p];
                        m.h_3_mod[p]          += w.h_3_mod[p];
                        m.n_3_mod[p]          += w.n_3_mod[p];
                        m.a_5_anc_all[p]      += w.a_5_anc_all[p];
                        m.c_5_anc_all[p]      += w.c_5_anc_all[p];
                        m.g_5_anc_all[p]      += w.g_5_anc_all[p];
                        m.t_5_anc_all[p]      += w.t_5_anc_all[p];
                    }
                    // Soft-EM fields (scalars, not per-position arrays).
                    m.sw_t5_anc  += w.sw_t5_anc;
                    m.sw_tc5_anc += w.sw_tc5_anc;
                    m.sw_h3_anc  += w.sw_h3_anc;
                    m.sw_n3_anc  += w.sw_n3_anc;
                    m.sw_sum     += w.sw_sum;
                }
            }

            lsd_prebuilt.edges        = lsd_fuse_edges;
            lsd_prebuilt.merged_stats = std::move(lsd_master);
            lsd_prebuilt.llr_bins     = std::move(merged_llr);
            lsd_prebuilt_ready        = true;

            // Compute ancient-fraction d_max using posterior-weighted (soft-EM)
            // accumulation. sw_* fields hold pos-0 counts weighted by
            // P(ancient|read,π) where π is estimated from bulk d_max/d_anc,
            // avoiding the selection-on-outcome bias of hard LLR>0 calls.
            // Hard-call counts (n_damaged/n_undamaged) are preserved for
            // informational use (n_reads in JSON).
            {
                int64_t hard_n_damaged = 0, hard_n_tot = 0;
                double sw_sum = 0.0;
                for (const auto& bin : lsd_prebuilt.llr_bins) {
                    hard_n_damaged += bin.n_damaged;
                    hard_n_tot     += bin.n_damaged + bin.n_undamaged;
                    sw_sum += bin.sw_sum;
                }
                if (hard_n_tot >= 10000 && sw_sum >= 10.0) {
                    dp.damaged_fraction_pi    = static_cast<float>(sw_sum / hard_n_tot);
                    dp.damaged_fraction_n     = hard_n_damaged;
                    dp.damaged_fraction_valid = true;
                    // bulk_d_max = π × d_anc  →  d_anc = bulk_d_max / π
                    // More robust than the sw_t5/sw_tc5 - bg5 approach, which
                    // collapses to 0 when posteriors are all ≈ prior (low π).
                    const double pi_est = dp.damaged_fraction_pi;
                    if (pi_est > 0.0f) {
                        dp.damaged_fraction_d5 = static_cast<float>(
                            std::clamp(static_cast<double>(dp.d_max_5prime) / pi_est, 0.0, 1.0));
                        dp.damaged_fraction_d3 = static_cast<float>(
                            std::clamp(static_cast<double>(dp.d_max_3prime) / pi_est, 0.0, 1.0));
                    }
                    // Modern d_max: hard-call counts from global.t_5_mod/tc_5_mod
                    // are still computed via lsd_accumulate; aggregate them here.
                    const double bg5 = lsd_cls_params.bg_5;
                    const double bg3 = lsd_cls_params.bg_3;
                    int64_t mod_t5 = 0, mod_tc5 = 0, mod_h3 = 0, mod_n3 = 0;
                    for (const auto& bin : lsd_prebuilt.llr_bins) {
                        for (int p = 0; p < LengthBinDamageProfile::N_POS; ++p) {
                            mod_t5  += bin.t_5_mod[p];
                            mod_tc5 += bin.tc_5_mod[p];
                            mod_h3  += bin.h_3_mod[p];
                            mod_n3  += bin.n_3_mod[p];
                            if (p == 0) break;  // only pos 0 needed for d_max peak
                        }
                    }
                    if (mod_tc5 >= 50)
                        dp.modern_fraction_d5 = static_cast<float>(
                            std::max(0.0, static_cast<double>(mod_t5) / mod_tc5 - bg5));
                    if (mod_n3 >= 50)
                        dp.modern_fraction_d3 = static_cast<float>(
                            std::max(0.0, static_cast<double>(mod_h3) / mod_n3 - bg3));
                }
            }
        }
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
    if (!dp.library_bic_winner_model.empty()) {
        std::cout << "  winner=" << dp.library_bic_winner_model
                  << "  second=" << dp.library_bic_second_model
                  << "  margin=" << std::fixed << std::setprecision(1) << dp.library_bic_margin;
        if (dp.library_artifact_contaminated && !dp.library_artifact_reasons.empty()) {
            std::cout << "  [artifact-contaminated:";
            for (size_t i = 0; i < dp.library_artifact_reasons.size(); ++i) {
                std::cout << (i ? "," : " ") << dp.library_artifact_reasons[i];
            }
            std::cout << "]";
        }
        std::cout << "\n";
    }
    if (dp.library_type_evaluable) {
        std::cout << "  post p_ds=" << std::fixed << std::setprecision(3) << dp.library_p_ds
                  << "  p_ss=" << dp.library_p_ss
                  << "  p_bias=" << dp.library_p_bias
                  << "  p_winner=" << dp.library_p_winner
                  << (dp.library_p_winner < taph::SampleDamageProfile::kLibraryTypeConfidenceThreshold
                          ? "  [low confidence]" : "")
                  << "\n";
    }
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
    if (n_stub_reads_checked > 0 && (!stubs.stubs5.empty() || !stubs.stubs3.empty())) {
        std::cout << "  adapter stubs:";
        if (!stubs.stubs5.empty())
            std::cout << " 5'=" << stubs.stubs5[0]
                      << " (" << std::fixed << std::setprecision(1)
                      << (static_cast<double>(n_stub5_hits) / n_stub_reads_checked * 100.0)
                      << "% of reads)";
        if (!stubs.stubs3.empty())
            std::cout << " 3'=" << stubs.stubs3[0]
                      << " (" << std::fixed << std::setprecision(1)
                      << (static_cast<double>(n_stub3_hits) / n_stub_reads_checked * 100.0)
                      << "% of reads)";
        std::cout << "\n";
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
    // bulk_dp was built before the oxoG pass; reuse it here.
    LengthStratifiedDamageProfile lsd;
    if (lb_opts.enabled() && paired_mode) {
        std::cerr << "Warning: --length-bins is not supported in paired mode; skipping.\n";
    } else if (lb_opts.enabled()) {
        lsd = estimate_damage_by_length(in_path, forced_lib, lb_opts,
                                        &merged_lsd_hist,
                                        static_cast<size_t>(n_threads), 0, &bulk_dp,
                                        lsd_prebuilt_ready ? &lsd_prebuilt : nullptr);
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
    if (!html_path.empty() && json_path.empty()) {
        auto otr = taph::compute_oxog_trinuc(dp);
        oxog_context_cosine = otr.cosine;
        oxog_gt_asymmetry_val   = otr.gt_asymmetry;
        oxog_gt_rate_val     = otr.gt_rate;
    }

    if (!json_path.empty()) {
        std::ofstream j(json_path);
        if (!j) {
            std::cerr << "Error: cannot write JSON: " << json_path << "\n";
            return 1;
        }
        taph::ProfileJsonInput pji;
        pji.sample_name          = in_path;
        pji.version              = std::string("fqdup/") + FQDUP_VERSION;
        pji.n_reads              = static_cast<uint64_t>(reads_scanned);
        pji.adapter_stubs_5prime = stubs.stubs5;
        pji.adapter_stubs_3prime = stubs.stubs3;
        if (n_stub_reads_checked > 0) {
            pji.adapter_stub5_read_fraction = static_cast<double>(n_stub5_hits) / n_stub_reads_checked;
            pji.adapter_stub3_read_fraction = static_cast<double>(n_stub3_hits) / n_stub_reads_checked;
            pji.adapter_stub_reads_checked  = n_stub_reads_checked;
        }
        pji.top_hex_enriched        = stubs.top_enriched;
        pji.top_hex_enriched_3prime = stubs.top_enriched_3prime;
        {
            auto hea = taph::compute_hex_end_asymmetry(
                dp, stubs.top_enriched, stubs.top_enriched_3prime);
            // Certify d_max_3: confounded when ends carry independent 6-mer families
            // (rc_overlap_topk==0) and 3' enriched hexamers are mostly non-damage-consistent.
            // When confounded on an SS library, override d_max_combined to 5prime_only so
            // the combined estimate reflects the reliable end rather than the artifact-inflated one.
            const auto& t3 = stubs.top_enriched_3prime;
            int k3 = std::min(5, (int)t3.size());
            double dmg_frac_3 = std::numeric_limits<double>::quiet_NaN();
            if (k3 > 0) {
                int m = 0;
                for (int i = 0; i < k3; ++i) if (t3[i].damage_consistent) ++m;
                dmg_frac_3 = static_cast<double>(m) / k3;
            }
            bool d3_confounded = (hea.rc_overlap_topk == 0)
                              && (!std::isnan(dmg_frac_3) && dmg_frac_3 < 0.5)
                              && (dp.fit_offset_3prime >= 1);
            if (d3_confounded && is_ss
                    && dp.d_max_source != taph::SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY) {
                dp.d_max_combined = dp.d_max_5prime;
                dp.d_max_source   = taph::SampleDamageProfile::DmaxSource::FIVE_PRIME_ONLY;
            }
            pji.hex_end_asymmetry = std::move(hea);
        }
        pji.adapter_clipped      = stubs.adapter_clipped;
        pji.adapter3_clipped     = stubs.adapter3_clipped;
        pji.flag_hex_artifact    = stubs.flag_hex_artifact;
        pji.short_read_frac      = short_read_frac;
        pji.s_oxog               = S_oxog;
        pji.se_s_oxog            = SE_oxog;
        pji.d_oriented           = D_oriented;
        pji.has_oxog_score       = has_oxog_score;
        pji.lsd                  = lsd.bins.empty() ? nullptr : &lsd;
        taph::profile_to_json(dp, j, pji);
        std::cout << "JSON written: " << json_path << "\n";
    }

    if (!html_path.empty()) {
        std::ofstream h(html_path);
        if (!h) {
            std::cerr << "Error: cannot write HTML: " << html_path << "\n";
            return 1;
        }

        // Compute summaries (may be re-computed if JSON was also written; cheap)
        auto pres_h  = taph::compute_preservation_summary(dp, is_ss,
            stubs.adapter_clipped, stubs.flag_hex_artifact,
            cpg_score_z, oxog_score_z, oxog_context_cosine, hex_shift_p);
        auto flags_h = taph::compute_library_qc_flags(dp, is_ss,
            stubs.flag_hex_artifact, hex_stats.jsd,
            hex_stats.entropy_terminal, short_read_frac);
        auto dcp_h   = taph::compute_damage_context_profile(
            dp, cpg_score_z, hex_shift_z,
            stubs.adapter_clipped, stubs.adapter3_clipped,
            stubs.flag_hex_artifact);

        // Adapter name lookup (mirrors JSON serialisation)
        static const std::pair<const char*, const char*> kAdapters[] = {
            {"ACACTC","TruSeq/P5"},{"AATGAT","TruSeq/Universal"},{"GATCGG","TruSeq/i7"},
            {"CTGTCT","Nextera/Tn5"},{"AGATCG","TruSeq/R1"},{"TGGAAT","TruSeq/R2"},{"GCGAAT","TruSeq/R2alt"},
        };
        std::string adapter_top_seq = stubs.stubs5.empty() ? "" : stubs.stubs5[0];
        std::string adapter_name_str = "unknown";
        if (!adapter_top_seq.empty())
            for (const auto& kv : kAdapters)
                if (adapter_top_seq == kv.first) { adapter_name_str = kv.second; break; }

        // Derive sample name from input path basename, strip extension
        auto sample_name = [&]() -> std::string {
            std::string base = in_path;
            auto sl = base.rfind('/');
            if (sl != std::string::npos) base = base.substr(sl + 1);
            for (const char* ext : {".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"}) {
                if (base.size() > strlen(ext) &&
                    base.compare(base.size() - strlen(ext), strlen(ext), ext) == 0)
                    base = base.substr(0, base.size() - strlen(ext));
            }
            return base;
        }();

        auto jv = [&](double v) -> std::string {
            if (std::isnan(v) || std::isinf(v)) return "null";
            std::ostringstream os;
            os << std::fixed << std::setprecision(6) << v;
            return os.str();
        };
        auto jb = [](bool v) { return v ? "true" : "false"; };

        h << fqdup_dmghtml::HTML_PRE;
        h << "const D = {\n";
        h << "  \"sample\": \"" << sample_name << "\",\n";
        h << "  \"n_reads\": " << reads_scanned << ",\n";
        h << "  \"library_type\": \"" << dp.library_type_str() << "\",\n";
        h << "  \"d5\": " << jv(dp.d_max_5prime) << ",\n";
        h << "  \"d3\": " << jv(dp.d_max_3prime) << ",\n";
        h << "  \"lam5\": " << jv(dp.lambda_5prime) << ",\n";
        h << "  \"lam3\": " << jv(dp.lambda_3prime) << ",\n";
        h << "  \"bg5\": " << jv(dp.fit_baseline_5prime) << ",\n";
        h << "  \"bg3\": " << jv(dp.fit_baseline_3prime) << ",\n";
        h << "  \"s_gt\": " << jv(dp.s_gt) << ",\n";
        h << "  \"oxog_gt_asymmetry\": " << jv(oxog_gt_asymmetry_val) << ",\n";
        h << "  \"oxog_gt_rate\": " << jv(oxog_gt_rate_val) << ",\n";
        h << "  \"oxog_cosine\": " << jv(oxog_context_cosine) << ",\n";

        // Per-position arrays (1-based positions)
        h << "  \"pos\": [";
        for (int p = 0; p < N_POS; ++p) { if (p) h << ","; h << (p + 1); }
        h << "],\n";

        h << "  \"ct5\": [";
        for (int p = 0; p < N_POS; ++p) {
            if (p) h << ",";
            double v = (dp.tc_total_5prime[p] >= MIN_COV) ? dp.t_freq_5prime[p] : -1.0;
            h << jv(v);
        }
        h << "],\n";

        if (!is_ss) {
            h << "  \"ga5\": [";
            for (int p = 0; p < N_POS; ++p) {
                if (p) h << ",";
                double v = (dp.ag_total_3prime[p] >= MIN_COV) ? dp.a_freq_3prime[p] : -1.0;
                h << jv(v);
            }
            h << "],\n";
        } else {
            h << "  \"ct3\": [";
            for (int p = 0; p < N_POS; ++p) {
                if (p) h << ",";
                double tc3 = dp.tc_total_3prime[p];
                double v = (tc3 >= MIN_COV) ? dp.t_freq_3prime[p] / tc3 : -1.0;
                h << jv(v);
            }
            h << "],\n";
        }

        h << "  \"gt5\": [";
        for (int p = 0; p < N_POS; ++p) {
            if (p) h << ",";
            double denom = dp.t_from_g_5prime[p] + dp.g_count_5prime[p];
            double v = (denom >= MIN_COV) ? dp.t_from_g_5prime[p] / denom : -1.0;
            h << jv(v);
        }
        h << "],\n";

        // Channels array (A=deamination, B=stop-codon, C=8-oxoG top, D=Chargaff, F=complement, G=hydantoin, H=AT-stop)
        {
            bool ch_f_det = dp.channel_f_valid && dp.channel_f_z > kOxChannelZDetect;
            bool ch_g_det = dp.channel_g_valid && dp.channel_g_z > kOxChannelZDetect;
            bool ch_h_det = dp.channel_h_valid && (dp.channel_h_z > kOxChannelZDetect || dp.channel_h_z_p2plus > kOxChannelZDetect);
            h << "  \"channels\": [\n";
            // A: deamination
            // Use presence (not strict validation) — validated=false fires when pos0 artifact
            // forces offset correction but damage is still clearly real.
            bool ch_a_det = !dp.damage_artifact &&
                            (dp.d_max_5prime > 0.01 || dp.d_max_3prime > 0.01);
            h << "    {\"id\":\"A\",\"name\":\"deamination\","
              << "\"detected\":" << jb(ch_a_det) << ","
              << "\"applicable\":true,"
              << "\"d5\":" << jv(dp.d_max_5prime) << "},\n";
            // B: stop-codon validator
            h << "    {\"id\":\"B\",\"name\":\"stop_codon_conversion\","
              << "\"detected\":" << jb(dp.channel_b_valid && !dp.channel_b_inverted) << ","
              << "\"applicable\":true},\n";
            // C: 8-oxoG top strand
            h << "    {\"id\":\"C\",\"name\":\"8_oxog_top_strand\","
              << "\"detected\":" << jb(dp.ox_damage_detected) << ","
              << "\"applicable\":true},\n";
            // D: Chargaff asymmetry
            h << "    {\"id\":\"D\",\"name\":\"chargaff_gt_asymmetry\","
              << "\"detected\":" << jb(std::abs(dp.ox_gt_asymmetry) > 0.01f) << ","
              << "\"applicable\":true,"
              << "\"d5\":" << jv(dp.s_gt) << "},\n";
            // F: complement 8-oxoG
            h << "    {\"id\":\"F\",\"name\":\"8_oxog_complement\","
              << "\"detected\":" << jb(ch_f_det) << ","
              << "\"applicable\":" << jb(!is_ss) << ","
              << "\"z_score\":" << jv(dp.channel_f_z) << "},\n";
            // G: hydantoin
            h << "    {\"id\":\"G\",\"name\":\"hydantoin_oxidation\","
              << "\"detected\":" << jb(ch_g_det) << ","
              << "\"applicable\":" << jb(!is_ss) << ","
              << "\"z_score\":" << jv(dp.channel_g_z) << "},\n";
            // H: AT-stop (depurination)
            h << "    {\"id\":\"H\",\"name\":\"depurination_at_stop\","
              << "\"detected\":" << jb(ch_h_det) << ","
              << "\"applicable\":true,"
              << "\"z_score\":" << jv(dp.channel_h_z) << "}\n";
            h << "  ],\n";
        }

        // Length bins
        if (!lsd.bins.empty()) {
            h << "  \"length_bins\": [\n";
            for (size_t b = 0; b < lsd.bins.size(); ++b) {
                const auto& lb = lsd.bins[b];
                if (b) h << ",\n";
                h << "    {\"lo\":" << lb.length_lo << ",\"hi\":" << lb.length_hi
                  << ",\"d5\":" << jv(lb.d_max_5prime) << ",\"d3\":" << jv(lb.d_max_3prime)
                  << ",\"lam5\":" << jv(lb.lambda_5prime) << ",\"lam3\":" << jv(lb.lambda_3prime)
                  << ",\"bg5\":" << jv(lb.bg_5prime) << ",\"bg3\":" << jv(lb.bg_3prime)
                  << ",\"n\":" << lb.n_reads << "}";
            }
            h << "\n  ]\n";
        } else {
            h << "  \"length_bins\": []\n";
        }

        // ── Preservation ─────────────────────────────────────────────────────
        h << "  ,\"anc_valid\": "      << jb(dp.damaged_fraction_valid)  << "\n";
        h << "  ,\"anc_d5\": "        << jv(dp.damaged_fraction_d5)     << "\n";
        h << "  ,\"anc_d3\": "        << jv(dp.damaged_fraction_d3)     << "\n";
        h << "  ,\"anc_pi\": "        << jv(dp.damaged_fraction_pi)     << "\n";
        h << "  ,\"anc_n\": "         << dp.damaged_fraction_n          << "\n";
        h << "  ,\"mod_d5\": "        << jv(dp.modern_fraction_d5)      << "\n";
        h << "  ,\"mod_d3\": "        << jv(dp.modern_fraction_d3)      << "\n";
        h << "  ,\"pres_score\": "    << jv(dp.preservation_score)    << "\n";
        h << "  ,\"pres_label\": \""  << pres_h.label                 << "\"\n";
        h << "  ,\"auth_eff\": "      << jv(pres_h.authenticity_eff)  << "\n";
        h << "  ,\"ox_eff\": "        << jv(pres_h.oxidation_eff)     << "\n";
        h << "  ,\"qcr_eff\": "       << jv(pres_h.qc_risk_eff)       << "\n";

        // ── Library typing ────────────────────────────────────────────────────
        h << "  ,\"lib_p_ds\": "      << jv(dp.library_p_ds_final)        << "\n";
        h << "  ,\"lib_p_ss\": "      << jv(dp.library_p_ss_final)        << "\n";
        h << "  ,\"lib_bic_model\": \"" << dp.library_bic_winner_model    << "\"\n";
        h << "  ,\"lib_bic_margin\": " << jv(dp.library_bic_margin)       << "\n";
        h << "  ,\"lib_artifact\": "   << jb(dp.library_artifact_contaminated) << "\n";

        // ── Damage context ────────────────────────────────────────────────────
        h << "  ,\"dom_process\": \""  << dcp_h.dominant_process_str << "\"\n";
        h << "  ,\"interpretation\": \"" << [&](std::string s) -> std::string {
            // minimal JSON string escaping for the interpretation sentence
            std::string out; out.reserve(s.size() + 8);
            for (char c : s) { if (c == '"') out += "\\\""; else if (c == '\\') out += "\\\\"; else out += c; }
            return out;
        }(dcp_h.interpretation) << "\"\n";
        h << "  ,\"dam_score\": "      << jv(dcp_h.terminal_deamination_score)  << "\n";
        h << "  ,\"ox_score\": "       << jv(dcp_h.oxidative_context_score)     << "\n";
        h << "  ,\"frag_score\": "     << jv(dcp_h.fragmentation_context_score) << "\n";
        h << "  ,\"art_score\": "      << jv(dcp_h.library_artifact_score)      << "\n";
        h << "  ,\"cpg_contrast\": "   << jv(dp.cpg_contrast)                   << "\n";
        h << "  ,\"cpg_z\": "          << jv(cpg_score_z)                       << "\n";

        // ── Depurination ──────────────────────────────────────────────────────
        h << "  ,\"depur_detected\": " << jb(dp.depurination_detected)     << "\n";
        h << "  ,\"depur_enrich5\": "  << jv(dp.purine_enrichment_5prime)  << "\n";
        h << "  ,\"depur_enrich3\": "  << jv(dp.purine_enrichment_3prime)  << "\n";
        h << "  ,\"depur_z\": "        << jv(depur_score_z)                << "\n";

        // ── Library QC ────────────────────────────────────────────────────────
        h << "  ,\"pos0_art5\": "      << jb(dp.position_0_artifact_5prime) << "\n";
        h << "  ,\"pos0_art3\": "      << jb(dp.position_0_artifact_3prime) << "\n";
        h << "  ,\"hex_ent5\": "       << jv(hex_stats.entropy_terminal)    << "\n";
        h << "  ,\"hex_ent_int\": "    << jv(hex_stats.entropy_interior)    << "\n";
        h << "  ,\"hex_jsd\": "        << jv(hex_stats.jsd)                 << "\n";
        h << "  ,\"hex_shift_z\": "    << jv(hex_shift_z)                   << "\n";
        h << "  ,\"short_frac\": "     << jv(short_read_frac < 0 ? 0.0 : short_read_frac) << "\n";
        h << "  ,\"stubs5\": [";
        for (size_t i = 0; i < stubs.stubs5.size(); ++i) { if (i) h << ","; h << "\"" << stubs.stubs5[i] << "\""; }
        h << "]\n";
        h << "  ,\"stubs3\": [";
        for (size_t i = 0; i < stubs.stubs3.size(); ++i) { if (i) h << ","; h << "\"" << stubs.stubs3[i] << "\""; }
        h << "]\n";
        h << "  ,\"adapter_seq\": "    << (adapter_top_seq.empty() ? "null" : "\"" + adapter_top_seq + "\"") << "\n";
        h << "  ,\"adapter_name\": "   << (adapter_top_seq.empty() ? "null" : "\"" + adapter_name_str + "\"") << "\n";
        h << "  ,\"top_hex\": [";
        {
            int n_out = 0;
            for (const auto& hr : stubs.top_enriched) {
                if (n_out >= 5) break;
                auto seq = taph::decode_hex(hr.idx);
                if (n_out > 0) h << ",";
                h << "{\"seq\":\"" << seq.data() << "\","
                  << "\"log2fc\":" << jv(hr.log2fc) << ","
                  << "\"dc\":" << jb(hr.damage_consistent) << "}";
                ++n_out;
            }
        }
        h << "]\n";
        h << "  ,\"qc_flags\": [";
        {
            bool first = true;
            auto ef = [&](const char* n) { if (!first) h << ","; h << "\"" << n << "\""; first = false; };
            if (flags_h.adapter_remnant_5prime)   ef("adapter_remnant_5prime");
            if (flags_h.adapter_remnant_3prime)   ef("adapter_remnant_3prime");
            if (flags_h.hexamer_composition_bias) ef("hexamer_composition_bias");
            if (flags_h.hexamer_terminal_shift)   ef("hexamer_terminal_shift");
            if (flags_h.short_read_spike)         ef("short_read_spike");
            if (flags_h.depurination)             ef("depurination");
            if (flags_h.ds_3prime_signal_absent)  ef("ds_3prime_signal_absent");
            if (flags_h.ga3_inward_displaced)     ef("ga3_inward_displaced");
            if (flags_h.hexamer_artifact_bias)    ef("hexamer_artifact_bias");
        }
        h << "]\n";

        h << "};\n";
        h << fqdup_dmghtml::HTML_POST;
        std::cout << "HTML written: " << html_path << "\n";
    }

    return 0;
}
