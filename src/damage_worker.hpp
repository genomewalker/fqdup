#pragma once
// Private worker queue and thread helpers for fqdup damage subcommand.
#include "fqdup/fastq_common.hpp"
#include "fqdup/damage_profile.hpp"
#include "taph/frame_selector_decl.hpp"
#include "taph/sample_damage_profile.hpp"
#include "taph/length_stratified_profile.hpp"
#include <array>
#include <climits>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

// LSD sig5 encoding:
//   5' ternary (T=1,C=2,other=0) — 5 positions = 3^5 = 243 values
//   3' quinary (A=1,G=2,T=3,C=4,other=0) — 5 positions = 5^5 = 3125 values
// ceiling: CpG context not encoded (needs next-base); upgrade via quinary 5' sig
static constexpr int N_LSD_SIG5_5P   = 243;   // 3^5 ternary
static constexpr int N_LSD_SIG5_3R   = 3125;  // 5^5 quinary
static constexpr int N_LSD_SIG5      = N_LSD_SIG5_5P * N_LSD_SIG5_3R;  // 759375
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

// ---- Single-pass OxoG + LSD accumulation structs -----------------------
// Ternary 10-position signature histogram cell. Second moments allow
// exact SE computation without re-reading the FASTQ.
struct OxSigCell {
    int64_t t=0, g=0, cv=0, av=0, n=0;  // first moments + count
    int64_t t2=0, tm=0, m2=0;            // Σ T^2, Σ T*M, Σ M^2 for SE
};

static constexpr int N_SIG10 = 59049;  // 3^10: ternary 10-position space


// ---- per-thread state --------------------------------------------------
struct WorkerState {
    taph::SampleDamageProfile           profile;
    std::array<uint64_t, LSD_HIST_BINS> lsd_hist{};
    int64_t reads_scanned = 0;
    int64_t reads_skipped = 0;
    int     len_min       = INT_MAX;
    int     len_max       = 0;
    int64_t len_sum       = 0;
    // OxoG ternary sig hists indexed [gc_bin * N_SIG10 + sig]
    std::vector<OxSigCell> ox_fwd;
    std::vector<OxSigCell> ox_rev;
    // LSD single-pass: per-(prov_bin, lsd_sig5) read counts.
    // Populated only when lsd_cnt is non-empty (set by damage.cpp before pass 1).
    taph::LengthBinStats    lbs;
    std::vector<int32_t>    lsd_cnt;        // size = n_prov_bins × N_LSD_SIG5
    std::vector<int>        lsd_prov_edges; // provisional bin boundaries (copy)
    WorkerState()
        : ox_fwd(taph::SampleDamageProfile::N_GC_BINS * N_SIG10)
        , ox_rev(taph::SampleDamageProfile::N_GC_BINS * N_SIG10) {}
};

// Accumulate one read into the OxoG sig histograms and LSD reservoir.
// Called after update_sample_profile inside worker loops.
inline static void worker_ox_accumulate(WorkerState& state, const std::string& seq, int L) {
    const int gc  = taph::SampleDamageProfile::get_gc_bin(seq);
    const int beg = L / 3, end_ = L - (L / 3);
    int T = 0, G = 0, Cv = 0, Av = 0;
    for (int j = beg; j < end_; ++j) {
        const char b = seq[j] & ~0x20u;
        if      (b == 'T') ++T;
        else if (b == 'G') ++G;
        else if (b == 'C') ++Cv;
        else if (b == 'A') ++Av;
    }
    // Compute both ternary sigs in one endpoint pass.
    // ALL reads go into BOTH ox_fwd and ox_rev — no proxy orientation decision.
    // is_ss is unknown during pass 1; oxog_from_sig_hist applies only merged_fwd
    // for SS (all reads, correct fwd counts) and both for DS (q-weighting selects
    // the dominant orientation per cell naturally post-fit).
    uint32_t fsig = 0, rsig = 0, fb = 1, rb = 1;
    const int k = (L < 10) ? L : 10;
    for (int i = 0; i < k; ++i) {
        const char fc = seq[i]     & ~0x20u;
        const char rc = seq[L-1-i] & ~0x20u;
        const int  fv = (fc == 'T') ? 1 : (fc == 'C') ? 2 : 0;
        const int  rv = (rc == 'A') ? 1 : (rc == 'G') ? 2 : 0;
        fsig += fv * fb; fb *= 3;
        rsig += rv * rb; rb *= 3;
    }
    // Fwd histogram: original orientation (correct for SS and DS fwd reads)
    {
        auto& cell = state.ox_fwd[gc * N_SIG10 + fsig];
        const int M = T + G;
        cell.t += T; cell.g += G; cell.cv += Cv; cell.av += Av; ++cell.n;
        cell.t2 += (int64_t)T * T; cell.tm += (int64_t)T * M; cell.m2 += (int64_t)M * M;
    }
    // Rev histogram: complement-mapped orientation (correct for DS rev reads)
    {
        const int Tr = Av, Gr = Cv, Cvr = G, Avr = T;
        auto& cell = state.ox_rev[gc * N_SIG10 + rsig];
        const int M = Tr + Gr;
        cell.t += Tr; cell.g += Gr; cell.cv += Cvr; cell.av += Avr; ++cell.n;
        cell.t2 += (int64_t)Tr*Tr; cell.tm += (int64_t)Tr*M; cell.m2 += (int64_t)M*M;
    }
    // LSD sig5: joint (5' ternary, 3' quinary) for single-pass reconstruction.
    if (!state.lsd_cnt.empty()) {
        // sig5f: first 5 trits of fsig (T=1, C=2, other=0)
        const uint32_t sig5f = fsig % N_LSD_SIG5_5P;
        // sig5r_4: 3' end with 5-way encoding (A=1, G=2, T=3, C=4, other=0)
        uint32_t sig5r_4 = 0, q5 = 1;
        const int k5 = (L < 5) ? L : 5;
        for (int i = 0; i < k5; ++i) {
            const char rc4 = seq[L-1-i] & ~0x20u;
            const int rv4 = (rc4=='A') ? 1 : (rc4=='G') ? 2 : (rc4=='T') ? 3 : (rc4=='C') ? 4 : 0;
            sig5r_4 += rv4 * q5; q5 *= 5;
        }
        const uint32_t lsd_sig = sig5f * N_LSD_SIG5_3R + sig5r_4;
        int cb = 0;
        const int Lcap = (L < LSD_L_MAX) ? L : LSD_L_MAX;
        for (int e : state.lsd_prov_edges) { if (Lcap >= e) ++cb; else break; }
        ++state.lsd_cnt[cb * N_LSD_SIG5 + lsd_sig];
    }
}

inline static void worker_fn(WorkQueue& queue, WorkerState& state) {
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
            if (!state.lsd_cnt.empty())
                state.lbs.update(seq, (L < LSD_L_MAX) ? L : LSD_L_MAX);
            worker_ox_accumulate(state, seq, L);
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

inline static void paired_worker_fn(PairedWorkQueue& queue, WorkerState& state) {
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
inline static void clip_worker_fn(WorkQueue& queue, WorkerState& state,
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
            if (!state.lsd_cnt.empty())
                state.lbs.update(seq, (L < LSD_L_MAX) ? L : LSD_L_MAX);
            worker_ox_accumulate(state, seq, L);
        }
    }
}
