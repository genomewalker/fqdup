#pragma once
// Private worker queue and thread helpers for fqdup damage subcommand.
#include "fqdup/fastq_common.hpp"
#include "fqdup/damage_profile.hpp"
#include "taph/frame_selector_decl.hpp"
#include "taph/sample_damage_profile.hpp"
#include <array>
#include <climits>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>
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
        }
    }
}
