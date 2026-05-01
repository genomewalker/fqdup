#pragma once
// T8.2 SyncmerIndex — compact hash→parent_id postings table for the opt-in
// indel-rescue path. Built once per length-window pass, queried by every
// child sketch in the window, then dropped.
//
// Layout:
//   spans[] is sorted by hash, deduped to one entry per unique hash. Each
//   span owns a contiguous slice of parent_ids[] via (begin, len). Hot
//   spans (whose postings list would have ≥ hash_hot parents) are dropped
//   inline at build time and counted in overflow_hashes.
//
// Build is single-pass (compute_sketch per eligible parent → push pairs →
// sort → group). All work uses 64-bit hashes and 32-bit parent_ids.
//
// Eligibility for indexing:
//   arena.is_eligible(id) AND L ∈ [min_len, max_len]
//   AND 2 ≤ bundle_occ < bundle_hot
// Singleton bundles (occ=1) cannot be rescued under within-bundle policy
// → skipped entirely. Hotspots (occ ≥ bundle_hot) are also skipped.

#include "fqdup/syncmer.hpp"
#include "arena.hpp"

#include "flat_hash_map.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace fqdup::derep_detail {

struct SyncmerIndex {
    struct Span {
        uint64_t hash;
        uint32_t begin;   // offset into parent_ids
        uint16_t len;     // # parents under this hash (after hot filter)
    };
    std::vector<Span>     spans;
    std::vector<uint32_t> parent_ids;

    uint64_t overflow_hashes = 0;
    uint64_t indexed_parents = 0;

    size_t bytes_used() const {
        return spans.size() * sizeof(Span) + parent_ids.size() * sizeof(uint32_t);
    }

    struct Hit { uint32_t parent_id; uint16_t shared_hits; };

    // Caller-owned scratch — one per worker thread, reused across all child
    // queries. Eliminates flat_hash_map allocation/rehash on the hot path.
    struct QueryScratch {
        ska::flat_hash_map<uint32_t, uint16_t> tally;
        std::vector<Hit> out;
        std::vector<Hit> hits_out;  // returned to caller; reused buffer
    };

    // Returns reference to scratch.hits_out, populated with (parent_id,
    // shared_hits) where shared_hits ≥ min_hits, capped to topk by
    // (shared_hits desc, parent_id asc). Caller MUST consume the result
    // before the next query() call on the same scratch.
    const std::vector<Hit>& query(QueryScratch& scratch,
                                  const uint64_t* child_hashes, int n_hashes,
                                  uint16_t min_hits, uint32_t topk) const {
        scratch.tally.clear();
        if (scratch.tally.bucket_count() < static_cast<size_t>(n_hashes) * 4)
            scratch.tally.reserve(static_cast<size_t>(n_hashes) * 4);

        for (int i = 0; i < n_hashes; ++i) {
            const uint64_t h = child_hashes[i];
            auto lo = std::lower_bound(
                spans.begin(), spans.end(), h,
                [](const Span& s, uint64_t v) { return s.hash < v; });
            if (lo == spans.end() || lo->hash != h) continue;
            const uint32_t* p = parent_ids.data() + lo->begin;
            const uint16_t  n = lo->len;
            for (uint16_t k = 0; k < n; ++k) {
                auto& slot = scratch.tally[p[k]];
                if (slot < UINT16_MAX) ++slot;
            }
        }

        scratch.out.clear();
        scratch.out.reserve(scratch.tally.size());
        for (auto& kv : scratch.tally) {
            if (kv.second >= min_hits)
                scratch.out.push_back(Hit{kv.first, kv.second});
        }
        // Use partial-sort when the survivor set exceeds topk; otherwise full sort.
        auto cmp = [](const Hit& a, const Hit& b) {
            if (a.shared_hits != b.shared_hits)
                return a.shared_hits > b.shared_hits;
            return a.parent_id < b.parent_id;
        };
        if (scratch.out.size() > topk) {
            std::partial_sort(scratch.out.begin(),
                              scratch.out.begin() + topk,
                              scratch.out.end(), cmp);
            scratch.out.resize(topk);
        } else {
            std::sort(scratch.out.begin(), scratch.out.end(), cmp);
        }
        scratch.hits_out = scratch.out;  // copy small bounded buffer
        return scratch.hits_out;
    }

    // T8 Step 6 — filtered query.
    //
    // Same as query() but the caller-supplied predicate `accept(parent_id)` is
    // applied INSIDE the per-hash posting loop, before tally accumulation.
    // Top-k is computed over filtered survivors, so length-mismatched or
    // id_count-rejected parents can never evict valid same-length parents
    // from the topk slots. Fixes a silent recall loss flagged by GPT-5.5
    // (post-topk filtering let wrong-length parents starve real ones).
    template <typename FilterFn>
    const std::vector<Hit>& query_filtered(QueryScratch& scratch,
                                           const uint64_t* child_hashes, int n_hashes,
                                           uint16_t min_hits, uint32_t topk,
                                           FilterFn&& accept) const {
        scratch.tally.clear();
        if (scratch.tally.bucket_count() < static_cast<size_t>(n_hashes) * 4)
            scratch.tally.reserve(static_cast<size_t>(n_hashes) * 4);

        for (int i = 0; i < n_hashes; ++i) {
            const uint64_t h = child_hashes[i];
            auto lo = std::lower_bound(
                spans.begin(), spans.end(), h,
                [](const Span& s, uint64_t v) { return s.hash < v; });
            if (lo == spans.end() || lo->hash != h) continue;
            const uint32_t* p = parent_ids.data() + lo->begin;
            const uint16_t  n = lo->len;
            for (uint16_t k = 0; k < n; ++k) {
                if (!accept(p[k])) continue;
                auto& slot = scratch.tally[p[k]];
                if (slot < UINT16_MAX) ++slot;
            }
        }

        scratch.out.clear();
        scratch.out.reserve(scratch.tally.size());
        for (auto& kv : scratch.tally) {
            if (kv.second >= min_hits)
                scratch.out.push_back(Hit{kv.first, kv.second});
        }
        auto cmp = [](const Hit& a, const Hit& b) {
            if (a.shared_hits != b.shared_hits)
                return a.shared_hits > b.shared_hits;
            return a.parent_id < b.parent_id;
        };
        if (scratch.out.size() > topk) {
            std::partial_sort(scratch.out.begin(),
                              scratch.out.begin() + topk,
                              scratch.out.end(), cmp);
            scratch.out.resize(topk);
        } else {
            std::sort(scratch.out.begin(), scratch.out.end(), cmp);
        }
        scratch.hits_out = scratch.out;
        return scratch.hits_out;
    }

    // T8 Step 7 — reset for reuse as a thread_local scratch index.
    void clear() {
        spans.clear();
        parent_ids.clear();
        overflow_hashes = 0;
        indexed_parents = 0;
    }
};

// Build a SyncmerIndex over eligible parents in length window [min_len, max_len].
// `bundle_occ_of` is a callable: uint32_t id → uint32_t occupancy.
template <typename BundleOccFn>
inline SyncmerIndex build_syncmer_index_impl(
    const SeqArena& arena,
    BundleOccFn&& bundle_occ_of,
    int min_len,
    int max_len,
    uint32_t bundle_hot,
    uint32_t hash_hot)
{
    SyncmerIndex idx;
    if (max_len < min_len) return idx;

    struct Pair { uint64_t hash; uint32_t pid; };
    std::vector<Pair> pairs;
    pairs.reserve(arena.size());

    std::vector<uint8_t> dec;
    std::vector<uint64_t> sketch_h;
    std::vector<uint16_t> sketch_pos;

    for (uint32_t id = 0; id < arena.size(); ++id) {
        if (!arena.is_eligible(id)) continue;
        int L = arena.length(id);
        if (L < min_len || L > max_len) continue;
        uint32_t occ = bundle_occ_of(id);
        if (occ < 2u || occ >= bundle_hot) continue;

        if ((int)dec.size() < L) dec.resize(L);
        arena.decode_range(id, 0, L, dec.data());
        fqdup::compute_sketch(dec.data(), L, sketch_h, sketch_pos);
        if (sketch_h.empty()) continue;
        ++idx.indexed_parents;
        for (uint64_t h : sketch_h) pairs.push_back(Pair{h, id});
    }

    if (pairs.empty()) return idx;

    std::sort(pairs.begin(), pairs.end(),
              [](const Pair& a, const Pair& b) {
                  if (a.hash != b.hash) return a.hash < b.hash;
                  return a.pid < b.pid;
              });

    // Group runs of equal hash into spans, dropping hot lists inline.
    idx.spans.reserve(pairs.size() / 4 + 16);
    idx.parent_ids.reserve(pairs.size());
    size_t i = 0;
    while (i < pairs.size()) {
        size_t j = i + 1;
        while (j < pairs.size() && pairs[j].hash == pairs[i].hash) ++j;
        size_t run = j - i;
        if (run >= hash_hot) {
            ++idx.overflow_hashes;
            i = j; continue;
        }
        // Dedup parent_ids within the run (a parent may emit the same hash
        // twice if its sketch contains repeats; collapse to one posting).
        uint32_t begin = static_cast<uint32_t>(idx.parent_ids.size());
        uint32_t last_pid = UINT32_MAX;
        uint16_t added = 0;
        for (size_t k = i; k < j; ++k) {
            if (pairs[k].pid == last_pid) continue;
            idx.parent_ids.push_back(pairs[k].pid);
            last_pid = pairs[k].pid;
            ++added;
        }
        idx.spans.push_back(SyncmerIndex::Span{pairs[i].hash, begin, added});
        i = j;
    }

    return idx;
}

// Same builder, but takes a pre-filtered list of parent ids. Skips the full
// arena scan + per-id eligibility/length/occ gating — caller has already
// done that. Used by T8 driver after pre-bucketing parents by length.
inline SyncmerIndex build_syncmer_index_from_pids(
    const SeqArena& arena,
    const uint32_t* pids, size_t n_pids,
    uint32_t hash_hot)
{
    SyncmerIndex idx;
    if (n_pids == 0) return idx;

    struct Pair { uint64_t hash; uint32_t pid; };
    std::vector<Pair> pairs;
    pairs.reserve(n_pids * 8);

    std::vector<uint8_t> dec;
    std::vector<uint64_t> sketch_h;
    std::vector<uint16_t> sketch_pos;

    for (size_t i = 0; i < n_pids; ++i) {
        uint32_t id = pids[i];
        int L = arena.length(id);
        if ((int)dec.size() < L) dec.resize(L);
        arena.decode_range(id, 0, L, dec.data());
        fqdup::compute_sketch(dec.data(), L, sketch_h, sketch_pos);
        if (sketch_h.empty()) continue;
        ++idx.indexed_parents;
        for (uint64_t h : sketch_h) pairs.push_back(Pair{h, id});
    }

    if (pairs.empty()) return idx;

    std::sort(pairs.begin(), pairs.end(),
              [](const Pair& a, const Pair& b) {
                  if (a.hash != b.hash) return a.hash < b.hash;
                  return a.pid < b.pid;
              });

    idx.spans.reserve(pairs.size() / 4 + 16);
    idx.parent_ids.reserve(pairs.size());
    size_t i = 0;
    while (i < pairs.size()) {
        size_t j = i + 1;
        while (j < pairs.size() && pairs[j].hash == pairs[i].hash) ++j;
        size_t run = j - i;
        if (run >= hash_hot) {
            ++idx.overflow_hashes;
            i = j; continue;
        }
        uint32_t begin = static_cast<uint32_t>(idx.parent_ids.size());
        uint32_t last_pid = UINT32_MAX;
        uint16_t added = 0;
        for (size_t k = i; k < j; ++k) {
            if (pairs[k].pid == last_pid) continue;
            idx.parent_ids.push_back(pairs[k].pid);
            last_pid = pairs[k].pid;
            ++added;
        }
        idx.spans.push_back(SyncmerIndex::Span{pairs[i].hash, begin, added});
        i = j;
    }
    return idx;
}

// T8 Step 7 — build into a caller-owned (typically thread_local) index from
// already-decoded sequences. Avoids re-decoding from the arena (the T8 worker
// has already decoded all bundle members into a flat buffer for direct-path
// banded-ed). Reuses pairs_scratch across calls to keep allocations to zero
// after warmup.
//
// `decoded[off[i] .. off[i]+L[i])` holds the i-th parent's bases (2-bit
// packed-decoded). `pids[i]` is the corresponding global parent id.
inline void build_syncmer_index_from_decoded(
    SyncmerIndex& idx,
    const uint8_t* decoded, const uint32_t* off, const int* L,
    const uint32_t* pids, size_t n_pids,
    std::vector<std::pair<uint64_t, uint32_t>>& pairs_scratch,
    uint32_t hash_hot)
{
    idx.clear();
    pairs_scratch.clear();
    if (n_pids == 0) return;

    thread_local std::vector<uint64_t> sketch_h;
    thread_local std::vector<uint16_t> sketch_pos;

    for (size_t i = 0; i < n_pids; ++i) {
        fqdup::compute_sketch(decoded + off[i], L[i], sketch_h, sketch_pos);
        if (sketch_h.empty()) continue;
        ++idx.indexed_parents;
        for (uint64_t h : sketch_h)
            pairs_scratch.emplace_back(h, pids[i]);
    }

    if (pairs_scratch.empty()) return;

    std::sort(pairs_scratch.begin(), pairs_scratch.end(),
              [](const std::pair<uint64_t, uint32_t>& a,
                 const std::pair<uint64_t, uint32_t>& b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });

    idx.spans.reserve(pairs_scratch.size() / 4 + 16);
    idx.parent_ids.reserve(pairs_scratch.size());
    size_t i = 0;
    while (i < pairs_scratch.size()) {
        size_t j = i + 1;
        while (j < pairs_scratch.size() && pairs_scratch[j].first == pairs_scratch[i].first) ++j;
        size_t run = j - i;
        if (run >= hash_hot) {
            ++idx.overflow_hashes;
            i = j; continue;
        }
        uint32_t begin = static_cast<uint32_t>(idx.parent_ids.size());
        uint32_t last_pid = UINT32_MAX;
        uint16_t added = 0;
        for (size_t k = i; k < j; ++k) {
            if (pairs_scratch[k].second == last_pid) continue;
            idx.parent_ids.push_back(pairs_scratch[k].second);
            last_pid = pairs_scratch[k].second;
            ++added;
        }
        idx.spans.push_back(SyncmerIndex::Span{pairs_scratch[i].first, begin, added});
        i = j;
    }
}

}  // namespace fqdup::derep_detail
