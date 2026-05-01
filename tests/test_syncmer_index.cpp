// T8.2 — SyncmerIndex build + query unit test.
// Build via add_executable in CMakeLists.txt; runs after build.

#include "derep_detail/syncmer_index.hpp"
#include "derep_detail/arena.hpp"
#include "fqdup/syncmer.hpp"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using fqdup::derep_detail::SeqArena;
using fqdup::derep_detail::SyncmerIndex;
using fqdup::derep_detail::build_syncmer_index_impl;

static std::string mutate_one_indel(const std::string& s, std::mt19937& rng) {
    std::uniform_int_distribution<int> pos(1, static_cast<int>(s.size()) - 2);
    std::uniform_int_distribution<int> base(0, 3);
    static const char A[] = {'A', 'C', 'G', 'T'};
    std::string out = s;
    int p = pos(rng);
    if (rng() & 1u) out.insert(out.begin() + p, A[base(rng)]);
    else            out.erase(out.begin() + p);
    return out;
}

static std::string random_seq(int L, std::mt19937& rng) {
    static const char A[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> base(0, 3);
    std::string s; s.reserve(L);
    for (int i = 0; i < L; ++i) s.push_back(A[base(rng)]);
    return s;
}

static std::string revcomp(const std::string& s) {
    std::string out; out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        switch (*it) {
            case 'A': out.push_back('T'); break;
            case 'T': out.push_back('A'); break;
            case 'C': out.push_back('G'); break;
            case 'G': out.push_back('C'); break;
            default:  out.push_back('N');
        }
    }
    return out;
}

int main() {
    std::mt19937 rng(0xC0FFEEu);

    // Test 1: RC-equivalence — sketch of a sequence equals sketch of its
    // reverse complement (canonical hashing).
    {
        std::string s = random_seq(100, rng);
        std::string r = revcomp(s);
        std::vector<uint8_t> a(s.size()), b(r.size());
        for (size_t i = 0; i < s.size(); ++i) {
            a[i] = (s[i] == 'A' ? 0 : s[i] == 'C' ? 1 : s[i] == 'G' ? 2 : 3);
            b[i] = (r[i] == 'A' ? 0 : r[i] == 'C' ? 1 : r[i] == 'G' ? 2 : 3);
        }
        std::vector<uint64_t> ha, hb;
        std::vector<uint16_t> pa, pb;
        fqdup::compute_sketch(a.data(), (int)a.size(), ha, pa);
        fqdup::compute_sketch(b.data(), (int)b.size(), hb, pb);
        assert(!ha.empty());
        // hashes are canonical; sketches should be set-equal.
        std::sort(ha.begin(), ha.end());
        std::sort(hb.begin(), hb.end());
        assert(ha == hb);
        std::printf("[ok] RC-equivalence: %zu hashes\n", ha.size());
    }

    // Test 2: a 1-indel pair shares ≥ τ syncmers; index round-trip finds it.
    {
        SeqArena arena;
        const int L = 100;
        std::string parent = random_seq(L, rng);
        std::string child  = mutate_one_indel(parent, rng);
        // Pad/extra unrelated parents to make the index non-trivial.
        std::vector<std::string> seqs;
        seqs.push_back(parent);
        for (int k = 0; k < 20; ++k) seqs.push_back(random_seq(L, rng));
        seqs.push_back(child);
        for (auto& s : seqs) arena.append(s);
        // Bundle occ: pretend everyone is in the same bundle of size > 2.
        auto bundle_occ_of = [&](uint32_t) -> uint32_t { return (uint32_t)seqs.size(); };

        SyncmerIndex idx = build_syncmer_index_impl(
            arena, bundle_occ_of, L - 1, L + 1, /*bundle_hot=*/1000, /*hash_hot=*/4096);
        assert(idx.indexed_parents > 0);

        // Query with the child sketch.
        std::vector<uint8_t> dec(arena.length(arena.size() - 1));
        arena.decode_range((uint32_t)arena.size() - 1, 0, (int)dec.size(), dec.data());
        std::vector<uint64_t> ch; std::vector<uint16_t> cp;
        fqdup::compute_sketch(dec.data(), (int)dec.size(), ch, cp);

        SyncmerIndex::QueryScratch qs;
        auto hits = idx.query(qs, ch.data(), (int)ch.size(), /*min_hits=*/2, /*topk=*/8);
        bool found_parent = false;
        uint16_t parent_shared = 0;
        for (auto& h : hits) {
            if (h.parent_id == 0) { found_parent = true; parent_shared = h.shared_hits; }
        }
        assert(found_parent);
        std::printf("[ok] 1-indel pair: parent shares %u syncmers (out of %zu child)\n",
                    parent_shared, ch.size());
    }

    // Test 3: hot-span overflow counter trips when hash_hot is small.
    {
        SeqArena arena;
        const int L = 100;
        // Give every read the same prefix to force shared syncmers → hot spans.
        std::string base = random_seq(L, rng);
        for (int k = 0; k < 30; ++k) {
            std::string s = base;
            // Mutate the last 5 bp only — shared region dominates the sketch.
            for (int j = 0; j < 5; ++j) {
                int p = L - 1 - j;
                s[p] = "ACGT"[(rng()) & 3u];
            }
            arena.append(s);
        }
        auto bundle_occ_of = [&](uint32_t) -> uint32_t { return 30u; };
        SyncmerIndex idx_lo = build_syncmer_index_impl(
            arena, bundle_occ_of, L, L, 1000, /*hash_hot=*/4);
        SyncmerIndex idx_hi = build_syncmer_index_impl(
            arena, bundle_occ_of, L, L, 1000, /*hash_hot=*/1024);
        assert(idx_lo.overflow_hashes > 0);
        assert(idx_hi.overflow_hashes == 0);
        std::printf("[ok] hot-span overflow: lo=%llu hi=%llu\n",
                    (unsigned long long)idx_lo.overflow_hashes,
                    (unsigned long long)idx_hi.overflow_hashes);
    }

    // Test 4 (T8 Step 6) — query_filtered preserves valid hits even when
    // wrong-length parents would otherwise consume top-k slots.
    //
    // Setup: 4 parents share a common syncmer-rich core. Two are flagged
    // "good" (filter accepts) and have distinctive variation; two are "bad"
    // (filter rejects) but emit just as many shared syncmers. With min_hits
    // small and topk=2, post-filter top-k could easily evict the goods if
    // the filter ran AFTER tally. The filtered overload runs filter BEFORE
    // tally, so only goods compete for the slots.
    {
        SeqArena arena;
        const int L = 100;
        std::string core = random_seq(L, rng);
        // 4 near-identical parents, each with one mismatch in different spots,
        // all sharing most syncmers with `core`.
        std::vector<std::string> seqs;
        for (int k = 0; k < 4; ++k) {
            std::string s = core;
            int p = 5 + k * 20;
            s[p] = "ACGT"[(rng()) & 3u];
            seqs.push_back(s);
        }
        for (auto& s : seqs) arena.append(s);
        auto bundle_occ_of = [&](uint32_t) -> uint32_t { return 4u; };
        SyncmerIndex idx = build_syncmer_index_impl(
            arena, bundle_occ_of, L, L, /*bundle_hot=*/1000, /*hash_hot=*/1024);
        assert(idx.indexed_parents == 4);

        // Query with `core` itself — every parent shares many syncmers.
        std::vector<uint8_t> dec(L);
        for (int i = 0; i < L; ++i)
            dec[i] = (core[i] == 'A' ? 0 : core[i] == 'C' ? 1 : core[i] == 'G' ? 2 : 3);
        std::vector<uint64_t> ch; std::vector<uint16_t> cp;
        fqdup::compute_sketch(dec.data(), L, ch, cp);

        SyncmerIndex::QueryScratch qs;
        // Filter: accept only parent ids 0 and 2 ("goods"). 1 and 3 = "bads".
        auto accept = [](uint32_t pid) { return pid == 0u || pid == 2u; };
        auto& hits = idx.query_filtered(qs, ch.data(), (int)ch.size(),
                                        /*min_hits=*/1, /*topk=*/2, accept);
        // top-k=2, min_hits=1 → both goods must survive; no bads.
        assert(hits.size() == 2);
        bool got_zero = false, got_two = false;
        for (auto& h : hits) {
            assert(h.parent_id != 1u && h.parent_id != 3u);
            if (h.parent_id == 0u) got_zero = true;
            if (h.parent_id == 2u) got_two  = true;
        }
        assert(got_zero && got_two);
        std::printf("[ok] query_filtered preserves valid hits under topk pressure\n");
    }

    std::printf("test_syncmer_index: all tests passed\n");
    return 0;
}
