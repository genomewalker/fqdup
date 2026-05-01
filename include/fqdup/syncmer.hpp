#pragma once
// Open-syncmer sketches for the opt-in indel-rescue path (T8).
// Canonical 64-bit hashes; deterministic, RC-equivalent.

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

namespace fqdup {

namespace syncmer_detail {

inline constexpr uint8_t kRC2bit[4] = {3, 2, 1, 0};

inline uint64_t splitmix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Canonical 2*k-bit packed k-mer hash (forward vs reverse-complement, take min).
// seq2bit[] holds bases in {0..3}; len is the read length.
inline uint64_t canonical_kmer_hash(const uint8_t* seq2bit, int pos, int k) {
    uint64_t fwd = 0, rev = 0;
    for (int i = 0; i < k; ++i) {
        const uint8_t b = seq2bit[pos + i];
        fwd = (fwd << 2) | b;
        rev = (rev << 2) | kRC2bit[seq2bit[pos + k - 1 - i]];
    }
    const uint64_t key = (fwd <= rev) ? fwd : rev;
    return splitmix64(key ^ (static_cast<uint64_t>(k) << 48));
}

}  // namespace syncmer_detail

// Sketch parameters.
struct SyncmerParams {
    int k;          // k-mer length
    int s;          // s-mer length (s < k)
    int t_offset;   // open-syncmer anchor offset (typically (k-s)/2)
    int cap;        // max sketch size after pick_topN_canonical
};

// Auto-select parameters by read length. RC-equivalence requires (k - s) even
// so the open-syncmer anchor offset t = (k-s)/2 maps to itself under window
// mirroring. We deviate from the GPT-5.5 spec's k=13/s=8 (k-s=5, asymmetric)
// to k=13/s=7 (k-s=6, symmetric) for that reason.
//   L < 75 → k=11, s=7, t=2, cap=8
//   L >= 75 → k=13, s=7, t=3, cap=16
inline SyncmerParams syncmer_params_for_length(int L) {
    if (L < 75) return SyncmerParams{11, 7, (11 - 7) / 2, 8};
    return SyncmerParams{13, 7, (13 - 7) / 2, 16};
}

// Compute open-syncmers for a 2-bit encoded sequence.
// A k-mer at position p is an open-syncmer iff the minimum s-mer hash
// inside [p, p+k-s] is located at offset t_offset (i.e. at p + t_offset).
// out_hashes / out_pos are appended (not cleared); caller may reuse buffers.
inline void compute_open_syncmers(const uint8_t* seq2bit, int L,
                                  const SyncmerParams& sp,
                                  std::vector<uint64_t>& out_hashes,
                                  std::vector<uint16_t>& out_pos) {
    const int k = sp.k, s = sp.s, t = sp.t_offset;
    if (L < k || s <= 0 || s >= k || t < 0 || t > k - s) return;

    const int n_smers_per_kmer = k - s + 1;
    // For each k-mer position p in [0, L-k], inspect the s-mers at p..p+k-s.
    // For each, compute s-mer hash on the fly (s is small, no rolling needed
    // for clarity; this is opt-in path, not the hot loop).
    for (int p = 0; p + k <= L; ++p) {
        uint64_t min_h = UINT64_MAX;
        int min_off = -1;
        for (int j = 0; j < n_smers_per_kmer; ++j) {
            uint64_t h = syncmer_detail::canonical_kmer_hash(seq2bit, p + j, s);
            if (h < min_h) { min_h = h; min_off = j; }
        }
        if (min_off == t) {
            out_hashes.push_back(syncmer_detail::canonical_kmer_hash(seq2bit, p, k));
            out_pos.push_back(static_cast<uint16_t>(p));
        }
    }
}

// Truncate a sketch to N entries deterministically: keep the N smallest
// canonical-hash values (ties broken by position). Stable across runs and
// thread counts. Operates in place on parallel hash/pos buffers.
inline void pick_topN_canonical(std::vector<uint64_t>& hashes,
                                std::vector<uint16_t>& pos,
                                int cap) {
    if (static_cast<int>(hashes.size()) <= cap) return;
    // Reuse a thread-local scratch buffer (no allocations on the hot path).
    thread_local std::vector<std::pair<uint64_t, uint16_t>> tmp;
    tmp.clear();
    tmp.reserve(hashes.size());
    for (size_t i = 0; i < hashes.size(); ++i) tmp.emplace_back(hashes[i], pos[i]);
    auto cmp = [](const std::pair<uint64_t, uint16_t>& a,
                  const std::pair<uint64_t, uint16_t>& b) {
        return a.first < b.first || (a.first == b.first && a.second < b.second);
    };
    std::nth_element(tmp.begin(), tmp.begin() + cap, tmp.end(), cmp);
    tmp.resize(cap);
    std::sort(tmp.begin(), tmp.end(), cmp);
    hashes.clear(); pos.clear();
    for (auto& kv : tmp) { hashes.push_back(kv.first); pos.push_back(kv.second); }
}

// Convenience: compute and cap in one call.
inline void compute_sketch(const uint8_t* seq2bit, int L,
                           std::vector<uint64_t>& out_hashes,
                           std::vector<uint16_t>& out_pos,
                           const SyncmerParams* override_sp = nullptr) {
    out_hashes.clear(); out_pos.clear();
    SyncmerParams sp = override_sp ? *override_sp : syncmer_params_for_length(L);
    compute_open_syncmers(seq2bit, L, sp, out_hashes, out_pos);
    pick_topN_canonical(out_hashes, out_pos, sp.cap);
}

}  // namespace fqdup
