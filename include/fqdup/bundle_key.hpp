#pragma once

#include "fqdup/cluster_format.hpp"

#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace fqdup::bundlekey {

inline std::uint64_t fnv1a(const char* p, std::size_t n) {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (std::size_t i = 0; i < n; ++i) {
        h ^= static_cast<std::uint8_t>(p[i]);
        h *= 0x100000001b3ULL;
    }
    return h;
}

// Default endpoint k for bundle-occupancy prior. k=12 caps collisions on
// capture hotspots; clamped to L/2 so short reads (≥30 bp) still split 5'/3'.
constexpr int kDefaultEndK = 12;

inline std::uint64_t from_decoded(const char* seq, std::uint32_t len, int end_k) {
    if (len == 0) return 0;
    int k = std::min<int>(end_k, static_cast<int>(len) / 2);
    if (k <= 0) k = static_cast<int>(len);
    std::uint64_t h5 = fnv1a(seq, k);
    std::uint64_t h3 = fnv1a(seq + len - k, k);
    return h5 ^ (h3 * 0x9E3779B97F4A7C15ULL) ^ (static_cast<std::uint64_t>(len) << 1);
}

// Length-agnostic variant — same endpoint hash, no length term. Used by the
// adjacent-length (L±1) probe so reads differing by one indel can co-locate
// in the same bundle. To stay collision-resistant against the length-aware
// key when both maps are kept side by side, we tag the high bit.
inline std::uint64_t from_decoded_no_len(const char* seq, std::uint32_t len, int end_k) {
    if (len == 0) return 0;
    int k = std::min<int>(end_k, static_cast<int>(len) / 2);
    if (k <= 0) k = static_cast<int>(len);
    std::uint64_t h5 = fnv1a(seq, k);
    std::uint64_t h3 = fnv1a(seq + len - k, k);
    return (h5 ^ (h3 * 0x9E3779B97F4A7C15ULL)) | 0x8000000000000000ULL;
}

inline std::uint64_t from_cluster(const fqdup::clusterfmt::ClusterRecord& r, int end_k) {
    if (r.parent_seq_len == 0) return 0;
    std::vector<char> buf(r.parent_seq_len);
    fqdup::clusterfmt::decode_2bit(r.parent_seq.data(), r.parent_seq_len, buf.data());
    return from_decoded(buf.data(), r.parent_seq_len, end_k);
}

} // namespace fqdup::bundlekey
