#pragma once
// 2-bit base encoding tables and small mixing helpers shared across derep.cpp.
// Internal to derep.cpp — not a public API.

#include <array>
#include <cstdint>
#include <utility>

namespace fqdup::derep_detail {

inline constexpr uint8_t kDec2bit[4] = {'A', 'C', 'G', 'T'};

inline const std::array<uint8_t, 256> kEnc2bit = [] {
    std::array<uint8_t, 256> lut{};
    lut.fill(0xFFu);
    lut[static_cast<uint8_t>('A')] = lut[static_cast<uint8_t>('a')] = 0u;
    lut[static_cast<uint8_t>('C')] = lut[static_cast<uint8_t>('c')] = 1u;
    lut[static_cast<uint8_t>('G')] = lut[static_cast<uint8_t>('g')] = 2u;
    lut[static_cast<uint8_t>('T')] = lut[static_cast<uint8_t>('t')] = 3u;
    return lut;
}();

inline const std::array<std::array<uint8_t, 4>, 256> kDec4 = [] {
    std::array<std::array<uint8_t, 4>, 256> lut{};
    for (int v = 0; v < 256; ++v) {
        lut[v][0] = kDec2bit[(v >> 6) & 0x3u];
        lut[v][1] = kDec2bit[(v >> 4) & 0x3u];
        lut[v][2] = kDec2bit[(v >> 2) & 0x3u];
        lut[v][3] = kDec2bit[(v >> 0) & 0x3u];
    }
    return lut;
}();

inline void split4_lens(int ilen, int& s0, int& s1, int& s2, int& s3) {
    s0 = ilen / 4;
    s1 = ilen * 2 / 4 - s0;
    s2 = ilen * 3 / 4 - s0 - s1;
    s3 = ilen - s0 - s1 - s2;
}

inline uint64_t splitmix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

inline uint64_t pair_key(uint64_t ha, uint64_t hb, int tag, int ilen) {
    if (ha > hb) std::swap(ha, hb);
    return splitmix64(ha ^ rotl64(hb, 23)
                      ^ (static_cast<uint64_t>(tag) << 56)
                      ^ static_cast<uint64_t>(ilen));
}

}  // namespace fqdup::derep_detail
