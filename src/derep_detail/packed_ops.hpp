#pragma once
// 2-bit packed sequence comparison primitives used by Phase 3 PCR error
// correction. Internal to derep.cpp.

#include <cstdint>
#include <cstring>
#include <limits>

namespace fqdup::derep_detail {

// Count Hamming mismatches between packed interiors [k5, k5+ilen).
// Returns the count, stopping at max_mm+1 (early exit if exceeded).
inline int packed_hamming(const uint8_t* pa, const uint8_t* pb,
                          int k5, int ilen, int max_mm) {
    int mismatches = 0;
    for (int i = 0; i < ilen && mismatches <= max_mm; ) {
        int pos    = k5 + i;
        int byte_i = pos >> 2;
        int bit_i  = pos & 3;
        if (bit_i == 0 && i + 4 <= ilen) {
            uint8_t diff = pa[byte_i] ^ pb[byte_i];
            if (diff) mismatches += __builtin_popcount((diff | (diff >> 1)) & 0x55u);
            i += 4;
        } else {
            uint8_t sh = 6 - 2 * bit_i;
            if (((pa[byte_i] >> sh) & 0x3u) != ((pb[byte_i] >> sh) & 0x3u)) ++mismatches;
            ++i;
        }
    }
    return mismatches;
}

// Returns true if the substitution a→b (or b→a) is damage-consistent
// using 2-bit encoded values. C↔T (xr=2) and G↔A/C↔A/G↔T (xr=1) are the
// defining aDNA damage signatures; A↔T / C↔G (xr=3) are transversions.
inline bool is_damage_sub_packed(uint8_t a, uint8_t b, bool /*protect_deamination*/) {
    uint8_t xr = a ^ b;
    return (xr == 1u || xr == 2u);
}

// True if the mismatch (pa_2b, pb_2b) at interior position ipos (0-based, length ilen)
// is damage-compatible: C→T within `margin` bases of the 5' end (DS+SS), or G→A within
// `margin` bases of the 3' end (DS only; both channels at both ends for SS).
// Uses the same terminal margin as kDamageEdgeMargin / kDamageEdgeMargin_T58 (5 bp).
inline bool is_terminal_damage_mismatch(uint8_t pa_2b, uint8_t pb_2b,
                                        int ipos, int ilen,
                                        bool ss_mode,
                                        int margin = 5) {
    uint8_t xr = pa_2b ^ pb_2b;
    if (xr != 2u) return false;             // not C↔T or G↔A
    bool is_ct = ((pa_2b & pb_2b) == 1u);   // C(1) & T(3) = 1
    bool is_ga = ((pa_2b & pb_2b) == 0u);   // G(2) & A(0) = 0
    bool at_5  = (ipos < margin);
    bool at_3  = (ipos >= ilen - margin);
    if (ss_mode) return (at_5 || at_3) && (is_ct || is_ga);
    return (at_5 && is_ct) || (at_3 && is_ga);
}

// Calibration-free binomial LRT for interior-transition absorption.
// Returns log(P(x|n,f0) / P(x|n,f1)) = x*log(f0/f1) + (n-x)*log((1-f0)/(1-f1)).
// Positive → H0 (PCR error) more likely; absorb if result > log(T).
// Typical: f0=0.10 (late-cycle PCR error fraction), f1=0.50 (real variant).
inline double interior_transition_lrt(uint64_t n_parent, uint64_t n_child,
                                       double f0, double f1) {
    uint64_t n = n_parent + n_child;
    if (n == 0) return 0.0;
    double x  = static_cast<double>(n_child);
    double nx = static_cast<double>(n - n_child);
    return x * std::log(f0 / f1) + nx * std::log((1.0 - f0) / (1.0 - f1));
}

// Wilson score lower 95% confidence bound for proportion x/n.
// Replaces raw alt_freq >= threshold SNP veto — avoids over-protection at
// low bundle depth (e.g. a 1/2 observation has point freq 50% but LB ≈ 1%).
inline double wilson_lower95(uint64_t x, uint64_t n) {
    if (n == 0) return 0.0;
    constexpr double z  = 1.96;
    constexpr double z2 = z * z;   // 3.8416
    double p  = static_cast<double>(x) / static_cast<double>(n);
    double nd = static_cast<double>(n);
    double num = p + z2 / (2.0 * nd)
                   - z * std::sqrt(p * (1.0 - p) / nd + z2 / (4.0 * nd * nd));
    return std::max(0.0, num / (1.0 + z2 / nd));
}

// Returns true if child should be absorbed by parent: exactly 1 interior
// mismatch AND the substitution is NOT damage-consistent.
inline bool packed_should_absorb(const uint8_t* pa, const uint8_t* pb,
                                 int k5, int ilen, bool protect_deamination) {
    int mismatches = 0;
    uint8_t mm_a = 0, mm_b = 0;

    for (int i = 0; i < ilen; ) {
        int pos   = k5 + i;
        int byte_i = pos >> 2;
        int bit_i  = pos & 3;

        if (bit_i == 0 && i + 4 <= ilen) {
            uint8_t diff = pa[byte_i] ^ pb[byte_i];
            if (diff) {
                uint8_t any = (diff | (diff >> 1)) & 0x55u;
                int byte_mismatches = __builtin_popcount(any);
                mismatches += byte_mismatches;
                if (mismatches > 1) return false;
                if (byte_mismatches == 1) {
                    uint8_t sh = static_cast<uint8_t>(__builtin_ctz(any));
                    mm_a = (pa[byte_i] >> sh) & 0x3u;
                    mm_b = (pb[byte_i] >> sh) & 0x3u;
                }
            }
            i += 4;
        } else {
            uint8_t sh = 6 - 2 * bit_i;
            uint8_t ba = (pa[byte_i] >> sh) & 0x3u;
            uint8_t bb = (pb[byte_i] >> sh) & 0x3u;
            if (ba != bb) {
                if (++mismatches > 1) return false;
                mm_a = ba; mm_b = bb;
            }
            ++i;
        }
    }
    return mismatches == 1 && !is_damage_sub_packed(mm_a, mm_b, protect_deamination);
}

struct MismatchInfo {
    bool     found;
    uint16_t position;  // interior-relative position (0-based) of the single mismatch
    uint8_t  base_a;    // parent's 2-bit base
    uint8_t  base_b;    // child's 2-bit base
};

// Like packed_should_absorb but also returns position and bases of the single mismatch.
inline MismatchInfo packed_find_mismatch(const uint8_t* pa, const uint8_t* pb,
                                         int k5, int ilen, bool protect_deamination) {
    int mismatches = 0;
    uint8_t mm_a = 0, mm_b = 0;
    uint16_t mm_pos = 0;

    for (int i = 0; i < ilen; ) {
        int pos    = k5 + i;
        int byte_i = pos >> 2;
        int bit_i  = pos & 3;

        if (bit_i == 0 && i + 4 <= ilen) {
            uint8_t diff = pa[byte_i] ^ pb[byte_i];
            if (diff) {
                uint8_t any = (diff | (diff >> 1)) & 0x55u;
                int nm = __builtin_popcount(any);
                mismatches += nm;
                if (mismatches > 1) return {false, 0, 0, 0};
                if (nm == 1) {
                    uint8_t sh = static_cast<uint8_t>(__builtin_ctz(any));
                    mm_a = (pa[byte_i] >> sh) & 0x3u;
                    mm_b = (pb[byte_i] >> sh) & 0x3u;
                    mm_pos = static_cast<uint16_t>(i + (3 - (sh >> 1)));
                }
            }
            i += 4;
        } else {
            uint8_t sh = static_cast<uint8_t>(6 - 2 * bit_i);
            uint8_t ba = (pa[byte_i] >> sh) & 0x3u;
            uint8_t bb = (pb[byte_i] >> sh) & 0x3u;
            if (ba != bb) {
                if (++mismatches > 1) return {false, 0, 0, 0};
                mm_a = ba; mm_b = bb;
                mm_pos = static_cast<uint16_t>(i);
            }
            ++i;
        }
    }
    if (mismatches != 1) return {false, 0, 0, 0};
    if (is_damage_sub_packed(mm_a, mm_b, protect_deamination)) return {false, 0, 0, 0};
    return {true, mm_pos, mm_a, mm_b};
}

// Returns up to 2 mismatches between packed interiors [k5, k5+ilen).
// No damage filtering — caller is responsible. count==3 means rejected.
struct MismatchInfo2 {
    int     count;
    uint16_t pos[2];
    uint8_t  base_a[2];
    uint8_t  base_b[2];
};

inline MismatchInfo2 packed_find_mismatches2(const uint8_t* pa, const uint8_t* pb,
                                             int k5, int ilen) {
    MismatchInfo2 r{};
    for (int i = 0; i < ilen; ) {
        int pos    = k5 + i;
        int byte_i = pos >> 2;
        int bit_i  = pos & 3;
        if (bit_i == 0 && i + 4 <= ilen) {
            uint8_t diff = pa[byte_i] ^ pb[byte_i];
            if (diff) {
                int nm = __builtin_popcount((diff | (diff >> 1)) & 0x55u);
                if (r.count + nm > 2) return {3, {}, {}, {}};
                for (int lane = 0; lane < 4; ++lane) {
                    int sh = 6 - 2 * lane;
                    uint8_t ba = (pa[byte_i] >> sh) & 0x3u;
                    uint8_t bb = (pb[byte_i] >> sh) & 0x3u;
                    if (ba != bb) {
                        r.pos[r.count]    = static_cast<uint16_t>(i + lane);
                        r.base_a[r.count] = ba;
                        r.base_b[r.count] = bb;
                        r.count++;
                    }
                }
            }
            i += 4;
        } else {
            uint8_t sh = static_cast<uint8_t>(6 - 2 * bit_i);
            uint8_t ba = (pa[byte_i] >> sh) & 0x3u;
            uint8_t bb = (pb[byte_i] >> sh) & 0x3u;
            if (ba != bb) {
                if (r.count >= 2) return {3, {}, {}, {}};
                r.pos[r.count]    = static_cast<uint16_t>(i);
                r.base_a[r.count] = ba;
                r.base_b[r.count] = bb;
                r.count++;
            }
            ++i;
        }
    }
    return r;
}

struct ChildMismatch {
    uint32_t parent_id;
    uint32_t child_id;
    uint16_t mismatch_pos;
    uint8_t  alt_base;
    uint8_t  parent_base;
    uint16_t mismatch_pos2;
    uint8_t  alt_base2;
    uint8_t  parent_base2;
    uint8_t  hamming;
    uint8_t  _pad[3];
    float    lr_score = std::numeric_limits<float>::quiet_NaN();
};
static_assert(sizeof(ChildMismatch) == 24, "ChildMismatch must be 24 bytes");

// Extract bases [start, start+count) from 2-bit packed data into dst[],
// byte-normalized (first base at MSB of dst[0], trailing bits zeroed).
// Returns number of bytes written.
inline int extract_packed_part(const uint8_t* packed, int start, int count, uint8_t* dst) {
    int nbytes = (count + 3) / 4;
    int src_byte = start >> 2;
    int src_shift = (start & 3) << 1;

    if (src_shift == 0) {
        std::memcpy(dst, packed + src_byte, nbytes);
    } else {
        int rshift = 8 - src_shift;
        for (int i = 0; i < nbytes; ++i) {
            dst[i] = static_cast<uint8_t>(
                (packed[src_byte + i] << src_shift) |
                (packed[src_byte + i + 1] >> rshift));
        }
    }
    int used = count & 3;
    if (used) dst[nbytes - 1] &= static_cast<uint8_t>(0xFFu << ((4 - used) << 1));
    return nbytes;
}

// Compute reverse-complement of a pre-extracted packed interior buffer.
inline void compute_interior_rc(const uint8_t* interior, int ilen, uint8_t* out) {
    int nbytes = (ilen + 3) / 4;
    for (int i = 0; i < nbytes; ++i) out[i] = 0;
    for (int i = 0; i < ilen; ++i) {
        int src = ilen - 1 - i;
        uint8_t base = (interior[src >> 2] >> (6 - 2 * (src & 3))) & 0x3u;
        uint8_t comp = base ^ 0x3u;
        out[i >> 2] |= static_cast<uint8_t>(comp << (6 - 2 * (i & 3)));
    }
}

// Count mismatches between two pre-extracted interior packed buffers (same ilen).
// Both buffers must be (ilen+3)/4 bytes with trailing bits zeroed (as produced
// by extract_packed_part). Returns actual count, or max_diffs+1 on early exit.
inline int hamming_packed_max(const uint8_t* pa, const uint8_t* pb,
                               int ilen, int max_diffs) {
    int count = 0;
    int nb    = (ilen + 3) / 4;
    for (int i = 0; i < nb; ++i) {
        uint8_t diff = pa[i] ^ pb[i];
        if (diff) {
            uint8_t odd  = static_cast<uint8_t>((diff >> 1) & 0x55u);
            uint8_t even = static_cast<uint8_t>( diff        & 0x55u);
            count += __builtin_popcount(static_cast<uint8_t>(odd | even));
            if (count > max_diffs) return max_diffs + 1;
        }
    }
    return count;
}

}  // namespace fqdup::derep_detail
