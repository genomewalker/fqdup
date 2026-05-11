#pragma once
// Phase B3: damage-aware H>2 deduplication helpers.
// Included only from derep.cpp (within its anonymous namespace) after
// packed_ops.hpp and damage_keys.hpp are already included.
//
// Rationale: for 2+ Myr aDNA with d_max ≥ 0.25, PCR copies of the same
// ancient molecule can accumulate 3–5 deamination differences in positions
// just outside the mask zone (p_damage 1–4%). Phase 3 pigeonhole (H≤2) never
// pairs these reads. K_deam hashing normalizes damage-probable positions so
// reads from the same molecule land in the same bucket; the diff filter and
// LRT then decide absorption.

#include "fqdup/damage_profile.hpp"
#include "packed_ops.hpp"
#include <xxhash.h>

namespace fqdup::b3 {

static constexpr int kMaxDiffs = 8;

struct DiffList {
    int      count;                  // actual count; kMaxDiffs+1 if exceeded
    uint16_t pos[kMaxDiffs];         // interior-relative positions (0-based)
    uint8_t  base_p[kMaxDiffs];      // parent 2-bit base
    uint8_t  base_c[kMaxDiffs];      // child  2-bit base
};

// Enumerate mismatches between two pre-extracted interior packed buffers.
// Stops filling at kMaxDiffs; sets count = kMaxDiffs+1 if there are more.
inline void packed_diff_list(const uint8_t* pa, const uint8_t* pb,
                              int ilen, DiffList& out) {
    out.count = 0;
    for (int i = 0; i < ilen; ++i) {
        uint8_t ba = (pa[i >> 2] >> (6 - 2 * (i & 3))) & 0x3u;
        uint8_t bb = (pb[i >> 2] >> (6 - 2 * (i & 3))) & 0x3u;
        if (ba != bb) {
            if (out.count >= kMaxDiffs) { out.count = kMaxDiffs + 1; return; }
            out.pos   [out.count] = static_cast<uint16_t>(i);
            out.base_p[out.count] = ba;
            out.base_c[out.count] = bb;
            ++out.count;
        }
    }
}

// Is this mismatch deamination-consistent at the given absolute read position?
// near_5 = pos_abs < L/2; only C→T allowed at 5' end, G→A at 3' end (DS),
// both at both ends (SS). p_damage must already exceed t_deam (caller's gate).
inline bool is_deam_direction(uint8_t ba, uint8_t bb, bool near_5, bool ss_mode) {
    uint8_t xr = ba ^ bb;
    if (xr != 2u) return false;                       // not C↔T or G↔A
    bool is_ct = ((ba & bb) == 1u);                   // C(01)↔T(11): low bit set
    bool is_ga = ((ba & bb) == 0u);                   // G(10)↔A(00): both bits clear
    if (ss_mode) return is_ct || is_ga;
    return near_5 ? is_ct : is_ga;
}

// True if ALL mismatches in the diff list are deamination-consistent:
// each must be C↔T / G↔A in the correct orientation AND at a position
// where P(damage | profile, L) > t_deam.
inline bool all_diffs_deam_consistent(const DiffList& diffs,
                                       int k5, int L,
                                       const DamageProfile& prof,
                                       float t_deam) {
    if (diffs.count < 1 || diffs.count > kMaxDiffs) return false;
    for (int i = 0; i < diffs.count; ++i) {
        int pos_abs = k5 + static_cast<int>(diffs.pos[i]);
        if (prof.p_damage_at(pos_abs, L) < static_cast<double>(t_deam)) return false;
        bool near_5 = (pos_abs < L / 2);
        if (!is_deam_direction(diffs.base_p[i], diffs.base_c[i], near_5, prof.ss_mode))
            return false;
    }
    return true;
}

// Normalize damage-probable bases in a pre-extracted interior buffer in-place.
// For each interior position ip (0-based), if P(damage at k5+ip) > t_deam:
//   DS: T→C at 5'-proximal positions; A→G at 3'-proximal positions.
//   SS: T→C and A→G at both ends.
// Reads from the same ancient molecule with different deamination states
// at these positions will yield the same normalized buffer.
inline void normalize_deam_bases(uint8_t* buf, int k5, int ilen, int L,
                                  const DamageProfile& prof, float t_deam) {
    for (int ip = 0; ip < ilen; ++ip) {
        int pos_abs = k5 + ip;
        if (prof.p_damage_at(pos_abs, L) < static_cast<double>(t_deam)) continue;
        int     byte_i    = ip >> 2;
        int     shift     = 6 - 2 * (ip & 3);
        uint8_t base      = (buf[byte_i] >> shift) & 0x3u;
        uint8_t normalized = base;
        bool near_5 = (pos_abs < L / 2);
        if (prof.ss_mode) {
            if (base == 3u) normalized = 1u;       // T→C (deaminated C)
            else if (base == 0u) normalized = 2u;  // A→G (deaminated G on RC)
        } else if (near_5) {
            if (base == 3u) normalized = 1u;       // T→C
        } else {
            if (base == 0u) normalized = 2u;       // A→G
        }
        if (normalized != base)
            buf[byte_i] = static_cast<uint8_t>(
                (buf[byte_i] & ~(0x3u << shift)) | (normalized << shift));
    }
}

// Compute K_deam bucket key: hash of the damage-normalized interior bytes.
// ilen is folded into the high 16 bits to prevent cross-length collisions.
// scratch must be at least (ilen+3)/4 bytes (caller-allocated, reused across calls).
inline uint64_t kdamage_hash(const uint8_t* arena_packed, int k5, int ilen, int L,
                              const DamageProfile& prof, float t_deam,
                              uint8_t* scratch) {
    using namespace fqdup::derep_detail;
    extract_packed_part(arena_packed, k5, ilen, scratch);
    normalize_deam_bases(scratch, k5, ilen, L, prof, t_deam);
    int nb = (ilen + 3) / 4;
    uint64_t h = XXH3_64bits(scratch, static_cast<size_t>(nb));
    return h ^ (static_cast<uint64_t>(static_cast<uint16_t>(ilen)) << 48);
}

// Count-based LRT identical to interior_transition_lrt — reused for B3.
inline double b3_count_lrt(uint64_t n_parent, uint64_t n_child,
                            double f0, double f1) {
    return fqdup::derep_detail::interior_transition_lrt(n_parent, n_child, f0, f1);
}

}  // namespace fqdup::b3
