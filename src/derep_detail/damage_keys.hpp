#pragma once
// Damage-aware canonical hashing helpers used by derep.cpp.
// All routines are heap-free — caller provides scratch buffers.

#include "fqdup/damage_profile.hpp"

#include <cctype>
#include <cstring>
#include <string>
#include <utility>

#include <xxhash.h>

namespace fqdup::derep_detail {

// Replace damage-prone terminal bases with neutral bytes before hashing.
// scratch must point to a buffer of at least seq.size() bytes.
inline void apply_damage_mask_inplace(const std::string& seq,
                                      const DamageProfile& prof,
                                      char* scratch) {
    int L = static_cast<int>(seq.size());
    std::memcpy(scratch, seq.data(), L);
    for (int i = 0; i < L; ++i) {
        char cu = static_cast<char>(
            std::toupper(static_cast<unsigned char>(scratch[i])));
        bool in_5zone = (i         < DamageProfile::MASK_POSITIONS) && prof.mask_pos[i];
        bool in_3zone = (L - 1 - i < DamageProfile::MASK_POSITIONS) && prof.mask_pos[L - 1 - i];
        if (prof.ss_mode) {
            if ((in_5zone || in_3zone) && (cu == 'C' || cu == 'T')) {
                scratch[i] = '\x01';
            } else if ((in_5zone || in_3zone) && (cu == 'G' || cu == 'A')) {
                scratch[i] = '\x02';
            }
        } else {
            if (in_5zone && (cu == 'C' || cu == 'T')) {
                scratch[i] = '\x01';
            } else if (in_3zone && (cu == 'G' || cu == 'A')) {
                scratch[i] = '\x02';
            }
        }
    }
}

inline std::pair<int,int> damage_zone_bounds(int L, const DamageProfile& prof) {
    if (!prof.enabled) return {0, 0};
    int k5 = 0;
    while (k5 < L && k5 < DamageProfile::MASK_POSITIONS && prof.mask_pos[k5]) k5++;
    int k3 = 0;
    while (k3 < L - k5 && k3 < DamageProfile::MASK_POSITIONS && prof.mask_pos[k3]) k3++;
    return {k5, k3};
}

// scratch1 and scratch2 must each be at least seq.size() bytes.
inline XXH128_hash_t damage_canonical_hash(const std::string& seq,
                                           const DamageProfile& prof,
                                           bool use_revcomp,
                                           char* scratch1,
                                           char* scratch2) {
    int L = static_cast<int>(seq.size());
    apply_damage_mask_inplace(seq, prof, scratch1);
    XXH128_hash_t h1 = XXH3_128bits(scratch1, L);
    if (!use_revcomp) return h1;

    for (int i = 0; i < L; ++i) {
        unsigned char c = static_cast<unsigned char>(seq[L - 1 - i]);
        switch (c) {
            case 'A': case 'a': scratch2[i] = (c == 'A') ? 'T' : 't'; break;
            case 'C': case 'c': scratch2[i] = (c == 'C') ? 'G' : 'g'; break;
            case 'G': case 'g': scratch2[i] = (c == 'G') ? 'C' : 'c'; break;
            case 'T': case 't': scratch2[i] = (c == 'T') ? 'A' : 'a'; break;
            default:            scratch2[i] = 'N'; break;
        }
    }
    for (int i = 0; i < L; ++i) {
        char cu = static_cast<char>(std::toupper(static_cast<unsigned char>(scratch2[i])));
        bool in_5zone = (i         < DamageProfile::MASK_POSITIONS) && prof.mask_pos[i];
        bool in_3zone = (L - 1 - i < DamageProfile::MASK_POSITIONS) && prof.mask_pos[L - 1 - i];
        if (prof.ss_mode) {
            if ((in_5zone || in_3zone) && (cu == 'C' || cu == 'T')) scratch2[i] = '\x01';
            else if ((in_5zone || in_3zone) && (cu == 'G' || cu == 'A')) scratch2[i] = '\x02';
        } else {
            if (in_5zone && (cu == 'C' || cu == 'T')) scratch2[i] = '\x01';
            else if (in_3zone && (cu == 'G' || cu == 'A')) scratch2[i] = '\x02';
        }
    }
    XXH128_hash_t h2 = XXH3_128bits(scratch2, L);
    if (h1.high64 < h2.high64 || (h1.high64 == h2.high64 && h1.low64 <= h2.low64))
        return h1;
    return h2;
}

// Canonical hash without heap allocation for reverse-complement construction.
inline XXH128_hash_t canonical_hash_noalloc(const std::string& seq,
                                            bool use_revcomp,
                                            char* scratch) {
    XXH128_hash_t h1 = XXH3_128bits(seq.data(), seq.size());
    if (!use_revcomp) return h1;

    int L = static_cast<int>(seq.size());
    for (int i = 0; i < L; ++i) {
        unsigned char c = static_cast<unsigned char>(seq[L - 1 - i]);
        switch (c) {
            case 'A': case 'a': scratch[i] = (c == 'A') ? 'T' : 't'; break;
            case 'C': case 'c': scratch[i] = (c == 'C') ? 'G' : 'g'; break;
            case 'G': case 'g': scratch[i] = (c == 'G') ? 'C' : 'c'; break;
            case 'T': case 't': scratch[i] = (c == 'T') ? 'A' : 'a'; break;
            default:            scratch[i] = 'N'; break;
        }
    }
    XXH128_hash_t h2 = XXH3_128bits(scratch, seq.size());
    if (h1.high64 < h2.high64 || (h1.high64 == h2.high64 && h1.low64 <= h2.low64))
        return h1;
    return h2;
}

}  // namespace fqdup::derep_detail
