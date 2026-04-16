// derep.cpp
// Single-file FASTQ deduplication with damage-aware hashing and PCR error correction.
// Operates on a single sorted FASTQ (typically the non-extended output of derep_pairs).
//
// REQUIRES: Input MUST be sorted by read ID (use 'fqdup sort' first).
//
// Strategy:
//   Pass 0 (optional): Estimate ancient DNA damage parameters from input
//   Pass 1: Stream input, build hash → position index with damage masking
//   Phase 3 (optional): PCR error correction via 4-way pigeonhole H≤2 Hamming search
//   Pass 2: Stream again, write representative unique records
//
// Memory: ~16 bytes per input record + seq arena if --error-correct is used.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include "dart/frame_selector_decl.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __linux__
#include <malloc.h>
#endif

#ifdef __GLIBC__
extern "C" {
    int mallctl(const char *name, void *oldp, size_t *oldlenp,
                void *newp, size_t newlen) __attribute__((weak));
}
#endif

// fastq_common.hpp provides: FastqRecord, FastqReader, FastqReaderIgzip,
// FastqWriter, SequenceFingerprint, SequenceFingerprintHash, canonical_hash,
// revcomp, trim_id, shell_escape_path, GZBUF_SIZE — all in anonymous namespace.
#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"
#include "flat_hash_map.hpp"

// All file-local types are in an anonymous namespace to give their member
// functions internal linkage, avoiding ODR violations with other TUs.
namespace {

// ============================================================================
// Configuration
// ============================================================================

struct ErrCorParams {
    bool     enabled          = false;
    double   snp_threshold    = 0.20; // sig_count_weighted/parent_count >= this → SNP protection
    uint32_t snp_min_count    = 1;    // absolute minimum sig_count for SNP veto on strict-dominance edges
    uint32_t snp_low_cov_cutoff = 10;   // parent count below which adaptive threshold applies
    double   snp_low_cov_factor = 1.75; // snp_threshold multiplier when parent_count < cutoff
    uint32_t bucket_cap       = 64;
    double   pcr_phi          = 5.3e-7;
    double   pcr_rate         = 0.0;  // still estimated for logging but no longer gates absorption
    uint32_t max_h2_count     = 2;    // H=2 absorption only when child count ≤ this
};

struct Phase3Stats {
    double decode_hash_parent_ms = 0;
    double insert_ms             = 0;
    double decode_hash_child_ms  = 0;
    double query_ms              = 0;
    double check_ms              = 0;
    uint64_t total_candidates    = 0;
    uint64_t cap_fires           = 0;
    uint64_t children_scanned    = 0;
    uint64_t parents_indexed     = 0;
    uint64_t children_found      = 0;  // H=1 child-parent pairs detected
    uint64_t snp_protected       = 0;  // children not absorbed due to SNP veto
    uint64_t absorbed            = 0;  // children absorbed as PCR errors
    void log() const {
        auto f = [](double ms){ return std::to_string(static_cast<int>(ms)) + " ms"; };
        log_info("Phase 3 timing:");
        log_info("  Parent decode+hash : " + f(decode_hash_parent_ms));
        log_info("  Parent insert      : " + f(insert_ms));
        log_info("  Child decode+hash  : " + f(decode_hash_child_ms));
        log_info("  Child query        : " + f(query_ms));
        log_info("  Candidate check    : " + f(check_ms));
        log_info("  Parents indexed    : " + std::to_string(parents_indexed));
        log_info("  Children scanned   : " + std::to_string(children_scanned));
        log_info("  Children found     : " + std::to_string(children_found));
        log_info("  SNP protected      : " + std::to_string(snp_protected));
        log_info("  Absorbed           : " + std::to_string(absorbed));
        log_info("  Total candidates   : " + std::to_string(total_candidates));
        if (children_scanned > 0)
            log_info("  Avg cand/child     : " +
                     std::to_string(static_cast<double>(total_candidates)/children_scanned));
        log_info("  Bucket cap fires   : " + std::to_string(cap_fires));
        if (cap_fires > 0)
            log_warn("Phase 3: " + std::to_string(cap_fires) +
                     " bucket cap fire(s) — some PCR error candidates were not evaluated. "
                     "Consider increasing --errcor-bucket-cap (current cap applied at index build).");
    }
};

// ============================================================================
// Index entry — single-file version (no non_len; seq_id for error correction)
// ============================================================================

struct IndexEntry {
    uint64_t record_index;  // sequential record number of the representative
    uint64_t count;         // total reads in this cluster
    uint32_t seq_id;        // index into SeqArena (populated when --error-correct)
    uint8_t  damage_score;  // terminal C→T (5') + G→A (3') count of representative

    IndexEntry() : record_index(0), count(1), seq_id(0), damage_score(0) {}
    explicit IndexEntry(uint64_t idx) : record_index(idx), count(1), seq_id(0), damage_score(0) {}
};

// ============================================================================
// Phase 3: PCR Error Correction data structures
//
// Count-stratified 4-way pigeonhole neighbor finding for H≤2 detection.
// Reads with count <= max_count that have a neighbour with count >= ratio ×
// their own count are classified as PCR error duplicates.
//
// 4-way pigeonhole: split interior [k5, L-k3) into four parts p0, p1, p2, p3.
// If Hamming(child, parent) ≤ 2 inside the interior, at least 2 of the 4
// parts match exactly.  We index all C(4,2)=6 pair combinations so any
// 1- or 2-error variant is found by at least one lookup.
// ============================================================================

// 2-bit base helpers — defined before SeqArena which uses them.
static constexpr uint8_t kDec2bit[4] = {'A', 'C', 'G', 'T'};
static const std::array<uint8_t, 256> kEnc2bit = [] {
    std::array<uint8_t, 256> lut{};
    lut.fill(0xFFu);
    lut[static_cast<uint8_t>('A')] = lut[static_cast<uint8_t>('a')] = 0u;
    lut[static_cast<uint8_t>('C')] = lut[static_cast<uint8_t>('c')] = 1u;
    lut[static_cast<uint8_t>('G')] = lut[static_cast<uint8_t>('g')] = 2u;
    lut[static_cast<uint8_t>('T')] = lut[static_cast<uint8_t>('t')] = 3u;
    return lut;
}();

static const std::array<std::array<uint8_t, 4>, 256> kDec4 = [] {
    std::array<std::array<uint8_t, 4>, 256> lut{};
    for (int v = 0; v < 256; ++v) {
        lut[v][0] = kDec2bit[(v >> 6) & 0x3u];
        lut[v][1] = kDec2bit[(v >> 4) & 0x3u];
        lut[v][2] = kDec2bit[(v >> 2) & 0x3u];
        lut[v][3] = kDec2bit[(v >> 0) & 0x3u];
    }
    return lut;
}();

// 2-bit packed sequence arena.  Reduces memory ~4× vs ASCII storage.
// Non-ACGT bases (N etc.) are encoded as A=0 and the sequence is flagged
// ineligible so Phase 3 skips it (avoids false PCR error absorptions).
struct SeqArena {
    std::vector<uint8_t>  packed;    // 2-bit data; (L+3)/4 bytes per sequence
    std::vector<uint64_t> offsets;   // byte offset of each sequence start
    std::vector<uint16_t> lengths;   // base lengths
    std::vector<bool>     eligible;  // false when sequence contains non-ACGT

    uint32_t append(const std::string& seq) {
        if (seq.size() > 65535u)
            throw std::runtime_error("Sequence too long for arena (>65535 bp)");
        if (offsets.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("SeqArena overflow: more than 4 G unique sequences");
        uint32_t id    = static_cast<uint32_t>(offsets.size());
        uint16_t L     = static_cast<uint16_t>(seq.size());
        int      nbytes = (L + 3) / 4;
        offsets.push_back(packed.size());
        lengths.push_back(L);
        packed.resize(packed.size() + nbytes, 0);
        uint8_t* dst = packed.data() + offsets[id];
        bool ok = true;
        int i = 0, byte = 0;
        while (i < L) {
            uint8_t v = 0;
            for (int lane = 0; lane < 4 && i < L; ++lane, ++i) {
                uint8_t b = kEnc2bit[static_cast<uint8_t>(seq[i])];
                if (b == 0xFFu) { ok = false; b = 0; }
                v |= static_cast<uint8_t>(b << (6 - 2 * lane));
            }
            dst[byte++] = v;
        }
        eligible.push_back(ok);
        return id;
    }

    uint32_t append_chars(const char* data_in, int L) {
        if (L > 65535)
            throw std::runtime_error("Sequence too long for arena (>65535 bp)");
        if (offsets.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("SeqArena overflow: more than 4 G unique sequences");
        uint32_t id    = static_cast<uint32_t>(offsets.size());
        int      nbytes = (L + 3) / 4;
        offsets.push_back(packed.size());
        lengths.push_back(static_cast<uint16_t>(L));
        packed.resize(packed.size() + nbytes, 0);
        uint8_t* dst = packed.data() + offsets[id];
        bool ok = true;
        int i = 0, byte = 0;
        while (i < L) {
            uint8_t v = 0;
            for (int lane = 0; lane < 4 && i < L; ++lane, ++i) {
                uint8_t b = kEnc2bit[static_cast<uint8_t>(data_in[i])];
                if (b == 0xFFu) { ok = false; b = 0; }
                v |= static_cast<uint8_t>(b << (6 - 2 * lane));
            }
            dst[byte++] = v;
        }
        eligible.push_back(ok);
        return id;
    }

    const uint8_t* data(uint32_t id)        const { return packed.data() + offsets[id]; }
    uint16_t       length(uint32_t id)      const { return lengths[id]; }
    bool           is_eligible(uint32_t id) const { return eligible[id]; }
    size_t         size()                   const { return offsets.size(); }

    // Decode bases [start, start+count) to ASCII scratch buffer.
    void decode_range(uint32_t id, int start, int count, uint8_t* out) const {
        const uint8_t* p = data(id);
        int i = 0;
        int pos = start;

        // Handle leading unaligned bases until we can decode full bytes.
        while (i < count && (pos & 3) != 0) {
            out[i++] = kDec2bit[(p[pos >> 2] >> (6 - 2 * (pos & 3))) & 0x3u];
            ++pos;
        }

        while (i + 4 <= count) {
            const auto& d = kDec4[p[pos >> 2]];
            std::memcpy(out + i, d.data(), 4);
            i += 4;
            pos += 4;
        }

        while (i < count) {
            out[i++] = kDec2bit[(p[pos >> 2] >> (6 - 2 * (pos & 3))) & 0x3u];
            ++pos;
        }
    }
};

struct FlatPairIndex {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> ids;

    template <typename F>
    void query(uint64_t key, F&& fn) const {
        auto it = std::lower_bound(keys.begin(), keys.end(), key);
        if (it == keys.end() || *it != key) return;
        size_t idx = static_cast<size_t>(it - keys.begin());
        for (uint32_t i = offsets[idx]; i < offsets[idx + 1]; ++i)
            fn(ids[i]);
    }
};

static inline void split4_lens(int ilen, int& s0, int& s1, int& s2, int& s3) {
    s0 = ilen / 4;
    s1 = ilen * 2 / 4 - s0;
    s2 = ilen * 3 / 4 - s0 - s1;
    s3 = ilen - s0 - s1 - s2;
}

static inline uint64_t splitmix64(uint64_t x) {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t pair_key(uint64_t ha, uint64_t hb, int tag, int ilen) {
    if (ha > hb) std::swap(ha, hb);
    return splitmix64(ha ^ rotl64(hb, 23)
                      ^ (static_cast<uint64_t>(tag) << 56)
                      ^ static_cast<uint64_t>(ilen));
}


// Count Hamming mismatches between packed interiors [k5, k5+ilen).
// Returns the count, stopping at max_mm+1 (early exit if exceeded).
static int packed_hamming(const uint8_t* pa, const uint8_t* pb,
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
// using 2-bit encoded values.
//
// XOR analysis (A=0,C=1,G=2,T=3):
//   xr=1 (01b): G↔T (10^11) or C↔A (01^00) — 8-oxoG, protect in damage mode
//   xr=2 (10b): G↔A (10^00) or C↔T (01^11) — deamination, protect in damage mode
//   xr=3 (11b): A↔T or C↔G — transversions, never damage
static bool is_damage_sub_packed(uint8_t a, uint8_t b, bool /*protect_deamination*/) {
    // C↔T (xr=2) and G↔A/C↔A/G↔T (xr=1) are the defining aDNA damage signatures.
    // They are ALWAYS protected from error-correction absorption regardless of whether
    // --collapse-damage is active: even without explicit damage mode these substitutions
    // are too suspicious to absorb (could be genuine variants or undetected damage).
    // A↔T / C↔G (xr=3): transversions — no known damage mechanism, always absorbable.
    uint8_t xr = a ^ b;
    return (xr == 1u || xr == 2u);
}

// Scan interior region [k5, k5+ilen) of two 2-bit packed sequences.
// Returns true if child should be absorbed by parent: exactly 1 interior
// mismatch AND the substitution is NOT damage-consistent.
// Fast path: full-byte XOR with 2-bit popcount via (x | x>>1) & 0x55.
static bool packed_should_absorb(const uint8_t* pa, const uint8_t* pb,
                                  int k5, int ilen, bool protect_deamination) {
    int mismatches = 0;
    uint8_t mm_a = 0, mm_b = 0;

    for (int i = 0; i < ilen; ) {
        int pos   = k5 + i;
        int byte_i = pos >> 2;
        int bit_i  = pos & 3;

        if (bit_i == 0 && i + 4 <= ilen) {
            // Full aligned byte: check all 4 bases at once
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
            // Partial byte: single base
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
// Returns {false,...} if mismatches != 1 or the substitution is damage-consistent.
static MismatchInfo packed_find_mismatch(const uint8_t* pa, const uint8_t* pb,
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
                    // sh is 0,2,4,6 bits from MSB side → offset within 4-base block
                    // Bit position from MSB: bit 7=base0, bit 5=base1, bit 3=base2, bit 1=base3
                    // sh=6→base0, sh=4→base1, sh=2→base2, sh=0→base3
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
// No damage filtering — caller is responsible.
// r.count == 3 means more than 2 mismatches (rejected).
struct MismatchInfo2 {
    int     count;
    uint16_t pos[2];
    uint8_t  base_a[2];  // parent's 2-bit base
    uint8_t  base_b[2];  // child's 2-bit base
};

static MismatchInfo2 packed_find_mismatches2(const uint8_t* pa, const uint8_t* pb,
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

// ── ChildMismatch record ────────────────────────────────────────────────────
struct ChildMismatch {
    uint32_t parent_id;
    uint32_t child_id;
    uint16_t mismatch_pos;   // first mismatch (interior-relative, 0-based)
    uint8_t  alt_base;       // child's 2-bit base (first mismatch)
    uint8_t  parent_base;    // parent's 2-bit base (first mismatch)
    uint16_t mismatch_pos2;  // second mismatch (only used when hamming==2)
    uint8_t  alt_base2;      // child's 2-bit base (second mismatch)
    uint8_t  parent_base2;   // parent's 2-bit base (second mismatch)
    uint8_t  hamming;        // 1 or 2
    uint8_t  _pad[3];
};  // 20 bytes
static_assert(sizeof(ChildMismatch) == 20, "ChildMismatch must be 20 bytes");

// Extract bases [start, start+count) from 2-bit packed data into dst[],
// byte-normalized (first base at MSB of dst[0], trailing bits zeroed).
// dst must have at least (count+3)/4 bytes.
// Returns number of bytes written.
static int extract_packed_part(const uint8_t* packed, int start, int count, uint8_t* dst) {
    int nbytes = (count + 3) / 4;
    int src_byte = start >> 2;
    int src_shift = (start & 3) << 1;  // 0,2,4,6 bits into first src byte

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
    // Zero trailing unused bits in last byte
    int used = count & 3;
    if (used) dst[nbytes - 1] &= static_cast<uint8_t>(0xFFu << ((4 - used) << 1));
    return nbytes;
}

// Compute reverse-complement of a pre-extracted packed interior buffer.
// interior[] holds `ilen` bases packed MSB-first starting at bit offset 0.
// Writes the RC into out[] (same size: (ilen+3)/4 bytes).
// Complement: each 2-bit base b → b^3 (A↔T=00↔11, C↔G=01↔10).
static void compute_interior_rc(const uint8_t* interior, int ilen, uint8_t* out) {
    int nbytes = (ilen + 3) / 4;
    for (int i = 0; i < nbytes; ++i) out[i] = 0;
    for (int i = 0; i < ilen; ++i) {
        int src = ilen - 1 - i;
        uint8_t base = (interior[src >> 2] >> (6 - 2 * (src & 3))) & 0x3u;
        uint8_t comp = base ^ 0x3u;
        out[i >> 2] |= static_cast<uint8_t>(comp << (6 - 2 * (i & 3)));
    }
}

// ============================================================================
// Ancient DNA Damage Model
// DamageProfile struct and estimate_damage() are in damage_profile.hpp/.cpp
// (shared with extend.cpp).  TU-specific helpers follow below.
// ============================================================================

// Replace damage-prone terminal bases with neutral bytes before hashing.
//
// Purpose: two reads from the same ancient DNA molecule can differ at terminal
// positions due to cytosine deamination (C→T at the 5' end, G→A at the 3' end).
// An exact-sequence hash would treat them as distinct — inflating unique cluster
// counts. By replacing any C or T at a 5'-damaged position with \x01, and any
// G or A at a 3'-damaged position with \x02, independently-deaminated copies of
// the same molecule hash to the same value and collapse into a single cluster.
//
// Which positions are masked is determined by DamageProfile::mask_pos[], an
// empirical array populated directly from observed per-position T/(T+C) and
// A/(A+G) frequencies (see estimate_damage()). No exponential model is evaluated
// at hash time — just an array lookup per position.
//
// Symmetry: mask_pos uses max(5'-signal, 3'-signal) semantics so that position p
// from the 5' end and position p from the 3' end are treated identically. This
// guarantees canonical_hash(seq) == canonical_hash(revcomp(seq)).
//
// scratch must point to a buffer of at least seq.size() bytes. Writes the
// masked sequence into scratch without heap allocation.
static void apply_damage_mask_inplace(const std::string& seq,
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
            // SS: C→T at both ends. Mask the full Watson-Crick transition pair at every
            // damage-zone position so that canonical_hash(seq)==canonical_hash(revcomp(seq)):
            // forward C/T → \x01; revcomp of those positions has G/A → \x02.
            if ((in_5zone || in_3zone) && (cu == 'C' || cu == 'T')) {
                scratch[i] = '\x01';
            } else if ((in_5zone || in_3zone) && (cu == 'G' || cu == 'A')) {
                scratch[i] = '\x02';
            }
        } else {
            // DS: C→T at 5', G→A at 3'.
            if (in_5zone && (cu == 'C' || cu == 'T')) {
                scratch[i] = '\x01';
            } else if (in_3zone && (cu == 'G' || cu == 'A')) {
                scratch[i] = '\x02';
            }
        }
    }
}

static std::pair<int,int> damage_zone_bounds(int L, const DamageProfile& prof) {
    if (!prof.enabled) return {0, 0};
    int k5 = 0;
    while (k5 < L && k5 < DamageProfile::MASK_POSITIONS && prof.mask_pos[k5]) k5++;
    int k3 = 0;
    while (k3 < L - k5 && k3 < DamageProfile::MASK_POSITIONS && prof.mask_pos[k3]) k3++;
    return {k5, k3};
}

// scratch1 and scratch2 must each be at least seq.size() bytes.
// No heap allocation — uses caller-provided scratch buffers.
static XXH128_hash_t damage_canonical_hash(const std::string& seq,
                                           const DamageProfile& prof,
                                           bool use_revcomp,
                                           char* scratch1,
                                           char* scratch2) {
    int L = static_cast<int>(seq.size());
    apply_damage_mask_inplace(seq, prof, scratch1);
    XXH128_hash_t h1 = XXH3_128bits(scratch1, L);
    if (!use_revcomp) return h1;

    // Build revcomp+mask into scratch2 without allocating
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
    // Canonical = lexicographically smaller of the two 128-bit hashes
    if (h1.high64 < h2.high64 || (h1.high64 == h2.high64 && h1.low64 <= h2.low64))
        return h1;
    return h2;
}

// Canonical hash without heap allocation for reverse-complement construction.
// Matches fastq_common.hpp canonical_hash() behavior.
static XXH128_hash_t canonical_hash_noalloc(const std::string& seq,
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


// ============================================================================
// DerepEngine — single-file two-pass deduplication
// ============================================================================

class DerepEngine {
public:
    DerepEngine(bool use_revcomp, bool write_clusters,
                const DamageProfile& profile = DamageProfile{},
                const ErrCorParams& errcor = ErrCorParams{},
                bool /*errcor_mode_explicit*/ = false,
                bool bucket_cap_explicit = false)
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          profile_(profile), errcor_(errcor),
          bucket_cap_explicit_(bucket_cap_explicit),
          total_reads_(0), errcor_absorbed_(0), n_unique_clusters_(0) {}

    void process(const std::string& in_path,
                 const std::string& out_path,
                 const std::string& cluster_path) {

        log_info("=== Two-pass single-file deduplication ===");
        log_info("Pass 1: Build lightweight index");
        log_info("Decompression: " + std::string(
#ifdef HAVE_RAPIDGZIP
            "rapidgzip (parallel multi-threaded)"
#elif defined(HAVE_ISAL)
            "ISA-L (hardware-accelerated)"
#else
            "zlib"
#endif
        ));

        pass1(in_path);

        size_t index_mb = (index_.size() *
                           sizeof(std::pair<SequenceFingerprint, IndexEntry>)) /
                          (1024 * 1024);
        log_info("Index size: " + std::to_string(index_mb) + " MB for " +
                 std::to_string(total_reads_) + " reads");

        if (errcor_.enabled) {
            // Compute D_eff from duplication ratio.  Used both for PCR rate estimation
            // and coverage regime auto-detection.
            // Under PCR kinetics: D_eff = log2(total_reads / unique_reads).
            // This is exact for uniform amplification; a lower bound when starting
            // copy number varies (conservative: under-estimates errors, avoids
            // false absorptions).
            double d_eff = 0.0;
            if (total_reads_ > index_.size() && index_.size() > 0) {
                d_eff = std::log2(static_cast<double>(total_reads_) /
                                  static_cast<double>(index_.size()));
            }
            if (errcor_.pcr_rate == 0.0 && errcor_.pcr_phi > 0.0 && d_eff > 0.0) {
                errcor_.pcr_rate = errcor_.pcr_phi * d_eff;
                log_info("Phase 3: D_eff=" +
                         std::to_string(d_eff).substr(0, 5) +
                         " estimated from duplication ratio " +
                         std::to_string(total_reads_) + "/" +
                         std::to_string(index_.size()) +
                         " (use --pcr-cycles for explicit value)");
            }

            log_info("Phase 3: parent-centric mismatch pattern detection");
            phase3_error_correct();
        }

        // Free arena_ — only needed for phase3, not pass2.
        { SeqArena empty; std::swap(arena_, empty); }

        // Return freed phase3 structures (build_map, shards, id_count, acc_count)
        // to the OS before pass2 allocates records_to_write.  Without this,
        // jemalloc holds freed pages in its arenas and they count toward RSS.
#ifdef __linux__
        if (mallctl != nullptr) {
            try {
                mallctl("thread.tcache.flush", nullptr, nullptr, nullptr, 0);
                mallctl("arena.0.purge", nullptr, nullptr, nullptr, 0);
            } catch (...) {}
        }
        malloc_trim(0);
#endif

        log_info("Pass 2: Write unique records");

        pass2(in_path, out_path, cluster_path);

        print_stats();
    }

private:
    XXH128_hash_t compute_hash(const std::string& seq, bool& is_forward) const {
        int L = static_cast<int>(seq.size());
        if (profile_.enabled) {
            if (scratch1_.size() < seq.size()) scratch1_.resize(seq.size());
            if (scratch2_.size() < seq.size()) scratch2_.resize(seq.size());
            apply_damage_mask_inplace(seq, profile_, scratch1_.data());
            XXH128_hash_t h1 = XXH3_128bits(scratch1_.data(), L);
            if (!use_revcomp_) { is_forward = true; return h1; }
            // Build RC into scratch2_
            for (int i = 0; i < L; ++i) {
                unsigned char c = static_cast<unsigned char>(seq[L - 1 - i]);
                switch (c) {
                    case 'A': case 'a': scratch2_[i] = (c == 'A') ? 'T' : 't'; break;
                    case 'C': case 'c': scratch2_[i] = (c == 'C') ? 'G' : 'g'; break;
                    case 'G': case 'g': scratch2_[i] = (c == 'G') ? 'C' : 'c'; break;
                    case 'T': case 't': scratch2_[i] = (c == 'T') ? 'A' : 'a'; break;
                    default:            scratch2_[i] = 'N'; break;
                }
            }
            // Apply damage mask to RC
            for (int i = 0; i < L; ++i) {
                char cu = static_cast<char>(std::toupper(static_cast<unsigned char>(scratch2_[i])));
                bool in_5zone = (i         < DamageProfile::MASK_POSITIONS) && profile_.mask_pos[i];
                bool in_3zone = (L - 1 - i < DamageProfile::MASK_POSITIONS) && profile_.mask_pos[L - 1 - i];
                if (profile_.ss_mode) {
                    if ((in_5zone || in_3zone) && (cu == 'C' || cu == 'T')) scratch2_[i] = '\x01';
                    else if ((in_5zone || in_3zone) && (cu == 'G' || cu == 'A')) scratch2_[i] = '\x02';
                } else {
                    if (in_5zone && (cu == 'C' || cu == 'T')) scratch2_[i] = '\x01';
                    else if (in_3zone && (cu == 'G' || cu == 'A')) scratch2_[i] = '\x02';
                }
            }
            XXH128_hash_t h2 = XXH3_128bits(scratch2_.data(), L);
            is_forward = (h1.high64 < h2.high64 ||
                          (h1.high64 == h2.high64 && h1.low64 <= h2.low64));
            return is_forward ? h1 : h2;
        }
        // No damage masking
        XXH128_hash_t h1 = XXH3_128bits(seq.data(), seq.size());
        if (!use_revcomp_) { is_forward = true; return h1; }
        if (scratch1_.size() < seq.size()) scratch1_.resize(seq.size());
        for (int i = 0; i < L; ++i) {
            unsigned char c = static_cast<unsigned char>(seq[L - 1 - i]);
            switch (c) {
                case 'A': case 'a': scratch1_[i] = (c == 'A') ? 'T' : 't'; break;
                case 'C': case 'c': scratch1_[i] = (c == 'C') ? 'G' : 'g'; break;
                case 'G': case 'g': scratch1_[i] = (c == 'G') ? 'C' : 'c'; break;
                case 'T': case 't': scratch1_[i] = (c == 'T') ? 'A' : 'a'; break;
                default:            scratch1_[i] = 'N'; break;
            }
        }
        XXH128_hash_t h2 = XXH3_128bits(scratch1_.data(), L);
        is_forward = (h1.high64 < h2.high64 ||
                      (h1.high64 == h2.high64 && h1.low64 <= h2.low64));
        return is_forward ? h1 : h2;
    }

    // Count terminal damage markers: T at masked 5' positions + A at masked 3' positions.
    // Higher score = more terminal deamination signal = more likely to be authentic ancient DNA.
    // Used to select the most-damaged read as cluster representative, maximising the
    // damage signal in the output.
    static uint8_t compute_damage_score(const std::string& seq,
                                        const DamageProfile& prof) {
        int score = 0;
        int L = static_cast<int>(seq.size());
        for (int p = 0; p < DamageProfile::MASK_POSITIONS && p < L; ++p) {
            if (!prof.mask_pos[p]) continue;
            if (std::toupper(static_cast<unsigned char>(seq[p]))       == 'T') ++score;
            if (std::toupper(static_cast<unsigned char>(seq[L-1-p]))   == 'A') ++score;
        }
        return static_cast<uint8_t>(std::min(score, 255));
    }

    void pass1(const std::string& in_path) {
        auto reader = make_fastq_reader(in_path);
        FastqRecord rec;
        uint64_t record_idx = 0;

        while (reader->read(rec)) {
            bool is_forward = true;
            XXH128_hash_t h = compute_hash(rec.seq, is_forward);
            SequenceFingerprint fp(h, rec.seq.size());

            auto [it, inserted] = index_.emplace(fp, IndexEntry(record_idx));
            if (inserted) {
                if (profile_.enabled)
                    it->second.damage_score = compute_damage_score(rec.seq, profile_);
                if (errcor_.enabled) {
                    if (is_forward) {
                        it->second.seq_id = arena_.append(rec.seq);
                    } else {
                        // Store in canonical orientation so Phase 3 comparisons are correct.
                        // Compute the unmasked revcomp into rc_buf_.
                        int L = static_cast<int>(rec.seq.size());
                        if (rc_buf_.size() < static_cast<size_t>(L)) rc_buf_.resize(L);
                        for (int i = 0; i < L; ++i) {
                            unsigned char c = static_cast<unsigned char>(rec.seq[L - 1 - i]);
                            switch (c) {
                                case 'A': case 'a': rc_buf_[i] = 'T'; break;
                                case 'C': case 'c': rc_buf_[i] = 'G'; break;
                                case 'G': case 'g': rc_buf_[i] = 'C'; break;
                                case 'T': case 't': rc_buf_[i] = 'A'; break;
                                default:            rc_buf_[i] = 'N'; break;
                            }
                        }
                        it->second.seq_id = arena_.append_chars(rc_buf_.data(), L);
                    }
                }
            } else {
                it->second.count++;
                // Maximize the ancient damage signal in the output: when a new read
                // in this cluster carries more terminal C→T / G→A than the current
                // representative, swap to it.  seq_id is intentionally kept pointing
                // to the first-occurrence sequence so Phase-3 EC is unaffected.
                if (profile_.enabled) {
                    uint8_t score = compute_damage_score(rec.seq, profile_);
                    if (score > it->second.damage_score) {
                        it->second.record_index = record_idx;
                        it->second.damage_score = score;
                    }
                }
            }

            record_idx++;
            total_reads_++;

            if ((total_reads_ % 1000000) == 0) {
                size_t unique = index_.size();
                double dup_pct = 100.0 * (1.0 - (double)unique / total_reads_);
                std::cerr << "\r[Pass 1] " << total_reads_ << " reads, "
                          << unique << " unique, "
                          << std::fixed << std::setprecision(1) << dup_pct << "% dedup"
                          << std::flush;
            }
        }

        std::cerr << "\r";
        log_info("Pass 1 complete: " + std::to_string(total_reads_) + " reads indexed");
    }

    // Phase 3: parent-centric mismatch pattern detection.
    // Indexes sequences with count > min_parent_count as parents, then for each
    // potential child (count <= min_parent_count) finds all H=1 parent neighbours.
    // A child is absorbed unless the mismatch pattern is recurrent (SNP veto):
    //   sig_count_weighted >= snp_min_count AND
    //   sig_count_weighted / parent_count >= snp_threshold.
    void phase3_error_correct() {
        if (arena_.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("Phase 3: arena size exceeds uint32_t range");
        const uint32_t N = static_cast<uint32_t>(arena_.size());
        if (N == 0) return;

        // Trailing padding for safe extraction in extract_packed_part
        arena_.packed.push_back(0);

        is_error_.assign(N, false);

        std::vector<uint64_t> id_count(N, 0);
        // Reverse map: seq_id → IndexEntry* for representative propagation during absorption.
        std::vector<IndexEntry*> seq_entry(N, nullptr);
        for (auto& [fp, entry] : index_) {
            id_count[entry.seq_id] = entry.count;
            seq_entry[entry.seq_id] = &entry;
        }

        // Accumulated mass per cluster: starts as Pass-1 count, grows as children
        // are absorbed.  Used for B3 boundary check only — the SNP veto still uses
        // id_count (original Pass-1 count) so the ratio threshold remains calibrated.
        std::vector<uint64_t> acc_count(id_count);

        // Scratch for packed interior hashing
        std::vector<uint8_t> scratch;
        scratch.reserve(65535);

        // Number of interior positions from each end considered "near-damage zone".
        // Positions just beyond the mask boundary retain residual deamination damage
        // (exponential tail, typically 5-10% at the first unmasked position).
        // C↔T / G↔A mismatches (xr=2) at these positions are treated as damage
        // variants and bypass the SNP veto, rather than being protected as SNPs.
        constexpr int kDamageEdgeMargin = 5;
        // Minimum interior length for a meaningful 3-way pigeonhole split.
        // Below this, splits degenerate (parts of 0-4 bp) and the hash keys
        // become unreliable, risking false absorptions.
        constexpr int kMinInteriorLen = 20;
        uint64_t short_interior_skipped = 0;

        struct InteriorLayout {
            int  k5 = 0, ilen = 0;
            int  s0 = 0, s1 = 0, s2 = 0, s3 = 0;
            int  nb0 = 0, nb1 = 0, nb2 = 0, nb3 = 0;
            int  nbytes = 0;
            bool ready = false;
        };
        std::vector<InteriorLayout> layout_cache(65536);
        auto get_layout = [&](int L) -> const InteriorLayout& {
            auto& lay = layout_cache[static_cast<uint16_t>(L)];
            if (!lay.ready) {
                auto [k5, k3] = damage_zone_bounds(L, profile_);
                lay.k5   = k5;
                lay.ilen = std::max(0, L - k5 - k3);
                if (lay.ilen >= kMinInteriorLen) {
                    split4_lens(lay.ilen, lay.s0, lay.s1, lay.s2, lay.s3);
                    lay.nb0    = (lay.s0 + 3) / 4;
                    lay.nb1    = (lay.s1 + 3) / 4;
                    lay.nb2    = (lay.s2 + 3) / 4;
                    lay.nb3    = (lay.s3 + 3) / 4;
                    lay.nbytes = lay.nb0 + lay.nb1 + lay.nb2 + lay.nb3;
                }
                lay.ready = true;
            }
            return lay;
        };

        Phase3Stats stats;
        using clk = std::chrono::steady_clock;

        // ── Phase A: index parents (count > min_parent_count) ──────────────
        struct BuildEntry { uint64_t key; uint32_t id; };
        // Tags 0-5: pair-hash for all C(4,2)=6 combinations of 4 parts.
        // 4-way pigeonhole: if H(child,parent)≤2, at least 2 of 4 parts match →
        // one of the 6 pair-keys fires.  Covers both H=1 and H=2.
        static constexpr int kPairA[6] = {0, 0, 0, 1, 1, 2};
        static constexpr int kPairB[6] = {1, 2, 3, 2, 3, 3};
        ska::flat_hash_map<int, std::array<std::vector<BuildEntry>, 6>> build_map;
        build_map.reserve(64);
        struct LenShard { std::array<FlatPairIndex, 6> pi; };

        uint64_t n_parents = 0;
        for (uint32_t id = 0; id < N; ++id) {
            if (!arena_.is_eligible(id)) continue;
            int L = arena_.length(id);
            const auto& lay = get_layout(L);
            if (lay.ilen < kMinInteriorLen) continue;

            auto t0 = clk::now();
            if (scratch.size() < static_cast<size_t>(lay.nbytes))
                scratch.resize(lay.nbytes);
            const uint8_t* psrc = arena_.data(id);
            int starts[4] = {lay.k5,
                             lay.k5 + lay.s0,
                             lay.k5 + lay.s0 + lay.s1,
                             lay.k5 + lay.s0 + lay.s1 + lay.s2};
            int sizes[4]  = {lay.s0, lay.s1, lay.s2, lay.s3};
            int nbs[4]    = {lay.nb0, lay.nb1, lay.nb2, lay.nb3};
            uint8_t* parts[4] = {scratch.data(),
                                 scratch.data() + lay.nb0,
                                 scratch.data() + lay.nb0 + lay.nb1,
                                 scratch.data() + lay.nb0 + lay.nb1 + lay.nb2};
            for (int p = 0; p < 4; ++p)
                extract_packed_part(psrc, starts[p], sizes[p], parts[p]);
            uint64_t h[4];
            for (int p = 0; p < 4; ++p)
                h[p] = XXH3_64bits(parts[p], nbs[p]);
            auto t1 = clk::now();

            auto& entries = build_map[lay.ilen];
            for (int t = 0; t < 6; ++t)
                entries[t].push_back({pair_key(h[kPairA[t]], h[kPairB[t]], t, lay.ilen), id});
            auto t2 = clk::now();

            stats.decode_hash_parent_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
            stats.insert_ms             += std::chrono::duration<double, std::milli>(t2 - t1).count();
            n_parents++;
        }
        stats.parents_indexed = n_parents;
        log_info("Phase 3: indexed " + std::to_string(n_parents) +
                 " sequences (directed monotone ascent, all eligible)");

        // Build CSR per shard per tag
        ska::flat_hash_map<int, LenShard> shards;
        shards.reserve(build_map.size());
        for (auto& [ilen, tag_entries] : build_map) {
            auto& sh = shards[ilen];
            for (int t = 0; t < 6; ++t) {
                auto& ev = tag_entries[t];
                std::sort(ev.begin(), ev.end(),
                          [](const BuildEntry& a, const BuildEntry& b){ return a.key < b.key; });
                auto& pi = sh.pi[t];
                size_t i = 0;
                pi.offsets.push_back(0);
                while (i < ev.size()) {
                    uint64_t k = ev[i].key;
                    pi.keys.push_back(k);
                    uint32_t cnt = 0;
                    while (i < ev.size() && ev[i].key == k) {
                        if (cnt < errcor_.bucket_cap) { pi.ids.push_back(ev[i].id); cnt++; }
                        else { stats.cap_fires++; }
                        i++;
                    }
                    pi.offsets.push_back(static_cast<uint32_t>(pi.ids.size()));
                }
            }
        }
        build_map.clear();  // free memory before allocating ChildMismatch vector

        // Bucket histogram
        std::array<uint64_t, 8> bhist{};
        for (const auto& [ilen, sh] : shards)
            for (const auto& pi : sh.pi)
                for (size_t ki = 0; ki < pi.keys.size(); ++ki) {
                    uint32_t blen = pi.offsets[ki + 1] - pi.offsets[ki];
                    unsigned b = blen == 0 ? 0 : 31 - __builtin_clz(blen);
                    bhist[std::min(b, 7u)]++;
                }
        std::string hstr;
        for (int i = 0; i < 8; ++i) hstr += (i ? "," : "") + std::to_string(bhist[i]);
        log_info("Phase 3 bucket histogram [1,2,3-4,5-8,9-16,17-32,33-64,65+]: " + hstr);

        // ── Phase B1: collect ChildMismatch records ─────────────────────────
        std::vector<ChildMismatch> mismatches;
        mismatches.reserve(std::min(static_cast<size_t>(N) * 2, static_cast<size_t>(1 << 24)));

        // Pre-allocated candidate buffer: 6 queries × bucket_cap (3 canonical + 3 RC pair tags).
        // Sized here so raising bucket_cap (e.g. in capture mode) does not silently overflow.
        std::vector<uint32_t> cand_storage(12u * errcor_.bucket_cap);

        for (uint32_t cid = 0; cid < N; ++cid) {
            if (is_error_[cid]) continue;
            if (!arena_.is_eligible(cid)) continue;

            int L = arena_.length(cid);
            const auto& lay = get_layout(L);
            if (lay.ilen < kMinInteriorLen) { short_interior_skipped++; continue; }

            auto shard_it = shards.find(lay.ilen);
            if (shard_it == shards.end()) continue;

            auto t0 = clk::now();
            // Scratch layout (per child):
            //   [0 .. nbytes)               : ci_parts  (4-part canonical, nb0+nb1+nb2+nb3 bytes)
            //   [nbytes .. nbytes+nf)       : ci_full   (full canonical interior)
            //   [nbytes+nf .. nbytes+2nf)   : crc_full  (RC of ci_full)
            //   [nbytes+2nf .. 2*nbytes+2nf): crc_parts (4-part RC, nb0+nb1+nb2+nb3 bytes)
            //   [2*nbytes+2nf .. 2*nbytes+3nf): pi_buf  (parent, per candidate)
            int nf = (lay.ilen + 3) / 4;
            size_t total_scratch = static_cast<size_t>(2 * lay.nbytes + 3 * nf);
            if (scratch.size() < total_scratch)
                scratch.resize(total_scratch);
            uint8_t* ci_parts  = scratch.data();
            uint8_t* ci_full   = ci_parts  + lay.nbytes;
            uint8_t* crc_full  = ci_full   + nf;
            uint8_t* crc_parts = crc_full  + nf;
            uint8_t* pi_buf    = crc_parts + lay.nbytes;

            const uint8_t* psrc_c = arena_.data(cid);
            int starts[4] = {lay.k5,
                             lay.k5 + lay.s0,
                             lay.k5 + lay.s0 + lay.s1,
                             lay.k5 + lay.s0 + lay.s1 + lay.s2};
            int sizes[4]  = {lay.s0, lay.s1, lay.s2, lay.s3};
            int nbs[4]    = {lay.nb0, lay.nb1, lay.nb2, lay.nb3};
            uint8_t* ci_p[4] = {ci_parts,
                                ci_parts + lay.nb0,
                                ci_parts + lay.nb0 + lay.nb1,
                                ci_parts + lay.nb0 + lay.nb1 + lay.nb2};
            for (int p = 0; p < 4; ++p)
                extract_packed_part(psrc_c, starts[p], sizes[p], ci_p[p]);
            uint64_t h[4], rh[4];
            for (int p = 0; p < 4; ++p) h[p] = XXH3_64bits(ci_p[p], nbs[p]);

            extract_packed_part(psrc_c, lay.k5, lay.ilen, ci_full);

            // RC of full interior for orientation-aware comparison and RC hash keys.
            // A 1-base change can flip which orientation produces the smaller canonical hash,
            // so parent and child may be stored in opposite orientations in the arena.
            compute_interior_rc(ci_full, lay.ilen, crc_full);

            int rc_starts[4] = {0, lay.s0, lay.s0+lay.s1, lay.s0+lay.s1+lay.s2};
            uint8_t* crc_p[4] = {crc_parts,
                                 crc_parts + lay.nb0,
                                 crc_parts + lay.nb0 + lay.nb1,
                                 crc_parts + lay.nb0 + lay.nb1 + lay.nb2};
            for (int p = 0; p < 4; ++p)
                extract_packed_part(crc_full, rc_starts[p], sizes[p], crc_p[p]);
            for (int p = 0; p < 4; ++p) rh[p] = XXH3_64bits(crc_p[p], nbs[p]);
            auto t1 = clk::now();
            stats.decode_hash_child_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

            stats.children_scanned++;

            // Collect unique parent candidates — query both canonical and RC hash keys.
            uint32_t* cand_buf = cand_storage.data();
            uint32_t n_cands = 0;
            auto collect = [&](uint32_t pid) {
                if (n_cands < static_cast<uint32_t>(cand_storage.size())) cand_buf[n_cands++] = pid;
            };

            auto tq0 = clk::now();
            auto& sh = shard_it->second;
            for (int t = 0; t < 6; ++t) {
                sh.pi[t].query(pair_key(h[kPairA[t]],  h[kPairB[t]],  t, lay.ilen), collect);
                sh.pi[t].query(pair_key(rh[kPairA[t]], rh[kPairB[t]], t, lay.ilen), collect);
            }
            auto tq1 = clk::now();
            stats.query_ms += std::chrono::duration<double, std::milli>(tq1 - tq0).count();

            std::sort(cand_buf, cand_buf + n_cands);
            n_cands = static_cast<uint32_t>(
                std::unique(cand_buf, cand_buf + n_cands) - cand_buf);
            stats.total_candidates += n_cands;

            auto tc0 = clk::now();
            for (uint32_t ci = 0; ci < n_cands; ++ci) {
                uint32_t pid = cand_buf[ci];
                if (is_error_[pid]) continue;
                if (!arena_.is_eligible(pid)) continue;
                if (arena_.length(pid) != static_cast<uint16_t>(L)) continue;
                if (pid == cid) continue;
                // Directed edge condition: absorb into higher-count sequences, OR into
                // equal-count sequences using seq_id as tiebreak (lower id = parent).
                // The tiebreak makes the equal-count DAG acyclic (edges always flow toward
                // lower seq_id) while still allowing singleton→singleton chains, which
                // enables H=2 absorption via intermediates: A(1)→B(1)→T(30).
                // Equal-count true SNP variants are still protected by the SNP veto
                // (both halves of the pair share the same mismatch → sig/parent = 100%).
                if (id_count[pid] < id_count[cid]) continue;
                if (id_count[pid] == id_count[cid] && pid >= cid) continue;

                // Extract parent's full interior for comparison.
                extract_packed_part(arena_.data(pid), lay.k5, lay.ilen, pi_buf);

                // Try direct comparison (H=1 — same canonical orientation)
                MismatchInfo mm = packed_find_mismatch(pi_buf, ci_full, 0, lay.ilen, profile_.enabled);
                if (mm.found) {
                    mismatches.push_back({pid, cid, mm.position, mm.base_b, mm.base_a,
                                          0, 0, 0, 1, {}});
                    stats.children_found++;
                    continue;
                }
                // Try RC comparison (H=1)
                mm = packed_find_mismatch(pi_buf, crc_full, 0, lay.ilen, profile_.enabled);
                if (mm.found) {
                    uint16_t cpos = static_cast<uint16_t>(lay.ilen - 1 - mm.position);
                    mismatches.push_back({pid, cid, cpos,
                                          static_cast<uint8_t>(mm.base_b ^ 0x3u),
                                          static_cast<uint8_t>(mm.base_a ^ 0x3u),
                                          0, 0, 0, 1, {}});
                    stats.children_found++;
                    continue;
                }
                // H=2 path: only for reads with count ≤ max_h2_count.
                // Both mismatches must be non-damage transversions (A↔T or C↔G, xr==3).
                if (id_count[cid] <= errcor_.max_h2_count) {
                    MismatchInfo2 mm2 = packed_find_mismatches2(pi_buf, ci_full, 0, lay.ilen);
                    if (mm2.count == 2 &&
                        !is_damage_sub_packed(mm2.base_a[0], mm2.base_b[0], false) &&
                        !is_damage_sub_packed(mm2.base_a[1], mm2.base_b[1], false)) {
                        mismatches.push_back({pid, cid,
                                              mm2.pos[0], mm2.base_b[0], mm2.base_a[0],
                                              mm2.pos[1], mm2.base_b[1], mm2.base_a[1],
                                              2, {}});
                        stats.children_found++;
                        continue;
                    }
                    mm2 = packed_find_mismatches2(pi_buf, crc_full, 0, lay.ilen);
                    if (mm2.count == 2 &&
                        !is_damage_sub_packed(mm2.base_a[0], mm2.base_b[0], false) &&
                        !is_damage_sub_packed(mm2.base_a[1], mm2.base_b[1], false)) {
                        // Convert RC positions/bases to child canonical frame
                        uint16_t p0 = static_cast<uint16_t>(lay.ilen - 1 - mm2.pos[0]);
                        uint16_t p1 = static_cast<uint16_t>(lay.ilen - 1 - mm2.pos[1]);
                        mismatches.push_back({pid, cid,
                                              p0, static_cast<uint8_t>(mm2.base_b[0] ^ 0x3u),
                                                  static_cast<uint8_t>(mm2.base_a[0] ^ 0x3u),
                                              p1, static_cast<uint8_t>(mm2.base_b[1] ^ 0x3u),
                                                  static_cast<uint8_t>(mm2.base_a[1] ^ 0x3u),
                                              2, {}});
                        stats.children_found++;
                    }
                }
            }
            stats.check_ms += std::chrono::duration<double, std::milli>(clk::now() - tc0).count();
        }

        log_info("Phase 3 B1: found " + std::to_string(stats.children_found) +
                 " directed edges from " + std::to_string(stats.children_scanned) + " sequences scanned");

        // ── Phase B2: directed ascent — sort high-count parents first, apply SNP veto ──
        // Processing highest-count parents first enables chain rerouting: when a parent
        // is itself absorbed (early PCR error), its children are redirected to the root's
        // count for the SNP veto.  This handles H=2 transitivity: A→B→C absorbs A even
        // when A is H=2 from C, because A→B uses B's count and B→C uses C's count.
        std::vector<uint32_t> parent_chain(N, UINT32_MAX);
        auto find_root_chain = [&](uint32_t id) -> uint32_t {
            while (parent_chain[id] != UINT32_MAX) id = parent_chain[id];
            return id;
        };

        std::sort(mismatches.begin(), mismatches.end(),
                  [&id_count](const ChildMismatch& a, const ChildMismatch& b) {
                      uint64_t ca = id_count[a.parent_id], cb = id_count[b.parent_id];
                      if (ca != cb) return ca > cb;  // descending count: absorb into highest count first
                      if (a.parent_id != b.parent_id) return a.parent_id < b.parent_id;
                      if (a.mismatch_pos != b.mismatch_pos) return a.mismatch_pos < b.mismatch_pos;
                      return a.alt_base < b.alt_base;
                  });

        size_t i = 0;
        while (i < mismatches.size()) {
            uint32_t pid = mismatches[i].parent_id;
            // Reroute to effective root if this parent was itself absorbed.
            uint32_t eff_pid = pid;
            if (is_error_[pid]) eff_pid = find_root_chain(pid);
            uint64_t parent_count = id_count[eff_pid];
            // SNP veto denominator uses pid's original count: sig was accumulated from
            // pid's children, so the correct ratio is sig/pid_count, not sig/eff_pid_count.
            // Using eff_pid's inflated count deflates the ratio and causes over-absorption
            // of real SNPs that cluster under an absorbed intermediate parent.
            uint64_t veto_count = id_count[pid];

            // Collect all children for this parent
            size_t parent_start = i;
            while (i < mismatches.size() && mismatches[i].parent_id == pid) ++i;
            size_t parent_end = i;

            // Accumulate sig_count per (pos, alt_base) for SNP detection.
            // Key: mismatch_pos * 4 + alt_base. Use unordered_map for O(1) lookup
            // per child instead of the previous O(m) linear scan.
            std::unordered_map<uint32_t, uint64_t> pos_alt_counts;
            pos_alt_counts.reserve(parent_end - parent_start);

            for (size_t j = parent_start; j < parent_end; ++j) {
                uint32_t key = static_cast<uint32_t>(mismatches[j].mismatch_pos) * 4
                               + mismatches[j].alt_base;
                pos_alt_counts[key] += id_count[mismatches[j].child_id];
                if (mismatches[j].hamming == 2) {
                    uint32_t key2 = static_cast<uint32_t>(mismatches[j].mismatch_pos2) * 4
                                    + mismatches[j].alt_base2;
                    pos_alt_counts[key2] += id_count[mismatches[j].child_id];
                }
            }

            // For each child, decide: SNP-protected or absorb
            const int pid_ilen = get_layout(static_cast<int>(arena_.length(pid))).ilen;
            for (size_t j = parent_start; j < parent_end; ++j) {
                const ChildMismatch& cm = mismatches[j];
                uint32_t key = static_cast<uint32_t>(cm.mismatch_pos) * 4 + cm.alt_base;

                // O(1) sig_count lookup
                uint64_t sig = 0;
                auto it = pos_alt_counts.find(key);
                if (it != pos_alt_counts.end()) sig = it->second;

                // H=2 path: both mismatches already verified non-damage in B1.
                // SNP veto: protect if either mismatch position has population support.
                if (cm.hamming == 2) {
                    uint32_t key2 = static_cast<uint32_t>(cm.mismatch_pos2) * 4 + cm.alt_base2;
                    uint64_t sig2 = 0;
                    auto it2 = pos_alt_counts.find(key2);
                    if (it2 != pos_alt_counts.end()) sig2 = it2->second;

                    double eff_snp = errcor_.snp_threshold *
                        (veto_count < errcor_.snp_low_cov_cutoff ? errcor_.snp_low_cov_factor : 1.0);
                    // On strict-dominance edges relax snp_min_count to 1; keep 2 on equal-count
                    // edges where singleton-vs-singleton false protection is a real risk.
                    uint32_t eff_min_count = (parent_count > id_count[cm.child_id]) ? 1u : errcor_.snp_min_count;
                    bool snp_veto_h2 = (sig  >= eff_min_count &&
                                        static_cast<double>(sig)  >= eff_snp * static_cast<double>(veto_count)) ||
                                       (sig2 >= eff_min_count &&
                                        static_cast<double>(sig2) >= eff_snp * static_cast<double>(veto_count));
                    if (snp_veto_h2) {
                        stats.snp_protected++;
                    } else {
                        if (!is_error_[cm.child_id]) {
                            is_error_[cm.child_id] = true;
                            parent_chain[cm.child_id] = eff_pid;
                            acc_count[eff_pid] += acc_count[cm.child_id];
                            acc_count[cm.child_id] = 0;
                            if (seq_entry[cm.child_id] && seq_entry[eff_pid] &&
                                seq_entry[cm.child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                                seq_entry[eff_pid]->record_index = seq_entry[cm.child_id]->record_index;
                                seq_entry[eff_pid]->damage_score = seq_entry[cm.child_id]->damage_score;
                            }
                            stats.absorbed++;
                            errcor_absorbed_++;
                        }
                    }
                    continue;  // skip the H=1 block below
                }

                // Damage-aware bypass: C↔T / G↔A (xr=2) mismatches at interior
                // positions adjacent to the damage zone are residual deamination,
                // not genuine SNPs — bypass the SNP veto for these.
                bool damage_bypass = profile_.enabled &&
                                     ((cm.alt_base ^ cm.parent_base) == 2u) &&
                                     (cm.mismatch_pos < kDamageEdgeMargin ||
                                      cm.mismatch_pos >=
                                          static_cast<uint16_t>(pid_ilen - kDamageEdgeMargin));

                double eff_snp_threshold = errcor_.snp_threshold *
                    (veto_count < errcor_.snp_low_cov_cutoff ? errcor_.snp_low_cov_factor : 1.0);
                // On strict-dominance edges relax snp_min_count to 1.
                // damage_bypass mismatches are still absorbed unless population support is strong
                // enough (sig >= 1 AND ratio >= threshold) indicating a real SNP at the margin.
                uint32_t eff_min_count = (parent_count > id_count[cm.child_id]) ? 1u : errcor_.snp_min_count;
                bool snp_veto = (sig >= (damage_bypass ? 1u : eff_min_count)) &&
                                (static_cast<double>(sig) >= eff_snp_threshold *
                                                              static_cast<double>(veto_count));

                if (snp_veto) {
                    stats.snp_protected++;
                } else {
                    if (!is_error_[cm.child_id]) {
                        is_error_[cm.child_id] = true;
                        parent_chain[cm.child_id] = eff_pid;  // track chain for rerouting descendants
                        acc_count[eff_pid] += acc_count[cm.child_id];  // propagate accumulated mass (incl. child's absorbed descendants)
                        acc_count[cm.child_id] = 0;  // zero to prevent double-counting if re-traversed
                        // Propagate representative: adopt child's record if it has stronger damage signal.
                        if (seq_entry[cm.child_id] && seq_entry[eff_pid] &&
                            seq_entry[cm.child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                            seq_entry[eff_pid]->record_index = seq_entry[cm.child_id]->record_index;
                            seq_entry[eff_pid]->damage_score = seq_entry[cm.child_id]->damage_score;
                        }
                        stats.absorbed++;
                        errcor_absorbed_++;
                    }
                }
            }
        }

        // Write accumulated cluster counts back into the index so pass2() and the
        // cluster TSV reflect absorbed PCR error reads, not just the original
        // representative count.  acc_count is a local — must be written back before
        // phase3_error_correct() returns.
        for (auto& [fp, entry] : index_) {
            if (entry.seq_id < static_cast<uint32_t>(acc_count.size()) &&
                !is_error_[entry.seq_id])
                entry.count = acc_count[entry.seq_id];
        }

        stats.log();
        if (short_interior_skipped > 0) {
            double pct = 100.0 * static_cast<double>(short_interior_skipped) /
                         static_cast<double>(N);
            log_info("Phase 3: " + std::to_string(short_interior_skipped) +
                     " reads skipped (interior < " + std::to_string(kMinInteriorLen) +
                     " bp after damage masking, " +
                     std::to_string(pct).substr(0, 4) + "%)");
            if (pct > 10.0)
                log_warn("Over 10% of reads have interiors too short for Phase 3 error correction "
                         "after damage masking. Consider reducing --mask-threshold or "
                         "adjusting damage zone parameters.");
        }
        log_info("Phase 3 complete: absorbed " + std::to_string(errcor_absorbed_) +
                 " sequences (H=1 directed-ascent + H=2 direct-detection, snp_threshold=" +
                 std::to_string(errcor_.snp_threshold) + ")");
    }

    void pass2(const std::string& in_path,
               const std::string& out_path,
               const std::string& cluster_path) {
        bool compress = (out_path.size() > 3 &&
                         out_path.substr(out_path.size() - 3) == ".gz");
        FastqWriter writer(out_path, compress);

        gzFile cluster_gz = nullptr;
        if (write_clusters_ && !cluster_path.empty()) {
            cluster_gz = gzopen(cluster_path.c_str(), "wb6");
            if (!cluster_gz)
                throw std::runtime_error("Cannot open cluster file: " + cluster_path);
            gzbuffer(cluster_gz, GZBUF_SIZE);
            if (gzprintf(cluster_gz, "hash\tseq_len\tcount\n") < 0)
                throw std::runtime_error("gzprintf failed while writing cluster header");
        }

        struct WriteEntry {
            uint64_t            record_index;
            SequenceFingerprint fingerprint;
            uint64_t            count;
        };
        std::vector<WriteEntry> records_to_write;
        records_to_write.reserve(index_.size());
        for (const auto& [fingerprint, entry] : index_) {
            if (!is_error_.empty() && is_error_[entry.seq_id]) continue;
            records_to_write.push_back({entry.record_index, fingerprint, entry.count});
        }
        std::sort(records_to_write.begin(), records_to_write.end(),
                  [](const WriteEntry& a, const WriteEntry& b) {
                      return a.record_index < b.record_index;
                  });
        const size_t n_to_write = records_to_write.size();
        n_unique_clusters_ = n_to_write;

        // Free index_ and is_error_ — no longer needed; records_to_write has everything.
        // Frees ~10 GB (index_) before the write loop begins.
        { decltype(index_) tmp; std::swap(index_, tmp); }
        { std::vector<bool> tmp; std::swap(is_error_, tmp); }
#ifdef __linux__
        if (mallctl != nullptr) {
            try {
                mallctl("thread.tcache.flush", nullptr, nullptr, nullptr, 0);
                mallctl("arena.0.purge", nullptr, nullptr, nullptr, 0);
            } catch (...) {}
        }
        malloc_trim(0);
#endif

        auto reader = make_fastq_reader(in_path);
        FastqRecord rec;
        uint64_t record_idx = 0;
        size_t written = 0;
        size_t next_write = 0;

        while (reader->read(rec)) {
            if (next_write < n_to_write &&
                records_to_write[next_write].record_index == record_idx) {
                const WriteEntry& entry = records_to_write[next_write];
                writer.write(rec);

                if (cluster_gz) {
                    if (gzprintf(cluster_gz, "%016lx%016lx\t%lu\t%lu\n",
                                 static_cast<unsigned long>(entry.fingerprint.hash_hi),
                                 static_cast<unsigned long>(entry.fingerprint.hash_lo),
                                 static_cast<unsigned long>(rec.seq.size()),
                                 static_cast<unsigned long>(entry.count)) < 0)
                        throw std::runtime_error(
                            "gzprintf failed while writing cluster record");
                }

                written++;
                next_write++;
                if ((written % 100000) == 0) {
                    std::cerr << "\r[Pass 2] " << written << " / " << n_to_write
                              << " unique records written" << std::flush;
                }
            }
            record_idx++;
        }

        if (cluster_gz && gzclose(cluster_gz) != Z_OK)
            throw std::runtime_error("gzclose failed writing cluster file");
        std::cerr << "\r";
        log_info("Pass 2 complete: " + std::to_string(written) + " unique reads written");
    }

    void print_stats() const {
        log_info("=== Final Statistics ===");
        log_info("Total reads processed: " + std::to_string(total_reads_));
        log_info("Unique clusters (dedup): " + std::to_string(n_unique_clusters_));

        if (total_reads_ > 0) {
            double dup_rate = 100.0 * (1.0 - (double)n_unique_clusters_ / total_reads_);
            log_info("Deduplication rate: " + std::to_string(dup_rate) + "%");
        }

        if (errcor_.enabled) {
            size_t output = n_unique_clusters_;
            log_info("PCR error sequences removed: " + std::to_string(errcor_absorbed_) +
                     " (snp_threshold=" + std::to_string(errcor_.snp_threshold) +
                     ", snp_min_count=" + std::to_string(errcor_.snp_min_count) + ")");
            log_info("Output sequences: " + std::to_string(output));
            if (n_unique_clusters_ + errcor_absorbed_ > 0) {
                double ecr = 100.0 * errcor_absorbed_ / (n_unique_clusters_ + errcor_absorbed_);
                log_info("Error correction rate: " + std::to_string(ecr) + "%");
            }
        }

#ifdef __linux__
        if (mallctl != nullptr) {
            try {
                mallctl("thread.tcache.flush", nullptr, nullptr, nullptr, 0);
                mallctl("arena.0.purge", nullptr, nullptr, nullptr, 0);
            } catch (...) {}
        }
        malloc_trim(0);
#endif
    }

    bool use_revcomp_;
    bool write_clusters_;
    DamageProfile profile_;
    ErrCorParams  errcor_;
    bool bucket_cap_explicit_;

    ska::flat_hash_map<SequenceFingerprint, IndexEntry, SequenceFingerprintHash> index_;

    // Phase 3 error correction state (populated only when errcor_.enabled)
    SeqArena          arena_;
    std::vector<bool> is_error_;

    // Reusable scratch buffers for allocation-free damage masking/hashing.
    // Grown lazily to max sequence length seen; never shrunk.
    mutable std::vector<char> scratch1_;
    mutable std::vector<char> scratch2_;
    mutable std::vector<char> rc_buf_;  // unmasked revcomp for canonical arena storage

    uint64_t total_reads_;
    uint64_t errcor_absorbed_;
    size_t   n_unique_clusters_;  // saved before index_ is freed in pass2
};

}  // anonymous namespace

// ============================================================================
// Public entry point
// ============================================================================

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: fqdup " << prog << " [OPTIONS]\n"
        << "\nSingle-file FASTQ deduplication with damage-aware hashing\n"
        << "and PCR error correction. Designed for sorted single-end FASTQ\n"
        << "(e.g. the non-extended output of 'fqdup derep_pairs').\n"
        << "\nInput MUST be sorted by read ID (use 'fqdup sort' first).\n"
        << "\nRequired:\n"
        << "  -i FILE      Input sorted FASTQ\n"
        << "  -o FILE      Output deduplicated FASTQ\n"
        << "\nOptional:\n"
        << "  -c FILE      Write cluster statistics to gzipped TSV\n"
        << "  --no-revcomp Disable reverse-complement matching (default: enabled)\n"
        << "  -h, --help   Show this help\n"
        << "\nPCR error correction (Phase 3 — ON by default):\n"
        << "  --no-error-correct           Disable PCR error duplicate removal\n"
        << "  --error-correct              Explicitly enable (already default)\n"
        << "  --errcor-snp-threshold FLOAT SNP veto: sig/parent_count threshold (default: 0.20)\n"
        << "  --errcor-snp-min-count INT   SNP veto: min absolute sig_count (default: 1)\n"
        << "                               On equal-count edges this is raised to 2 automatically.\n"
        << "  --errcor-bucket-cap INT      Max pair-key bucket size (default: 64)\n"
        << "  --errcor-snp-cutoff INT      Parent count below which SNP threshold is tightened (default: 10)\n"
        << "  --errcor-snp-factor FLOAT    SNP threshold multiplier at low coverage (default: 1.75)\n"
        << "\nAncient DNA damage variant collapsing (OFF by default):\n"
        << "  --collapse-damage        Collapse reads that differ only by terminal\n"
        << "                           deamination into one cluster (Pass 0 estimation)\n"
        << "                           WARNING: distorts per-position damage frequencies.\n"
        << "                           Do NOT use if running metaDMG or mapDamage downstream.\n"
        << "  --library-type TYPE      Library prep: auto (default), ds, ss\n"
        << "                           auto = DART infers from data; ds/ss = override\n"
        << "  --damage-dmax  FLOAT     Set d_max for both 5' and 3' ends manually\n"
        << "  --damage-dmax5 FLOAT     Set d_max for 5' end only\n"
        << "  --damage-dmax3 FLOAT     Set d_max for 3' end only\n"
        << "  --damage-lambda FLOAT    Set lambda (decay) for both ends\n"
        << "  --damage-lambda5 FLOAT   Set lambda for 5' end only\n"
        << "  --damage-lambda3 FLOAT   Set lambda for 3' end only\n"
        << "  --damage-bg FLOAT        Background deamination rate (default: 0.02)\n"
        << "  --mask-threshold FLOAT   Mask positions with excess damage > T (default: 0.05)\n"
        << "  --pcr-cycles INT         Number of PCR cycles (informational; for log only)\n"
        << "  --pcr-efficiency FLOAT   PCR efficiency per cycle, 0-1 (default: 1.0)\n"
        << "  --pcr-error-rate FLOAT   Error rate sub/base/doubling (default: 5.3e-7 = Q5)\n"
        << "\nMemory: ~16 bytes per input read + 2-bit seq arena (~1 byte/4 bp) for error correction.\n";
}

int derep_main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string in_path, out_path, cluster_path;
    bool use_revcomp = true;

    ErrCorParams errcor;
    errcor.enabled = true;
    bool bucket_cap_explicit  = false;   // default on: aDNA primary use case

    bool   damage_auto    = false;  // default off: damage-aware hashing distorts
                                    // downstream damage analysis (e.g. DART) when
                                    // run on fqdup output. Use --damage-auto explicitly.
    dart::SampleDamageProfile::LibraryType forced_library_type =
        dart::SampleDamageProfile::LibraryType::UNKNOWN;  // auto-detect by default
    double damage_dmax5   = -1.0;
    double damage_dmax3   = -1.0;
    double damage_lambda5 = 0.5;
    double damage_lambda3 = 0.5;
    double damage_bg      = 0.02;
    double mask_threshold = 0.05;
    int    pcr_cycles     = 0;
    double pcr_efficiency = 1.0;
    double pcr_phi        = 5.3e-7;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-i" && i + 1 < argc) {
            in_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg == "-c" && i + 1 < argc) {
            cluster_path = argv[++i];
        } else if (arg == "--no-revcomp") {
            use_revcomp = false;
        } else if (arg == "--error-correct") {
            errcor.enabled = true;
        } else if (arg == "--no-error-correct") {
            errcor.enabled = false;
        } else if (arg == "--errcor-snp-threshold" && i + 1 < argc) {
            errcor.snp_threshold = std::stod(argv[++i]);
        } else if (arg == "--errcor-snp-min-count" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 0) { std::cerr << "Error: --errcor-snp-min-count must be >= 0, got " << v << "\n"; return 1; }
            errcor.snp_min_count = static_cast<uint32_t>(v);
        } else if (arg == "--errcor-bucket-cap" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 1) { std::cerr << "Error: --errcor-bucket-cap must be >= 1, got " << v << "\n"; return 1; }
            errcor.bucket_cap = static_cast<uint32_t>(v);
            bucket_cap_explicit = true;
        } else if (arg == "--errcor-snp-cutoff" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 1) { std::cerr << "Error: --errcor-snp-cutoff must be >= 1, got " << v << "\n"; return 1; }
            errcor.snp_low_cov_cutoff = static_cast<uint32_t>(v);
        } else if (arg == "--errcor-snp-factor" && i + 1 < argc) {
            errcor.snp_low_cov_factor = std::stod(argv[++i]);
        } else if (arg == "--library-type" && i + 1 < argc) {
            std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                forced_library_type = dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                forced_library_type = dart::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt == "auto")
                forced_library_type = dart::SampleDamageProfile::LibraryType::UNKNOWN;
            else { std::cerr << "Error: Unknown --library-type: " << lt << " (use auto, ds, ss)\n"; return 1; }
        } else if (arg == "--collapse-damage") {
            damage_auto = true;
        } else if (arg == "--damage-dmax" && i + 1 < argc) {
            damage_dmax5 = damage_dmax3 = std::stod(argv[++i]);
        } else if (arg == "--damage-dmax5" && i + 1 < argc) {
            damage_dmax5 = std::stod(argv[++i]);
        } else if (arg == "--damage-dmax3" && i + 1 < argc) {
            damage_dmax3 = std::stod(argv[++i]);
        } else if (arg == "--damage-lambda" && i + 1 < argc) {
            damage_lambda5 = damage_lambda3 = std::stod(argv[++i]);
        } else if (arg == "--damage-lambda5" && i + 1 < argc) {
            damage_lambda5 = std::stod(argv[++i]);
        } else if (arg == "--damage-lambda3" && i + 1 < argc) {
            damage_lambda3 = std::stod(argv[++i]);
        } else if (arg == "--damage-bg" && i + 1 < argc) {
            damage_bg = std::stod(argv[++i]);
        } else if (arg == "--mask-threshold" && i + 1 < argc) {
            mask_threshold = std::stod(argv[++i]);
        } else if (arg == "--pcr-cycles" && i + 1 < argc) {
            pcr_cycles = std::stoi(argv[++i]);
        } else if (arg == "--pcr-efficiency" && i + 1 < argc) {
            pcr_efficiency = std::stod(argv[++i]);
        } else if (arg == "--pcr-error-rate" && i + 1 < argc) {
            pcr_phi = std::stod(argv[++i]);
        }
    }

    if (in_path.empty() || out_path.empty()) {
        std::cerr << "Error: Missing required arguments (-i and -o)\n\n";
        print_usage(argv[0]);
        return 1;
    }

    // Validate parameter ranges
    if (pcr_efficiency < 0.0 || pcr_efficiency > 1.0) {
        std::cerr << "Error: --pcr-efficiency must be in [0, 1], got "
                  << pcr_efficiency << "\n";
        return 1;
    }
    if (pcr_phi <= 0.0) {
        std::cerr << "Error: --pcr-error-rate must be > 0, got " << pcr_phi << "\n";
        return 1;
    }
    if (pcr_cycles < 0) {
        std::cerr << "Error: --pcr-cycles must be >= 0, got " << pcr_cycles << "\n";
        return 1;
    }
    if (mask_threshold <= 0.0 || mask_threshold >= 1.0) {
        std::cerr << "Error: --mask-threshold must be in (0, 1), got "
                  << mask_threshold << "\n";
        return 1;
    }
    if (errcor.snp_threshold < 0.0 || errcor.snp_threshold > 1.0) {
        std::cerr << "Error: --errcor-snp-threshold must be in [0, 1], got "
                  << errcor.snp_threshold << "\n";
        return 1;
    }
    // damage_dmax sentinel is -1.0 (not set); only validate explicitly-set values (> -0.5)
    if (damage_dmax5 > -0.5 && (damage_dmax5 < 0.0 || damage_dmax5 > 1.0)) {
        std::cerr << "Error: --damage-dmax5 must be in [0, 1], got " << damage_dmax5 << "\n";
        return 1;
    }
    if (damage_dmax3 > -0.5 && (damage_dmax3 < 0.0 || damage_dmax3 > 1.0)) {
        std::cerr << "Error: --damage-dmax3 must be in [0, 1], got " << damage_dmax3 << "\n";
        return 1;
    }
    if (damage_dmax5 > 1.0 || damage_dmax3 > 1.0) {
        std::cerr << "Error: --damage-dmax values must be <= 1.0\n";
        return 1;
    }
    if (damage_lambda5 <= 0.0 || damage_lambda3 <= 0.0) {
        std::cerr << "Error: --damage-lambda values must be > 0\n";
        return 1;
    }

    init_logger("fqdup-derep.log");
    log_info("=== fqdup derep: Single-file two-pass deduplication ===");
    log_info("Input (sorted): " + in_path);
    log_info("Output: " + out_path);
    if (!cluster_path.empty())
        log_info("Cluster output: " + cluster_path);
    log_info("Reverse-complement: " + std::string(use_revcomp ? "enabled" : "disabled"));
    log_info("Damage variant collapsing: " + std::string(damage_auto ? "on (--collapse-damage)" :
             (damage_dmax5 >= 0.0 ? "manual" : "disabled (default)")));
    log_info("PCR error correction: " + std::string(errcor.enabled ? "enabled (default)" : "disabled (--no-error-correct)"));

    // Thread polymerase error rate into Phase 3 adaptive threshold.
    // When --pcr-cycles is given, D_eff is known and pcr_rate is set directly.
    // When --pcr-cycles is not given, D_eff is estimated after Pass 1 from the
    // observed duplication ratio: D_eff = log2(total_reads / unique_reads).
    errcor.pcr_phi = pcr_phi;
    double pcr_total_rate = 0.0;
    if (pcr_cycles > 0) {
        double D_eff = pcr_cycles * std::log2(1.0 + pcr_efficiency);
        pcr_total_rate = pcr_phi * D_eff;
        errcor.pcr_rate = pcr_total_rate;
        log_info("PCR model: phi=" + std::to_string(pcr_phi) +
                 " D_eff=" + std::to_string(D_eff).substr(0, 5) +
                 " rate=" + std::to_string(pcr_total_rate));
    }

    DamageProfile profile;
    profile.mask_threshold = mask_threshold;
    profile.pcr_error_rate = pcr_total_rate;

    try {
        if (damage_auto) {
            profile = estimate_damage(in_path, mask_threshold, forced_library_type);
            profile.pcr_error_rate = pcr_total_rate;
        } else if (damage_dmax5 >= 0.0) {
            profile.d_max_5prime  = damage_dmax5;
            profile.d_max_3prime  = (damage_dmax3 >= 0.0) ? damage_dmax3 : damage_dmax5;
            profile.lambda_5prime = damage_lambda5;
            profile.lambda_3prime = damage_lambda3;
            profile.background    = damage_bg;
            profile.mask_threshold = mask_threshold;
            profile.pcr_error_rate = pcr_total_rate;
            profile.enabled       = (damage_dmax5 > 0.0);
            profile.ss_mode       = (forced_library_type ==
                                     dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            profile.populate_mask_from_model();
            profile.print_info(/*typical_read_length=*/0);
        }
        if (damage_auto && errcor.enabled && !profile.enabled) {
            log_warn("WARNING: --error-correct is ON but no ancient DNA damage detected. "
                     "On modern DNA this may absorb genuine low-frequency variants. "
                     "Use --no-error-correct to disable.");
        }

        DerepEngine engine(use_revcomp, !cluster_path.empty(), profile, errcor,
                           false, bucket_cap_explicit);
        engine.process(in_path, out_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
