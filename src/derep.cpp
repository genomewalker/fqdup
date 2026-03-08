// derep.cpp
// Single-file FASTQ deduplication with damage-aware hashing and PCR error correction.
// Operates on a single sorted FASTQ (typically the non-extended output of derep_pairs).
//
// REQUIRES: Input MUST be sorted by read ID (use 'fqdup sort' first).
//
// Strategy:
//   Pass 0 (optional): Estimate ancient DNA damage parameters from input
//   Pass 1: Stream input, build hash → position index with damage masking
//   Phase 3 (optional): PCR error correction via 3-way pigeonhole Hamming search
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
    uint32_t min_parent_count = 3;    // sequences with count > this are indexed as parents
    double   snp_threshold    = 0.05; // sig_count_weighted/parent_count >= this → SNP protection
    uint32_t snp_min_count    = 2;    // absolute minimum sig_count for SNP veto
    uint32_t bucket_cap       = 64;
    double   pcr_phi          = 5.3e-7;
    double   pcr_rate         = 0.0;  // still estimated for logging but no longer gates absorption
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
// Count-stratified 3-way pigeonhole neighbor finding for edit-distance-1
// detection.  Reads with count <= max_count that have a neighbour with
// count >= ratio × their own count are classified as PCR error duplicates.
//
// 3-way pigeonhole: split interior [k5, L-k3) into three parts p0, p1, p2.
// If Hamming(child, parent) == 1 inside the interior, at least two of the
// three parts match exactly.  We index all three (matched-pair) combinations
// so any 1-error variant is found by exactly one lookup.
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

static inline void split3_lens(int ilen, int& s0, int& s1, int& s2) {
    s0 = ilen / 3;
    s1 = (ilen * 2) / 3 - s0;
    s2 = ilen - s0 - s1;
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

// Returns true if the substitution a→b (or b→a) is damage-consistent
// using 2-bit encoded values.
//
// XOR analysis (A=0,C=1,G=2,T=3):
//   xr=1 (01b): G↔T (10^11) or C↔A (01^00) — 8-oxoG, always protect
//   xr=2 (10b): G↔A (10^00) or C↔T (01^11) — deamination, protect when enabled
//   xr=3 (11b): A↔T or C↔G — transversions, never damage
static bool is_damage_sub_packed(uint8_t a, uint8_t b, bool protect_deamination) {
    uint8_t xr = a ^ b;
    if (xr == 1u) return true;
    if (xr == 2u && protect_deamination) return true;
    return false;
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

// ── ChildMismatch record ────────────────────────────────────────────────────
struct ChildMismatch {
    uint32_t parent_id;
    uint32_t child_id;
    uint16_t mismatch_pos;  // interior-relative position (0-based)
    uint8_t  alt_base;      // child's 2-bit base
    uint8_t  parent_base;   // parent's 2-bit base
};  // 12 bytes
static_assert(sizeof(ChildMismatch) == 12, "ChildMismatch must be 12 bytes");

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

// ============================================================================
// Ancient DNA Damage Model
// ============================================================================

struct DamageProfile {
    // Fitted exponential model parameters (for logging and expected_mismatches only).
    double d_max_5prime   = 0.0;
    double d_max_3prime   = 0.0;
    double lambda_5prime  = 0.5;
    double lambda_3prime  = 0.5;
    double background     = 0.02;
    double mask_threshold = 0.05;
    double pcr_error_rate = 0.0;
    bool   enabled        = false;
    bool   ss_mode        = false;  // true = single-stranded library (C→T at both ends)

    // Empirical per-position mask (primary masking authority).
    //
    // mask_pos[p] = true means position p (measured from either terminus) sits
    // in the damage zone and will be replaced with a neutral byte before hashing.
    //
    // Populated by:
    //   - estimate_damage(): directly from observed T/(T+C) and A/(A+G) excesses
    //     at each position (no model assumption, no OLS fit error).
    //   - populate_mask_from_model(): from the fitted exponential when damage
    //     parameters are supplied manually (--damage-dmax / --damage-lambda).
    //
    // Why empirical masking is more precise than recomputing the exponential at
    // hash time:
    //   The exponential is a convenient approximation, but real damage profiles
    //   can deviate — UDG-treated libraries, mixed fragment lengths, or unusual
    //   protocols can produce non-exponential patterns. By recording which
    //   positions actually exceed the threshold in the observed frequency data,
    //   we mask exactly the positions that are damaged, not the positions the
    //   model predicts to be damaged. This eliminates both OLS fit error and
    //   model misspecification as sources of masking inaccuracy.
    //
    // Symmetry invariant: canonical_hash(seq) == canonical_hash(revcomp(seq)).
    //   apply_damage_mask() uses max(mask_pos[i], mask_pos[j]) semantics
    //   (matching the original max(dmg_5, dmg_3) logic), so the same set of
    //   position indices is masked regardless of which strand is processed.
    static constexpr int MASK_POSITIONS = 15;
    bool mask_pos[MASK_POSITIONS] = {};

    // Populate mask_pos from the fitted exponential model.
    // Used when damage parameters are given manually rather than estimated.
    void populate_mask_from_model() {
        for (int p = 0; p < MASK_POSITIONS; ++p) {
            double exc5 = d_max_5prime * std::exp(-lambda_5prime * p);
            double exc3 = d_max_3prime * std::exp(-lambda_3prime * p);
            mask_pos[p] = std::max(exc5, exc3) > mask_threshold;
        }
    }

    // Helpers retained for expected_mismatches() and informational logging only.
    double dmg_5(int p) const { return d_max_5prime * std::exp(-lambda_5prime * p); }
    double dmg_3(int p) const { return d_max_3prime * std::exp(-lambda_3prime * p); }

    double expected_mismatches(int L) const {
        double e = 0.0;
        for (int p = 0; p < L; ++p) {
            double d5 = dmg_5(p);
            double d3 = dmg_3(L - 1 - p);
            e += 2.0 * d5 * (1.0 - d5);
            e += 2.0 * d3 * (1.0 - d3);
        }
        e += 2.0 * L * pcr_error_rate;
        return e;
    }

    int mismatch_tolerance(int L) const {
        double lam = expected_mismatches(L);
        if (lam < 1e-9) return 0;
        double cumP = 0.0, pois = std::exp(-lam);
        int k = 0;
        while (cumP + pois < 0.99 && k < 100) {
            cumP += pois;
            ++k;
            pois *= lam / k;
        }
        return k;
    }

    void print_info(int typical_read_length) const {
        log_info("--- Damage-Aware Deduplication ---");
        log_info("  5'-end d_max:  " + std::to_string(d_max_5prime));
        log_info("  3'-end d_max:  " + std::to_string(d_max_3prime));
        log_info("  5'-end lambda: " + std::to_string(lambda_5prime));
        log_info("  3'-end lambda: " + std::to_string(lambda_3prime));
        log_info("  Background:    " + std::to_string(background));
        log_info("  Mask threshold:" + std::to_string(mask_threshold));
        if (pcr_error_rate > 0.0)
            log_info("  PCR error rate:" + std::to_string(pcr_error_rate) + " per bp");
        // Log which positions are actually masked.
        int n_masked = 0;
        std::string pos_str;
        for (int p = 0; p < MASK_POSITIONS; ++p) {
            if (mask_pos[p]) {
                if (n_masked > 0) pos_str += ',';
                pos_str += std::to_string(p);
                ++n_masked;
            }
        }
        if (n_masked > 0)
            log_info("  Masked positions: " + pos_str + " (" +
                     std::to_string(n_masked) + " bp each end)");
        else
            log_info("  Masked positions: none");
        if (typical_read_length > 0) {
            double e   = expected_mismatches(typical_read_length);
            int    tol = mismatch_tolerance(typical_read_length);
            log_info("  Expected mismatches (L=" + std::to_string(typical_read_length) +
                     "): " + std::to_string(e) +
                     ", 99th-pct tolerance: " + std::to_string(tol));
        }
    }
};

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
// Damage estimation — DART full pipeline (Pass 0)
// Uses libdart-damage: dart::FrameSelector + SampleDamageProfile
// Channels: A (T/C nucleotide freq) + B (stop codon conversion)
// GC-stratified d_max, Durbin mixture model, joint BIC, pos-0 artifact detection
// ============================================================================

static DamageProfile estimate_damage_impl(FastqReaderBase& reader,
                                          const std::string& path,
                                          double mask_threshold,
                                          dart::SampleDamageProfile::LibraryType forced_lib =
                                              dart::SampleDamageProfile::LibraryType::UNKNOWN) {
    FastqRecord rec;
    dart::SampleDamageProfile dart_profile;
    dart_profile.forced_library_type = forced_lib;

    int      reads_scanned = 0;
    int      typical_len   = 0;
    uint64_t record_pos    = 0;

    while (reader.read(rec)) {
        int L = static_cast<int>(rec.seq.size());
        record_pos++;
        if (L < 30) continue;
        if (L > typical_len) typical_len = L;
        dart::FrameSelector::update_sample_profile(dart_profile, rec.seq);
        reads_scanned++;
    }

    dart::FrameSelector::finalize_sample_profile(dart_profile);

    // Use DART's calibrated per-end values (damage_rate[0] = raw pos0 excess) —
    // these match DART's JSON output and are on the same scale as metaDMG's Δ.
    double d_max_5   = dart_profile.d_max_5prime;
    double d_max_3   = dart_profile.d_max_3prime;
    double lambda_5  = dart_profile.lambda_5prime;
    double lambda_3  = dart_profile.lambda_3prime;
    double bg_5      = dart_profile.fit_baseline_5prime;
    double bg_3      = dart_profile.fit_baseline_3prime;
    double background = (bg_5 + bg_3) / 2.0;

    // d_max_combined is the primary damage signal — mixture model d_reference
    // (metaDMG proxy) when converged; falls back to GC-weighted or joint model.
    double d_max_combined = dart_profile.d_max_combined;

    const bool is_ss = (dart_profile.library_type ==
                        dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED);

    DamageProfile profile;
    profile.d_max_5prime   = d_max_5;
    profile.d_max_3prime   = d_max_3;
    profile.lambda_5prime  = lambda_5;
    profile.lambda_3prime  = lambda_3;
    profile.background     = background;
    profile.mask_threshold = mask_threshold;
    profile.ss_mode        = is_ss;
    // Enable if combined or per-end d_max exceeds threshold, OR if DART validates
    // damage even when d_max_5prime/3prime = 0 (e.g. pos0-artifact samples).
    profile.enabled        = (d_max_combined > 0.02 || d_max_5 > 0.02 || d_max_3 > 0.02
                              || dart_profile.damage_validated);

    // Empirical per-position mask from DART's normalized frequencies.
    // DS: max(5' C→T excess, 3' G→A excess) > threshold.
    // SS: max(5' C→T excess, 3' T/(T+C) excess) > threshold — both ends have C→T.
    constexpr double MIN_COV = 100.0;
    double bg_tc = dart_profile.baseline_t_freq /
                   (dart_profile.baseline_t_freq + dart_profile.baseline_c_freq + 1e-9);
    for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dart_profile.tc_total_5prime[p] >= MIN_COV)
            excess_5 = dart_profile.t_freq_5prime[p] - bg_5;
        if (is_ss) {
            // SS: 3' damage signal is T/(T+C) excess (same type as 5')
            if (dart_profile.tc_total_3prime[p] >= MIN_COV) {
                double tc3 = dart_profile.tc_total_3prime[p];
                double ct3_freq = (tc3 > 0)
                    ? dart_profile.t_freq_3prime[p] / tc3
                    : bg_tc;
                excess_3 = ct3_freq - bg_tc;
            }
        } else {
            // DS: 3' damage signal is A/(A+G) excess
            if (dart_profile.ag_total_3prime[p] >= MIN_COV)
                excess_3 = dart_profile.a_freq_3prime[p] - bg_3;
        }
        profile.mask_pos[p] = (excess_5 > mask_threshold) || (excess_3 > mask_threshold);
    }

    log_info("Pass 0: DART damage estimation — " +
             std::to_string(reads_scanned) + " reads in " + path);
    log_info("  Library type: " + std::string(dart_profile.library_type_str()) +
             (dart_profile.library_type_auto_detected ? " (auto-detected)" : " (forced)"));
    log_info("  5'-end: d_max=" + std::to_string(d_max_5) +
             " lambda=" + std::to_string(lambda_5) +
             " bg=" + std::to_string(bg_5));
    log_info("  3'-end: d_max=" + std::to_string(d_max_3) +
             " lambda=" + std::to_string(lambda_3) +
             " bg=" + std::to_string(bg_3));
    log_info("  d_max_combined=" + std::to_string(d_max_combined) +
             " (source=" + dart_profile.d_max_source_str() + ")" +
             (dart_profile.mixture_converged ? " [mixture]" : "") +
             (dart_profile.damage_validated  ? " [validated]" : "") +
             (dart_profile.damage_artifact   ? " [ARTIFACT]" : ""));
    if (profile.enabled && typical_len > 0) {
        profile.print_info(typical_len);
    } else {
        log_info("  Damage below threshold — standard exact hashing will be used");
    }
    return profile;
}

static DamageProfile estimate_damage(const std::string& path, double mask_threshold,
                                     dart::SampleDamageProfile::LibraryType forced_lib =
                                         dart::SampleDamageProfile::LibraryType::UNKNOWN) {
    auto reader = make_fastq_reader(path);
    return estimate_damage_impl(*reader, path, mask_threshold, forced_lib);
}

// ============================================================================
// DerepEngine — single-file two-pass deduplication
// ============================================================================

class DerepEngine {
public:
    DerepEngine(bool use_revcomp, bool write_clusters,
                const DamageProfile& profile = DamageProfile{},
                const ErrCorParams& errcor = ErrCorParams{})
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          profile_(profile), errcor_(errcor),
          total_reads_(0), errcor_absorbed_(0) {}

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
            // Auto-estimate D_eff from duplication ratio when --pcr-cycles not given.
            // Under PCR kinetics: mean copies per molecule = (1+E)^n, so
            //   D_eff = log2((1+E)^n) = log2(total_reads / unique_reads).
            // This is exact for uniform amplification; a lower bound when starting
            // copy number varies (conservative: under-estimates errors, avoids
            // false absorptions).
            if (errcor_.pcr_rate == 0.0 && errcor_.pcr_phi > 0.0 &&
                total_reads_ > index_.size() && index_.size() > 0) {
                double dr     = static_cast<double>(total_reads_) /
                                static_cast<double>(index_.size());
                double d_eff  = std::log2(dr);
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
        for (const auto& [fp, entry] : index_)
            id_count[entry.seq_id] = entry.count;

        // Scratch for packed interior hashing
        std::vector<uint8_t> scratch;
        scratch.reserve(65535);

        struct InteriorLayout {
            int  k5 = 0, ilen = 0, s0 = 0, s1 = 0, s2 = 0;
            int  nb0 = 0, nb1 = 0, nb2 = 0, nbytes = 0;
            bool ready = false;
        };
        std::vector<InteriorLayout> layout_cache(65536);
        auto get_layout = [&](int L) -> const InteriorLayout& {
            auto& lay = layout_cache[static_cast<uint16_t>(L)];
            if (!lay.ready) {
                auto [k5, k3] = damage_zone_bounds(L, profile_);
                lay.k5   = k5;
                lay.ilen = std::max(0, L - k5 - k3);
                if (lay.ilen >= 3) {
                    split3_lens(lay.ilen, lay.s0, lay.s1, lay.s2);
                    lay.nb0    = (lay.s0 + 3) / 4;
                    lay.nb1    = (lay.s1 + 3) / 4;
                    lay.nb2    = (lay.s2 + 3) / 4;
                    lay.nbytes = lay.nb0 + lay.nb1 + lay.nb2;
                }
                lay.ready = true;
            }
            return lay;
        };

        Phase3Stats stats;
        using clk = std::chrono::steady_clock;

        // ── Phase A: index parents (count > min_parent_count) ──────────────
        struct BuildEntry { uint64_t key; uint32_t id; };
        ska::flat_hash_map<int, std::array<std::vector<BuildEntry>, 3>> build_map;
        build_map.reserve(64);
        struct LenShard { std::array<FlatPairIndex, 3> pi; };

        uint64_t n_parents = 0;
        for (uint32_t id = 0; id < N; ++id) {
            if (id_count[id] <= errcor_.min_parent_count) continue;
            if (!arena_.is_eligible(id)) continue;
            int L = arena_.length(id);
            const auto& lay = get_layout(L);
            if (lay.ilen < 3) continue;

            auto t0 = clk::now();
            if (scratch.size() < static_cast<size_t>(lay.nbytes))
                scratch.resize(lay.nbytes);
            const uint8_t* psrc = arena_.data(id);
            extract_packed_part(psrc, lay.k5,                    lay.s0, scratch.data());
            extract_packed_part(psrc, lay.k5 + lay.s0,           lay.s1, scratch.data() + lay.nb0);
            extract_packed_part(psrc, lay.k5 + lay.s0 + lay.s1, lay.s2, scratch.data() + lay.nb0 + lay.nb1);
            uint64_t h0 = XXH3_64bits(scratch.data(),                     lay.nb0);
            uint64_t h1 = XXH3_64bits(scratch.data() + lay.nb0,           lay.nb1);
            uint64_t h2 = XXH3_64bits(scratch.data() + lay.nb0 + lay.nb1, lay.nb2);
            auto t1 = clk::now();

            auto& entries = build_map[lay.ilen];
            entries[0].push_back({pair_key(h0, h1, 0, lay.ilen), id});
            entries[1].push_back({pair_key(h0, h2, 1, lay.ilen), id});
            entries[2].push_back({pair_key(h1, h2, 2, lay.ilen), id});
            auto t2 = clk::now();

            stats.decode_hash_parent_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
            stats.insert_ms             += std::chrono::duration<double, std::milli>(t2 - t1).count();
            n_parents++;
        }
        stats.parents_indexed = n_parents;
        log_info("Phase 3: indexed " + std::to_string(n_parents) +
                 " parent sequences (count>" + std::to_string(errcor_.min_parent_count) + ")");

        // Build CSR per shard per tag
        ska::flat_hash_map<int, LenShard> shards;
        shards.reserve(build_map.size());
        for (auto& [ilen, tag_entries] : build_map) {
            auto& sh = shards[ilen];
            for (int t = 0; t < 3; ++t) {
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

        for (uint32_t cid = 0; cid < N; ++cid) {
            if (id_count[cid] > errcor_.min_parent_count) continue;  // not a child
            if (is_error_[cid]) continue;
            if (!arena_.is_eligible(cid)) continue;

            int L = arena_.length(cid);
            const auto& lay = get_layout(L);
            if (lay.ilen < 3) continue;

            auto shard_it = shards.find(lay.ilen);
            if (shard_it == shards.end()) continue;

            auto t0 = clk::now();
            if (scratch.size() < static_cast<size_t>(lay.nbytes))
                scratch.resize(lay.nbytes);
            const uint8_t* psrc_c = arena_.data(cid);
            extract_packed_part(psrc_c, lay.k5,                    lay.s0, scratch.data());
            extract_packed_part(psrc_c, lay.k5 + lay.s0,           lay.s1, scratch.data() + lay.nb0);
            extract_packed_part(psrc_c, lay.k5 + lay.s0 + lay.s1, lay.s2, scratch.data() + lay.nb0 + lay.nb1);
            uint64_t h0 = XXH3_64bits(scratch.data(),                     lay.nb0);
            uint64_t h1 = XXH3_64bits(scratch.data() + lay.nb0,           lay.nb1);
            uint64_t h2 = XXH3_64bits(scratch.data() + lay.nb0 + lay.nb1, lay.nb2);
            auto t1 = clk::now();
            stats.decode_hash_child_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

            stats.children_scanned++;

            // Collect unique parent candidates
            uint32_t cand_buf[192];
            uint32_t n_cands = 0;
            auto collect = [&](uint32_t pid) {
                if (n_cands < 192) cand_buf[n_cands++] = pid;
            };

            auto tq0 = clk::now();
            auto& sh = shard_it->second;
            sh.pi[0].query(pair_key(h0, h1, 0, lay.ilen), collect);
            sh.pi[1].query(pair_key(h0, h2, 1, lay.ilen), collect);
            sh.pi[2].query(pair_key(h1, h2, 2, lay.ilen), collect);
            auto tq1 = clk::now();
            stats.query_ms += std::chrono::duration<double, std::milli>(tq1 - tq0).count();

            std::sort(cand_buf, cand_buf + n_cands);
            n_cands = static_cast<uint32_t>(
                std::unique(cand_buf, cand_buf + n_cands) - cand_buf);
            stats.total_candidates += n_cands;

            const uint8_t* child_packed = arena_.data(cid);
            auto tc0 = clk::now();
            for (uint32_t ci = 0; ci < n_cands; ++ci) {
                uint32_t pid = cand_buf[ci];
                if (is_error_[pid]) continue;
                if (!arena_.is_eligible(pid)) continue;
                if (arena_.length(pid) != static_cast<uint16_t>(L)) continue;

                MismatchInfo mm = packed_find_mismatch(arena_.data(pid), child_packed,
                                                       lay.k5, lay.ilen, profile_.enabled);
                if (!mm.found) continue;

                mismatches.push_back({pid, cid, mm.position, mm.base_b, mm.base_a});
                stats.children_found++;
            }
            stats.check_ms += std::chrono::duration<double, std::milli>(clk::now() - tc0).count();
        }

        log_info("Phase 3 B1: found " + std::to_string(stats.children_found) +
                 " child-parent pairs from " + std::to_string(stats.children_scanned) + " children scanned");

        // ── Phase B2: sort and apply SNP veto ──────────────────────────────
        // Sort by parent_id (primary), then mismatch_pos, then alt_base for grouping
        std::sort(mismatches.begin(), mismatches.end(),
                  [](const ChildMismatch& a, const ChildMismatch& b) {
                      if (a.parent_id != b.parent_id) return a.parent_id < b.parent_id;
                      if (a.mismatch_pos != b.mismatch_pos) return a.mismatch_pos < b.mismatch_pos;
                      return a.alt_base < b.alt_base;
                  });

        size_t i = 0;
        while (i < mismatches.size()) {
            uint32_t pid = mismatches[i].parent_id;
            uint64_t parent_count = id_count[pid];

            // Collect all children for this parent
            size_t parent_start = i;
            while (i < mismatches.size() && mismatches[i].parent_id == pid) ++i;
            size_t parent_end = i;

            // Accumulate sig_count per (pos, alt_base) for SNP detection.
            // Key: mismatch_pos * 4 + alt_base. Entries already sorted by (pos, alt).
            std::vector<std::pair<uint32_t, uint64_t>> pos_alt_counts;
            pos_alt_counts.reserve(parent_end - parent_start);

            {
                size_t j = parent_start;
                while (j < parent_end) {
                    uint16_t cur_pos = mismatches[j].mismatch_pos;
                    uint8_t  cur_alt = mismatches[j].alt_base;
                    uint32_t key     = static_cast<uint32_t>(cur_pos) * 4 + cur_alt;
                    uint64_t weight  = 0;
                    while (j < parent_end && mismatches[j].mismatch_pos == cur_pos &&
                           mismatches[j].alt_base == cur_alt) {
                        weight += id_count[mismatches[j].child_id];
                        ++j;
                    }
                    pos_alt_counts.push_back({key, weight});
                }
            }

            // For each child, decide: SNP-protected or absorb
            for (size_t j = parent_start; j < parent_end; ++j) {
                const ChildMismatch& cm = mismatches[j];
                uint32_t key = static_cast<uint32_t>(cm.mismatch_pos) * 4 + cm.alt_base;

                // Find sig_count for this (pos, alt)
                uint64_t sig = 0;
                for (const auto& [k, w] : pos_alt_counts) {
                    if (k == key) { sig = w; break; }
                }

                bool snp_veto = (sig >= errcor_.snp_min_count) &&
                                (static_cast<double>(sig) >= errcor_.snp_threshold *
                                                              static_cast<double>(parent_count));

                if (snp_veto) {
                    stats.snp_protected++;
                } else {
                    if (!is_error_[cm.child_id]) {
                        is_error_[cm.child_id] = true;
                        stats.absorbed++;
                        errcor_absorbed_++;
                    }
                }
            }
        }

        stats.log();
        log_info("Phase 3 complete: absorbed " + std::to_string(errcor_absorbed_) +
                 " PCR error sequences (min_parent=" + std::to_string(errcor_.min_parent_count) +
                 ", snp_threshold=" + std::to_string(errcor_.snp_threshold) +
                 ", snp_min_count=" + std::to_string(errcor_.snp_min_count) + ")");
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
            if (cluster_gz) {
                gzbuffer(cluster_gz, GZBUF_SIZE);
                if (gzprintf(cluster_gz, "hash\tseq_len\tcount\n") < 0)
                    throw std::runtime_error(
                        "gzprintf failed while writing cluster header");
            }
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

        if (cluster_gz) gzclose(cluster_gz);
        std::cerr << "\r";
        log_info("Pass 2 complete: " + std::to_string(written) + " unique reads written");
    }

    void print_stats() const {
        log_info("=== Final Statistics ===");
        log_info("Total reads processed: " + std::to_string(total_reads_));
        log_info("Unique clusters (dedup): " + std::to_string(index_.size()));

        if (total_reads_ > 0) {
            double dup_rate = 100.0 * (1.0 - (double)index_.size() / total_reads_);
            log_info("Deduplication rate: " + std::to_string(dup_rate) + "%");
        }

        if (errcor_.enabled) {
            size_t output = index_.size() - errcor_absorbed_;
            log_info("PCR error sequences removed: " + std::to_string(errcor_absorbed_) +
                     " (min_parent=" + std::to_string(errcor_.min_parent_count) +
                     ", snp_threshold=" + std::to_string(errcor_.snp_threshold) + ")");
            log_info("Output sequences: " + std::to_string(output));
            if (index_.size() > 0) {
                double ecr = 100.0 * errcor_absorbed_ / index_.size();
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
        << "  --errcor-min-parent INT      Min count to be indexed as parent (default: 3)\n"
        << "  --errcor-snp-threshold FLOAT SNP veto: sig/parent_count threshold (default: 0.05)\n"
        << "  --errcor-snp-min-count INT   SNP veto: minimum absolute sig_count (default: 2)\n"
        << "  --errcor-bucket-cap INT      Max pair-key bucket size (default: 64)\n"
        << "\nAncient DNA damage-aware deduplication (OFF by default):\n"
        << "  --damage-auto            Enable damage estimation and masking (Pass 0)\n"
        << "  --no-damage              Explicitly disable (already default)\n"
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
    errcor.enabled = true;   // default on: aDNA primary use case

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
        } else if (arg == "--errcor-min-parent" && i + 1 < argc) {
            errcor.min_parent_count = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--errcor-snp-threshold" && i + 1 < argc) {
            errcor.snp_threshold = std::stod(argv[++i]);
        } else if (arg == "--errcor-snp-min-count" && i + 1 < argc) {
            errcor.snp_min_count = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--errcor-bucket-cap" && i + 1 < argc) {
            errcor.bucket_cap = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--library-type" && i + 1 < argc) {
            std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                forced_library_type = dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                forced_library_type = dart::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt == "auto")
                forced_library_type = dart::SampleDamageProfile::LibraryType::UNKNOWN;
            else { std::cerr << "Error: Unknown --library-type: " << lt << " (use auto, ds, ss)\n"; return 1; }
        } else if (arg == "--damage-auto") {
            damage_auto = true;
        } else if (arg == "--no-damage") {
            damage_auto = false;
            damage_dmax5 = damage_dmax3 = -1.0;  // clear any earlier manual values
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
    log_info("Damage-aware mode: " + std::string(damage_auto ? "auto (--damage-auto)" :
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

        DerepEngine engine(use_revcomp, !cluster_path.empty(), profile, errcor);
        engine.process(in_path, out_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
