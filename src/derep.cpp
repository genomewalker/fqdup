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
    bool     enabled    = false;
    double   ratio      = 50.0;   // parent must have count >= ratio × child count
    uint32_t max_count  = 5;      // static child ceiling (hard cap; also minimum when adaptive)
    uint32_t bucket_cap = 64;     // max pair-key bucket size (limits low-complexity bloat)
    // Adaptive child ceiling from PCR model (Potapov & Ong 2017):
    //   pcr_rate = phi × D_eff  (sub/base/full-PCR-run)
    // Per parent P: effective_max = max(max_count, ceil(P × pcr_rate / 3)).
    // The /3 is for the specific-substitution spectrum under the equal-rates assumption.
    //
    // pcr_phi: polymerase error rate (sub/base/doubling). Q5=5.3e-7 (default).
    //   Always set from --pcr-error-rate (or default). Used to auto-estimate pcr_rate
    //   from the observed duplication ratio when --pcr-cycles is not given.
    //
    // pcr_rate: phi × D_eff. Set explicitly when pcr_cycles>0; otherwise estimated
    //   after Pass 1 from D_eff = log2(total_reads / unique_reads). 0 = not yet set.
    double   pcr_phi    = 5.3e-7;  // polymerase error rate (sub/base/doubling)
    double   pcr_rate   = 0.0;     // phi × D_eff; 0 = use static max_count only
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
static uint8_t enc2bit(uint8_t c) {
    c &= 0xDFu;
    if (c == 'A') return 0u;
    if (c == 'C') return 1u;
    if (c == 'G') return 2u;
    if (c == 'T') return 3u;
    return 0xFFu;
}

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
        for (int i = 0; i < (int)L; ++i) {
            uint8_t b = enc2bit(static_cast<uint8_t>(seq[i]));
            if (b == 0xFFu) { ok = false; b = 0; }
            dst[i >> 2] |= b << (6 - 2 * (i & 3));
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
        for (int i = 0; i < count; ++i) {
            int pos = start + i;
            out[i] = kDec2bit[(p[pos >> 2] >> (6 - 2 * (pos & 3))) & 0x3u];
        }
    }
};

struct SpillNode { uint32_t id; uint32_t next; };

struct BucketVal {
    uint32_t first      = 0;
    uint32_t spill_head = 0;
    uint32_t len        = 0;
};

struct PairIndex {
    ska::flat_hash_map<uint64_t, BucketVal> map;
    std::vector<SpillNode> pool;

    PairIndex() { pool.push_back({0, 0}); }  // index 0 = sentinel

    void insert(uint64_t key, uint32_t id, uint32_t cap) {
        auto& bv = map[key];
        if (bv.len >= cap) return;
        if (bv.len == 0) {
            bv.first = id;
        } else {
            if (pool.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
                throw std::runtime_error("PairIndex pool overflow (> 4 G spill nodes)");
            uint32_t node = static_cast<uint32_t>(pool.size());
            pool.push_back({id, bv.spill_head});
            bv.spill_head = node;
        }
        bv.len++;
    }

    template <typename F>
    void query(uint64_t key, F&& fn) const {
        auto it = map.find(key);
        if (it == map.end() || it->second.len == 0) return;
        fn(it->second.first);
        uint32_t idx = it->second.spill_head;
        while (idx != 0) { fn(pool[idx].id); idx = pool[idx].next; }
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
                mismatches += __builtin_popcount(any);
                if (mismatches > 1) return false;
                // Find the single mismatching 2-bit lane
                for (int lane = 0; lane < 4; ++lane) {
                    uint8_t sh = 6 - 2 * lane;
                    uint8_t ba = (pa[byte_i] >> sh) & 0x3u;
                    uint8_t bb = (pb[byte_i] >> sh) & 0x3u;
                    if (ba != bb) { mm_a = ba; mm_b = bb; break; }
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
        if (in_5zone && (cu == 'C' || cu == 'T')) {
            scratch[i] = '\x01';
        } else if (in_3zone && (cu == 'G' || cu == 'A')) {
            scratch[i] = '\x02';
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
        if (in_5zone && (cu == 'C' || cu == 'T')) {
            scratch2[i] = '\x01';
        } else if (in_3zone && (cu == 'G' || cu == 'A')) {
            scratch2[i] = '\x02';
        }
    }
    XXH128_hash_t h2 = XXH3_128bits(scratch2, L);
    // Canonical = lexicographically smaller of the two 128-bit hashes
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
                                          double mask_threshold) {
    FastqRecord rec;
    dart::SampleDamageProfile dart_profile;

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

    DamageProfile profile;
    profile.d_max_5prime   = d_max_5;
    profile.d_max_3prime   = d_max_3;
    profile.lambda_5prime  = lambda_5;
    profile.lambda_3prime  = lambda_3;
    profile.background     = background;
    profile.mask_threshold = mask_threshold;
    // Enable if combined or per-end d_max exceeds threshold, OR if DART validates
    // damage even when d_max_5prime/3prime = 0 (e.g. pos0-artifact samples).
    profile.enabled        = (d_max_combined > 0.02 || d_max_5 > 0.02 || d_max_3 > 0.02
                              || dart_profile.damage_validated);

    // Empirical per-position mask from DART's normalized frequencies.
    // Uses max(5'-excess, 3'-excess) > mask_threshold; coverage-gated by tc_total.
    constexpr double MIN_COV = 100.0;
    for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dart_profile.tc_total_5prime[p] >= MIN_COV)
            excess_5 = dart_profile.t_freq_5prime[p] - bg_5;
        if (dart_profile.ag_total_3prime[p] >= MIN_COV)
            excess_3 = dart_profile.a_freq_3prime[p] - bg_3;
        profile.mask_pos[p] = (excess_5 > mask_threshold) || (excess_3 > mask_threshold);
    }

    log_info("Pass 0: DART damage estimation — " +
             std::to_string(reads_scanned) + " reads in " + path);
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

static DamageProfile estimate_damage(const std::string& path, double mask_threshold) {
    auto reader = make_fastq_reader(path);
    return estimate_damage_impl(*reader, path, mask_threshold);
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

            log_info("Phase 3: PCR error correction");
            phase3_error_correct();
        }

        log_info("Pass 2: Write unique records");

        pass2(in_path, out_path, cluster_path);

        print_stats();
    }

private:
    XXH128_hash_t compute_hash(const std::string& seq) const {
        if (profile_.enabled) {
            if (scratch1_.size() < seq.size()) scratch1_.resize(seq.size());
            if (scratch2_.size() < seq.size()) scratch2_.resize(seq.size());
            return damage_canonical_hash(seq, profile_, use_revcomp_,
                                         scratch1_.data(), scratch2_.data());
        }
        return canonical_hash(seq, use_revcomp_);
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
            XXH128_hash_t h = compute_hash(rec.seq);
            SequenceFingerprint fp(h, rec.seq.size());

            auto it = index_.find(fp);
            if (it == index_.end()) {
                IndexEntry entry(record_idx);
                if (profile_.enabled)
                    entry.damage_score = compute_damage_score(rec.seq, profile_);
                if (errcor_.enabled)
                    entry.seq_id = arena_.append(rec.seq);
                index_.emplace(fp, entry);
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

    // Phase 3: absorb PCR error duplicates.
    // A cluster C (count <= max_count) is absorbed into parent P
    // (count >= ratio*C.count) if their sequences differ by at most 1
    // substitution OUTSIDE the damage zone.
    void phase3_error_correct() {
        if (arena_.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("Phase 3: arena size exceeds uint32_t range");
        const uint32_t N = static_cast<uint32_t>(arena_.size());
        if (N == 0) return;

        is_error_.assign(N, false);

        std::vector<uint64_t> id_count(N, 0);
        for (const auto& [fp, entry] : index_)
            id_count[entry.seq_id] = entry.count;

        // Adaptive child ceiling: when PCR model is active, raise the pre-filter
        // so that sequences whose count is explained by a large parent's error rate
        // are also evaluated as potential children.
        // child_prefilter = max(static max_count, ceil(max_parent_count × phi×D / 3))
        uint64_t child_prefilter = errcor_.max_count;
        if (errcor_.pcr_rate > 0.0) {
            uint64_t max_parent = 0;
            for (uint32_t id = 0; id < N; ++id)
                if (id_count[id] > errcor_.max_count)
                    max_parent = std::max(max_parent, id_count[id]);
            if (max_parent > 0) {
                uint64_t adaptive = static_cast<uint64_t>(
                    std::ceil(static_cast<double>(max_parent) * errcor_.pcr_rate / 3.0));
                child_prefilter = std::max(child_prefilter, adaptive);
                if (adaptive > errcor_.max_count)
                    log_info("Phase 3: adaptive child ceiling = " +
                             std::to_string(child_prefilter) +
                             " (PCR model: max_parent=" + std::to_string(max_parent) +
                             " × " + std::to_string(errcor_.pcr_rate) + "/3)");
            }
        }

        struct LenShard { std::array<PairIndex, 3> pi; };
        ska::flat_hash_map<int, LenShard> shards;
        shards.reserve(64);

        // Scratch buffer for decoding packed interior regions before hashing.
        std::vector<uint8_t> scratch;
        scratch.reserve(65535);

        uint64_t n_parents = 0;
        for (uint32_t id = 0; id < N; ++id) {
            if (id_count[id] <= errcor_.max_count) continue;
            if (!arena_.is_eligible(id)) continue;
            int L = arena_.length(id);
            auto [k5, k3] = damage_zone_bounds(L, profile_);
            int ilen = std::max(0, L - k5 - k3);
            if (ilen < 3) continue;

            // Decode interior to ASCII scratch for hashing; packed data for Hamming.
            scratch.resize(ilen);
            arena_.decode_range(id, k5, ilen, scratch.data());
            int s0, s1, s2; split3_lens(ilen, s0, s1, s2);
            uint64_t h0 = XXH3_64bits(scratch.data(),          s0);
            uint64_t h1 = XXH3_64bits(scratch.data() + s0,     s1);
            uint64_t h2 = XXH3_64bits(scratch.data() + s0 + s1, s2);

            auto& sh = shards[ilen];
            sh.pi[0].insert(pair_key(h0, h1, 0, ilen), id, errcor_.bucket_cap);
            sh.pi[1].insert(pair_key(h0, h2, 1, ilen), id, errcor_.bucket_cap);
            sh.pi[2].insert(pair_key(h1, h2, 2, ilen), id, errcor_.bucket_cap);
            n_parents++;
        }
        log_info("Phase 3: indexed " + std::to_string(n_parents) +
                 " parent sequences (count>" + std::to_string(errcor_.max_count) + ")");

        std::vector<uint64_t> seen_epoch(N, 0);
        uint64_t epoch = 0;

        for (uint32_t cid = 0; cid < N; ++cid) {
            if (id_count[cid] > child_prefilter) continue;
            if (is_error_[cid]) continue;
            if (!arena_.is_eligible(cid)) continue;

            int L = arena_.length(cid);
            auto [k5, k3] = damage_zone_bounds(L, profile_);
            int ilen = std::max(0, L - k5 - k3);
            if (ilen < 3) continue;

            auto shard_it = shards.find(ilen);
            if (shard_it == shards.end()) continue;

            scratch.resize(ilen);
            arena_.decode_range(cid, k5, ilen, scratch.data());
            int s0, s1, s2; split3_lens(ilen, s0, s1, s2);
            uint64_t h0 = XXH3_64bits(scratch.data(),           s0);
            uint64_t h1 = XXH3_64bits(scratch.data() + s0,      s1);
            uint64_t h2 = XXH3_64bits(scratch.data() + s0 + s1, s2);

            epoch++;
            bool absorbed = false;
            const uint8_t* child_packed = arena_.data(cid);

            auto check = [&](uint32_t pid) {
                if (absorbed) return;
                if (seen_epoch[pid] == epoch) return;
                seen_epoch[pid] = epoch;
                if (is_error_[pid]) return;
                if (!arena_.is_eligible(pid)) return;
                if (arena_.length(pid) != static_cast<uint16_t>(L)) return;
                if (static_cast<double>(id_count[pid]) <
                    errcor_.ratio * id_count[cid]) return;
                // Adaptive child ceiling: per-parent effective max count.
                // With PCR model active: ceil(parent_count × phi×D / 3).
                uint64_t eff_max = errcor_.max_count;
                if (errcor_.pcr_rate > 0.0) {
                    uint64_t adaptive = static_cast<uint64_t>(
                        std::ceil(static_cast<double>(id_count[pid]) *
                                  errcor_.pcr_rate / 3.0));
                    eff_max = std::max(eff_max, adaptive);
                }
                if (id_count[cid] > eff_max) return;
                // packed_should_absorb: Hamming==1 in interior AND not damage-consistent.
                // G↔T/C↔A (8-oxoG) always protected; C↔T/G↔A protected when damage active.
                if (!packed_should_absorb(arena_.data(pid), child_packed,
                                          k5, ilen, profile_.enabled))
                    return;
                absorbed = true;
            };

            auto& sh = shard_it->second;
            sh.pi[0].query(pair_key(h0, h1, 0, ilen), check);
            if (!absorbed) sh.pi[1].query(pair_key(h0, h2, 1, ilen), check);
            if (!absorbed) sh.pi[2].query(pair_key(h1, h2, 2, ilen), check);

            if (absorbed) {
                is_error_[cid] = true;
                errcor_absorbed_++;
            }
        }

        log_info("Phase 3 complete: absorbed " + std::to_string(errcor_absorbed_) +
                 " PCR error sequences (ratio=" + std::to_string(errcor_.ratio) +
                 ", max-count=" + std::to_string(errcor_.max_count) + ")");
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

        ska::flat_hash_map<uint64_t, SequenceFingerprint> records_to_write;
        records_to_write.reserve(index_.size());
        for (const auto& [fingerprint, entry] : index_) {
            if (!is_error_.empty() && is_error_[entry.seq_id]) continue;
            records_to_write[entry.record_index] = fingerprint;
        }
        const size_t n_to_write = records_to_write.size();

        auto reader = make_fastq_reader(in_path);
        FastqRecord rec;
        uint64_t record_idx = 0;
        size_t written = 0;

        while (reader->read(rec)) {
            auto it = records_to_write.find(record_idx);
            if (it != records_to_write.end()) {
                writer.write(rec);

                if (cluster_gz) {
                    const SequenceFingerprint& fp = it->second;
                    auto ei = index_.find(fp);
                    if (ei == index_.end())
                        throw std::runtime_error(
                            "Internal error: missing fingerprint for cluster output");
                    const IndexEntry& entry = ei->second;
                    if (gzprintf(cluster_gz, "%016lx%016lx\t%lu\t%lu\n",
                                 static_cast<unsigned long>(fp.hash_hi),
                                 static_cast<unsigned long>(fp.hash_lo),
                                 static_cast<unsigned long>(rec.seq.size()),
                                 static_cast<unsigned long>(entry.count)) < 0)
                        throw std::runtime_error(
                            "gzprintf failed while writing cluster record");
                }

                written++;
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
            log_info("PCR error sequences removed: " + std::to_string(errcor_absorbed_));
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
        << "  --no-error-correct       Disable PCR error duplicate removal\n"
        << "  --error-correct          Explicitly enable (already default)\n"
        << "  --errcor-ratio FLOAT     Min count ratio parent/child (default: 50.0)\n"
        << "  --errcor-max-count INT   Only absorb clusters with count <= N (default: 5)\n"
        << "  --errcor-bucket-cap INT  Max pair-key bucket size (default: 64)\n"
        << "\nAncient DNA damage-aware deduplication (ON by default):\n"
        << "  --no-damage              Disable damage estimation and masking\n"
        << "  --damage-auto            Explicitly enable (already default)\n"
        << "  --damage-dmax  FLOAT     Set d_max for both 5' and 3' ends manually\n"
        << "  --damage-dmax5 FLOAT     Set d_max for 5' end only\n"
        << "  --damage-dmax3 FLOAT     Set d_max for 3' end only\n"
        << "  --damage-lambda FLOAT    Set lambda (decay) for both ends\n"
        << "  --damage-lambda5 FLOAT   Set lambda for 5' end only\n"
        << "  --damage-lambda3 FLOAT   Set lambda for 3' end only\n"
        << "  --damage-bg FLOAT        Background deamination rate (default: 0.02)\n"
        << "  --mask-threshold FLOAT   Mask positions with excess damage > T (default: 0.05)\n"
        << "  --pcr-cycles INT         Number of PCR cycles\n"
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

    bool   damage_auto    = true;   // default on: aDNA primary use case
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
        } else if (arg == "--errcor-ratio" && i + 1 < argc) {
            errcor.ratio = std::stod(argv[++i]);
        } else if (arg == "--errcor-max-count" && i + 1 < argc) {
            errcor.max_count = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--errcor-bucket-cap" && i + 1 < argc) {
            errcor.bucket_cap = static_cast<uint32_t>(std::stoul(argv[++i]));
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
    if (errcor.enabled && errcor.ratio <= 0.0) {
        std::cerr << "Error: --errcor-ratio must be > 0, got " << errcor.ratio << "\n";
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
    log_info("Damage-aware mode: " + std::string(damage_auto ? "auto (default)" :
             (damage_dmax5 >= 0.0 ? "manual" : "disabled (--no-damage)")));
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
            profile = estimate_damage(in_path, mask_threshold);
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
            profile.populate_mask_from_model();
            profile.print_info(/*typical_read_length=*/0);
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
