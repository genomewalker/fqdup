#define fqdup fqdup_r8qkeep_native_headers
#define fqalign fqdup_r8qkeep_native_alignment
#define fqquality fqdup_r8qkeep_native_quality
#define fqrepair fqdup_r8qkeep_native_repair
#include "fqdup/terminal_construct.hpp"
#include "fqdup/quality_compatibility.hpp"
#include "fqdup/adaptive_gap_evidence.hpp"
#include "fqdup/expected_error_qc.hpp"
#include "fqdup/terminal_learning.hpp"
#include <memory>
#include "fqdup/adapter_catalogue.hpp"
// fqdup merge — ultrafast paired-end read overlap detection and merging
//
// Geometry (aDNA dominant case, insert < read length):
//   R1   = [P5_adapter][insert_fwd][P7_RC_tail...]
//   rcR2 = [P5_adapter][insert_fwd][P7_tail...]      (RC of R2)
//   → R1[0:ov] == rcR2[0:ov] for ov = len(P5) + len(insert)
//   → Mismatch rate jumps at ov > ov_best (P7_RC ≠ P7 adapter sequences differ)
//
// Algorithm: build incremental mismatch prefix mm[0..L], find largest ov in
//   [min_ov, min(L1,L2)] where mm[ov]/ov ≤ max_mm_rate.
//   Falls back to 3'-of-R1 vs 5'-of-rcR2 scan for long inserts (no read-through).
//
// Thread layout:
//   reader thread → PairQueue (bounded) → N-1 worker threads
//                                       → MergeOutQueue (ordered) → writer thread

#include "fqdup/fastq_common.hpp"
#include "taph/frame_selector_decl.hpp"
#include "taph/profile_json.hpp"

#include <array>
#include <condition_variable>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring>
#include <deque>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

#ifdef __AVX2__
#include <immintrin.h>
#endif

static constexpr int PAIR_BATCH_SZ = 2048;
static constexpr int MAX_READ_LEN  = 600;  // packed buffer ceiling

// ============================================================================
// Structs
// ============================================================================

namespace fqdup_native_xc_internal {
struct ReadPair {
    FastqRecord r1, r2;
};

struct PairBatch {
    uint64_t id = 0;
    std::vector<ReadPair> pairs;
};

struct MergeRecord {
    FastqRecord merged;    // merged read (is_merged), trimmed R1 (unmerged), or orphan read (is_orphan)
    FastqRecord unmerged2; // trimmed R2 when is_merged=false and !is_orphan
    bool is_merged = false;
    bool is_orphan = false; // one mate discarded; surviving mate is in merged
    bool orphan_r1 = false; // when is_orphan: true=survivor was R1 (molecule 5'), false=R2 (RC 3')
};

struct MergeBatch {
    uint64_t id = 0;
    std::vector<MergeRecord> records;
    std::string outcomes_jsonl; // optional ordered per-input disposition audit
    int64_t n_merged     = 0;
    int64_t n_unmerged   = 0;
    int64_t n_orphan     = 0;
    int64_t n_orphan_r1  = 0; // survivor = R1 (molecule 5', C->T frame)
    int64_t n_orphan_r2  = 0; // survivor = R2 (RC molecule 3'; ds biological, ss dA-tail)
    int64_t n_dropped    = 0; // both mates failed QC after clipping
};

// Append a JSON string without allowing header whitespace to change its framing.
static void append_outcome_json_string(std::string& dst, const std::string& value) {
    static constexpr char hex[] = "0123456789abcdef";
    dst.push_back('"');
    for (unsigned char c : value) {
        switch (c) {
            case '"': dst += "\\\""; break;
            case '\\': dst += "\\\\"; break;
            case '\b': dst += "\\b"; break;
            case '\f': dst += "\\f"; break;
            case '\n': dst += "\\n"; break;
            case '\r': dst += "\\r"; break;
            case '\t': dst += "\\t"; break;
            default:
                if (c < 32) {
                    dst += "\\u00";
                    dst.push_back(hex[c >> 4]); dst.push_back(hex[c & 15]);
                } else dst.push_back(static_cast<char>(c));
        }
    }
    dst.push_back('"');
}

// A scoped disposition guard covers every early return/continue in the worker.
// Header copies and JSON accumulation are disabled unless outcomes were requested.
struct PairOutcomeGuard {
    MergeBatch* batch = nullptr;
    std::string r1_header, r2_header;
    int64_t merged = 0, unmerged = 0, orphan_r1 = 0, orphan_r2 = 0, dropped = 0;
    PairOutcomeGuard(MergeBatch& out, const ReadPair& pair, bool enabled) {
        if (!enabled) return;
        batch = &out;
        r1_header = pair.r1.header; r2_header = pair.r2.header;
        merged = out.n_merged; unmerged = out.n_unmerged;
        orphan_r1 = out.n_orphan_r1; orphan_r2 = out.n_orphan_r2;
        dropped = out.n_dropped;
    }
    ~PairOutcomeGuard() {
        if (!batch) return;
        const char* outcome = "unknown";
        const int64_t dm = batch->n_merged - merged;
        const int64_t du = batch->n_unmerged - unmerged;
        const int64_t d1 = batch->n_orphan_r1 - orphan_r1;
        const int64_t d2 = batch->n_orphan_r2 - orphan_r2;
        const int64_t dd = batch->n_dropped - dropped;
        if (dm + du + d1 + d2 + dd == 1) {
            if (dm == 1) outcome = "merged";
            else if (du == 1) outcome = "unmerged";
            else if (d1 == 1) outcome = "orphan_r1";
            else if (d2 == 1) outcome = "orphan_r2";
            else if (dd == 1) outcome = "dropped";
        }
        auto& text = batch->outcomes_jsonl;
        text += "{\"r1_header\":"; append_outcome_json_string(text, r1_header);
        text += ",\"r2_header\":"; append_outcome_json_string(text, r2_header);
        text += ",\"outcome\":\""; text += outcome; text += "\"}\n";
    }
};

// ============================================================================
// Queues (same bounded-FIFO / ordered-output pattern as trim.cpp)
// ============================================================================

struct PairQueue {
    std::mutex              mtx;
    std::condition_variable cv_ne, cv_nf;
    std::deque<PairBatch>   q;
    bool                    done     = false;
    int                     max_depth;
    explicit PairQueue(int d) : max_depth(d) {}

    void push(PairBatch&& b) {
        std::unique_lock lk(mtx);
        cv_nf.wait(lk, [&]{ return (int)q.size() < max_depth || done; });
        q.push_back(std::move(b));
        cv_ne.notify_one();
    }
    bool pop(PairBatch& b) {
        std::unique_lock lk(mtx);
        cv_ne.wait(lk, [&]{ return !q.empty() || done; });
        if (q.empty()) return false;
        b = std::move(q.front()); q.pop_front();
        cv_nf.notify_one();
        return true;
    }
    void set_done() {
        std::unique_lock lk(mtx);
        done = true;
        cv_ne.notify_all(); cv_nf.notify_all();
    }
};

struct MergeOutQueue {
    std::mutex mtx;
    std::condition_variable cv, cv_space;
    std::map<uint64_t, MergeBatch> pending;
    bool done = false;
    uint64_t next_expected = 0;
    size_t max_depth;
    explicit MergeOutQueue(size_t depth) : max_depth(std::max<size_t>(1, depth)) {}

    void push(MergeBatch&& b) {
        std::unique_lock lk(mtx);
        // Reserve admission for the missing next batch even when later batches
        // fill the queue. Blocking it would deadlock the ordered consumer.
        cv_space.wait(lk, [&] {
            return pending.size() < max_depth || b.id == next_expected || done;
        });
        pending.emplace(b.id, std::move(b));
        cv.notify_one();
    }
    bool pop_ordered(uint64_t expected, MergeBatch& out) {
        std::unique_lock lk(mtx);
        cv.wait(lk, [&]{ return pending.count(expected) || done; });
        auto it = pending.find(expected);
        if (it == pending.end()) return false;
        out = std::move(it->second);
        pending.erase(it);
        next_expected = expected + 1;
        cv_space.notify_all();
        return true;
    }
    void set_done() {
        std::unique_lock lk(mtx);
        done = true;
        cv.notify_all();
        cv_space.notify_all();
    }
};

// Per-position overlap substitution counts.
// r1_base vs rc2_base at each insert position, before consensus.
// A=0 C=1 G=2 T=3 — same as BASE2BIT encoding.
// Signals:
//   fwd[p][T][C] = C→T deamination at p from 5' end (top strand)
//   rev[p][G][A] = C→T deamination at p from 3' end (bottom strand, seen as G→A)
//   all [T][G]   = OxoG (8-oxoG→T misread, position-independent)
// Length bins for per-bin damage estimation. Edges: [0,35) [35,50) [50,75) [75,100) [100,150) [150+)
struct LenBins {
    static constexpr int N = 6;
    static constexpr int32_t EDGES[N + 1] = {0, 35, 50, 75, 100, 150, 32767};
    static int bin(int len) {
        for (int i = 0; i < N; ++i) if (len < EDGES[i+1]) return i;
        return N - 1;
    }
};

struct OverlapSubstCounts {
    static constexpr int MAX_POS = 30;
    int64_t fwd[MAX_POS][4][4]; // position from 5' end (all lengths combined)
    int64_t rev[MAX_POS][4][4]; // position from 3' end (all lengths combined)
    int64_t all[4][4];          // all positions combined
    int64_t n_pairs  = 0;
    int64_t n_bases  = 0;
    // v2 extension: per-length-bin matrices
    int64_t fwd_len[LenBins::N][MAX_POS][4][4];
    int64_t rev_len[LenBins::N][MAX_POS][4][4];
    int64_t bin_n_pairs[LenBins::N];
    int64_t bin_n_bases[LenBins::N];

    OverlapSubstCounts() { std::memset(this, 0, sizeof(*this)); }

    void merge(const OverlapSubstCounts& o) {
        for (int p = 0; p < MAX_POS; ++p)
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    { fwd[p][a][b] += o.fwd[p][a][b]; rev[p][a][b] += o.rev[p][a][b]; }
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) all[a][b] += o.all[a][b];
        n_pairs += o.n_pairs; n_bases += o.n_bases;
        for (int bi = 0; bi < LenBins::N; ++bi) {
            bin_n_pairs[bi] += o.bin_n_pairs[bi];
            bin_n_bases[bi] += o.bin_n_bases[bi];
            for (int p = 0; p < MAX_POS; ++p)
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b) {
                        fwd_len[bi][p][a][b] += o.fwd_len[bi][p][a][b];
                        rev_len[bi][p][a][b] += o.rev_len[bi][p][a][b];
                    }
        }
    }
};

// ============================================================================
// 2-bit packing: A=0, C=1, G=2, T=3, N→0
// 4 bases per byte, base i at bits 2*(i%4)
// ============================================================================

static const std::array<uint8_t, 256> BASE2BIT = []() {
    std::array<uint8_t, 256> t{};
    t[(uint8_t)'A'] = t[(uint8_t)'a'] = 0;
    t[(uint8_t)'C'] = t[(uint8_t)'c'] = 1;
    t[(uint8_t)'G'] = t[(uint8_t)'g'] = 2;
    t[(uint8_t)'T'] = t[(uint8_t)'t'] = 3;
    return t;
}();

static void pack_2bit(const char* seq, int L, uint8_t* out) {
    int n_bytes = (L + 3) / 4;
    std::memset(out, 0, n_bytes);
    for (int i = 0; i < L; ++i)
        out[i / 4] |= (uint8_t)(BASE2BIT[(uint8_t)seq[i]] << (2 * (i % 4)));
}

// ============================================================================
// Mismatch prefix: mm[i] = number of base mismatches in positions [0, i)
// mm has length L+1; mm[0]=0.
// ============================================================================

// skip_start: treat positions [0, skip_start) as always matching (terminal damage zone).
static void mismatch_prefix(const uint8_t* pa, const uint8_t* pb,
                            int* mm, int L, int skip_start = 0) {
    int n_bytes = (L + 3) / 4;

    // XOR packed arrays. For each packed byte:
    //   xor_byte ^ 0 → bits differ wherever bases differ
    //   (xor | (xor >> 1)) & 0x55 → one set bit per mismatched 2-bit pair
    int acc = 0;
    mm[0] = 0;

#ifdef __AVX2__
    alignas(32) uint8_t xbuf[MAX_READ_LEN / 4 + 32] = {};
    int i = 0;
    for (; i + 32 <= n_bytes; i += 32) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(pa + i));
        __m256i b = _mm256_loadu_si256((const __m256i*)(pb + i));
        _mm256_storeu_si256((__m256i*)(xbuf + i), _mm256_xor_si256(a, b));
    }
    for (; i < n_bytes; ++i) xbuf[i] = pa[i] ^ pb[i];

    for (int byte = 0; byte < n_bytes; ++byte) {
        uint8_t x  = xbuf[byte];
        uint8_t mb = (x | (x >> 1)) & 0x55u;
        int base   = byte * 4;
        for (int k = 0; k < 4 && base + k < L; ++k) {
            if (base + k >= skip_start) acc += (mb >> (k * 2)) & 1;
            mm[base + k + 1] = acc;
        }
    }
#else
    for (int byte = 0; byte < n_bytes; ++byte) {
        uint8_t x  = pa[byte] ^ pb[byte];
        uint8_t mb = (x | (x >> 1)) & 0x55u;
        int base   = byte * 4;
        for (int k = 0; k < 4 && base + k < L; ++k) {
            if (base + k >= skip_start) acc += (mb >> (k * 2)) & 1;
            mm[base + k + 1] = acc;
        }
    }
#endif
}

// ============================================================================
// Find best overlap: scan ov in [min_ov, L] for largest ov where mm[ov]/ov
// ≤ max_mm_rate. Stops early once 3+ consecutive positions exceed 2×threshold.
// Returns -1 if not found.
// ============================================================================

// skip_start: positions [0,skip_start) were excluded from mm[] — divide by effective length.
static int find_best_ov(const int* mm, int L, int min_ov, float max_mm_rate,
                        int skip_start = 0) {
    int best  = -1;
    int grace = 0;
    for (int ov = min_ov; ov <= L; ++ov) {
        int eff = ov - skip_start;
        if (eff <= 0) { best = ov; continue; }  // all positions skipped → accept
        float rate = (float)mm[ov] / (float)eff;
        if (rate <= max_mm_rate) {
            best  = ov;
            grace = 0;
        } else {
            if (best >= 0 && ++grace >= 4)
                break;
        }
    }
    return best;
}

// Quality-weighted overlap finder: mismatches at position i are penalized by
// min(Q1[i], Q2[i]) / 30 so low-quality disagreements (terminal damage) don't
// block merging. Drop-in replacement for find_best_ov in the hot path.
static int find_best_ov_qwt(
        const char* r1s, const char* r1q,
        const char* r2s, const char* r2q,
        int L, int min_ov, float max_mm_rate, int skip_start = 0) {
    float cumqmm = 0.f;
    int   best = -1, grace = 0;
    for (int ov = 1; ov <= L; ++ov) {
        int i = ov - 1;
        if (i >= skip_start && r1s[i] != r2s[i] &&
            r1s[i] != 'N' && r2s[i] != 'N') {
            int q1 = (uint8_t)r1q[i] - 33;
            int q2 = (uint8_t)r2q[i] - 33;
            cumqmm += std::min(q1, q2) / 30.f;
        }
        if (ov < min_ov) continue;
        int eff = ov - skip_start;
        if (eff <= 0) { best = ov; continue; }
        if (cumqmm / (float)eff <= max_mm_rate) {
            best = ov; grace = 0;
        } else if (best >= 0 && ++grace >= 4) {
            break;
        }
    }
    return best;
}

// ============================================================================
// Long-insert fallback: compare 3' tail of R1 vs 5' of rcR2 at each overlap.
// Returns best overlap length (in terms of suffix of R1 / prefix of rcR2),
// or -1 if none passes.
// ============================================================================

static inline int count_informative_mismatches(const char* a, const char* b, int n, int stop_after) {
    int mm = 0;
    int i = 0;
#ifdef __AVX2__
    const __m256i n_base = _mm256_set1_epi8('N');
    for (; i + 32 <= n; i += 32) {
        const __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
        const __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i));
        const __m256i ignored = _mm256_or_si256(
            _mm256_cmpeq_epi8(va, vb),
            _mm256_or_si256(_mm256_cmpeq_epi8(va, n_base),
                            _mm256_cmpeq_epi8(vb, n_base)));
        const unsigned mask = static_cast<unsigned>(_mm256_movemask_epi8(ignored));
        mm += 32 - __builtin_popcount(mask);
        if (mm > stop_after) return mm;
    }
#endif
    for (; i < n; ++i) {
        mm += (a[i] != b[i] && a[i] != 'N' && b[i] != 'N') ? 1 : 0;
        if (mm > stop_after) return mm;
    }
    return mm;
}

static bool overlap_is_low_complexity(const char* seq, int len);

static bool has_exact_long_self_overlap(const char* seq, int len) {
    // An exact repeated terminal segment gives a paired overlap a second
    // plausible placement in a longer insert. Without a reference or insert
    // prior, preserve the mates instead of choosing either geometry.
    for (int shift = 20; shift <= len - 20; ++shift) {
        if (std::memcmp(seq, seq + shift, len - shift) == 0) return true;
    }
    return false;
}

static int find_tail_head_ov(const char* r1, int L1,
                             const char* rc2, int L2,
                             int min_ov, float max_mm_rate,
                             bool require_unique = false) {
    int best = -1;
    const int max_ov = std::min(L1, L2) - (require_unique ? 1 : 0);
    for (int ov = min_ov; ov <= max_ov; ++ov) {
        const char* a = r1  + (L1 - ov);
        const char* b = rc2;
        // Early rejection avoids scanning the rest of an impossible overlap.
        // Ceil keeps the bound conservative; the original float comparison below
        // remains authoritative at exact threshold boundaries.
        const int stop_after = (require_unique && std::isfinite(max_mm_rate) &&
                                max_mm_rate >= 0.f && max_mm_rate <= 1.f) ?
            static_cast<int>(std::ceil(max_mm_rate * ov)) : ov;
        const int mm = count_informative_mismatches(a, b, ov, stop_after);
        if ((float)mm / (float)ov <= max_mm_rate) {
            // A short homopolymer/repeat match is not molecule-identity evidence.
            // Accuracy-first mode retains the mates rather than constructing a longer
            // molecule from a low-complexity tail/head coincidence.
            if (require_unique && overlap_is_low_complexity(a, ov)) continue;
            // Accuracy-first mode requires one identifiable insert geometry.
            // A second passing offset means sequence alone cannot choose safely.
            if (require_unique && best >= 0) return -2;
            best = ov;
        } else if (!require_unique && best >= 0 &&
                   mm - (int)((float)best * max_mm_rate) >= 3) {
            break;
        }
    }
    return best;
}

// ============================================================================
// Bayesian consensus at overlap.
// agree:    merged_q = min(q1 + q2, 60)  — Bayesian posterior; leeHom-equivalent
// disagree: call higher-qual base, merged_q = |q1 - q2| (min 1)
// ============================================================================

// r2_start: offset into rc2_seq for the start of the insert.
//   Short insert (Phase 0): r2_start = L2 - ov  (insert is at the END of RC(R2))
//   Long insert  (Phase 1): r2_start = 0         (insert fills all of RC(R2))
static void consensus_merge(const FastqRecord& r1,
                            const std::string& rc2_seq,
                            const std::string& rc2_qual,
                            int ov, int r2_start,
                            FastqRecord& out) {
    out.header = r1.header;
    auto& h = out.header;
    if (h.size() >= 2 && h[h.size()-2] == '/' &&
        (h.back() == '1' || h.back() == '2'))
        h.resize(h.size() - 2);
    out.plus = "+";
    // The overlap can never exceed either read. A spurious longer ov (seen on degenerate
    // poly-A/poly-G pairs whose overlap search mis-locks) would index past r1.seq / rc2_seq,
    // reading the string's NUL terminator and adjacent heap as bases+quals — emitting 0x00
    // bases and impossible qualities. Clamp to the valid span of both reads.
    if (r2_start < 0) r2_start = 0;
    int max_ov = std::min((int)r1.seq.size(), (int)rc2_seq.size() - r2_start);
    if (ov > max_ov) ov = std::max(max_ov, 0);
    out.seq.resize(ov);
    out.qual.resize(ov);

    for (int i = 0; i < ov; ++i) {
        char b1 = r1.seq[i];
        char b2 = rc2_seq[r2_start + i];
        int  q1 = (uint8_t)r1.qual[i]              - 33;
        int  q2 = (uint8_t)rc2_qual[r2_start + i]  - 33;

        // An unknown call provides no evidence for a canonical mate base.
        // Handle both orientations before combining agreeing observations.
        if (b1 == 'N' && b2 == 'N') {
            out.seq[i] = 'N';
            out.qual[i] = '!';
        } else if (b2 == 'N') {
            out.seq[i] = b1;
            out.qual[i] = (char)(q1 + 33);
        } else if (b1 == b2) {
            out.seq[i]  = b1;
            out.qual[i] = (char)(std::min(q1 + q2, 60) + 33);
        } else if (b1 == 'N') {
            out.seq[i]  = b2;
            out.qual[i] = (char)(q2 + 33);
        } else {
            if (q1 >= q2) {
                out.seq[i]  = b1;
                out.qual[i] = (char)(std::max(q1 - q2, 1) + 33);
            } else {
                out.seq[i]  = b2;
                out.qual[i] = (char)(std::max(q2 - q1, 1) + 33);
            }
        }
    }
}

// ============================================================================
// Adapter trimming for unmerged reads: semi-global search for adapter prefix
// in read suffix. Trims at first position where adapter matches with ≤2 mm.
// No-op if adapter is empty.
// ============================================================================

// Trim adapter sequence from the 5' end of a read (adapter complement bleedthrough artifact).
// Only trims if first alen (≤12) bases match adapter[0:alen] with ≤1 mismatch.
static void trim_adapter_5p(FastqRecord& rec, const std::string& adapter, int min_len) {
    if (adapter.empty() || (int)rec.seq.size() < min_len) return;
    int alen = std::min((int)adapter.size(), 12);
    if ((int)rec.seq.size() < alen) return;
    int mm = 0;
    for (int k = 0; k < alen; ++k)
        mm += (rec.seq[k] != adapter[k]) ? 1 : 0;
    if (mm > 1) return;
    int remaining = (int)rec.seq.size() - alen;
    if (remaining < min_len) { rec.seq.clear(); rec.qual.clear(); return; }
    fqdup::apply_retained_interval(rec.seq, rec.qual, alen, rec.seq.size());
}

// Trim a learned 5' construct (variable-length ss linker) from the read 5'.
// Anchors on the first 12 bp (<=1 mm), then confirms the FULL construct matches within
// max_mm_rate before removing exactly its length. Unlike trim_adapter_5p (12 bp cap), this
// removes the whole var-length linker so no 8 bp tail survives on the biological insert (I1).
static void trim_construct_5p(FastqRecord& rec, const std::string& construct,
                              float max_mm_rate, int min_len) {
    int clen = (int)construct.size();
    if (clen < 12 || (int)rec.seq.size() < clen) return;
    int anchor_mm = 0;
    for (int k = 0; k < 12; ++k)
        anchor_mm += (rec.seq[k] != construct[k]) ? 1 : 0;
    if (anchor_mm > 1) return;
    int mm = 0;
    for (int k = 0; k < clen; ++k)
        mm += (rec.seq[k] != construct[k]) ? 1 : 0;
    if ((float)mm > max_mm_rate * clen) return;
    if ((int)rec.seq.size() - clen < min_len) { rec.seq.clear(); rec.qual.clear(); return; }
    fqdup::apply_retained_interval(rec.seq, rec.qual, clen, rec.seq.size());
}

static void trim_adapter(FastqRecord& rec, const std::string& adapter, int min_len) {
    if (adapter.empty() || (int)rec.seq.size() < min_len) return;
    int L    = (int)rec.seq.size();
    int alen = (int)adapter.size();
    for (int start = std::max(0, L - alen - 5); start < L; ++start) {
        int cmp_len = std::min(alen, L - start);
        if (cmp_len < 6) break;
        int mm = 0;
        for (int k = 0; k < cmp_len; ++k)
            mm += (rec.seq[start + k] != adapter[k]) ? 1 : 0;
        float rate = (float)mm / (float)cmp_len;
        if (rate <= 0.10f) {
            if (start < min_len) {
                rec.seq.clear(); rec.qual.clear();
            } else {
                rec.seq.resize(start);
                rec.qual.resize(start);
            }
            return;
        }
    }
}

// Trim homopolymer run from 3' end (poly-G = NextSeq/NovaSeq dark cycle artifact).
static void trim_polybase(FastqRecord& rec, char base, int min_run, int min_length) {
    int L = (int)rec.seq.size();
    int run = 0, trim_to = L;
    for (int i = L - 1; i >= 0; --i) {
        if (rec.seq[i] == base) { ++run; if (run >= min_run) trim_to = i; }
        else break;
    }
    if (trim_to < L && trim_to >= min_length) {
        rec.seq.resize(trim_to);
        rec.qual.resize(trim_to);
    }
}

// Data-derived 3' poly-G run-length threshold. A trailing G-run is contamination
// (NovaSeq/NextSeq two-color dark-cycle run-through) only when longer than the
// interior G-composition would produce by chance: k = smallest run with
// P(run) < 1e-3 under interior G fraction f, i.e. k = ceil(log(1e-3)/log(f)).
// Floored at 4 so a G-rich genome can't push k so high the real tail survives.
// Safe by construction: on 4-channel/HiSeq data (no run-through) a random G-run
// reaching length k occurs with prob <1e-3, so trimming is a no-op.
static int derive_polyg_k(const std::vector<ReadPair>& scan_buf) {
    int64_t g = 0, n = 0;
    for (const auto& pr : scan_buf) {
        for (const std::string* s : {&pr.r1.seq, &pr.r2.seq}) {
            int L = (int)s->size(), lo = L / 3, hi = 2 * L / 3;  // interior third, unbiased by 3' tail
            for (int i = lo; i < hi; ++i) { ++n; if ((*s)[i] == 'G') ++g; }
        }
    }
    if (n == 0 || g == 0) return 4;
    double f = (double)g / (double)n;
    int k = (int)std::ceil(std::log(1e-3) / std::log(f));
    return k < 4 ? 4 : k;
}

// A two-color failure can begin at the same sequencing cycle across many pairs.
// Learn that cycle only from concordant long G tails in BOTH mates. A lone
// biological G tail cannot establish an instrument-wide boundary.
static int derive_paired_polyg_cycle(const std::vector<ReadPair>& scan_buf, int k) {
    if (k <= 0) return -1;
    std::unordered_map<int, int> starts;
    int eligible = 0;
    for (const auto& pair : scan_buf) {
        auto start = [k](const std::string& s) {
            int i = static_cast<int>(s.size());
            while (i > 0 && s[i - 1] == 'G') --i;
            return static_cast<int>(s.size()) - i >= k ? i : -1;
        };
        const int a = start(pair.r1.seq), b = start(pair.r2.seq);
        if (a < 0 || b < 0 || std::abs(a - b) > k - 1) continue;
        ++starts[std::max(a, b)];
        ++eligible;
    }
    int cycle = -1, support = 0;
    for (const auto& [position, n] : starts) {
        if (n > support) { cycle = position; support = n; }
    }
    const int floor = std::max(100, static_cast<int>(scan_buf.size() / 100));
    return support >= floor && support * 5 >= eligible * 4 ? cycle : -1;
}

static float shannon_entropy(const std::string& seq) {
    if (seq.empty()) return 0.f;
    int cnt[4] = {};
    for (char c : seq) {
        if      (c == 'A') cnt[0]++;
        else if (c == 'C') cnt[1]++;
        else if (c == 'G') cnt[2]++;
        else if (c == 'T') cnt[3]++;
    }
    float L = (float)seq.size(), h = 0.f;
    for (int i = 0; i < 4; ++i) {
        if (!cnt[i]) continue;
        float p = cnt[i] / L;
        h -= p * std::log2(p);
    }
    return h;
}

// Base-composition complexity of seq[lo,hi): Shannon entropy (bits, max 2) and the
// single-base dominance fraction (max over A/C/G/T). N excluded from the denominator.
// Both are the raw statistics whose merged-read distribution the unmerged QC gate is
// derived from (no hardcoded threshold — see detect_merge_params / passes_qc).
struct SeqComplexity { float entropy; float max_frac; };
static SeqComplexity seq_complexity(const std::string& seq, int lo, int hi) {
    int cnt[4] = {}, tot = 0;
    for (int i = lo; i < hi; ++i) {
        switch (seq[i]) {
            case 'A': cnt[0]++; tot++; break;  case 'C': cnt[1]++; tot++; break;
            case 'G': cnt[2]++; tot++; break;  case 'T': cnt[3]++; tot++; break;
        }
    }
    if (!tot) return {0.f, 1.f};
    float h = 0.f, mx = 0.f;
    for (int i = 0; i < 4; ++i) {
        if (!cnt[i]) continue;
        float p = (float)cnt[i] / (float)tot;
        h -= p * std::log2(p);
        if (p > mx) mx = p;
    }
    return {h, mx};
}

// Self-contained SDUST-style low-complexity score (Morgulis 2006): windowed triplet
// redundancy. In each window score = Σ_t c_t(c_t-1)/2 over the 64 triplets, normalised by
// the window's maximum (all-identical) score, i.e. the fraction of the maximum triplet
// redundancy. Random DNA → ~1/64 (0.016); a period-p tandem repeat → ~1/p (period-1=1.0,
// -2=0.5, -3=0.33, …). The read is low-complexity if ANY window exceeds DUST_MAX — the
// windowing catches PARTIAL low-complexity tracts (clean prefix + poly/repeat tail) that a
// whole-read statistic misses, and unlike base-composition entropy it sees k-mer-level
// repeats (GTCGTC…, GACGAC…). Self-contained: threshold is a structural constant, not tuned.
static constexpr int   DUST_W   = 64;      // analysis window (bp): ~whole-read for aDNA, so
                                           // it flags reads that are LARGELY low-complexity,
                                           // not longer reads with a short clean+repeat mix
                                           // (measured: W=64 drops the same 65% of promiscuous
                                           // hits as W=32 but a third the false positives).
static constexpr float DUST_MAX = 0.10f;   // random ~0.016; period-≤~10 repeat ≥0.10
static float dust_worst(const std::string& s) {
    const int L = (int)s.size();
    if (L < 5) return 0.f;
    auto b2 = [](char c) -> int { switch (c) { case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3; } return -1; };
    // triplet id (0..63) at each start, -1 if it spans an N
    std::vector<int> tri; tri.reserve(L);
    for (int i = 0; i + 3 <= L; ++i) {
        int a = b2(s[i]), b = b2(s[i+1]), c = b2(s[i+2]);
        tri.push_back((a<0||b<0||c<0) ? -1 : (a*16 + b*4 + c));
    }
    const int nt = (int)tri.size();
    int wt = std::min(DUST_W - 2, nt);         // triplets per window
    if (nt < 1 || wt < 2) return 0.f;
    int cnt[64] = {}; long score = 0; int valid = 0;
    auto add = [&](int t){ if (t >= 0) { score += cnt[t]; cnt[t]++; valid++; } };
    auto del = [&](int t){ if (t >= 0) { cnt[t]--; score -= cnt[t]; valid--; } };
    for (int i = 0; i < wt; ++i) add(tri[i]);
    float worst = 0.f;
    auto eval = [&]{ if (valid >= 2) { float mx = (float)valid * (valid - 1) / 2.f;
        float d = mx > 0.f ? (float)score / mx : 0.f; if (d > worst) worst = d; } };
    eval();
    for (int i = wt; i < nt; ++i) { del(tri[i - wt]); add(tri[i]); eval(); }
    return worst;
}

// Candidate overlaps are evaluated only after passing the mismatch gate. Compute the
// same validated SDUST statistic directly on the span: overlap search is a hot path, and
// constructing a temporary string/vector here costs more than the merge itself.
static bool overlap_is_low_complexity(const char* seq, int len) {
    if (len < 5) return false;
    auto b2 = [](char c) -> int { switch (c) { case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3; } return -1; };
    std::array<int, MAX_READ_LEN> tri{};
    const int nt = len - 2;
    for (int i = 0; i < nt; ++i) {
        const int a = b2(seq[i]), b = b2(seq[i + 1]), c = b2(seq[i + 2]);
        tri[i] = (a < 0 || b < 0 || c < 0) ? -1 : (a * 16 + b * 4 + c);
    }
    const int wt = std::min(DUST_W - 2, nt);
    if (wt < 2) return false;
    int cnt[64] = {}, valid = 0;
    long score = 0;
    auto add = [&](int t) { if (t >= 0) { score += cnt[t]; ++cnt[t]; ++valid; } };
    auto del = [&](int t) { if (t >= 0) { --cnt[t]; score -= cnt[t]; --valid; } };
    auto exceeds = [&]() {
        if (valid < 2) return false;
        const float maximum = static_cast<float>(valid) * (valid - 1) / 2.f;
        return static_cast<float>(score) / maximum > DUST_MAX;
    };
    for (int i = 0; i < wt; ++i) add(tri[i]);
    if (exceeds()) return true;
    for (int i = wt; i < nt; ++i) {
        del(tri[i - wt]); add(tri[i]);
        if (exceeds()) return true;
    }
    return false;
}

// Structural low-complexity window. A run-through/low-complexity tract is scored on its
// worst fixed-width window so the statistic is length-invariant: the merged-insert reference
// and the full-length unmerged mate are compared on the SAME 20bp-window scale, and a
// "clean prefix + poly-G tail" partial run-through is caught even when whole-read entropy
// stays moderate. 20bp is a structural analysis width (cf. min_ov), not a data threshold.
// ceiling: a low-complexity tract shorter than 20bp hides; upgrade: tie W to derived poly-G k.
static constexpr int     COMPLEXITY_W         = 20;
static constexpr int64_t COMPLEXITY_REF_FLOOR = 500;  // min merged inserts to enable the gate

// Worst window over seq[lo,hi): MIN Shannon entropy and MAX single-base dominance across all
// COMPLEXITY_W-wide windows (rolling counts, O(L)). N excluded from each window denominator.
static SeqComplexity worst_window(const std::string& seq, int lo, int hi, int W = COMPLEXITY_W) {
    int L = hi - lo;
    if (L <= W) return seq_complexity(seq, lo, hi);
    auto idx = [](char c) -> int {
        switch (c) { case 'A': return 0; case 'C': return 1;
                     case 'G': return 2; case 'T': return 3; } return -1; };
    int cnt[4] = {}, valid = 0;
    for (int i = lo; i < lo + W; ++i) { int b = idx(seq[i]); if (b >= 0) { cnt[b]++; valid++; } }
    float worst_ent = 2.f, worst_dom = 0.f;
    auto eval = [&]() {
        if (!valid) { worst_dom = 1.f; worst_ent = 0.f; return; }
        int mx = 0; for (int j = 0; j < 4; ++j) if (cnt[j] > mx) mx = cnt[j];
        float dom = (float)mx / (float)valid;
        if (dom > worst_dom) worst_dom = dom;
        float h = 0.f;
        for (int j = 0; j < 4; ++j) { if (!cnt[j]) continue;
            float p = (float)cnt[j] / (float)valid; h -= p * std::log2(p); }
        if (h < worst_ent) worst_ent = h;
    };
    eval();
    for (int i = lo + W; i < hi; ++i) {
        int bo = idx(seq[i - W]); if (bo >= 0) { cnt[bo]--; valid--; }
        int bn = idx(seq[i]);     if (bn >= 0) { cnt[bn]++; valid++; }
        eval();
    }
    return {worst_ent, worst_dom};
}

// Conservative homopolymer/degenerate test for the pre-merge fast reject. Fixed hard
// threshold (NOT the learned gate): a read is obvious junk only if its worst window is
// near-degenerate (a single base ≥90% or entropy <0.4 bits — poly-A/poly-G tracts). This
// catches only the extreme tail cheaply, before spending overlap-detection work; the
// data-derived merged-read gate downstream handles the borderline low-complexity cases.
static constexpr float HOMOPOLYMER_DOM = 0.90f;
static constexpr float HOMOPOLYMER_ENT = 0.40f;
static bool is_homopolymer_read(const std::string& seq) {
    if ((int)seq.size() < COMPLEXITY_W) {
        auto c = seq_complexity(seq, 0, (int)seq.size());
        return c.max_frac >= HOMOPOLYMER_DOM || c.entropy < HOMOPOLYMER_ENT;
    }
    auto w = worst_window(seq, 0, (int)seq.size());
    return w.max_frac >= HOMOPOLYMER_DOM || w.entropy < HOMOPOLYMER_ENT;
}

// Adapter-dimer test for an unmerged mate's 5' end. A zero-insert dimer reads straight
// into the read-through adapter, but a 1-base indel (e.g. missing leading A) shifts the
// match so clip_unmerged's exact anchored find_adapter_in misses it and the full adapter
// survives at high entropy. Slide both the read start and the adapter start over [0,max_off)
// and accept a 12-mer window with <=2 mismatches: reject when the learned adapter matches
// within the first few bases. Uses the ALREADY-LEARNED per-mate read-through adapter; empty
// adapter (undetected) → no dimer gate for that mate.
static bool is_adapter_dimer_5p(const std::string& seq, const std::string& adapter,
                                int max_off = 4, int anchor = 12, int max_mm = 2) {
    int L = (int)seq.size(), A = (int)adapter.size();
    if (A < anchor || L < anchor) return false;
    for (int rs = 0; rs < max_off && rs + anchor <= L; ++rs) {
        for (int as = 0; as < max_off && as + anchor <= A; ++as) {
            int mm = 0;
            for (int k = 0; k < anchor && mm <= max_mm; ++k)
                mm += (seq[rs + k] != adapter[as + k]);
            if (mm <= max_mm) return true;
        }
    }
    return false;
}

// Myers (1999) bit-parallel semi-global approximate matcher: minimum edit distance of a
// pattern (<=64bp -> one machine word) against any substring of a text, in O(text) word
// ops. Semi-global (free start/end in text) => shift-in 0 (edlib "HW"). This is the
// mapper-grade primitive (edlib/minimap2 use it); adapter trimmers mostly use heavier
// edit-distance DP (cutadapt) or k-mers (bbduk). Here it detects adapter-fragment merged
// reads (dimers / mismerges) that position-anchored 5'/3' checks miss.
struct BitMatcher {
    uint64_t Peq[4] = {0, 0, 0, 0};
    int m = 0;
    static int b2(char c) { switch (c) { case 'A': return 0; case 'C': return 1;
                                         case 'G': return 2; case 'T': return 3; } return -1; }
    void build(const std::string& p) {
        m = (int)std::min<size_t>(p.size(), 64);
        Peq[0] = Peq[1] = Peq[2] = Peq[3] = 0;
        for (int i = 0; i < m; ++i) { int b = b2(p[i]); if (b >= 0) Peq[b] |= (1ULL << i); }
    }
    struct Hit { int dist = INT_MAX; int end = -1; };   // min edit distance + text end index
    Hit search_pos(const char* T, int n) const {
        Hit h;
        if (m == 0) return h;
        uint64_t VP = ~0ULL, VN = 0; int score = m; uint64_t hb = 1ULL << (m - 1);
        for (int j = 0; j < n; ++j) {
            int b = b2(T[j]); uint64_t Eq = (b >= 0) ? Peq[b] : 0;
            uint64_t Xv = Eq | VN;
            uint64_t Xh = (((Eq & VP) + VP) ^ VP) | Eq;
            uint64_t Ph = VN | ~(Xh | VP);
            uint64_t Mh = VP & Xh;
            if (Ph & hb) score++; else if (Mh & hb) score--;
            Ph <<= 1; Mh <<= 1;
            VP = Mh | ~(Xv | Ph); VN = Ph & Xv;
            if (score < h.dist) { h.dist = score; h.end = j; }
        }
        return h;
    }
    int search(const char* T, int n) const { return search_pos(T, n).dist; }
};

// A merged read is an adapter fragment (dimer / adapter-only mismerge) when the WHOLE read
// aligns as a substring of the learned adapter (either orientation) within a small edit
// budget. Read-as-pattern (<=64bp) searched in the adapter text — this catches the SHORT
// (15-24bp) adapter fragments that a full-adapter match can't (read shorter than adapter).
// Tiered: an exact 11-mer of the read must appear in the adapter reference first (fast
// reject: real inserts share no adapter 11-mer, so they skip the bit-vector entirely).
static bool is_adapter_fragment(const std::string& read,
                                const std::string& adapter1, const std::string& adapter2) {
    const int L = (int)read.size();
    if (L < 12 || (int)adapter1.size() < 12) return false;
    // Fast reject: real inserts share no exact adapter 11-mer. Probe the read 5', middle and
    // 3' so a single INTERIOR indel (which frameshifts the sequence past its site) cannot hide
    // every seed — a merged interior-indel dimer keeps a clean 11-mer on one side of the indel.
    // A 5'-only probe missed these (the AGATCGGAAGAG G-run indel sits in the first 11bp). This
    // only decides whether Myers RUNS; the k-budget below still gates the drop, so broadening
    // the seed cannot false-drop a real insert (which won't align full-length within ~12%).
    bool seed_hit = false;
    for (int off : {0, (L - 11) / 2, L - 11}) {
        const std::string kmer = read.substr(off, 11);
        if (adapter1.find(kmer) != std::string::npos || adapter2.find(kmer) != std::string::npos) {
            seed_hit = true; break;
        }
    }
    if (!seed_hit) return false;
    BitMatcher bm; bm.build(read);
    const int k = 1 + L / 8;   // ~12% edit budget over the read
    return bm.search(adapter1.data(), (int)adapter1.size()) <= k ||
           bm.search(adapter2.data(), (int)adapter2.size()) <= k;
}

// Indel-tolerant 5' adapter-dimer test for an UNMERGED mate. is_adapter_dimer_5p uses a
// rigid Hamming frame (read-start/adapter-start slide, no indels) and misses zero-insert
// dimers carrying a single INTERIOR indel in the read-through adapter — the AGATCGGAAGAG
// G-run is the common site (measured ss unmerged leak ~1.6%: the full adapter survives at
// high entropy). Fast-reject on an exact adapter-5' seed sitting at the read 5' (real
// inserts almost never carry it), then confirm with a Myers ±2 semi-global anchor of the
// adapter 5' segment against the read 5' window: a dimer when the whole segment aligns
// within edit distance max_ed and ends near its expected read position. Same edit-budget
// specificity as is_adapter_fragment, so it does not eat real inserts.
// ceiling: an indel inside the first `seedn` bases breaks the exact seed and escapes;
// upgrade: a second interior seed. Not observed in the ss leak (indels sit at pos >=6).
static bool is_adapter_dimer_5p_indel(const std::string& seq, const std::string& adapter,
                                      int max_off = 4, int seedn = 6, int seg = 20, int max_ed = 2) {
    const int L = (int)seq.size(), A = (int)adapter.size();
    if (A < seedn + max_ed || L < seg) return false;
    size_t sp = seq.find(adapter.substr(0, seedn));
    if (sp == std::string::npos || (int)sp >= max_off) return false;  // adapter 5' must sit at the read 5'
    const int plen = std::min(A, seg);
    BitMatcher bm; bm.build(adapter.substr(0, plen));
    const int win = std::min(L, plen + max_ed + max_off);
    BitMatcher::Hit h = bm.search_pos(seq.data(), win);
    // whole adapter 5' segment aligns within budget, ending near where a 5' dimer would put it
    return h.dist <= max_ed
        && h.end >= plen - 1 - max_ed
        && h.end <= plen - 1 + max_ed + max_off;
}

// Verification instrumentation for the Myers adapter-fragment drop (highest-risk path):
// a false positive eats short damaged inserts and biases damage estimation. These count
// the drops and their length histogram so a run can confirm the drops are dimer-shaped
// (peaked at adapter length ~15-24bp), not shaving the real short-fragment tail.
static std::atomic<int64_t> g_frag_drop{0};
static std::atomic<int64_t> g_frag_len_hist[64];
// Audit hook: FQDUP_TECH_DUMP=<path> dumps every dropped-as-technical sequence so a run
// can confirm the drops are genuine adapter/construct, not real damaged inserts.
static std::mutex g_tech_dump_mx;
static std::ofstream* tech_dump_stream() {
    static std::ofstream* s = [] {
        const char* p = std::getenv("FQDUP_TECH_DUMP");
        return p ? new std::ofstream(p) : nullptr;
    }();
    return s;
}
static inline void note_frag_drop(const std::string& seq) {
    g_frag_drop.fetch_add(1, std::memory_order_relaxed);
    g_frag_len_hist[std::min((int)seq.size(), 63)].fetch_add(1, std::memory_order_relaxed);
    if (auto* s = tech_dump_stream()) {
        std::lock_guard<std::mutex> lk(g_tech_dump_mx);
        (*s) << seq << '\n';
    }
}
// Merged-read technical-sequence QC across ALL learned constructs (both orientations
// pre-expanded in `techs`). A MERGED read is [insert][3'adapter], so its 5' is the insert
// and it is technical only when the WHOLE read is adapter/construct — a zero-insert dimer
// that overlap-merged into pure adapter (is_adapter_fragment, indel-tolerant, exact-k-mer
// fast-reject). The 5'-anchor dimer check (is_adapter_dimer_5p) is intentionally NOT used
// here: it is calibrated for UNMERGED mates that read into adapter at the 5', and on merged
// reads it false-drops real inserts whose 5' coincidentally resembles the adapter start
// (measured: 98% of its drops were real DNA, <70% adapter identity).
static std::atomic<int64_t> g_reason_frag{0};
static bool is_technical_read(const std::string& read, const std::vector<std::string>& techs) {
    for (size_t ti = 0; ti < techs.size(); ++ti) {
        const auto& t = techs[ti];
        if (is_adapter_fragment(read, t, t)) { g_reason_frag.fetch_add(1, std::memory_order_relaxed); return true; }
    }
    return false;
}

// Scan from the 3' end of seq for an adapter prefix tail of length t in [1,max_tail].
// Returns insert_len (= L - t) for the longest match, or -1.
static int find_adapter_tail(const std::string& seq, const std::string& adapter,
                             int max_tail, int min_cand) {
    int L    = (int)seq.size();
    int alen = (int)adapter.size();
    for (int t = std::min(max_tail, alen); t >= 1; --t) {
        int cand = L - t;
        if (cand < min_cand) break;
        bool ok = true;
        for (int k = 0; k < t && ok; ++k)
            ok = (seq[cand + k] == adapter[k]);
        if (ok) return cand;
    }
    return -1;
}

// ============================================================================
// RC with reversed quality
// ============================================================================

static void rc_record(const FastqRecord& r2, std::string& seq, std::string& qual) {
    int L = (int)r2.seq.size();
    seq.resize(L); qual.resize(L);
    static const std::array<char, 256> COMP = []() {
        std::array<char, 256> t;
        t.fill('N');
        t[(uint8_t)'A'] = 'T'; t[(uint8_t)'a'] = 'T';
        t[(uint8_t)'C'] = 'G'; t[(uint8_t)'c'] = 'G';
        t[(uint8_t)'G'] = 'C'; t[(uint8_t)'g'] = 'C';
        t[(uint8_t)'T'] = 'A'; t[(uint8_t)'t'] = 'A';
        t[(uint8_t)'N'] = 'N'; t[(uint8_t)'n'] = 'N';
        return t;
    }();
    for (int i = 0; i < L; ++i) {
        seq[i]  = COMP[(uint8_t)r2.seq[L-1-i]];
        qual[i] = r2.qual[L-1-i];
    }
}

// ============================================================================
// Worker: processes PairBatches, emits MergeBatches
// ============================================================================

// Optional development quality stage; defaults retain native behavior.
} // namespace fqdup_native_xc_internal
#include "fqdup/pair_quality/quality_guard_xc.hpp"
namespace fqdup_native_xc_internal {

struct MergeOpts {
    double experimental_gapped_mott = -1.0;
    bool experimental_disagreement_risk = false;
    int experimental_selective_disagreement = 0;
    std::string experimental_pair_mixture_path;
    int experimental_pair_mixture_mode = 1;
    std::string experimental_full_pair_mixture_path;
    std::shared_ptr<const fqquality::PairMixture> experimental_full_pair_mixture;
    bool experimental_merged_disagreement_floor = false;
    std::shared_ptr<const fqquality::PairMixture> experimental_pair_mixture;
    int   min_ov        = 11;
    float max_mm_rate   = 0.08f;
    // Default 16 bp; 0=disabled. Enforced on ALL emit paths via passes_qc (merged +
    // unmerged mates). Raise to 30nt for DART: predict cannot use shorter reads, so
    // carrying them through sort/derep is wasted work and a low-complexity
    // candidate-explosion source.
    int   min_length    = 16;
    int   max_length    = 0;   // 0=disabled (no upper cap); drop reads longer than this on all emit paths
    int   skip_terminal = 0;
    int   clip_5p       = 0;   // 0=disabled; hard-clip N bases from R1 5' end before merge
    int   poly_g_min_run = 0;   // 0=disabled; trim 3' poly-G runs >= this length
    int   paired_poly_g_cycle = -1; // strong cohort evidence for a shared two-color failure cycle
    float min_entropy   = 0.0f; // 0=disabled; Shannon entropy floor (bits, max=2.0)
    float max_n_rate    = 1.0f; // 1=disabled; max fraction of N bases
    // Data-derived low-complexity gate for unmerged mates (from merged-mate reference), per mate.
    float complexity_entropy_lo_r1 = 0.f,  complexity_dom_hi_r1 = 1.0f;
    float complexity_entropy_lo_r2 = 0.f,  complexity_dom_hi_r2 = 1.0f;
    std::string adapter1;
    std::string adapter2;
    std::string adapter_5p_linker;              // library 5' construct clipped from unmerged reads
    std::string forced_library_type;            // "ss"/"ds" declared override; ""=auto-detect
    std::string terminal_model_path;
    std::shared_ptr<const fqdup::CompiledTerminalCatalogue> terminal_model;
    bool  known_terminal_constructs = false;      // explicit exact terminal residual cleaning on merged output
    bool  strict_library_type = false;          // abort if declared type contradicts confident detection
    std::string json_out;                       // --json: comprehensive lossless merge-QC report
    std::string outcomes_out;                   // optional per-input JSONL disposition audit
    std::vector<std::string> extra_adapters1;   // additional R1 adapters to try (from --adapter-fasta)
    std::vector<std::string> tech_seqs;          // ALL learned technical constructs (multi-adapter QC)
    bool  use_internal_panel = true;            // aDNA construct read-through table; default ON, --no-internal-panel disables
    bool  pair_supported_boundaries = true;     // accuracy default: require mate support, preserve ambiguity
    int   trim_quality = -1;                    // explicit quality-tail cutoff; disabled by default
    double qc_max_expected_errors = -1;         // opt-in development expected-error budget
    double qc_max_error_rate = -1;              // opt-in development retained mean-error budget
    bool  experimental_gapped_clipping = false; // separate opt-in clipping experiment
    bool  experimental_adaptive_overlap = false; // development-only quality model
    bool  experimental_short_tail = false;     // opt-in candidate, not independently validated
    bool  complexity_filter = false;            // preserve biological low-complexity reads by default
    bool  disable_polyg = false;
    std::string damage_out;     // path for paired damage profile JSON; empty=disabled
    std::string subst_out;      // path for overlap substitution matrix TSV; empty=disabled
    std::string subst_binary;   // path for binary .bsubst format; empty=disabled
};

// ============================================================================
// Pre-scan: auto-detect adapter sequences, library geometry, UDG status.
// Reads the first n_scan pairs; reader is positioned just after scan_buf on return.
// ============================================================================

struct DetectedParams {
    std::string adapter1;      // dominant learned construct on R1 past overlap
    std::string adapter2;      // its revcomp (RC(R2) side)
    std::vector<std::string> tech_seqs;   // ALL learned technical constructs (P7/P5/splint/index/...)
    // Per-construct provenance for the QC JSON (parallel to tech_seqs, before revcomp expansion).
    struct TechInfo { std::string seq; int64_t support = 0; double ic_body = 0.0; };
    std::vector<TechInfo> tech_info;
    bool   is_ss           = false;
    bool   type_from_panel = false;  // is_ss set by construct-panel stem vote (overrides geometry)
    float  type_confidence = 0.f;    // |prefix_agree_rate - 0.5| * 2, in [0,1]
    bool   is_udg          = false;
    bool   is_half_udg     = false;
    int    skip_terminal   = 0;
    float  damage_5p       = 0.f;
    float  damage_3p       = 0.f;
    std::string adapter_5p_linker;   // library-specific fixed 5' construct (adapter-dimer prefix)
    // Low-complexity QC gate derived from the merged-MATE insert reference (worst-window
    // entropy/dominance), kept separate per mate because the run-through artifact is R2-heavy.
    // A mate is junk if its worst-window entropy < entropy_lo OR dominance > dom_hi.
    // Defaults (0 / 1.0) = gate disabled (the small-sample fallback below the ref floor).
    float   complexity_entropy_lo_r1 = 0.f,  complexity_dom_hi_r1 = 1.0f;
    float   complexity_entropy_lo_r2 = 0.f,  complexity_dom_hi_r2 = 1.0f;
    int64_t complexity_ref_n_r1 = 0,         complexity_ref_n_r2 = 0;
};

// Forward declaration (defined below)
static int find_adapter_in(const std::string& seq, const std::string& adapter,
                           int min_pos, int max_mm, int alen_req);

static const std::array<int,256>& make_base_idx() {
    static std::array<int,256> t;
    t.fill(-1);
    t[(uint8_t)'A'] = t[(uint8_t)'a'] = 0;
    t[(uint8_t)'C'] = t[(uint8_t)'c'] = 1;
    t[(uint8_t)'G'] = t[(uint8_t)'g'] = 2;
    t[(uint8_t)'T'] = t[(uint8_t)'t'] = 3;
    return t;
}
static const std::array<int,256>& BASE_IDX = make_base_idx();

static DetectedParams detect_merge_params(
        std::vector<ReadPair>& scan_buf,
        int min_ov, float mm_loose = 0.15f) {

    // Per-position (0..29) mismatch tallies across all pairs with a found overlap
    std::vector<int64_t> pos_mm(30, 0), pos_total(30, 0);
    // 5' prefix agreement (ds geometry check): R1[0:8] vs rcR2[0:8]
    int64_t prefix_agree = 0, prefix_total = 0;
    // Adapter suffix k-mer frequency
    std::unordered_map<std::string, int> a1_freq, a2_freq;
    // Full read-through adapter learned by overlap-consensus: r1[p1:] (the region past the
    // insert boundary) IS adapter1, read straight from the data. Accumulate per-position
    // ACGT counts across the sample; a column-majority walk to the coverage/agreement cliff
    // yields the full-length adapter (30+bp) instead of a 12bp seed — making read-through
    // trimming and adapter-dimer detection shift-robust with no magic anchor. O(MAXADAPT)/pair.
    static constexpr int MAXADAPT   = 48;
    static constexpr int ADAPT_MINN = 100;   // min supporting pairs before we trust a learned construct
    static constexpr int MAX_TECH   = 6;     // max distinct technical constructs learned (P7/P5/splint/index/...)
    // Min beyond-seed consensus information (bits) to accept a construct as a real fixed
    // adapter rather than a phantom cluster of real inserts sharing a leading 10-mer.
    // 0 during calibration (accept all + print); set to the data-derived cut after.
    static constexpr double MIN_IC_BODY = 0.0;
    // Multi-construct learning: the region past the insert is technical sequence, read from
    // the data. Cluster those post-insert regions by their leading 10-mer (different constructs
    // start differently); each cluster's per-position majority consensus is one full technical
    // sequence. This is adapter-AGNOSTIC (harvested at the overlap boundary, not via a known
    // adapter), so it learns EVERY construct present, not just the top-1.
    struct TechCluster { std::array<std::array<int64_t, 4>, MAXADAPT> cons{}; int64_t n = 0; };
    std::unordered_map<uint32_t, TechCluster> tech_clusters;
    // R1 5' 12-mer frequency: a fixed library construct (adapter-dimer prefix) spikes here;
    // real inserts spread thin (<1%). Used to detect a 5' linker to clip from unmerged reads.
    std::unordered_map<std::string, int> r1_5p_freq;
    int64_t n_5p = 0;
    // ss/ds construct-panel stem vote (Ellesmere supp §4.3.2): the SCR splint core GGAAGAGCGTCG
    // reads through R1 3' with NO upstream AGATC in ss libraries, but sits behind the TruSeq stem
    // (AGATCGGAAGAGCGTCG) in ds. Presence/absence of the AGATC prefix is a positive ss/ds signal.
    static const std::string SS_SPLINT_CORE = "GGAAGAGCGTCG";
    int64_t ss_stem = 0, ds_stem = 0;

    // Merged-mate complexity reference: worst-window entropy + base-dominance of each
    // overlap-verified insert, kept per mate (R1[0:ov], R2[0:ov]). These trusted real
    // fragments define the distribution the unmerged low-complexity gate is derived from —
    // no hardcoded entropy/dominance constant. Same worst-window statistic as the gate target.
    std::vector<float> ref_ent_r1, ref_dom_r1, ref_ent_r2, ref_dom_r2;

    alignas(32) uint8_t pa[MAX_READ_LEN / 4 + 4] = {};
    alignas(32) uint8_t pb[MAX_READ_LEN / 4 + 4] = {};
    int mm[MAX_READ_LEN + 1];
    std::string rc2_seq, rc2_qual;

    for (auto& pr : scan_buf) {
        int L1 = (int)pr.r1.seq.size();
        int L2 = (int)pr.r2.seq.size();
        if (L1 < min_ov || L2 < min_ov || L1 > MAX_READ_LEN || L2 > MAX_READ_LEN) continue;

        rc_record(pr.r2, rc2_seq, rc2_qual);
        int L = std::min(L1, L2);
        pack_2bit(pr.r1.seq.data(), L, pa);
        pack_2bit(rc2_seq.data(),   L, pb);
        // Detect overlap via adapter search in R1 (correct geometry).
        // For short inserts: R1=[insert][adapter1]; adapter1 at position=ov.
        // For long inserts (no adapter): use loose d=0 Hamming as fallback.
        int ov = -1;
        int r2s = 0;  // start of insert in rc2_seq

        // Phase 0 adapter search
        static const std::string TRUSEQ_R1  = fqdup::adapters::truseq_3p;
        static const std::string TRUSEQ_RC2 = fqdup::adapters::truseq_5p;
        int p1 = find_adapter_in(pr.r1.seq, TRUSEQ_R1, min_ov, 2, 12);
        if (p1 >= min_ov && p1 <= L1) {
            ov  = p1;
            r2s = L2 - p1;
        }
        // Fallback: loose d=0 Hamming (works for long inserts where RC(R2) starts with insert)
        if (ov < min_ov) {
            mismatch_prefix(pa, pb, mm, L, 0);
            int ov_d0 = find_best_ov(mm, L, min_ov, mm_loose, 0);
            if (ov_d0 >= min_ov) { ov = ov_d0; r2s = 0; }
        }
        if (ov < min_ov) continue;

        // Overlap-verified insert → per-mate worst-window complexity reference. Only inserts
        // at least one window wide contribute, so short-fragment sampling noise never widens
        // the reference tails. R2's insert sits at its 5' too (R2=[insert_RC][adapter2]).
        if (ov >= COMPLEXITY_W) {
            auto w1 = worst_window(pr.r1.seq, 0, ov);
            ref_ent_r1.push_back(w1.entropy); ref_dom_r1.push_back(w1.max_frac);
            if (ov <= L2) {
                auto w2 = worst_window(pr.r2.seq, 0, ov);
                ref_ent_r2.push_back(w2.entropy); ref_dom_r2.push_back(w2.max_frac);
            }
        }

        // ---- ds geometry check: R1[0:4] vs insert portion of RC(R2) ----
        {
            int chk = std::min(4, ov);
            bool agree = true;
            for (int i = 0; i < chk && agree; ++i)
                agree = (pr.r1.seq[i] == rc2_seq[r2s + i]);
            ++prefix_total;
            if (agree) ++prefix_agree;
        }

        // ---- per-position mismatch (damage profile) using correct RC(R2) slice ----
        for (int i = 0; i < std::min(ov, 30); ++i) {
            int j = r2s + i;
            if (j >= L2) break;
            bool is_mm = (pr.r1.seq[i] != rc2_seq[j] &&
                          pr.r1.seq[i] != 'N' && rc2_seq[j] != 'N');
            pos_mm[i]    += is_mm ? 1 : 0;
            pos_total[i] += 1;
        }
        // 3' terminal: last 5 positions of overlap (positions ov-5..ov-1 from R1)
        for (int k = 0; k < 5 && ov - 1 - k >= 0; ++k) {
            int i = ov - 1 - k;
            int j = r2s + i;
            if (j >= L2) continue;
            bool is_mm = (pr.r1.seq[i] != rc2_seq[j] &&
                          pr.r1.seq[i] != 'N' && rc2_seq[j] != 'N');
            if (k + 20 < (int)pos_mm.size()) {
                pos_mm[20 + k]    += is_mm ? 1 : 0;
                pos_total[20 + k] += 1;
            }
        }

        // ---- adapter k-mer collection ----
        if (p1 >= min_ov && p1 + 12 <= L1)
            a1_freq[pr.r1.seq.substr(p1, 12)]++;
        // multi-construct consensus: the region past the overlap boundary `ov` is technical
        // sequence, whatever the construct (harvested at the ADAPTER-AGNOSTIC overlap boundary,
        // not via a known adapter). Cluster by its leading 10-mer so different constructs
        // (P7 / P5 / splint / index) accumulate into separate consensuses.
        if (ov >= min_ov && ov + 10 <= L1) {
            uint32_t key = 0; bool ok = true;
            for (int j = 0; j < 10 && ok; ++j) {
                switch (pr.r1.seq[ov + j]) {
                    case 'A': key = (key << 2);          break;
                    case 'C': key = (key << 2) | 1u;     break;
                    case 'G': key = (key << 2) | 2u;     break;
                    case 'T': key = (key << 2) | 3u;     break;
                    default:  ok = false;                break;
                }
            }
            if (ok) {
                TechCluster& tc = tech_clusters[key];
                int amax = std::min(MAXADAPT, L1 - ov);
                for (int j = 0; j < amax; ++j) {
                    switch (pr.r1.seq[ov + j]) {
                        case 'A': tc.cons[j][0]++; break; case 'C': tc.cons[j][1]++; break;
                        case 'G': tc.cons[j][2]++; break; case 'T': tc.cons[j][3]++; break;
                    }
                }
                tc.n++;
            }
        }
        // adapter2 in RC(R2) at position r2s - 12 (just before the insert)
        if (r2s >= 12)
            a2_freq[rc2_seq.substr(r2s - 12, 12)]++;
        if (L1 >= 12) { r1_5p_freq[pr.r1.seq.substr(0, 12)]++; ++n_5p; }
        // ss/ds vote: the splint core reads through R1 3'. In ds it sits behind the TruSeq stem
        // (AGATC upstream), in ss it does not. Score ONLY reads where the stem is in-frame and
        // decisive: undeterminable (core within 5bp of read start) and ambiguous (2 mismatches to
        // AGATC, near the random-sequence expectation) reads cast no vote — an absent/degraded stem
        // is not positive ss evidence. Symmetric tolerance on both arms avoids a ds->ss bias that
        // would drop the ds 3' G->A terminus.
        int sp = find_adapter_in(pr.r1.seq, SS_SPLINT_CORE, 0, 1, 12);
        if (sp >= 5) {
            static const char STEM[5] = {'A','G','A','T','C'};
            int d0 = 0;
            for (int k = 0; k < 5; ++k) d0 += (pr.r1.seq[sp - 5 + k] != STEM[k]);
            if      (d0 <= 1) ++ds_stem;   // TruSeq stem present -> ds
            else if (d0 >= 3) ++ss_stem;   // clearly no stem -> ss
        }
    }

    DetectedParams d;

    // ---- derive the unmerged low-complexity gate from the merged-mate reference ----
    // Gate an unmerged mate when its worst window is MORE extreme than the merged mates'
    // own worst-window tails: entropy below their 1st percentile OR base-dominance above
    // their 99th percentile. The 1%/99% choice keeps the gate off the body of the real
    // distribution (≤1% of genuine biased windows sit past each tail) while poly-G
    // (entropy≈0, dominance≈1) falls far outside it. Below COMPLEXITY_REF_FLOOR trusted
    // inserts the gate stays disabled (defaults 0 / 1.0) — only clip + min-length + dimer QC
    // apply — so a tiny/blank library cannot derive a degenerate gate.
    auto quantile = [](std::vector<float>& v, double q) {
        std::sort(v.begin(), v.end());
        return v[(size_t)std::llround(q * (double)(v.size() - 1))];
    };
    auto derive = [&](std::vector<float>& ent, std::vector<float>& dom,
                      float& ent_lo, float& dom_hi, int64_t& n) {
        n = (int64_t)ent.size();
        if (n >= COMPLEXITY_REF_FLOOR) { ent_lo = quantile(ent, 0.01); dom_hi = quantile(dom, 0.99); }
    };
    derive(ref_ent_r1, ref_dom_r1, d.complexity_entropy_lo_r1, d.complexity_dom_hi_r1, d.complexity_ref_n_r1);
    derive(ref_ent_r2, ref_dom_r2, d.complexity_entropy_lo_r2, d.complexity_dom_hi_r2, d.complexity_ref_n_r2);

    if (prefix_total == 0) return d;

    // ---- library geometry ----
    float prefix_agree_rate = (float)prefix_agree / (float)prefix_total;
    d.is_ss = (prefix_agree_rate < 0.3f);
    d.type_confidence = std::min(1.f, std::fabs(prefix_agree_rate - 0.5f) * 2.f);
    // Construct-panel stem vote overrides geometry when it carries strong positive sequence
    // evidence: the geometry check alone miscalls SCR-splint ss libraries as ds (memory: Med11,
    // FLB01mAss4 both called ds conf~0.95 despite being ss). The splint stem is decisive.
    int64_t stem_tot = ss_stem + ds_stem;
    if (stem_tot >= 100) {
        float ss_frac = (float)ss_stem / (float)stem_tot;
        if (ss_frac >= 0.6f) {
            d.is_ss = true;  d.type_from_panel = true;
            d.type_confidence = std::max(d.type_confidence, ss_frac);
        } else if (ss_frac <= 0.4f) {
            d.is_ss = false; d.type_from_panel = true;
            d.type_confidence = std::max(d.type_confidence, 1.f - ss_frac);
        }
    }

    // ---- damage rates ----
    auto rate = [&](int i) -> float {
        return pos_total[i] > 0 ? (float)pos_mm[i] / (float)pos_total[i] : 0.f;
    };
    // 5' damage: average over positions 0-3
    float r5 = 0.f;
    for (int i = 0; i < 4; ++i) r5 += rate(i);
    d.damage_5p = r5 / 4.f;
    // 3' damage: positions stored in slots 20-23
    float r3 = 0.f;
    for (int i = 20; i < 24; ++i) r3 += rate(i);
    d.damage_3p = r3 / 4.f;

    // Library typing against THIS library's own interior mismatch null (positions 8-19: past the
    // 5' damage decay, before the 3' terminal slots 20-23). No hardcoded damage-rate cut: an end is
    // "damaged" iff its terminal C→T excess over the interior clears zero at 95% (pooled binomial SE
    // for the terminal and interior rates). This self-calibrates to each library's own overlap
    // seq-error + composition floor instead of a panel-tuned 0.02/0.005/0.01/2.5.
    auto pooled = [&](int lo, int hi) {
        int64_t mm = 0, tot = 0;
        for (int i = lo; i < hi && i < (int)pos_mm.size(); ++i) { mm += pos_mm[i]; tot += pos_total[i]; }
        double p = tot > 0 ? (double)mm / (double)tot : 0.0;
        double se = tot > 0 ? std::sqrt(p * (1.0 - p) / (double)tot) : 1.0;
        return std::pair<double, double>{p, se};
    };
    auto [p5, se5]     = pooled(0, 4);    // 5' terminal
    auto [p3, se3]     = pooled(20, 24);  // 3' terminal (slots 20-23)
    auto [pint, seint] = pooled(8, 20);   // interior baseline
    const bool dmg5 = (p5 - pint) - 1.96 * std::sqrt(se5 * se5 + seint * seint) > 0.0;
    const bool dmg3 = (p3 - pint) - 1.96 * std::sqrt(se3 * se3 + seint * seint) > 0.0;
    if (dmg5 && !dmg3) {
        d.is_half_udg   = true;    // 5' above interior, 3' not → half-UDG
        d.skip_terminal = 4;
    } else if (dmg5 && dmg3) {
        d.skip_terminal = 4;       // both ends above interior → untreated
    } else if (!dmg5 && !dmg3) {
        d.is_udg        = true;    // neither end above interior → UDG (or no signal)
        d.skip_terminal = 0;
    } else {
        d.skip_terminal = 4;       // only 3' above interior (rare) → still trim terminus
    }

    // ---- adapter sequences ----
    // adapter2 = what appears in RC(R2) at insert boundary = RC(TruSeq_P7) = CTCTTCCGATCT
    // adapter1 = what appears in R1 at insert boundary = RC(adapter2)
    // Fall back to standard TruSeq if detection is unclear.
    static const std::string TRUSEQ_R1  = fqdup::adapters::truseq_3p;   // appears in R1 at ov
    static const std::string TRUSEQ_RC2 = fqdup::adapters::truseq_5p;   // appears in RC(R2) at ov

    auto best_kmer = [](const std::unordered_map<std::string,int>& freq) -> std::string {
        std::string best; int best_n = 0;
        for (auto& [k,v] : freq) if (v > best_n) { best_n = v; best = k; }
        return best;
    };
    // Full read-through adapter1 from the R1 post-insert consensus. Walk columns while
    // coverage stays >= half the peak AND the majority base clears 60% — that is the fixed
    // adapter body; where either drops we have run into the (random) insert continuation, so
    // stop. Below ADAPT_MINN supporting pairs the consensus is untrusted → fall back to the
    // 12bp seed / TruSeq default (a safe adapter, never a degenerate one).
    // Consensus of one cluster: walk columns while coverage >= half the peak AND the majority
    // base clears 60% — the fixed construct body; where either drops we've hit the random insert
    // continuation, so stop.
    auto consensus_of = [&](const TechCluster& tc) -> std::string {
        int64_t peak = tc.n; std::string a;
        for (int j = 0; j < MAXADAPT; ++j) {
            int64_t tot = 0, mx = 0; int mb = 0;
            for (int b = 0; b < 4; ++b) { tot += tc.cons[j][b]; if (tc.cons[j][b] > mx) { mx = tc.cons[j][b]; mb = b; } }
            if (tot == 0 || tot < peak / 2) break;
            if ((double)mx / (double)tot < 0.60) break;
            a += "ACGT"[mb];
        }
        return a;
    };
    // Learn ALL technical constructs: rank clusters by support, consensus the top MAX_TECH with
    // a real body (>=12bp) and enough support. tech_seqs[0] is the dominant read-through adapter.
    std::vector<std::pair<int64_t, const TechCluster*>> ranked;
    ranked.reserve(tech_clusters.size());
    for (auto& kv : tech_clusters) ranked.push_back({kv.second.n, &kv.second});
    std::sort(ranked.begin(), ranked.end(),
              [](const std::pair<int64_t,const TechCluster*>& a,
                 const std::pair<int64_t,const TechCluster*>& b){ return a.first > b.first; });
    // Likelihood-ratio construct gate (fixed-adapter hypothesis vs random-DNA null).
    // Clusters are seeded on a shared leading 10-mer, so the first SEEDLEN columns carry
    // ~2 bits each by construction — for BOTH a real adapter and a phantom cluster of real
    // inserts that merely share a 10-mer. The discriminator is the information content of
    // the consensus BEYOND the seed: a fixed construct stays determined (~1.7 bits/col) for
    // 20-40 more columns; a phantom decays to background (~0 bits) at once. IC_body is that
    // beyond-seed information in bits; only constructs carrying real fixed sequence pass.
    constexpr int SEEDLEN = 10;
    auto ic_body_of = [&](const TechCluster& tc, int len) -> double {
        double ic = 0.0;
        for (int j = SEEDLEN; j < len && j < MAXADAPT; ++j) {
            int64_t tot = 0;
            for (int b = 0; b < 4; ++b) tot += tc.cons[j][b];
            if (tot == 0) continue;
            double H = 0.0;
            for (int b = 0; b < 4; ++b) {
                double p = (double)tc.cons[j][b] / (double)tot;
                if (p > 0) H -= p * std::log2(p);
            }
            ic += (2.0 - H);
        }
        return ic;
    };
    const bool tech_dbg = std::getenv("FQDUP_TECH_LEARN") != nullptr;
    for (auto& r : ranked) {
        if ((int)d.tech_seqs.size() >= MAX_TECH || r.first < ADAPT_MINN) break;
        std::string a = consensus_of(*r.second);
        double icb = ic_body_of(*r.second, (int)a.size());
        if (tech_dbg)
            std::cerr << "[tech-learn] support=" << r.first << " len=" << a.size()
                      << " ic_body=" << icb << "  " << a << "\n";
        if ((int)a.size() >= 12 && icb >= MIN_IC_BODY) {
            d.tech_seqs.push_back(a);
            d.tech_info.push_back({a, r.first, icb});
        }
    }
    if (!d.tech_seqs.empty()) {
        d.adapter1 = d.tech_seqs[0];          // dominant learned construct (30+bp typical)
        d.adapter2 = revcomp(d.adapter1);
    } else {
        std::string a2 = best_kmer(a2_freq);
        if (a2.empty() || a2_freq[a2] < 100) a2 = TRUSEQ_RC2;
        d.adapter2 = a2;
        d.adapter1 = revcomp(d.adapter2);
        d.tech_seqs.push_back(d.adapter1);    // never leave the list empty (safe fallback)
    }

    // ---- 5' library linker (adapter-dimer prefix) ----
    // A single 12-mer dominating R1 5' ends (>=5% of reads) that is NOT the read-through
    // adapter is a fixed library construct fused to insert-less dimers (observed:
    // CGCAATGCTCAT...GGACTCAA + P7). Register it so unmerged reads get it clipped and the
    // residual zero-insert body falls below min-length. Real inserts spread <1% per 12-mer.
    if (n_5p > 0) {
        std::string top = best_kmer(r1_5p_freq);
        auto is_known = [&](const std::string& s) {
            return s.compare(0, 8, TRUSEQ_R1, 0, 8) == 0
                || s.compare(0, 8, TRUSEQ_RC2, 0, 8) == 0
                || (!d.adapter1.empty() && s.compare(0, 8, d.adapter1, 0, 8) == 0)
                || (!d.adapter2.empty() && s.compare(0, 8, d.adapter2, 0, 8) == 0);
        };
        if (!top.empty() && (float)r1_5p_freq[top] / (float)n_5p >= 0.05f && !is_known(top)) {
            // Greedy-extend the 12-mer seed to the full var-length construct. Second pass over
            // reads whose 5' matches the seed: build a per-position majority consensus. In dimers
            // the construct is linker+read-through adapter, so the consensus reads linker then
            // adapter — the linker length is where the universal adapter begins (protocol constant,
            // not a calibration). Fallback: interior-purity boundary (self-calibrated, no magic).
            static const int LKMAX = 45;
            int64_t cnt[LKMAX][4] = {};
            for (auto& pr : scan_buf) {
                const std::string& s = pr.r1.seq;
                if ((int)s.size() < 12) continue;
                int amm = 0;
                for (int k = 0; k < 12; ++k) amm += (s[k] != top[k]) ? 1 : 0;
                if (amm > 1) continue;
                int n = std::min((int)s.size(), LKMAX);
                for (int j = 0; j < n; ++j) {
                    int b = BASE_IDX[(uint8_t)s[j]];
                    if (b >= 0) ++cnt[j][b];
                }
            }
            std::string cons(LKMAX, 'N');
            float maj[LKMAX] = {};
            for (int j = 0; j < LKMAX; ++j) {
                int64_t tot = cnt[j][0]+cnt[j][1]+cnt[j][2]+cnt[j][3];
                if (tot == 0) { maj[j] = 0.f; continue; }
                int bb = 0; for (int b = 1; b < 4; ++b) if (cnt[j][b] > cnt[j][bb]) bb = b;
                cons[j] = "ACGT"[bb];
                maj[j]  = (float)cnt[j][bb] / (float)tot;
            }
            int linker_len = 12;
            int adp = find_adapter_in(cons, TRUSEQ_R1, 12, 2, 12);
            if (adp >= 12) {
                linker_len = adp;
            } else {
                float m_seed = 0.f; for (int j = 0; j < 12; ++j) m_seed += maj[j]; m_seed /= 12.f;
                float m_int = 0.f; int ni = 0;
                for (int j = 30; j < LKMAX; ++j) if (maj[j] > 0.f) { m_int += maj[j]; ++ni; }
                m_int = ni ? m_int / ni : 0.30f;
                float thresh = 0.5f * (m_seed + m_int);
                int j = 12;
                while (j < LKMAX && maj[j] >= thresh) ++j;
                linker_len = j;
            }
            d.adapter_5p_linker = cons.substr(0, linker_len);
        }
    }

    return d;
}

// Exact bytewise mismatch count with early exit. Adapter probes are short, but
// they run at every candidate boundary; vectorizing the common 8/16/32-byte
// blocks removes that inner scalar loop without changing its mismatch contract.
static inline int count_mismatches_limited(const char* lhs, const char* rhs,
                                           int len, int limit) {
    int mm = 0, i = 0;
#ifdef __AVX2__
    for (; i + 32 <= len; i += 32) {
        const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(lhs + i));
        const __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(rhs + i));
        const unsigned eq = static_cast<unsigned>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)));
        mm += 32 - __builtin_popcount(eq);
        if (mm > limit) return mm;
    }
    if (i + 16 <= len) {
        const __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(lhs + i));
        const __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i*>(rhs + i));
        const unsigned eq = static_cast<unsigned>(_mm_movemask_epi8(_mm_cmpeq_epi8(a, b)));
        mm += 16 - __builtin_popcount(eq);
        if (mm > limit) return mm;
        i += 16;
    }
    if (i + 8 <= len) {
        const __m128i a = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(lhs + i));
        const __m128i b = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(rhs + i));
        const unsigned eq = static_cast<unsigned>(_mm_movemask_epi8(_mm_cmpeq_epi8(a, b))) & 0xffu;
        mm += 8 - __builtin_popcount(eq);
        if (mm > limit) return mm;
        i += 8;
    }
#endif
    for (; i < len; ++i)
        if (lhs[i] != rhs[i] && ++mm > limit) break;
    return mm;
}

// Find first occurrence of adapter (first 12 bp) in seq at position ≥ min_pos, ≤ max_mm mismatches.
static int find_adapter_in(const std::string& seq, const std::string& adapter,
                           int min_pos, int max_mm = 2, int alen_req = 12) {
    if (adapter.empty()) return -1;
    int alen = (int)std::min(adapter.size(), (size_t)alen_req);
    int L    = (int)seq.size();
    for (int i = min_pos; i <= L - alen; ++i) {
        const int mm = count_mismatches_limited(seq.data() + i, adapter.data(), alen, max_mm);
        if (mm <= max_mm) return i;
    }
    return -1;
}

static void accum_overlap_subst(OverlapSubstCounts& cnt,
                                const std::string& r1_seq,
                                const std::string& rc2_seq,
                                int r2_offset, int best_ov, int skip_terminal) {
    ++cnt.n_pairs;
    const int lb = LenBins::bin(best_ov);
    ++cnt.bin_n_pairs[lb];
    for (int i = skip_terminal; i < best_ov; ++i) {
        int b1 = BASE_IDX[(uint8_t)r1_seq[i]];
        int b2 = BASE_IDX[(uint8_t)rc2_seq[r2_offset + i]];
        if (b1 < 0 || b2 < 0) continue;
        ++cnt.n_bases;
        ++cnt.bin_n_bases[lb];
        cnt.all[b1][b2]++;
        if (i < OverlapSubstCounts::MAX_POS) {
            cnt.fwd[i][b1][b2]++;
            cnt.fwd_len[lb][i][b1][b2]++;
        }
        int rp = best_ov - 1 - i;
        if (rp < OverlapSubstCounts::MAX_POS) {
            cnt.rev[rp][b1][b2]++;
            cnt.rev_len[lb][rp][b1][b2]++;
        }
    }
}

// Apply explicit tail quality cutoff after mate consensus has been considered.
static void trim_quality_tail(FastqRecord& rec, int cutoff,
                              double max_expected_errors = -1, double max_error_rate = -1) {
    if (max_expected_errors >= 0 || max_error_rate >= 0) {
        const auto keep = fqdup::longest_error_budget_prefix(
            rec.seq, rec.qual, max_expected_errors, max_error_rate);
        rec.seq.resize(keep); rec.qual.resize(keep);
        return;
    }
    if (cutoff < 0) return;
    // A high-quality terminal base anchors a biological tail even when the
    // preceding scores are poor. Without mate support, quality alone cannot
    // establish that those preceding bases are errors. Remove only the
    // consecutive low-quality terminal run.
    size_t keep = rec.qual.size();
    while (keep) {
        // An unknown base cannot anchor a quality tail, regardless of its
        // reported score. Match consensus N-neutrality; cutoff zero stays a no-op.
        const int quality = rec.seq[keep - 1] == 'N' ? 0 :
            static_cast<unsigned char>(rec.qual[keep - 1]) - 33;
        if (quality >= cutoff) break;
        --keep;
    }
    rec.seq.resize(keep);
    rec.qual.resize(keep);
}

static bool passes_qc(const FastqRecord& rec, const MergeOpts& opts,
                      float complexity_entropy_lo, float complexity_dom_hi) {
    if ((int)rec.seq.size() < opts.min_length) return false;
    if (opts.max_length > 0 && (int)rec.seq.size() > opts.max_length) return false;
    if (opts.max_n_rate < 1.0f) {
        int ns = (int)std::count(rec.seq.begin(), rec.seq.end(), 'N');
        if ((float)ns / (float)rec.seq.size() > opts.max_n_rate) return false;
    }
    if (opts.min_entropy > 0.f && shannon_entropy(rec.seq) < opts.min_entropy) return false;
    if (!opts.complexity_filter) return true;
    // Self-contained SDUST-style low-complexity gate (windowed triplet redundancy, fixed
    // threshold, no calibration). Catches homopolymer + di-/tri-nucleotide tandem repeats
    // (incl. partial tracts) that base composition misses. These reads map across the whole
    // tree with zero taxonomic signal and dominate downstream mapping cost, so drop them at QC.
    if (dust_worst(rec.seq) > DUST_MAX) return false;
    // Data-derived low-complexity gate (merged-mate reference). One-sided: reject only the
    // low-entropy / high-dominance extreme (poly-G run-through, low-complexity tracts).
    if (complexity_entropy_lo > 0.f || complexity_dom_hi < 1.0f) {
        auto w = worst_window(rec.seq, 0, (int)rec.seq.size());
        if (w.entropy < complexity_entropy_lo) return false;
        if (w.max_frac > complexity_dom_hi)   return false;
    }
    return true;
}

// Full adapter/artifact clip for an unmerged read (no overlap geometry to lean on).
// Order matters: the whole-read read-through scan runs first so an all-adapter dimer
// (adapter at pos 0) collapses to empty, and a linker+adapter dimer collapses to just
// the linker, which the 5' clips then strip below min-length.
//   adapter_3p = read-through adapter in this read's orientation (R1:adapter1, R2:adapter2)
//   adapter_5p = adapter-complement bleedthrough at 5'      (R1:adapter2, R2:adapter1)
//   linker5p   = detected library 5' construct (may be empty)
// Diagnostic (behind FQDUP_CLIP_DEBUG): per-mate attribution of which clip step took a
// read below min_length. Step legend printed at exit. Zero cost when dbg==nullptr.
//   0 wholeread-cut@detected-adapter  1 wholeread-cut@universal  2 linker-5p
//   3 adapter5p-5p  4 adapter3p-5p  5 universal-5p  6 polyG-3p  7 already-short (pre-clip)
static std::atomic<int64_t> g_clip_death[2][8];
static const bool g_clip_dbg = std::getenv("FQDUP_CLIP_DEBUG") != nullptr;
struct ClipDebug { int mate; int min_len; bool dead=false; };
static inline void clip_track(ClipDebug* d, int step, size_t before, size_t after) {
    if (!d || d->dead) return;
    if ((int)before >= d->min_len && (int)after < d->min_len) {
        g_clip_death[d->mate][step].fetch_add(1, std::memory_order_relaxed);
        d->dead = true;
    }
}

static void clip_unmerged(FastqRecord& rec, const std::string& adapter_3p,
                          const std::string& adapter_5p, const std::string& linker5p,
                          int poly_g_run, int min_len, ClipDebug* dbg = nullptr) {
    // Universal Illumina read-through prefix — present at the insert boundary of EVERY
    // read-through/adapter-dimer regardless of what auto-detection picked for adapter1/2
    // (a library can carry TruSeq read-through while detection locks onto a Nextera stub).
    // Protocol constant, not a calibration.
    static const std::string UNIVERSAL_3P = fqdup::adapters::truseq_3p;
    // Whole-read scan: cut at the earliest read-through hit (detected adapter OR universal).
    // A pure adapter dimer hits at pos 0 -> empty; linker+adapter hits right after the
    // linker -> only the linker survives, which the 5' clips then strip below min-length.
    if (dbg && (int)rec.seq.size() < min_len) {   // already short before any clip
        g_clip_death[dbg->mate][7].fetch_add(1, std::memory_order_relaxed); dbg->dead = true;
    }
    int cut = -1, cut_src = -1;
    int pd = adapter_3p.empty() ? -1 : find_adapter_in(rec.seq, adapter_3p,   0, 2, 12);
    int pu =                          find_adapter_in(rec.seq, UNIVERSAL_3P, 0, 2, 12);
    if (pd >= 0)                    { cut = pd; cut_src = 0; }
    if (pu >= 0 && (cut < 0 || pu < cut)) { cut = pu; cut_src = 1; }
    size_t sz0 = rec.seq.size();
    if (cut >= 0) { rec.seq.resize(cut); rec.qual.resize(cut); }
    clip_track(dbg, cut_src < 0 ? 1 : cut_src, sz0, rec.seq.size());
    size_t s;
    s = rec.seq.size(); if (!linker5p.empty())   trim_construct_5p(rec, linker5p, 0.1f, min_len); clip_track(dbg, 2, s, rec.seq.size());
    s = rec.seq.size(); if (!adapter_5p.empty()) trim_adapter_5p(rec, adapter_5p, min_len);       clip_track(dbg, 3, s, rec.seq.size());
    s = rec.seq.size(); if (!adapter_3p.empty()) trim_adapter_5p(rec, adapter_3p, min_len);       clip_track(dbg, 4, s, rec.seq.size());  // dimer: read-through at 5'
    s = rec.seq.size(); trim_adapter_5p(rec, UNIVERSAL_3P, min_len);                              clip_track(dbg, 5, s, rec.seq.size());  // dimer: universal at 5'
    s = rec.seq.size(); if (poly_g_run > 0)      trim_polybase(rec, 'G', poly_g_run, min_len);    clip_track(dbg, 6, s, rec.seq.size());
}

// Apply the same conservative terminal policy as trim before output QC on either
// consensus path. Pair geometry/consensus are already established here.
static void clean_merged_terminals(FastqRecord& rec,const MergeOpts& opts) {
    static const fqdup::CompiledTerminalCatalogue fallback({
        {"truseq_5p_12", fqdup::adapters::truseq_5p, true},
        {"truseq_3p_12", fqdup::adapters::truseq_3p, false}
    });
    const auto& catalogue=opts.terminal_model?*opts.terminal_model:fallback;
    auto five = catalogue.evaluate(rec.seq,true);
    auto three = catalogue.evaluate(rec.seq,false);
    if (five.bases + three.bases > rec.seq.size()) return;
    if (five.bases || three.bases)
        fqdup::apply_retained_interval(rec.seq, rec.qual, five.bases, rec.seq.size()-three.bases);
}

// Separately supported mate boundaries; this never constructs a gapped consensus.
static std::pair<int, int> find_gapped_clip_boundaries(
        const FastqRecord& r1, const FastqRecord& r2, const MergeOpts& opts) {
    auto anchors = [&](const std::string& read, const std::string& adapter) {
        std::vector<int> positions;
        const int required = std::min(opts.min_ov, static_cast<int>(adapter.size()));
        if (required <= 0) return positions;
        for (int n = opts.min_ov; n + required <= static_cast<int>(read.size()); ++n) {
            const int observed = std::min(static_cast<int>(adapter.size()),
                                          static_cast<int>(read.size()) - n);
            const int allowed = static_cast<int>(std::floor(opts.max_mm_rate * observed));
            if (count_mismatches_limited(read.data() + n, adapter.data(), observed, allowed) <= allowed)
                positions.push_back(n);
        }
        return positions;
    };
    const auto one = anchors(r1.seq, opts.adapter1), two = anchors(r2.seq, opts.adapter2);
    if (one.empty() || two.empty()) return {-1, -1};
    std::string rc2, rq2;
    rc_record(r2, rc2, rq2);
    auto compatible = [&](const char* a, int m, const char* b, int n) {
        const int limit = static_cast<int>(std::floor(opts.max_mm_rate * std::min(m, n)));
        if (std::abs(m - n) > limit) return false;
        // Unknown bases are not evidence for a separately positioned clipping action.
        for (int i = 0; i < m; ++i) if (a[i] != 'A' && a[i] != 'C' && a[i] != 'G' && a[i] != 'T') return false;
        for (int i = 0; i < n; ++i) if (b[i] != 'A' && b[i] != 'C' && b[i] != 'G' && b[i] != 'T') return false;
        const int inf = limit + 1;
        std::vector<int> prev(n + 1, inf), curr(n + 1, inf);
        for (int j = 0; j <= std::min(n, limit); ++j) prev[j] = j;
        for (int i = 1; i <= m; ++i) {
            std::fill(curr.begin(), curr.end(), inf);
            if (i <= limit) curr[0] = i;
            int row_min = curr[0];
            for (int j = std::max(1, i - limit); j <= std::min(n, i + limit); ++j) {
                curr[j] = std::min({prev[j] + 1, curr[j - 1] + 1,
                                    prev[j - 1] + (a[i - 1] != b[j - 1])});
                row_min = std::min(row_min, curr[j]);
            }
            if (row_min > limit) return false;
            prev.swap(curr);
        }
        return prev[n] <= limit;
    };
    // Equal full reads can be biological adapter mimics; never replace that
    // unresolved full-molecule interpretation with a clipping hypothesis.
    if (compatible(r1.seq.data(), static_cast<int>(r1.seq.size()),
                   rc2.data(), static_cast<int>(rc2.size()))) return {-2, -2};
    std::pair<int, int> selected{-1, -1};
    for (int n1 : one) for (int n2 : two) {
        if (!compatible(r1.seq.data(), n1, rc2.data() + rc2.size() - n2, n2)) continue;
        if (selected.first >= 0) return {-1, -1};
        selected = {n1, n2};
    }
    return selected;
}

// Adaptive, one-contiguous-gap geometry for separately clipped mates. No
// gapped consensus is built. Adapter hypotheses and gap placements spend one
// conservative sequence-evidence budget; only a unique boundary pair acts.
static std::pair<int, int> find_adaptive_gapped_clip_boundaries(
        const FastqRecord& r1, const FastqRecord& r2, const MergeOpts& opts) {
    const int l1 = static_cast<int>(r1.seq.size()), l2 = static_cast<int>(r2.seq.size());
    const std::size_t extent = static_cast<std::size_t>(std::max(l1, l2)) + 1;
    const std::size_t trials = 4 * extent * extent * extent;
    auto anchors = [&](const FastqRecord& read, const std::string& adapter) {
        std::vector<int> positions;
        if (adapter.empty()) return positions;
        for (int n = opts.min_ov; n + opts.min_ov <= static_cast<int>(read.seq.size()); ++n) {
            const int observed = std::min(static_cast<int>(adapter.size()),
                                          static_cast<int>(read.seq.size()) - n);
            if (observed < opts.min_ov) continue;
            if (fqdup::adaptive_alignment_supported(
                    read.seq.data() + n, read.qual.data() + n,
                    adapter.data(), nullptr, observed, observed, trials, false))
                positions.push_back(n);
        }
        return positions;
    };
    const auto one = anchors(r1, opts.adapter1), two = anchors(r2, opts.adapter2);
    if (one.empty() || two.empty()) return {-1, -1};
    std::string rc2, rq2;
    rc_record(r2, rc2, rq2);
    // A credible full-read biological interpretation defeats a shorter
    // adapter-like internal boundary, including one mate disagreement.
    if (l1 == l2 && (r1.seq == rc2 || fqdup::adaptive_alignment_supported(
            r1.seq.data(), r1.qual.data(), rc2.data(), rq2.data(),
            l1, opts.min_ov, trials, true))) return {-2, -2};
    const int last_anchor = std::max(*std::max_element(one.begin(), one.end()),
                                     *std::max_element(two.begin(), two.end()));
    if (fqdup::full_read_gap_alternative(r1.seq, r1.qual, rc2, rq2,
                                         opts.min_ov, trials, last_anchor)) return {-2, -2};
    std::pair<int, int> selected{-1, -1};
    for (int n1 : one) for (int n2 : two) {
        if (n1 == n2) continue;  // Existing adaptive path handles no gap.
        if (!fqdup::one_run_gap_supported(
                std::string_view(r1.seq.data(), n1),
                std::string_view(r1.qual.data(), n1),
                std::string_view(rc2.data() + l2 - n2, n2),
                std::string_view(rq2.data() + l2 - n2, n2),
                opts.min_ov, trials)) continue;
        if (selected.first >= 0) return {-1, -1};
        selected = {n1, n2};
    }
    return selected;
}

static void merge_worker(PairQueue& in_q, MergeOutQueue& out_q,
                         const MergeOpts& opts,
                         taph::SampleDamageProfile* prof_out,
                         OverlapSubstCounts* subst_out) {
    fqalign::QualityProcessor quality_processor;
    std::string rc2_seq, rc2_qual;
    taph::SampleDamageProfile local_prof;
    OverlapSubstCounts local_subst;
    const bool do_profile = (prof_out != nullptr);
    const bool do_subst   = (subst_out != nullptr);

    // Emit an unmerged pair with full adapter/linker/poly-G clipping and QC.
    // Zero-insert dimers collapse below min-length and are dropped; one survivor → orphan.
    auto emit_unmerged = [&](FastqRecord& u1, FastqRecord& u2, MergeBatch& out,
                             bool trim = true, int paired_boundary = -1, bool allow_gapped = false) {
        ClipDebug d1{0, opts.min_length}, d2{1, opts.min_length};
        if (opts.pair_supported_boundaries || paired_boundary >= 0) {
            std::pair<int, int> separate_bounds{-1, -1};
            // When both declared adapters and both observed mates are fully
            // periodic, a one-run-gap scan has many equivalent placements.
            // No single gap geometry can justify a clip, so avoid the costly
            // search and retain the pair-supported ungapped decision.
            auto periodic = [](const std::string& sequence) {
                if (sequence.size() < 8) return false;
                for (std::size_t period = 1; period <= 4; ++period) {
                    bool repeats = true;
                    for (std::size_t i = period; i < sequence.size(); ++i)
                        if (sequence[i] != sequence[i % period]) { repeats = false; break; }
                    if (repeats) return true;
                }
                return false;
            };
            auto two_periodic_runs = [](const std::string& sequence) {
                if (sequence.size() < 24) return false;
                for (std::size_t first_period = 1; first_period <= 4; ++first_period) {
                    std::size_t prefix = first_period;
                    while (prefix < sequence.size() &&
                           sequence[prefix] == sequence[prefix % first_period]) ++prefix;
                    if (prefix < 16) continue;
                    for (std::size_t last_period = 1; last_period <= 4; ++last_period) {
                        std::size_t start = sequence.size() - last_period;
                        while (start > 0 &&
                               sequence[start - 1] == sequence[start - 1 + last_period]) --start;
                        if (sequence.size() - start >= 8 && start <= prefix) return true;
                    }
                }
                return false;
            };
            const bool periodic_ambiguity = opts.experimental_adaptive_overlap &&
                periodic(opts.adapter1) && periodic(opts.adapter2) &&
                two_periodic_runs(u1.seq) && two_periodic_runs(u2.seq);
            if (opts.experimental_gapped_clipping && !periodic_ambiguity) {
                separate_bounds = opts.experimental_adaptive_overlap
                    ? find_adaptive_gapped_clip_boundaries(u1, u2, opts)
                    : find_gapped_clip_boundaries(u1, u2, opts);
                if (separate_bounds.first == -2) paired_boundary = -1;
                if (opts.experimental_adaptive_overlap &&
                    separate_bounds.first >= 0 && paired_boundary >= 0 &&
                    (separate_bounds.first != paired_boundary ||
                     separate_bounds.second != paired_boundary)) {
                    // Two different supported boundaries are unresolved.
                    paired_boundary = -1;
                    separate_bounds = {-1, -1};
                }
            }
            // Never fall back to an independent single-read adapter-motif decision.
            // Clip adapters only when both mates found one identical boundary.
            trim = false;
            if (paired_boundary >= 0 &&
                paired_boundary <= static_cast<int>(u1.seq.size()) &&
                paired_boundary <= static_cast<int>(u2.seq.size())) {
                u1.seq.resize(paired_boundary); u1.qual.resize(paired_boundary);
                u2.seq.resize(paired_boundary); u2.qual.resize(paired_boundary);
            }
            if (allow_gapped && paired_boundary == -1) {
                const auto bounds = separate_bounds;
                if (bounds.first >= 0) {
                    u1.seq.resize(bounds.first); u1.qual.resize(bounds.first);
                    u2.seq.resize(bounds.second); u2.qual.resize(bounds.second);
                }
            }
            // A terminal G run is ambiguous in one read: it may be real biology.
            // Treat it as a two-color instrument artifact only when both mates show
            // the same cycle boundary, then retain each biological prefix.
            auto polyg_boundary = [&](const FastqRecord& rec) {
                if (opts.poly_g_min_run <= 0) return -1;
                int i = static_cast<int>(rec.seq.size());
                while (i > 0 && rec.seq[i - 1] == 'G') --i;
                const int run = static_cast<int>(rec.seq.size()) - i;
                return run >= opts.poly_g_min_run && i >= opts.min_length ? i : -1;
            };
            const int g1 = polyg_boundary(u1), g2 = polyg_boundary(u2);
            // A real G immediately before the failing cycle extends the apparent
            // run. Allow only the uncertainty implied by half the learned run
            // threshold. Retain a sub-threshold guard band at the later boundary:
            // this removes the technical run while never claiming ambiguous Gs as
            // artifact evidence.
            const int guard = opts.poly_g_min_run - 1;
            // The retained guard is also the maximum biologically plausible shift in
            // the observed run start. Using a narrower slop creates an unstable gap:
            // a genuine G immediately before the same technical failure cycle makes
            // one mate abstain even though the retained guard already protects it.
            const int boundary_slop = std::max(1, guard);
            if (g1 >= 0 && g2 >= 0 && std::abs(g1 - g2) <= boundary_slop) {
                const bool learned = opts.paired_poly_g_cycle >= std::max(g1, g2) &&
                    opts.paired_poly_g_cycle < static_cast<int>(u1.seq.size()) &&
                    opts.paired_poly_g_cycle < static_cast<int>(u2.seq.size());
                const int boundary = learned ? opts.paired_poly_g_cycle : std::max(g1, g2) + guard;
                if (boundary <= static_cast<int>(u1.seq.size()) &&
                    boundary <= static_cast<int>(u2.seq.size())) {
                    size_t s1 = u1.seq.size(), s2 = u2.seq.size();
                    u1.seq.resize(boundary); u1.qual.resize(boundary);
                    u2.seq.resize(boundary); u2.qual.resize(boundary);
                    clip_track(g_clip_dbg ? &d1 : nullptr, 6, s1, u1.seq.size());
                    clip_track(g_clip_dbg ? &d2 : nullptr, 6, s2, u2.seq.size());
                }
            }
        }
        if (trim) {
            clip_unmerged(u1, opts.adapter1, opts.adapter2, opts.adapter_5p_linker,
                          opts.poly_g_min_run, opts.min_length, g_clip_dbg ? &d1 : nullptr);
            clip_unmerged(u2, opts.adapter2, opts.adapter1, opts.adapter_5p_linker,
                          opts.poly_g_min_run, opts.min_length, g_clip_dbg ? &d2 : nullptr);
        }
        trim_quality_tail(u1, opts.trim_quality, opts.qc_max_expected_errors, opts.qc_max_error_rate);
        trim_quality_tail(u2, opts.trim_quality, opts.qc_max_expected_errors, opts.qc_max_error_rate);
        // Per-mate QC: low-complexity gate against that mate's merged reference, plus a
        // 5' adapter-dimer reject (learned read-through adapter shifted match clip missed).
        bool ok1 = passes_qc(u1, opts, opts.complexity_entropy_lo_r1, opts.complexity_dom_hi_r1)
                   && (!trim || (!is_adapter_dimer_5p(u1.seq, opts.adapter1)
                   && !is_adapter_dimer_5p_indel(u1.seq, opts.adapter1)));
        bool ok2 = passes_qc(u2, opts, opts.complexity_entropy_lo_r2, opts.complexity_dom_hi_r2)
                   && (!trim || (!is_adapter_dimer_5p(u2.seq, opts.adapter2)
                   && !is_adapter_dimer_5p_indel(u2.seq, opts.adapter2)));
        MergeRecord mr;
        mr.is_merged = false;
        if (ok1 && ok2) {
            mr.merged = std::move(u1); mr.unmerged2 = std::move(u2); ++out.n_unmerged;
        } else if (ok1 || ok2) {
            mr.is_orphan = true; mr.orphan_r1 = ok1;
            mr.merged = ok1 ? std::move(u1) : std::move(u2); ++out.n_orphan;
            if (ok1) ++out.n_orphan_r1; else ++out.n_orphan_r2;
        } else {
            ++out.n_dropped;
            return;  // both fail QC → drop
        }
        out.records.push_back(std::move(mr));
    };

    PairBatch batch;
    while (in_q.pop(batch)) {
        MergeBatch out;
        out.id = batch.id;
        out.records.reserve(batch.pairs.size());

        for (auto& pr : batch.pairs) {
            PairOutcomeGuard outcome_guard(out, pr, !opts.outcomes_out.empty());
            NativeQualityGuard<MergeBatch,FastqRecord> quality_guard(out,pr.r1,pr.r2,quality_processor,opts.experimental_gapped_mott,opts.experimental_disagreement_risk,opts.experimental_selective_disagreement,opts.experimental_pair_mixture.get(),opts.experimental_pair_mixture_mode,opts.experimental_full_pair_mixture.get(),opts.experimental_merged_disagreement_floor);
            FastqRecord& r1 = pr.r1;
            FastqRecord& r2 = pr.r2;

            // Poly-G run-through is trimmed only off the unmerged mate outputs
            // (in clip_unmerged); merged reads are left as-is — the insert-only
            // consensus already excludes the 3' adapter/run-through region.

            if (opts.clip_5p > 0 && (int)r1.seq.size() > opts.clip_5p) {
                r1.seq  = r1.seq.substr(opts.clip_5p);
                r1.qual = r1.qual.substr(opts.clip_5p);
            }

            int L1 = (int)r1.seq.size();
            int L2 = (int)r2.seq.size();

            // Conservative pre-merge homopolymer reject: when BOTH raw mates are
            // near-degenerate the insert is junk regardless of merge outcome, so drop
            // it before spending overlap-detection work (speed + smaller output). The
            // learned merged-read gate downstream handles borderline low-complexity.
            if (opts.complexity_filter && is_homopolymer_read(r1.seq) && is_homopolymer_read(r2.seq)) {
                ++out.n_dropped;
                continue;
            }

            if (L1 < opts.min_ov || L2 < opts.min_ov ||
                L1 > MAX_READ_LEN  || L2 > MAX_READ_LEN) {
                // Too short or too long to overlap: clip adapters and emit as unmerged
                emit_unmerged(r1, r2, out);
                continue;
            }

            rc_record(r2, rc2_seq, rc2_qual);

            int L_cmp = std::min(L1, L2);
            // Below the normal overlap minimum, two complete exact adapter
            // anchors and an exact canonical paired body can still locate a
            // short insert. Retain the two trimmed mates; a nine-base body
            // alone is too little sequence evidence for a consensus merge.
            if ((opts.experimental_adaptive_overlap ||
                 (opts.experimental_short_tail && opts.experimental_gapped_clipping)) &&
                !opts.adapter1.empty() && !opts.adapter2.empty()) {
                int short_boundary = -1;
                for (int n = std::max(1, opts.min_length); n < opts.min_ov && n <= L_cmp; ++n) {
                    if (L1 - n < static_cast<int>(opts.adapter1.size()) ||
                        L2 - n < static_cast<int>(opts.adapter2.size()) ||
                        r1.seq.compare(n, opts.adapter1.size(), opts.adapter1) != 0 ||
                        r2.seq.compare(n, opts.adapter2.size(), opts.adapter2) != 0)
                        continue;
                    bool body_exact = true;
                    for (int j = 0; j < n && body_exact; ++j) {
                        const char a = r1.seq[j];
                        body_exact = (a == 'A' || a == 'C' || a == 'G' || a == 'T') &&
                                     a == rc2_seq[L2 - n + j];
                    }
                    if (!body_exact) continue;
                    if (short_boundary == -1) short_boundary = n;
                    else { short_boundary = -2; break; }
                }
                if (short_boundary >= 0 && !(L1 == L2 && r1.seq == rc2_seq)) {
                    emit_unmerged(r1, r2, out, false, short_boundary);
                    continue;
                }
            }

            // Exact paired geometry with two short terminal anchors rescues
            // near-read-length inserts. A four-base tail on each mate is too
            // short to accept independently; the full exact paired body and
            // uniqueness make it useful as joint clipping evidence. Do not
            // turn this into a consensus merge or override a full biological
            // mate alignment.
            const bool combined_damaged_short = !opts.experimental_adaptive_overlap &&
                opts.experimental_short_tail && opts.experimental_gapped_clipping;
            int deferred_exact_short = -1;
            if (((opts.experimental_adaptive_overlap && opts.experimental_short_tail) ||
                 combined_damaged_short) &&
                !opts.adapter1.empty() && !opts.adapter2.empty() &&
                !(L1 == L2 && r1.seq == rc2_seq)) {
                const int near_lo = std::max(opts.min_ov, L_cmp - opts.min_ov + 1);
                int exact_short = -1;
                bool exact_anchors = false;
                for (int n = near_lo; n < L_cmp; ++n) {
                    const int tail1 = std::min(L1 - n, static_cast<int>(opts.adapter1.size()));
                    const int tail2 = std::min(L2 - n, static_cast<int>(opts.adapter2.size()));
                    if (combined_damaged_short) {
                        // One damaged base per visible 5-8 bp tail requires the
                        // full exact paired body and an unambiguous boundary.
                        if (tail1 < 5 || tail1 > 8 || tail2 < 5 || tail2 > 8) continue;
                    } else if (tail1 < 4 || tail2 < 4) continue;
                    bool exact = true;
                    int errors1 = 0, errors2 = 0;
                    for (int j = 0; j < tail1; ++j)
                        errors1 += r1.seq[n + j] != opts.adapter1[j];
                    for (int j = 0; j < tail2; ++j)
                        errors2 += r2.seq[n + j] != opts.adapter2[j];
                    const bool both_exact = errors1 == 0 && errors2 == 0;
                    exact = combined_damaged_short ?
                        ((errors1 == 1 && errors2 == 1) || both_exact) : both_exact;
                    for (int j = 0; j < n && exact; ++j) {
                        const char a = r1.seq[j];
                        exact = (a == 65 || a == 67 || a == 71 || a == 84) &&
                                a == rc2_seq[L2 - n + j];
                    }
                    if (!exact) continue;
                    if (exact_short == -1) {
                        exact_short = n;
                        exact_anchors = both_exact;
                    }
                    else { exact_short = -2; break; }
                }
                bool full_exact_alternative = false;
                if (exact_short >= 0)
                    for (int n = near_lo; n <= L_cmp && !full_exact_alternative; ++n) {
                        bool exact = true;
                        for (int j = 0; j < n && exact; ++j) {
                            const char a = r1.seq[L1 - n + j];
                            exact = (a == 65 || a == 67 || a == 71 || a == 84) &&
                                    a == rc2_seq[j];
                        }
                        full_exact_alternative = exact;
                    }
                if (exact_short >= 0 && !full_exact_alternative) {
                    if (combined_damaged_short && exact_anchors)
                        deferred_exact_short = exact_short;
                    else {
                        emit_unmerged(r1, r2, out, false, exact_short);
                        continue;
                    }
                }
            }

            // Phase 0: adapter-supported boundary search.
            // Geometry: R1=[insert][adapter1], RC(R2)=[RC(adapter2)][insert].
            // When adapter found at p1 in R1: insert_len=p1, insert in RC(R2) is at [L2-p1:L2].
            // Cross-validate with adapter2 in RC(R2): appears at L2-p1-alen2 (near start of RC(R2)).
            int best_ov   = -1;
            int r2_offset = 0;
            int paired_boundary = -1;
            int adaptive_tail = -1;
            if (opts.experimental_adaptive_overlap) {
                // Exhaustive proposal and acceptance use one model. Fixed-rate
                // legacy searches below are bypassed for this opt-in path.
                const std::size_t trials = static_cast<std::size_t>(4 * L_cmp + 1);
                auto supports = [&](const FastqRecord& read,
                                    const std::string& adapter, int n, bool partial) {
                    const int observed = std::min(static_cast<int>(adapter.size()),
                                                  static_cast<int>(read.seq.size()) - n);
                    const int required = std::min(static_cast<int>(adapter.size()), opts.min_ov);
                    if (adapter.empty() || observed <= 0 || (!partial && observed < required))
                        return false;
                    return fqdup::adaptive_alignment_supported(
                        read.seq.data() + n, read.qual.data() + n,
                        adapter.data(), nullptr, observed, observed, trials, false);
                };
                auto geometry = [&](int first, int second, int length) {
                    return fqdup::adaptive_alignment_supported(
                        r1.seq.data() + first, r1.qual.data() + first,
                        rc2_seq.data() + second, rc2_qual.data() + second,
                        length, opts.min_ov, trials, true);
                };
                int prefix = -1;
                for (int n = opts.min_ov; n <= L_cmp; ++n) {
                    const bool full1 = supports(r1, opts.adapter1, n, false);
                    const bool full2 = supports(r2, opts.adapter2, n, false);
                    // Two concordant adapter observations can locate tails
                    // without a consensus merge only with positive body
                    // evidence. Exact canonical bodies can still be retained
                    // separately when their repetition makes a merge unsafe.
                    // Mere quality compatibility at a low-Q accidental
                    // alignment is insufficient to remove biological bases.
                    bool body_exact = full1 && full2;
                    if (body_exact)
                        for (int j = 0; j < n && body_exact; ++j) {
                            const char a = r1.seq[j], b = rc2_seq[L2 - n + j];
                            body_exact = (a == 'A' || a == 'C' || a == 'G' || a == 'T') && a == b;
                        }
                    const bool paired_body = full1 && full2 &&
                        (body_exact || geometry(0, L2 - n, n));
                    if (paired_body) {
                        if (paired_boundary == -1) paired_boundary = n;
                        else paired_boundary = -2;
                    }
                    const bool short_both = opts.experimental_short_tail &&
                        supports(r1, opts.adapter1, n, true) &&
                        supports(r2, opts.adapter2, n, true);
                    if (!(full1 || full2 || short_both) || !geometry(0, L2 - n, n)) continue;
                    if (prefix >= 0) { prefix = -2; break; }
                    prefix = n;
                }
                for (int n = opts.min_ov; n <= L_cmp; ++n) {
                    if (!geometry(L1 - n, 0, n)) continue;
                    if (adaptive_tail >= 0) { adaptive_tail = -2; break; }
                    adaptive_tail = n;
                }
                // A credible full-read mate alignment is an alternative to
                // treating a shorter internal adapter-like motif as technical.
                // Preserve both reads under that ambiguity, including a full
                // biological overlap with a sequencing disagreement.
                if (paired_boundary >= 0 && L1 == L2 &&
                    (r1.seq == rc2_seq || geometry(0, 0, L1)))
                    paired_boundary = -1;
                // An exact near-full suffix/prefix mate overlap is an
                // alternative biological geometry even in periodic sequence.
                // The sequence-diversity floor correctly refuses to *merge*
                // such reads, but must not be used to justify clipping them.
                bool full_exact_tail_alternative = false;
                if (paired_boundary >= 0) {
                    const int near_lo = std::max(opts.min_ov, L_cmp - opts.min_ov + 1);
                    for (int n = near_lo; n <= L_cmp; ++n) {
                        const int offset = L1 - n;
                        bool exact = true;
                        for (int j = 0; j < n; ++j) {
                            const char a = r1.seq[offset + j], b = rc2_seq[j];
                            if ((a != 'A' && a != 'C' && a != 'G' && a != 'T') || a != b) {
                                exact = false;
                                break;
                            }
                        }
                        if (exact) { full_exact_tail_alternative = true; break; }
                    }
                }
                if (full_exact_tail_alternative) paired_boundary = -1;
                // A near-full insert with two exact short adapter tails and
                // an exact paired body competes with a shorter, low-quality
                // full-anchor boundary. Preserve both mates under ambiguity.
                bool exact_short_boundary_alternative = false;
                if (paired_boundary >= 0) {
                    const int near_lo = std::max(paired_boundary + 1,
                        std::max(opts.min_ov, L_cmp - opts.min_ov + 1));
                    for (int n = near_lo; n < L_cmp; ++n) {
                        const int tail1 = std::min(L1 - n, static_cast<int>(opts.adapter1.size()));
                        const int tail2 = std::min(L2 - n, static_cast<int>(opts.adapter2.size()));
                        if (tail1 < 4 || tail2 < 4) continue;
                        bool exact = true;
                        for (int j = 0; j < tail1 && exact; ++j)
                            exact = r1.seq[n + j] == opts.adapter1[j];
                        for (int j = 0; j < tail2 && exact; ++j)
                            exact = r2.seq[n + j] == opts.adapter2[j];
                        for (int j = 0; j < n && exact; ++j)
                            exact = r1.seq[j] == rc2_seq[L2 - n + j];
                        if (exact) { exact_short_boundary_alternative = true; break; }
                    }
                }
                if (exact_short_boundary_alternative) paired_boundary = -1;
                // A sequencing indel followed by terminal truncation can make
                // full biological mates look like shorter, adapter-bearing
                // inserts. Preserve the reads under that competing geometry.
                const std::size_t full_gap_trials =
                    4 * static_cast<std::size_t>(L_cmp + 1) *
                    static_cast<std::size_t>(L_cmp + 1) *
                    static_cast<std::size_t>(L_cmp + 1);
                const bool full_gap_alternative = (prefix >= 0 || adaptive_tail >= 0 ||
                    paired_boundary >= 0) &&
                    fqdup::full_read_gap_alternative(r1.seq, r1.qual, rc2_seq, rc2_qual,
                                                     opts.min_ov, full_gap_trials,
                                                     std::max(prefix, paired_boundary) >= 0
                                                         ? std::max(prefix, paired_boundary)
                                                         : L_cmp);
                if (full_gap_alternative) paired_boundary = -1;
                bool conflict = prefix == -2 || adaptive_tail == -2 ||
                                (prefix >= 0 && adaptive_tail >= 0) ||
                                (paired_boundary >= 0 && adaptive_tail >= 0 &&
                                 paired_boundary != adaptive_tail) ||
                                (paired_boundary >= 0 && prefix >= 0 &&
                                 paired_boundary != prefix) ||
                                full_gap_alternative || full_exact_tail_alternative ||
                                exact_short_boundary_alternative;
                // A near-full unanchored prefix remains an alternative to a
                // longer-insert tail/head interpretation or partial-anchor cut.
                if (!conflict && (adaptive_tail >= 0 || (prefix >= 0 && opts.experimental_short_tail))) {
                    const int near_lo = std::max(opts.min_ov, L_cmp - opts.min_ov + 1);
                    for (int n = near_lo; n <= L_cmp; ++n)
                        if (n != prefix && geometry(0, L2 - n, n)) { conflict = true; break; }
                }
                if (conflict) {
                    emit_unmerged(r1, r2, out, false, paired_boundary);
                    continue;
                }
                if (prefix >= 0) { best_ov = prefix; r2_offset = L2 - prefix; }
            }
            if (!opts.experimental_adaptive_overlap && opts.pair_supported_boundaries) {
                auto adapter_supports = [&](const std::string& read,
                                            const std::string& adapter, int n) {
                    if (n < 0 || n >= static_cast<int>(read.size()) || adapter.empty()) return false;
                    const int observed = std::min(static_cast<int>(adapter.size()),
                                                  static_cast<int>(read.size()) - n);
                    const int required = std::min(static_cast<int>(adapter.size()), opts.min_ov);
                    if (observed < required) return false;
                    const int allowed = static_cast<int>(std::ceil(opts.max_mm_rate * observed));
                    const int mm = count_mismatches_limited(
                        read.data() + n, adapter.data(), observed, allowed);
                    return mm <= allowed;
                };
                for (int n = opts.min_ov; n <= L_cmp; ++n) {
                    const bool supports1 = adapter_supports(r1.seq, opts.adapter1, n);
                    const bool supports2 = adapter_supports(r2.seq, opts.adapter2, n);
                    if (supports1 && supports2) {
                        if (paired_boundary == -1)
                            paired_boundary = n;
                        else
                            paired_boundary = -2;
                    }
                    if (!supports1 && !supports2) continue;
                    int offset = L2 - n, informative = 0;
                    float mismatches = 0;
                    for (int j = opts.skip_terminal; j < n; ++j) {
                        if (r1.seq[j] == 'N' || rc2_seq[offset + j] == 'N') continue;
                        ++informative;
                        if (r1.seq[j] != rc2_seq[offset + j])
                            mismatches += std::min(static_cast<int>(r1.qual[j]) - 33,
                                                   static_cast<int>(rc2_qual[offset + j]) - 33) / 30.f;
                    }
                    if (informative < opts.min_ov ||
                        mismatches > opts.max_mm_rate * informative) continue;
                    // More than one mate-supported candidate is unresolved.
                    // Preserve original mates rather than choosing a shorter cut.
                    if (best_ov >= opts.min_ov || best_ov == -2) {
                        best_ov = -2;
                    } else {
                        best_ov = n;
                        r2_offset = offset;
                    }
                }
                if (best_ov == -1 && paired_boundary >= opts.min_ov) {
                    // Two complete exact adapters identify a short-insert boundary
                    // even when two biological read errors exceed the fixed mismatch
                    // rate. Require independent sequence evidence for the body and
                    // reject a competing longer-insert tail/head geometry.
                    const int n = paired_boundary;
                    const bool exact_anchors = n <= L_cmp &&
                        L1 - n >= static_cast<int>(opts.adapter1.size()) &&
                        L2 - n >= static_cast<int>(opts.adapter2.size()) &&
                        r1.seq.compare(n, opts.adapter1.size(), opts.adapter1) == 0 &&
                        r2.seq.compare(n, opts.adapter2.size(), opts.adapter2) == 0;
                    if (exact_anchors && !(L1 == L2 && r1.seq == rc2_seq)) {
                        int informative = 0, mismatches = 0;
                        for (int j = opts.skip_terminal; j < n; ++j) {
                            const char x = r1.seq[j], y = rc2_seq[L2 - n + j];
                            if (x == 'N' || y == 'N') continue;
                            ++informative;
                            mismatches += x != y;
                        }
                        if (informative >= opts.min_ov) {
                            const double fraction = static_cast<double>(informative - mismatches) / informative;
                            double divergence = 0.0;
                            if (fraction > 0.25) {
                                divergence = fraction * std::log(fraction / 0.25);
                                if (fraction < 1.0)
                                    divergence += (1.0 - fraction) *
                                        std::log((1.0 - fraction) / 0.75);
                            }
                            if (informative * divergence >= std::log(1e6) &&
                                find_tail_head_ov(r1.seq.data(), L1, rc2_seq.data(), L2,
                                                  opts.min_ov, opts.max_mm_rate, true) == -1) {
                                best_ov = n;
                                r2_offset = L2 - n;
                            }
                        }
                    }
                }
                if (best_ov == -2) {
                    emit_unmerged(r1, r2, out, false,
                                  paired_boundary >= 0 ? paired_boundary : deferred_exact_short);
                    continue;
                }
            } else if (!opts.experimental_adaptive_overlap) {
                int p1    = find_adapter_in(r1.seq,  opts.adapter1, opts.min_ov);
                // adapter2 in RC(R2) is at position L2 - insert_len - 12; search from 0 toward min_ov
                int p2raw = find_adapter_in(rc2_seq, opts.adapter2, 0);
                int ov2   = (p2raw >= 0) ? (L2 - p2raw - 12) : -1;

                if (p1 >= opts.min_ov && ov2 >= opts.min_ov) {
                    // Both found: they should agree. Use the average or the consistent one.
                    best_ov   = (std::abs(p1 - ov2) <= 2) ? p1 : std::min(p1, ov2);
                    r2_offset = L2 - best_ov;
                } else if (p1 >= opts.min_ov) {
                    best_ov   = p1;
                    r2_offset = L2 - p1;
                } else if (ov2 >= opts.min_ov) {
                    best_ov   = ov2;
                    r2_offset = L2 - ov2;
                }
                // Bounds check: r2_offset must be non-negative
                if (best_ov >= opts.min_ov && r2_offset < 0) best_ov = -1;
            }

            // Phase 0-extra: try additional adapter sequences (--adapter-fasta with multiple pairs).
            if (!opts.experimental_adaptive_overlap && best_ov < opts.min_ov && !opts.extra_adapters1.empty()) {
                for (const auto& ea1 : opts.extra_adapters1) {
                    int pe = find_adapter_in(r1.seq, ea1, opts.min_ov, 2, 12);
                    if (pe >= opts.min_ov && pe <= L_cmp) {
                        int r2s = L2 - pe;
                        if (r2s >= 0) {
                            float wmm = 0.f;
                            int   eff = std::max(1, pe - opts.skip_terminal);
                            for (int i = opts.skip_terminal; i < pe && r2s+i < L2; ++i) {
                                if (r1.seq[i] != rc2_seq[r2s+i] &&
                                    r1.seq[i] != 'N' && rc2_seq[r2s+i] != 'N') {
                                    int q1 = (uint8_t)r1.qual[i]       - 33;
                                    int q2 = (uint8_t)rc2_qual[r2s+i]  - 33;
                                    wmm += std::min(q1, q2) / 30.f;
                                }
                            }
                            if (wmm / (float)eff <= opts.max_mm_rate) {
                                best_ov = pe; r2_offset = r2s; break;
                            }
                        }
                    }
                }
            }

            // Phase 0b: progressively shorter adapter anchors for near-read-length inserts.
            // 8bp ≤1mm: catches inserts 90-93bp; 6bp exact: catches 94-95bp.
            // Validation uses the CORRECT RC(R2) slice [L2-cand:L2].
            if (!opts.experimental_adaptive_overlap && best_ov < opts.min_ov && !opts.adapter1.empty()) {
                const std::pair<int,int> anchors[] = {{8,1},{6,0}};
                for (auto [alen, max_mm_b] : anchors) {
                    int p1b = find_adapter_in(r1.seq, opts.adapter1, opts.min_ov, max_mm_b, alen);
                    if (p1b >= opts.min_ov && p1b <= L_cmp) {
                        int r2s = L2 - p1b;
                        if (r2s >= 0) {
                            float wmm = 0.f;
                            int   eff = std::max(1, p1b - opts.skip_terminal);
                            for (int i = opts.skip_terminal; i < p1b; ++i) {
                                if (r1.seq[i] != rc2_seq[r2s + i] &&
                                    r1.seq[i] != 'N' && rc2_seq[r2s + i] != 'N') {
                                    int q1 = (uint8_t)r1.qual[i]        - 33;
                                    int q2 = (uint8_t)rc2_qual[r2s + i] - 33;
                                    wmm += std::min(q1, q2) / 30.f;
                                }
                            }
                            float thresh = opts.max_mm_rate * (alen >= 8 ? 1.0f : 0.6f);
                            if (wmm / (float)eff <= thresh) {
                                best_ov   = p1b;
                                r2_offset = r2s;
                                break;
                            }
                        }
                    }
                }
            }

            // Phase 0c: adapter tail scan — catches insert 96-100bp (only 1-5bp of
            // adapter visible at 3' end of R1). Scans from the longest matching tail
            // downward; strict validation threshold (half max_mm_rate) suppresses FP.
            if (!opts.experimental_adaptive_overlap && best_ov < opts.min_ov && !opts.adapter1.empty()) {
                int cand = find_adapter_tail(r1.seq, opts.adapter1, 5, opts.min_ov);
                if (opts.pair_supported_boundaries && cand >= opts.min_ov) {
                    const int observed = std::min(L2 - cand,
                                                  static_cast<int>(opts.adapter2.size()));
                    const int required = std::min(static_cast<int>(opts.adapter2.size()),
                                                  opts.min_ov);
                    const int allowed = static_cast<int>(std::ceil(opts.max_mm_rate * observed));
                    int mm = 0;
                    for (int j = 0; j < observed && mm <= allowed; ++j)
                        mm += r2.seq[cand + j] != opts.adapter2[j];
                    if (observed < required || mm > allowed) cand = -1;
                }
                if (cand >= opts.min_ov) {
                    int r2s = L2 - cand;
                    if (r2s >= 0) {
                        float wmm = 0.f;
                        int   eff = std::max(1, cand - opts.skip_terminal);
                        for (int i = opts.skip_terminal; i < cand && r2s + i < L2; ++i) {
                            if (r1.seq[i] != rc2_seq[r2s + i] &&
                                r1.seq[i] != 'N' && rc2_seq[r2s + i] != 'N') {
                                int q1 = (uint8_t)r1.qual[i]        - 33;
                                int q2 = (uint8_t)rc2_qual[r2s + i] - 33;
                                wmm += std::min(q1, q2) / 30.f;
                            }
                        }
                        if (wmm / (float)eff <= opts.max_mm_rate * 0.5f) {
                            best_ov   = cand;
                            r2_offset = r2s;
                        }
                    }
                }
            }

            // Phase 1: quality-weighted Hamming at d=0, retained for legacy
            // mode only. Prefix-prefix equality cannot distinguish an insert exactly
            // one read long from two periodic, non-overlapping insert ends.
            if (best_ov < opts.min_ov && !opts.pair_supported_boundaries) {
                // Use 0.10 threshold here (fastp default) — the adapter-free path
                // has no anchor bias; quality weighting already suppresses noise.
                float ph1_rate = std::max(opts.max_mm_rate, 0.10f);
                best_ov = find_best_ov_qwt(
                    r1.seq.data(), r1.qual.data(),
                    rc2_seq.data(), rc2_qual.data(),
                    L_cmp, opts.min_ov, ph1_rate, opts.skip_terminal);
                // r2_offset = 0 (already set above)
            }

            // Experimental partial paired anchors. A short adapter tail alone is
            // weak evidence: require both mates and reject every competing
            // prefix/suffix or tail/head geometry, including full-read equality.
            if (!opts.experimental_adaptive_overlap && opts.experimental_short_tail && opts.pair_supported_boundaries &&
                best_ov < opts.min_ov) {
                int candidate = -1;
                auto partial_anchor = [&](const std::string& read,
                                          const std::string& adapter, int n) {
                    const int observed = std::min(static_cast<int>(adapter.size()),
                                                  static_cast<int>(read.size()) - n);
                    if (observed <= 0 || observed >= opts.min_ov) return false;
                    const int allowed = static_cast<int>(std::floor(opts.max_mm_rate * observed));
                    return count_mismatches_limited(read.data() + n, adapter.data(),
                                                    observed, allowed) <= allowed;
                };
                auto geometry_passes = [&](int n) {
                    int informative = 0, mm = 0;
                    const int offset = L2 - n;
                    for (int j = opts.skip_terminal; j < n; ++j) {
                        const char a = r1.seq[j], b = rc2_seq[offset + j];
                        if (a == 'N' || b == 'N') continue;
                        ++informative;
                        mm += a != b;
                    }
                    return informative >= opts.min_ov &&
                           static_cast<float>(mm) <= opts.max_mm_rate * informative;
                };
                const int near_lo = std::max(opts.min_ov, L_cmp - opts.min_ov + 1);
                for (int n = near_lo; n < L_cmp; ++n) {
                    if (!partial_anchor(r1.seq, opts.adapter1, n) ||
                        !partial_anchor(r2.seq, opts.adapter2, n) || !geometry_passes(n)) continue;
                    if (candidate >= 0) { candidate = -2; break; }
                    candidate = n;
                }
                if (candidate >= opts.min_ov) {
                    bool conflict = false;
                    for (int n = opts.min_ov; n <= L_cmp && !conflict; ++n)
                        if (n != candidate && geometry_passes(n)) conflict = true;
                    if (!conflict && find_tail_head_ov(r1.seq.data(), L1,
                                                      rc2_seq.data(), L2,
                                                      opts.min_ov, opts.max_mm_rate, false) < 0) {
                        best_ov = candidate;
                        r2_offset = L2 - candidate;
                    }
                }
            }

            // Phase 2: tail-head scan for inserts > read_length.
            // R1[L1-ov : L1] overlaps RC(R2)[0 : ov].
            // Full merged = R1[0:L1-ov] + consensus(overlap) + RC(R2)[ov:L2].
            //
            // A repeat-rich molecule at or just below one read length can also produce a
            // unique-looking internal tail/head overlap. Before accepting that longer-insert
            // interpretation, check the otherwise-unanchorable near-boundary geometries
            // (adapter tail shorter than min_ov). If one agrees, sequence alone admits two
            // insert lengths; accuracy-first mode must retain the original mates.
            bool near_boundary_conflict = false;
            if (!opts.experimental_adaptive_overlap && best_ov < opts.min_ov && opts.pair_supported_boundaries) {
                const int near_lo = std::max(opts.min_ov, L_cmp - opts.min_ov + 1);
                for (int n = near_lo; n <= L_cmp && !near_boundary_conflict; ++n) {
                    const int offset = L2 - n;
                    int informative = 0, mm = 0;
                    for (int j = opts.skip_terminal; j < n; ++j) {
                        const char a = r1.seq[j], b = rc2_seq[offset + j];
                        if (a == "N"[0] || b == "N"[0]) continue;
                        ++informative;
                        mm += a != b;
                    }
                    near_boundary_conflict = informative >= opts.min_ov &&
                        static_cast<float>(mm) <= opts.max_mm_rate * informative;
                }
            }
            if (best_ov < opts.min_ov) {
                // Scan from the declared minimum so a shorter passing offset
                // still vetoes a competing longer geometry. A unique unanchored
                // overlap below 20 bases is too weak to construct a molecule.
                best_ov = opts.experimental_adaptive_overlap ? adaptive_tail :
                    find_tail_head_ov(r1.seq.data(), L1,
                                            rc2_seq.data(), L2,
                                            opts.min_ov, opts.max_mm_rate,
                                            opts.pair_supported_boundaries);
                if (!opts.experimental_adaptive_overlap && opts.pair_supported_boundaries &&
                    best_ov >= opts.min_ov && best_ov < 20)
                    best_ov = -1;
                if (best_ov >= opts.min_ov && near_boundary_conflict) {
                    emit_unmerged(r1, r2, out, false, paired_boundary);
                    continue;
                }
                if (best_ov == -2) {
                    emit_unmerged(r1, r2, out, false,
                                  paired_boundary >= 0 ? paired_boundary : deferred_exact_short);
                    continue;
                }
                if (best_ov >= opts.min_ov) {
                    // A short unanchored overlap with multiple confident disagreements
                    // is weak evidence for a shared molecule. Longer overlaps can
                    // tolerate the same number of errors because they contain more
                    // concordant bases. Keep this test local to the tail/head path.
                    int informative = 0, confident_mm = 0;
                    for (int j = 0; j < best_ov; ++j) {
                        const int i1 = L1 - best_ov + j;
                        const char a = r1.seq[i1], b = rc2_seq[j];
                        if (a == 'N' || b == 'N') continue;
                        ++informative;
                        if (a != b && r1.qual[i1] - 33 >= 25 &&
                            rc2_qual[j] - 33 >= 25) ++confident_mm;
                    }
                    if (opts.pair_supported_boundaries && confident_mm >= 2 &&
                        25 * confident_mm > informative) {
                        emit_unmerged(r1, r2, out, false, paired_boundary);
                        continue;
                    }
                    if (opts.pair_supported_boundaries &&
                        (has_exact_long_self_overlap(r1.seq.data(), L1) ||
                         has_exact_long_self_overlap(rc2_seq.data(), L2))) {
                        emit_unmerged(r1, r2, out, false, paired_boundary);
                        continue;
                    }
                    if (do_profile)
                        taph::FrameSelector::update_sample_profile_paired(local_prof, r1.seq, r2.seq);
                    // Phase 2 overlap is at R1[L1-best_ov:] vs RC(R2)[0:best_ov] — interior
                    // of the insert, not the 5'/3' termini. Per-position fwd/rev semantics
                    // would be wrong here, so skip subst accumulation for Phase 2 reads.
                    FastqRecord fake_r1;
                    fake_r1.header = r1.header;
                    fake_r1.plus   = "+";
                    fake_r1.seq    = r1.seq.substr(L1 - best_ov, best_ov);
                    fake_r1.qual   = r1.qual.substr(L1 - best_ov, best_ov);
                    FastqRecord ov_cons;
                    consensus_merge(fake_r1,
                                    rc2_seq.substr(0, best_ov),
                                    rc2_qual.substr(0, best_ov),
                                    best_ov, 0, ov_cons);
                    MergeRecord mr;
                    mr.is_merged = true;
                    // Build full merged: prefix + consensus overlap + suffix
                    mr.merged.header = ov_cons.header;
                    mr.merged.plus   = "+";
                    mr.merged.seq    = r1.seq.substr(0, L1 - best_ov)
                                     + ov_cons.seq
                                     + rc2_seq.substr(best_ov);
                    mr.merged.qual   = r1.qual.substr(0, L1 - best_ov)
                                     + ov_cons.qual
                                     + rc2_qual.substr(best_ov);
                    if (!opts.pair_supported_boundaries && !opts.adapter2.empty()) trim_adapter_5p(mr.merged, opts.adapter2, opts.min_length);
                    if (!opts.pair_supported_boundaries && !opts.adapter1.empty()) trim_adapter(mr.merged, opts.adapter1, opts.min_length);
                    // Adapter dimer on a MERGED read (real aDNA inserts are [insert][3'adapter],
                    // never 5'adapter): the whole read is adapter/construct that overlap-merged
                    // into pure adapter — the shifted/mismatch match the strict 5' trim missed.
                    // is_technical_read (Myers, indel-tolerant) drops it. NO low-complexity gate
                    // on merged reads: the gate is derived FROM the merged inserts it would gate
                    // (circular — rejects ~1-2% genuine inserts by construction); overlap
                    // verification IS the complexity check. passes_qc(0.f,1.0f) keeps only the
                    // absolute floors (min_length / max_n_rate / min_entropy).
                    if (opts.known_terminal_constructs) clean_merged_terminals(mr.merged,opts);
                    trim_quality_tail(mr.merged, opts.trim_quality,
                                      opts.qc_max_expected_errors, opts.qc_max_error_rate);
                    if (is_technical_read(mr.merged.seq, opts.tech_seqs) ||
                        (opts.complexity_filter && is_homopolymer_read(mr.merged.seq))) {
                        ++out.n_dropped; note_frag_drop(mr.merged.seq);
                    } else if (passes_qc(mr.merged, opts, 0.f, 1.0f)) {
                        ++out.n_merged;
                        out.records.push_back(std::move(mr));
                    } else {
                        emit_unmerged(r1, r2, out, true, paired_boundary);
                    }
                    continue;
                }
            }

            if (best_ov < opts.min_ov) {
                // No overlap found: clip adapters/linker/poly-G, QC, handle orphans/dimers
                if (deferred_exact_short >= 0 && paired_boundary < 0)
                    emit_unmerged(r1, r2, out, false, deferred_exact_short);
                else
                    emit_unmerged(r1, r2, out, true, paired_boundary,
                                  opts.experimental_gapped_clipping);
                continue;
            }

            // Quality weights may rank consensus evidence, but Q0 cannot make
            // disagreeing sequences establish a molecule boundary. Apply the
            // sequence evidence floor to every anchored prefix merge.
            // The tail/head path above already uses sequence-only geometry.
            if (!opts.experimental_adaptive_overlap && opts.pair_supported_boundaries) {
                int informative = 0, mismatches = 0;
                for (int j = opts.skip_terminal; j < best_ov; ++j) {
                    const char a = r1.seq[j], b = rc2_seq[r2_offset + j];
                    if (a == "N"[0] || b == "N"[0]) continue;
                    ++informative;
                    mismatches += a != b;
                }
                bool sequence_supported = informative >= opts.min_ov &&
                    static_cast<float>(mismatches) <= opts.max_mm_rate * informative;
                if (!sequence_supported && informative >= opts.min_ov &&
                    paired_boundary == best_ov) {
                    // Unique full anchors on BOTH mates can rescue a genuinely
                    // erroneous overlap. Still require quality-independent sequence
                    // evidence: a conservative Chernoff score against independent
                    // uniform random bases (match probability 1/4). This is an
                    // evidence floor, not a genomic false-discovery guarantee.
                    const double matches = static_cast<double>(informative - mismatches);
                    const double fraction = matches / informative;
                    if (fraction > 0.25) {
                        double divergence = fraction * std::log(fraction / 0.25);
                        if (fraction < 1.0)
                            divergence += (1.0 - fraction) *
                                std::log((1.0 - fraction) / 0.75);
                        sequence_supported = informative * divergence >= std::log(1e6);
                    }
                }
                if (!sequence_supported) {
                    // Rejected geometry supplies no valid shared clipping boundary.
                    emit_unmerged(r1, r2, out, false, -1);
                    continue;
                }
            }

            // Merge via Bayesian consensus
            if (do_profile)
                taph::FrameSelector::update_sample_profile_paired(local_prof, r1.seq, r2.seq);
            MergeRecord mr;
            mr.is_merged = true;
            consensus_merge(r1, rc2_seq, rc2_qual, best_ov, r2_offset, mr.merged);
            // Trim adapter artifacts from merged output:
            // 5': adapter complement (CTCTTCCGATCT) from library-prep bleedthrough
            // 3': adapter sequence surviving in Phase-2 suffix or mismerges
            // ss library 5' linker precedes the insert on the merged 5' too — clip the full
            // learned construct before the adapter trims, else its tail corrupts the insert 5'.
            if (!opts.pair_supported_boundaries && !opts.adapter_5p_linker.empty())
                trim_construct_5p(mr.merged, opts.adapter_5p_linker, 0.1f, opts.min_length);
            if (!opts.pair_supported_boundaries && !opts.adapter2.empty()) trim_adapter_5p(mr.merged, opts.adapter2, opts.min_length);
            if (!opts.pair_supported_boundaries && !opts.adapter1.empty()) trim_adapter(mr.merged, opts.adapter1, opts.min_length);

            if (opts.known_terminal_constructs) clean_merged_terminals(mr.merged,opts);
            trim_quality_tail(mr.merged, opts.trim_quality,
                              opts.qc_max_expected_errors, opts.qc_max_error_rate);

            // Adapter dimer on the merged read (5' adapter, shifted/mismatch match the
            // strict trim missed) → no real insert; drop. Dominant bucket-explosion lineage.
            if (is_technical_read(mr.merged.seq, opts.tech_seqs) ||
                (opts.complexity_filter && is_homopolymer_read(mr.merged.seq))) {
                ++out.n_dropped; note_frag_drop(mr.merged.seq);
            // Homopolymer inserts already dropped above; the data-derived complexity gate stays
            // off here (circular — see above), so passes_qc applies absolute floors only.
            } else if (!passes_qc(mr.merged, opts, 0.f, 1.0f)) {
                emit_unmerged(r1, r2, out, true, paired_boundary);
            } else {
                if (do_subst)
                    accum_overlap_subst(local_subst, r1.seq, rc2_seq, r2_offset, best_ov, 0);
                ++out.n_merged;
                out.records.push_back(std::move(mr));
            }
        }

        out_q.push(std::move(out));
    }
    if (do_profile || do_subst) {
        static std::mutex out_mtx;
        std::lock_guard<std::mutex> lk(out_mtx);
        if (do_profile) taph::FrameSelector::merge_sample_profiles(*prof_out, local_prof);
        if (do_subst)   subst_out->merge(local_subst);
    }
}


} // namespace fqdup_native_xc_internal
#undef fqdup
#undef fqalign
#undef fqquality
#undef fqrepair
#include "native.hpp"
#include "memory_status.hpp"
namespace fqdup::r8qkeep {
static NativeResult merge_native_chunk(const std::vector<FamilyPair>& pairs, const NativeSettings& settings,bool trace,bool audit) {
    using namespace ::fqdup_native_xc_internal;
    const auto memory=[&](const char* name){if(trace)std::cerr<<"MERGE_MEMORY native_"<<name<<memory_status()<<'\n';};
    memory("start");
    MergeOpts opts;
    opts.adapter1 = settings.adapter1;
    opts.adapter2 = settings.adapter2;
    opts.min_ov = settings.min_overlap;
    opts.min_length = settings.min_length;
    opts.max_mm_rate = settings.max_mismatch_rate;
    opts.forced_library_type = settings.library_type;
    opts.use_internal_panel = false;
    opts.disable_polyg = true;
    opts.experimental_short_tail = true;
    opts.experimental_gapped_clipping = true;
    opts.outcomes_out = "in-memory";
    PairBatch input; input.id = 0; input.pairs.reserve(pairs.size());
    for (const auto& pair : pairs) {
        ReadPair native;
        native.r1.header = pair.r[0].h; native.r1.seq = pair.r[0].s;
        native.r1.plus = "+"; native.r1.qual = pair.r[0].q;
        native.r2.header = pair.r[1].h; native.r2.seq = pair.r[1].s;
        native.r2.plus = "+"; native.r2.qual = pair.r[1].q;
        input.pairs.push_back(std::move(native));
    }
    memory("input_copied");
    PairQueue in(1); MergeOutQueue out(1);
    in.push(std::move(input)); in.set_done();
    merge_worker(in, out, opts, nullptr, nullptr);
    MergeBatch batch;
    if (!out.pop_ordered(0, batch)) throw std::runtime_error("native batch missing");
    memory("batch_ready");
    NativeResult result;result.records.reserve(pairs.size());if(audit)result.outcomes_jsonl = batch.outcomes_jsonl;
    auto append = [](std::string& text, const FastqRecord& r) {
        text += r.header + "\n" + r.seq + "\n" + r.plus + "\n" + r.qual + "\n";
    };
    if(audit)for (const auto& record : batch.records) {
        if (record.is_merged) { if (!record.merged.seq.empty()) append(result.merged_fq, record.merged); }
        else if (record.is_orphan) { if (!record.merged.seq.empty()) append(result.orphan_fq, record.merged); }
        else {
            if (!record.merged.seq.empty()) append(result.r1_fq, record.merged);
            if (!record.unmerged2.seq.empty()) append(result.r2_fq, record.unmerged2);
        }
    }
    std::istringstream outcomes(batch.outcomes_jsonl); std::string line;
    memory("parent_strings");
    std::size_t row = 0, emitted = 0;
    if(audit)result.native_outcomes_json = "[";
    while (std::getline(outcomes, line)) {
        if (row >= pairs.size()) throw std::runtime_error("extra native disposition");
        const auto marker = line.rfind(",\"outcome\":\"");
        if (marker == std::string::npos || line.size() < marker + 14) throw std::runtime_error("native disposition format");
        const auto begin = marker + 12;
        const std::string kind = line.substr(begin, line.find('"', begin) - begin);
        const auto& pair = pairs[row];
        NativeRecord native;native.raw=&pair;native.outcome=kind;
        if (kind != "dropped") {
            if (emitted >= batch.records.size()) throw std::runtime_error("missing native output");
            auto& record = batch.records[emitted++];
            native.first={std::move(record.merged.header),std::move(record.merged.seq),std::move(record.merged.qual)};
            if(kind=="unmerged")native.second={std::move(record.unmerged2.header),std::move(record.unmerged2.seq),std::move(record.unmerged2.qual)};
        }
        if(audit)result.native_tsv+=native_record_tsv(native);
        result.records.push_back(std::move(native));
        if(audit){if(row)result.native_outcomes_json += ',';result.native_outcomes_json += line;}++row;
    }
    if(row != pairs.size() || emitted != batch.records.size()) throw std::runtime_error("native disposition cardinality");
    if(audit)result.native_outcomes_json += "]\n";
    memory("serialized");
    return result;
}
NativeResult merge_native(const std::vector<FamilyPair>& pairs,const NativeSettings& settings,bool trace,bool audit){
    NativeResult result;result.records.reserve(pairs.size());if(audit)result.native_outcomes_json="[";
    bool first=true;
    for(std::size_t begin=0;begin<pairs.size();begin+=PAIR_BATCH_SZ){const auto end=std::min(pairs.size(),begin+PAIR_BATCH_SZ);
        std::vector<FamilyPair> chunk(pairs.begin()+begin,pairs.begin()+end);
        auto part=merge_native_chunk(chunk,settings,trace&&begin==0,audit);
        if(audit)result.native_tsv+=part.native_tsv;
        for(std::size_t i=0;i<part.records.size();++i){part.records[i].raw=&pairs[begin+i];result.records.push_back(std::move(part.records[i]));}
        if(audit){result.merged_fq+=part.merged_fq;result.r1_fq+=part.r1_fq;result.r2_fq+=part.r2_fq;result.orphan_fq+=part.orphan_fq;result.outcomes_jsonl+=part.outcomes_jsonl;
            if(!first)result.native_outcomes_json+=',';result.native_outcomes_json+=part.native_outcomes_json.substr(1,part.native_outcomes_json.size()-3);first=false;}
    }
    if(audit)result.native_outcomes_json+="]\n";
    if(trace)std::cerr<<"MERGE_MEMORY native_chunks_done"<<memory_status()<<'\n';
    return result;
}
}
