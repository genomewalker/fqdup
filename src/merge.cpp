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
    int64_t n_merged     = 0;
    int64_t n_unmerged   = 0;
    int64_t n_orphan     = 0;
    int64_t n_orphan_r1  = 0; // survivor = R1 (molecule 5', C->T frame)
    int64_t n_orphan_r2  = 0; // survivor = R2 (RC molecule 3'; ds biological, ss dA-tail)
    int64_t n_dropped    = 0; // both mates failed QC after clipping
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
    std::mutex                       mtx;
    std::condition_variable          cv;
    std::map<uint64_t, MergeBatch>   pending;
    bool                             done = false;

    void push(MergeBatch&& b) {
        std::unique_lock lk(mtx);
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
        return true;
    }
    void set_done() {
        std::unique_lock lk(mtx);
        done = true;
        cv.notify_all();
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

static int find_tail_head_ov(const char* r1, int L1,
                             const char* rc2, int L2,
                             int min_ov, float max_mm_rate) {
    int best = -1;
    for (int ov = min_ov; ov <= std::min(L1, L2); ++ov) {
        int mm = 0;
        const char* a = r1  + (L1 - ov);
        const char* b = rc2;
        for (int i = 0; i < ov; ++i)
            mm += (a[i] != b[i] && a[i] != 'N' && b[i] != 'N') ? 1 : 0;
        if ((float)mm / (float)ov <= max_mm_rate)
            best = ov;
        else if (best >= 0 && mm - (int)((float)best * max_mm_rate) >= 3)
            break;
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
    out.seq.resize(ov);
    out.qual.resize(ov);

    for (int i = 0; i < ov; ++i) {
        char b1 = r1.seq[i];
        char b2 = rc2_seq[r2_start + i];
        int  q1 = (uint8_t)r1.qual[i]              - 33;
        int  q2 = (uint8_t)rc2_qual[r2_start + i]  - 33;

        if (b1 == b2 || b2 == 'N') {
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
    rec.seq  = rec.seq.substr(alen);
    rec.qual = rec.qual.substr(alen);
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
    rec.seq  = rec.seq.substr(clen);
    rec.qual = rec.qual.substr(clen);
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

struct MergeOpts {
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
    float min_entropy   = 0.0f; // 0=disabled; Shannon entropy floor (bits, max=2.0)
    float max_n_rate    = 1.0f; // 1=disabled; max fraction of N bases
    // Data-derived low-complexity gate for unmerged mates (from merged-mate reference), per mate.
    float complexity_entropy_lo_r1 = 0.f,  complexity_dom_hi_r1 = 1.0f;
    float complexity_entropy_lo_r2 = 0.f,  complexity_dom_hi_r2 = 1.0f;
    std::string adapter1;
    std::string adapter2;
    std::string adapter_5p_linker;              // library 5' construct clipped from unmerged reads
    std::string forced_library_type;            // "ss"/"ds" declared override; ""=auto-detect
    std::string json_out;                       // --json: comprehensive lossless merge-QC report
    std::vector<std::string> extra_adapters1;   // additional R1 adapters to try (from --adapter-fasta)
    std::vector<std::string> tech_seqs;          // ALL learned technical constructs (multi-adapter QC)
    bool  use_internal_panel = true;            // aDNA construct read-through table; default ON, --no-internal-panel disables
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
        static const std::string TRUSEQ_R1  = "AGATCGGAAGAG";
        static const std::string TRUSEQ_RC2 = "CTCTTCCGATCT";
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
        // would drop the ds 3' G->A terminus (fable review, HIGH).
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
    // seq-error + composition floor instead of a panel-tuned 0.02/0.005/0.01/2.5 (fable review #8).
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
    static const std::string TRUSEQ_R1  = "AGATCGGAAGAG";   // appears in R1 at ov
    static const std::string TRUSEQ_RC2 = "CTCTTCCGATCT";   // appears in RC(R2) at ov

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

// Find first occurrence of adapter (first 12 bp) in seq at position ≥ min_pos, ≤ max_mm mismatches.
static int find_adapter_in(const std::string& seq, const std::string& adapter,
                           int min_pos, int max_mm = 2, int alen_req = 12) {
    if (adapter.empty()) return -1;
    int alen = (int)std::min(adapter.size(), (size_t)alen_req);
    int L    = (int)seq.size();
    for (int i = min_pos; i <= L - alen; ++i) {
        int mm = 0;
        for (int k = 0; k < alen; ++k)
            if (seq[i + k] != adapter[k] && ++mm > max_mm) break;
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

static bool passes_qc(const FastqRecord& rec, const MergeOpts& opts,
                      float complexity_entropy_lo, float complexity_dom_hi) {
    if ((int)rec.seq.size() < opts.min_length) return false;
    if (opts.max_length > 0 && (int)rec.seq.size() > opts.max_length) return false;
    if (opts.max_n_rate < 1.0f) {
        int ns = (int)std::count(rec.seq.begin(), rec.seq.end(), 'N');
        if ((float)ns / (float)rec.seq.size() > opts.max_n_rate) return false;
    }
    if (opts.min_entropy > 0.f && shannon_entropy(rec.seq) < opts.min_entropy) return false;
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
    static const std::string UNIVERSAL_3P = "AGATCGGAAGAG";
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

static void merge_worker(PairQueue& in_q, MergeOutQueue& out_q,
                         const MergeOpts& opts,
                         taph::SampleDamageProfile* prof_out,
                         OverlapSubstCounts* subst_out) {
    std::string rc2_seq, rc2_qual;
    taph::SampleDamageProfile local_prof;
    OverlapSubstCounts local_subst;
    const bool do_profile = (prof_out != nullptr);
    const bool do_subst   = (subst_out != nullptr);

    // Emit an unmerged pair with full adapter/linker/poly-G clipping and QC.
    // Zero-insert dimers collapse below min-length and are dropped; one survivor → orphan.
    auto emit_unmerged = [&](FastqRecord& u1, FastqRecord& u2, MergeBatch& out) {
        ClipDebug d1{0, opts.min_length}, d2{1, opts.min_length};
        clip_unmerged(u1, opts.adapter1, opts.adapter2, opts.adapter_5p_linker,
                      opts.poly_g_min_run, opts.min_length, g_clip_dbg ? &d1 : nullptr);
        clip_unmerged(u2, opts.adapter2, opts.adapter1, opts.adapter_5p_linker,
                      opts.poly_g_min_run, opts.min_length, g_clip_dbg ? &d2 : nullptr);
        // Per-mate QC: low-complexity gate against that mate's merged reference, plus a
        // 5' adapter-dimer reject (learned read-through adapter shifted match clip missed).
        bool ok1 = passes_qc(u1, opts, opts.complexity_entropy_lo_r1, opts.complexity_dom_hi_r1)
                   && !is_adapter_dimer_5p(u1.seq, opts.adapter1)
                   && !is_adapter_dimer_5p_indel(u1.seq, opts.adapter1);
        bool ok2 = passes_qc(u2, opts, opts.complexity_entropy_lo_r2, opts.complexity_dom_hi_r2)
                   && !is_adapter_dimer_5p(u2.seq, opts.adapter2)
                   && !is_adapter_dimer_5p_indel(u2.seq, opts.adapter2);
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
            if (is_homopolymer_read(r1.seq) && is_homopolymer_read(r2.seq)) {
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

            // Phase 0: adapter search in R1 (12bp anchor, ≤2 mm).
            // Geometry: R1=[insert][adapter1], RC(R2)=[RC(adapter2)][insert].
            // When adapter found at p1 in R1: insert_len=p1, insert in RC(R2) is at [L2-p1:L2].
            // Cross-validate with adapter2 in RC(R2): appears at L2-p1-alen2 (near start of RC(R2)).
            int best_ov   = -1;
            int r2_offset = 0;  // 0 = long-insert/d=0 path; L2-best_ov = short-insert path
            {
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
            if (best_ov < opts.min_ov && !opts.extra_adapters1.empty()) {
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
            if (best_ov < opts.min_ov && !opts.adapter1.empty()) {
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
            if (best_ov < opts.min_ov && !opts.adapter1.empty()) {
                int cand = find_adapter_tail(r1.seq, opts.adapter1, 5, opts.min_ov);
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

            // Phase 1: quality-weighted Hamming at d=0 (long-insert path).
            // For long inserts (insert ≥ read length), RC(R2) starts with the insert,
            // so d=0 comparison R1[0:] vs RC(R2)[0:] finds the overlap correctly.
            // r2_offset stays 0 for this path.
            if (best_ov < opts.min_ov) {
                // Use 0.10 threshold here (fastp default) — the adapter-free path
                // has no anchor bias; quality weighting already suppresses noise.
                float ph1_rate = std::max(opts.max_mm_rate, 0.10f);
                best_ov = find_best_ov_qwt(
                    r1.seq.data(), r1.qual.data(),
                    rc2_seq.data(), rc2_qual.data(),
                    L_cmp, opts.min_ov, ph1_rate, opts.skip_terminal);
                // r2_offset = 0 (already set above)
            }

            // Phase 2: tail-head scan for inserts > read_length.
            // R1[L1-ov : L1] overlaps RC(R2)[0 : ov].
            // Full merged = R1[0:L1-ov] + consensus(overlap) + RC(R2)[ov:L2].
            if (best_ov < opts.min_ov) {
                best_ov = find_tail_head_ov(r1.seq.data(), L1,
                                            rc2_seq.data(), L2,
                                            opts.min_ov, opts.max_mm_rate);
                if (best_ov >= opts.min_ov) {
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
                    if (!opts.adapter2.empty()) trim_adapter_5p(mr.merged, opts.adapter2, opts.min_length);
                    if (!opts.adapter1.empty()) trim_adapter(mr.merged, opts.adapter1, opts.min_length);
                    // Adapter dimer on a MERGED read (real aDNA inserts are [insert][3'adapter],
                    // never 5'adapter): the whole read is adapter/construct that overlap-merged
                    // into pure adapter — the shifted/mismatch match the strict 5' trim missed.
                    // is_technical_read (Myers, indel-tolerant) drops it. NO low-complexity gate
                    // on merged reads: the gate is derived FROM the merged inserts it would gate
                    // (circular — rejects ~1-2% genuine inserts by construction); overlap
                    // verification IS the complexity check. passes_qc(0.f,1.0f) keeps only the
                    // absolute floors (min_length / max_n_rate / min_entropy).
                    if (is_technical_read(mr.merged.seq, opts.tech_seqs)) {
                        ++out.n_dropped; note_frag_drop(mr.merged.seq);
                    } else if (passes_qc(mr.merged, opts, 0.f, 1.0f)) {
                        ++out.n_merged;
                        out.records.push_back(std::move(mr));
                    } else {
                        emit_unmerged(r1, r2, out);
                    }
                    continue;
                }
            }

            if (best_ov < opts.min_ov) {
                // No overlap found: clip adapters/linker/poly-G, QC, handle orphans/dimers
                emit_unmerged(r1, r2, out);
                continue;
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
            if (!opts.adapter_5p_linker.empty())
                trim_construct_5p(mr.merged, opts.adapter_5p_linker, 0.1f, opts.min_length);
            if (!opts.adapter2.empty()) trim_adapter_5p(mr.merged, opts.adapter2, opts.min_length);
            if (!opts.adapter1.empty()) trim_adapter(mr.merged, opts.adapter1, opts.min_length);

            // Adapter dimer on the merged read (5' adapter, shifted/mismatch match the
            // strict trim missed) → no real insert; drop. Dominant bucket-explosion lineage.
            if (is_technical_read(mr.merged.seq, opts.tech_seqs)) {
                ++out.n_dropped; note_frag_drop(mr.merged.seq);
            // No low-complexity gate on merged reads (circular — see above); absolute floors only.
            } else if (!passes_qc(mr.merged, opts, 0.f, 1.0f)) {
                emit_unmerged(r1, r2, out);
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

// ============================================================================
// Usage
// ============================================================================

static void usage() {
    std::cerr <<
        "Usage: fqdup merge -1 R1.fq.gz -2 R2.fq.gz -o merged.fq.gz [options]\n\n"
        "Detect read-pair overlap and merge into single collapsed reads.\n"
        "Optimized for ancient DNA: high damage-tolerance (--max-mm-rate 0.08 default),\n"
        "Bayesian quality consensus at overlap, incremental-Hamming overlap detection.\n\n"
        "Required:\n"
        "  -1 FILE            R1 (forward) reads (.fastq.gz or plain)\n"
        "  -2 FILE            R2 (reverse) reads (.fastq.gz or plain)\n"
        "  -o FILE            Output: merged reads (.fastq.gz)\n\n"
        "Optional output:\n"
        "  --r1-out FILE      Unmerged R1 reads (adapter-trimmed if --adapter1 given)\n"
        "  --r2-out FILE      Unmerged R2 reads (adapter-trimmed if --adapter2 given)\n"
        "  --orphan-out FILE  Reads where one mate was discarded by QC filters\n"
        "  --damage-out FILE  Strand-resolved paired damage profile JSON\n"
        "  --subst-out FILE   Per-position overlap substitution matrix TSV\n"
        "                     (encodes deamination, G->A, OxoG directly from R1 vs RC(R2))\n\n"
        "Overlap:\n"
        "  --min-overlap N    Minimum overlap length (default: 11)\n"
        "  --max-mm-rate F    Max mismatch rate in overlap (default: 0.08)\n"
        "  --min-length N     Discard reads shorter than N bp, all emit paths (default: 16; 0 = off)\n"
        "  --max-length N     Discard reads longer than N bp, all emit paths (default: off)\n"
        "  --json FILE        Write comprehensive lossless merge-QC report (JSON)\n"
        "  --clip-r1-5p N     Hard-clip N bases from R1 5' end before overlap (removes adapter stubs)\n"
        "  --min-entropy F    Discard low-complexity merged reads; Shannon entropy floor in bits\n"
        "                     (0=disabled; poly-G≈0, random≈2.0; default: 0)\n"
        "  --max-n-rate F     Discard merged reads with N fraction above F (default: 1.0=off)\n\n"
        "Adapter trimming (for unmerged pairs):\n"
        "  --adapter1 SEQ     R1 adapter sequence (Illumina P7 RC)\n"
        "  --adapter2 SEQ     R2 adapter sequence (Illumina P7)\n"
        "  --adapter-fasta F  FASTA with adapter pairs (odd=R1, even=R2); multiple pairs supported\n"
        "  --no-internal-panel  Disable built-in aDNA construct read-through trimming (default: ON;\n"
        "                       --internal-panel kept as a no-op for back-compat)\n\n"
        "Performance:\n"
        "  -p N               Threads (default: all cores)\n"
        "  -h, --help         Show this help\n\n"
        "Notes:\n"
        "  Overlap detection uses prefix Hamming (d=0 alignment), which is optimal\n"
        "  for aDNA where insert < read length (full read-through). Falls back to\n"
        "  3'/5' tail-head alignment for long inserts.\n"
        "  Output read name: R1 name with /1 or /2 suffix stripped.\n";
}

// ============================================================================
// Load adapters from FASTA file
// ============================================================================

// Loads adapters from FASTA. First sequence → a1 (R1 adapter), second → a2 (R2 adapter).
// Additional sequences at odd positions (3rd, 5th...) go into extra_adapters1 for Phase 0-extra.
// Pairs in FASTA: (R1_seq, R2_seq, R1_seq, R2_seq, ...) — only R1 seqs are used for searching.
static bool load_adapter_fasta(const std::string& path,
                               std::string& a1, std::string& a2,
                               std::vector<std::string>* extra) {
    std::ifstream f(path);
    if (!f.good()) return false;
    std::vector<std::string> seqs;
    std::string line;
    int which = -1;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') { ++which; seqs.emplace_back(); continue; }
        if (which >= 0) seqs.back() += line;
    }
    if (seqs.empty()) return false;
    a1 = seqs[0];
    // adapter2 in Phase 0 is searched in RC(R2) — it must be the RC of the R2 adapter fwd.
    // FASTA stores R2 adapter fwd; convert to the form that appears in RC(R2).
    if (seqs.size() >= 2 && !seqs[1].empty())
        a2 = revcomp(seqs[1].substr(0, std::min((int)seqs[1].size(), 12)));
    // Subsequent R1 adapters: seqs[2], seqs[4], ...
    if (extra) {
        for (size_t i = 2; i < seqs.size(); i += 2)
            if (!seqs[i].empty()) extra->push_back(seqs[i]);
    }
    return true;
}

// ============================================================================
// Entry point
// ============================================================================

int merge_main(int argc, char** argv) {
    std::vector<std::string> r1_paths, r2_paths;
    std::string out_path, r1_out_path, r2_out_path, orphan_out_path;
    std::string orphan_r1_out_path, orphan_r2_out_path;
    std::string adapter_fasta;
    MergeOpts opts;
    int n_threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      ((a == "-1")               && i+1 < argc) r1_paths.push_back(argv[++i]);
        else if ((a == "-2")               && i+1 < argc) r2_paths.push_back(argv[++i]);
        else if ((a == "-o")               && i+1 < argc) out_path     = argv[++i];
        else if (a == "--r1-out"           && i+1 < argc) r1_out_path      = argv[++i];
        else if (a == "--r2-out"           && i+1 < argc) r2_out_path      = argv[++i];
        else if (a == "--orphan-out"       && i+1 < argc) orphan_out_path  = argv[++i];
        else if (a == "--orphanr1-out"     && i+1 < argc) orphan_r1_out_path = argv[++i];
        else if (a == "--orphanr2-out"     && i+1 < argc) orphan_r2_out_path = argv[++i];
        else if (a == "--library-type"     && i+1 < argc) opts.forced_library_type = argv[++i];
        else if (a == "--damage-out"       && i+1 < argc) opts.damage_out  = argv[++i];
        else if (a == "--subst-out"        && i+1 < argc) opts.subst_out   = argv[++i];
        else if (a == "--subst-binary"     && i+1 < argc) opts.subst_binary = argv[++i];
        else if ((a == "-p" || a == "--threads") && i+1 < argc) n_threads = std::stoi(argv[++i]);
        else if (a == "--min-overlap"      && i+1 < argc) opts.min_ov       = std::stoi(argv[++i]);
        else if (a == "--max-mm-rate"      && i+1 < argc) opts.max_mm_rate  = std::stof(argv[++i]);
        else if (a == "--min-length"       && i+1 < argc) opts.min_length   = std::stoi(argv[++i]);
        else if (a == "--max-length"       && i+1 < argc) opts.max_length   = std::stoi(argv[++i]);
        else if (a == "--json"             && i+1 < argc) opts.json_out     = argv[++i];
        else if (a == "--adapter1"         && i+1 < argc) opts.adapter1       = argv[++i];
        else if (a == "--adapter2"         && i+1 < argc) opts.adapter2       = argv[++i];
        else if (a == "--adapter-fasta"    && i+1 < argc) adapter_fasta       = argv[++i];
        else if (a == "--internal-panel")                 opts.use_internal_panel = true;   // no-op back-compat (now default)
        else if (a == "--no-internal-panel")              opts.use_internal_panel = false;
        else if (a == "--clip-r1-5p"       && i+1 < argc) opts.clip_5p         = std::stoi(argv[++i]);
        else if (a == "--poly-g")                         opts.poly_g_min_run = 10;
        else if (a == "--poly-g-min-run"   && i+1 < argc) opts.poly_g_min_run = std::stoi(argv[++i]);
        else if (a == "--min-entropy"      && i+1 < argc) opts.min_entropy    = std::stof(argv[++i]);
        else if (a == "--max-n-rate"       && i+1 < argc) opts.max_n_rate     = std::stof(argv[++i]);
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else { std::cerr << "Error: unknown option '" << a << "'\n"; usage(); return 1; }
    }

    if (r1_paths.empty() || r2_paths.empty() || out_path.empty()) {
        std::cerr << "Error: -1, -2, and -o are required\n";
        usage();
        return 1;
    }

    if (opts.min_length < 0) {
        std::cerr << "merge: --min-length must be >= 0 (0 = disabled)\n";
        return 1;
    }
    if (opts.max_length != 0 && opts.max_length < opts.min_length) {
        std::cerr << "merge: --max-length must be >= --min-length (or 0 to disable)\n";
        return 1;
    }

    if (!adapter_fasta.empty()) {
        if (!load_adapter_fasta(adapter_fasta, opts.adapter1, opts.adapter2,
                                &opts.extra_adapters1)) {
            std::cerr << "Error: cannot read adapter FASTA: " << adapter_fasta << "\n";
            return 1;
        }
        std::cerr << "Adapters from " << adapter_fasta << ":\n"
                  << "  adapter1: " << opts.adapter1.substr(0, 20) << "...\n"
                  << "  adapter2: " << opts.adapter2.substr(0, 20) << "...\n";
        if (!opts.extra_adapters1.empty())
            std::cerr << "  extra adapters: " << opts.extra_adapters1.size() << " additional\n";
    }

    // Distribute decompression threads: half to R1, half to R2
    int io_threads  = std::max(1, n_threads / 2);
    int wrk_threads = std::max(1, n_threads - 2);  // 1 reader, 1 writer
    int wrt_threads = std::min(n_threads, 16);

    auto r1_rdr = make_chained_fastq_reader(r1_paths, static_cast<size_t>(io_threads));
    auto r2_rdr = make_chained_fastq_reader(r2_paths, static_cast<size_t>(io_threads));

    // ---- Pre-scan: auto-detect library type, adapters, UDG status ----
    static constexpr int64_t SCAN_READS = 200'000;
    std::vector<ReadPair> scan_buf;
    scan_buf.reserve(SCAN_READS);
    {
        ReadPair pr;
        while ((int64_t)scan_buf.size() < SCAN_READS &&
               r1_rdr->read(pr.r1) && r2_rdr->read(pr.r2))
            scan_buf.push_back(pr);
    }

    // Only auto-detect when user didn't override adapter/skip manually
    bool user_override = !opts.adapter1.empty() || !adapter_fasta.empty();
    DetectedParams det = detect_merge_params(scan_buf, opts.min_ov);

    if (!user_override) {
        opts.adapter1 = det.adapter1;
        opts.adapter2 = det.adapter2;
        opts.adapter_5p_linker = det.adapter_5p_linker;
        opts.skip_terminal = det.skip_terminal;
        // All learned technical constructs, in BOTH orientations, for the merged-read QC
        // (dimer / fragment / embedded checks scan against every construct, not just the top one).
        opts.tech_seqs = det.tech_seqs;
        for (const auto& t : det.tech_seqs) opts.tech_seqs.push_back(revcomp(t));
        // Internal aDNA construct panel (protocol constants; Ellesmere supp §4.3.2 / Kapp 2021 /
        // Illumina TruSeq): extra R1-side read-through constructs tried in phase-0-extra. Replicates
        // the published fastp --adapter_fasta panel. Gated OFF by default (--internal-panel);
        // validated 3'-terminus-neutral before enabling.
        if (opts.use_internal_panel && opts.extra_adapters1.empty()) {
            opts.extra_adapters1 = {
                "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",   // Illumina universal read-through
                "ATCTCGTATGCCGTCTTCTGCTTG",             // P7 index read-through
                "GTGTAGATCTCGGTGGTCGCCGTATCATT",        // P5 read-2 segment
            };
            opts.extra_adapters1.push_back(det.is_ss ? "GGAAGAGCGTCGTGTAGGGAAAGAGTGT"
                                                     : "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
        }
    } else if (opts.adapter1.empty() && !opts.adapter2.empty()) {
        opts.adapter1 = revcomp(opts.adapter2);
    } else if (opts.adapter2.empty() && !opts.adapter1.empty()) {
        opts.adapter2 = revcomp(opts.adapter1);
    }

    // Canonical QC (on by default, independent of adapter override): clean 3'
    // poly-G run-through off the unmerged mate outputs (merged reads are already
    // clean via overlap-clip). Data-derived run-length; a no-op on data without a
    // real G-run. CLI --poly-g[-min-run] overrides (leaves poly_g_min_run != 0).
    if (opts.poly_g_min_run == 0)
        opts.poly_g_min_run = derive_polyg_k(scan_buf);

    // Canonical QC: data-derived low-complexity gate on unmerged mates, learned per mate from
    // the merged-mate worst-window reference (no literal). Disabled (0/1.0) when the reference
    // is below floor. Independent of adapter override — always applied to the unmerged stream.
    opts.complexity_entropy_lo_r1 = det.complexity_entropy_lo_r1;
    opts.complexity_dom_hi_r1     = det.complexity_dom_hi_r1;
    opts.complexity_entropy_lo_r2 = det.complexity_entropy_lo_r2;
    opts.complexity_dom_hi_r2     = det.complexity_dom_hi_r2;

    // Resolve library type: an authoritative --library-type declaration always wins over
    // geometry auto-detection (memory: a named-"ss" library was proven ds by geometry+BIC).
    bool is_ss = det.is_ss;
    std::string lib_type_source = det.type_from_panel ? "panel" : "detected";
    float lib_type_conf = det.type_confidence;
    if (opts.forced_library_type == "ss")      { is_ss = true;  lib_type_source = "declared"; lib_type_conf = 1.f; }
    else if (opts.forced_library_type == "ds") { is_ss = false; lib_type_source = "declared"; lib_type_conf = 1.f; }
    // ds libraries carry no extra 5' linker; a linker learned on a ds library would risk
    // clipping a real biological 5', so gate the linker on ss (I1 / design §2.2).
    if (!is_ss) opts.adapter_5p_linker.clear();

    std::string lib_type = is_ss ? "ss" : "ds";
    std::string udg_type = det.is_udg ? "UDG" : (det.is_half_udg ? "half-UDG" : "untreated");
    std::cerr << "fqdup merge: " << n_threads << " threads (workers=" << wrk_threads << ")\n"
              << "  library=" << lib_type << " (" << lib_type_source
              << " conf=" << lib_type_conf << ") damage=" << udg_type
              << " skip-terminal=" << opts.skip_terminal << "\n"
              << "  damage_5p=" << det.damage_5p << " damage_3p=" << det.damage_3p << "\n"
              << "  adapter1=" << (opts.adapter1.empty() ? "(none)" : opts.adapter1) << "\n"
              << "  adapter2=" << (opts.adapter2.empty() ? "(none)" : opts.adapter2) << "\n"
              << "  5p-linker=" << (opts.adapter_5p_linker.empty() ? "(none)" : opts.adapter_5p_linker) << "\n"
              << "  min-overlap=" << opts.min_ov
              << " max-mm-rate=" << opts.max_mm_rate
              << " min-length=" << opts.min_length;
    if (opts.max_length > 0)
        std::cerr << " max-length=" << opts.max_length;
    if (opts.clip_5p > 0)
        std::cerr << " clip-r1-5p=" << opts.clip_5p;
    if (opts.poly_g_min_run > 0)
        std::cerr << " poly-g=" << opts.poly_g_min_run;
    if (opts.min_entropy > 0.f)
        std::cerr << " min-entropy=" << opts.min_entropy;
    std::cerr << "\n  complexity-gate R1[n=" << det.complexity_ref_n_r1
              << " ent>=" << opts.complexity_entropy_lo_r1 << " dom<=" << opts.complexity_dom_hi_r1
              << "] R2[n=" << det.complexity_ref_n_r2
              << " ent>=" << opts.complexity_entropy_lo_r2 << " dom<=" << opts.complexity_dom_hi_r2 << "]";
    if (opts.max_n_rate < 1.0f)
        std::cerr << " max-n-rate=" << opts.max_n_rate;
    std::cerr << "\n";

    bool compress_out = out_path.size() >= 3 &&
                        out_path.compare(out_path.size() - 3, 3, ".gz") == 0;
    FastqWriter merged_writer(out_path, compress_out, wrt_threads);

    std::unique_ptr<FastqWriter> r1_writer, r2_writer, orphan_writer;
    if (!r1_out_path.empty()) {
        bool c = r1_out_path.size() >= 3 && r1_out_path.compare(r1_out_path.size()-3,3,".gz")==0;
        r1_writer = std::make_unique<FastqWriter>(r1_out_path, c, wrt_threads);
    }
    if (!r2_out_path.empty()) {
        bool c = r2_out_path.size() >= 3 && r2_out_path.compare(r2_out_path.size()-3,3,".gz")==0;
        r2_writer = std::make_unique<FastqWriter>(r2_out_path, c, wrt_threads);
    }
    if (!orphan_out_path.empty()) {
        bool c = orphan_out_path.size() >= 3 && orphan_out_path.compare(orphan_out_path.size()-3,3,".gz")==0;
        orphan_writer = std::make_unique<FastqWriter>(orphan_out_path, c, wrt_threads);
    }
    // Mate-preserved orphan streams. When given, orphan-R1 (molecule 5', C->T frame) and
    // orphan-R2 (RC molecule 3'; ds biological, ss dA-tail) are routed separately so the
    // estimator never conflates the two channels. Fallback: combined --orphan-out.
    std::unique_ptr<FastqWriter> orphan_r1_writer, orphan_r2_writer;
    if (!orphan_r1_out_path.empty()) {
        bool c = orphan_r1_out_path.size() >= 3 && orphan_r1_out_path.compare(orphan_r1_out_path.size()-3,3,".gz")==0;
        orphan_r1_writer = std::make_unique<FastqWriter>(orphan_r1_out_path, c, wrt_threads);
    }
    if (!orphan_r2_out_path.empty()) {
        bool c = orphan_r2_out_path.size() >= 3 && orphan_r2_out_path.compare(orphan_r2_out_path.size()-3,3,".gz")==0;
        orphan_r2_writer = std::make_unique<FastqWriter>(orphan_r2_out_path, c, wrt_threads);
    }

    PairQueue    pair_q(2 * wrk_threads);
    MergeOutQueue out_q;

    int64_t n_pairs      = 0;
    int64_t n_merged     = 0;
    int64_t n_unmerged   = 0;
    int64_t n_orphan     = 0;
    int64_t n_orphan_r1  = 0;
    int64_t n_orphan_r2  = 0;
    int64_t n_dropped    = 0;

    // Writer thread
    std::thread writer_thr([&] {
        uint64_t expected = 0;
        MergeBatch batch;
        while (out_q.pop_ordered(expected, batch)) {
            n_merged    += batch.n_merged;
            n_unmerged  += batch.n_unmerged;
            n_orphan    += batch.n_orphan;
            n_orphan_r1 += batch.n_orphan_r1;
            n_orphan_r2 += batch.n_orphan_r2;
            n_dropped   += batch.n_dropped;
            for (auto& mr : batch.records) {
                ++n_pairs;
                if (mr.is_merged) {
                    if (!mr.merged.seq.empty())
                        merged_writer.write(mr.merged);
                } else if (mr.is_orphan) {
                    // Prefer the mate-specific stream; fall back to the combined orphan stream.
                    FastqWriter* ow = mr.orphan_r1 ? orphan_r1_writer.get() : orphan_r2_writer.get();
                    if (!ow) ow = orphan_writer.get();
                    if (ow && !mr.merged.seq.empty())
                        ow->write(mr.merged);
                } else {
                    if (r1_writer && !mr.merged.seq.empty())
                        r1_writer->write(mr.merged);
                    if (r2_writer && !mr.unmerged2.seq.empty())
                        r2_writer->write(mr.unmerged2);
                }
            }
            ++expected;
        }
    });

    taph::SampleDamageProfile damage_prof;
    taph::SampleDamageProfile* prof_ptr  = opts.damage_out.empty() ? nullptr : &damage_prof;
    OverlapSubstCounts subst_counts;
    OverlapSubstCounts* subst_ptr = (opts.subst_out.empty() && opts.subst_binary.empty()) ? nullptr : &subst_counts;

    // Worker threads
    std::vector<std::thread> workers;
    workers.reserve(wrk_threads);
    for (int t = 0; t < wrk_threads; ++t)
        workers.emplace_back(merge_worker,
                             std::ref(pair_q), std::ref(out_q), std::cref(opts),
                             prof_ptr, subst_ptr);

    // Reader (main thread): drain scan_buf first, then continue from open readers
    {
        ReadPair pair;
        std::vector<ReadPair> buf;
        buf.reserve(PAIR_BATCH_SZ);
        uint64_t batch_id = 0;

        for (auto& sp : scan_buf) {
            buf.push_back(std::move(sp));
            if ((int)buf.size() == PAIR_BATCH_SZ) {
                PairBatch pb; pb.id = batch_id++; pb.pairs = std::move(buf);
                pair_q.push(std::move(pb));
                buf.clear(); buf.reserve(PAIR_BATCH_SZ);
            }
        }
        scan_buf.clear(); scan_buf.shrink_to_fit();

        while (r1_rdr->read(pair.r1) && r2_rdr->read(pair.r2)) {
            buf.push_back(std::move(pair));
            if ((int)buf.size() == PAIR_BATCH_SZ) {
                PairBatch pb;
                pb.id    = batch_id++;
                pb.pairs = std::move(buf);
                pair_q.push(std::move(pb));
                buf.clear();
                buf.reserve(PAIR_BATCH_SZ);
            }
        }
        if (!buf.empty()) {
            PairBatch pb;
            pb.id    = batch_id++;
            pb.pairs = std::move(buf);
            pair_q.push(std::move(pb));
        }
    }

    pair_q.set_done();
    for (auto& w : workers) w.join();
    out_q.set_done();
    writer_thr.join();

    if (prof_ptr) {
        taph::FrameSelector::finalize_sample_profile(damage_prof);
        std::ofstream jf(opts.damage_out);
        if (!jf.good()) {
            std::cerr << "Error: cannot write damage profile: " << opts.damage_out << "\n";
        } else {
            taph::ProfileJsonInput pji;
            pji.sample_name = r1_paths[0];
            pji.n_reads     = static_cast<uint64_t>(n_merged);
            taph::profile_to_json(damage_prof, jf, pji);
            std::cerr << "Damage profile:   " << opts.damage_out << "\n";
        }
    }

    if (subst_ptr && !opts.subst_out.empty()) {
        std::ofstream sf(opts.subst_out);
        if (!sf.good()) {
            std::cerr << "Error: cannot write subst matrix: " << opts.subst_out << "\n";
        } else {
            static const char* BASES = "ACGT";
            sf << "# overlap substitution matrix from fqdup merge\n"
               << "# r1=R1 base, rc2=RC(R2) base at same insert position before consensus\n"
               << "# Deamination 5': r1=T rc2=C at small fwd_pos\n"
               << "# Deamination 3': r1=G rc2=A at small rev_pos\n"
               << "# OxoG:           r1=T rc2=G (all positions)\n"
               << "# n_pairs=" << subst_ptr->n_pairs
               << " n_bases=" << subst_ptr->n_bases << "\n";
            sf << "strand\tpos\tr1\trc2\tcount\n";
            for (int p = 0; p < OverlapSubstCounts::MAX_POS; ++p)
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b) {
                        if (subst_ptr->fwd[p][a][b])
                            sf << "fwd\t" << p << "\t" << BASES[a] << "\t" << BASES[b]
                               << "\t" << subst_ptr->fwd[p][a][b] << "\n";
                        if (subst_ptr->rev[p][a][b])
                            sf << "rev\t" << p << "\t" << BASES[a] << "\t" << BASES[b]
                               << "\t" << subst_ptr->rev[p][a][b] << "\n";
                    }
            // all-positions matrix
            sf << "# all positions combined\nall\t.\t";
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    sf << BASES[a] << BASES[b] << "=" << subst_ptr->all[a][b] << "\t";
            sf << "\n";
            std::cerr << "Subst matrix:     " << opts.subst_out << "\n";
        }
    }

    // Binary .bsubst writer
    // v2 layout: magic(8) n_pairs(8) n_bases(8) n_pos(4) fwd[30][4][4] rev[30][4][4] all[4][4]
    //            n_len_bins(4) bin_lo[N](4*N) bin_hi[N](4*N) bin_n_pairs[N](8*N) bin_n_bases[N](8*N)
    //            fwd_len[N][30][4][4] rev_len[N][30][4][4]
    // All integers little-endian. N = LenBins::N = 6.
    if (subst_ptr && !opts.subst_binary.empty()) {
        FILE* bf = fopen(opts.subst_binary.c_str(), "wb");
        if (!bf) {
            std::cerr << "Error: cannot write binary subst: " << opts.subst_binary << "\n";
        } else {
            static const uint8_t MAGIC[8] = {'B','S','U','B','S','T',0x03,0x00};
            fwrite(MAGIC, 1, 8, bf);
            fwrite(&subst_ptr->n_pairs, 8, 1, bf);
            fwrite(&subst_ptr->n_bases, 8, 1, bf);
            int32_t npos = OverlapSubstCounts::MAX_POS;
            fwrite(&npos, 4, 1, bf);
            fwrite(subst_ptr->fwd, sizeof(subst_ptr->fwd), 1, bf);
            fwrite(subst_ptr->rev, sizeof(subst_ptr->rev), 1, bf);
            fwrite(subst_ptr->all, sizeof(subst_ptr->all), 1, bf);
            // v3 library trailer: sits between `all` and the v2 len-bin block. Carries the
            // library-type contract + merge counts so `fqdup damage` can splice a "library"
            // block into the profile JSON (no sidecar, no separate pass). Read by V3 readers;
            // V1/V2 readers stop after `all` and never see it.
            {
                uint8_t  ss8  = is_ss ? 1 : 0;
                uint8_t  src8 = (lib_type_source == "declared") ? 1 : 0;
                uint8_t  bal8 = ((n_orphan_r1 + n_orphan_r2 == n_orphan) &&
                                 (n_merged + n_unmerged + n_orphan_r1 + n_orphan_r2 + n_dropped
                                  == n_pairs + n_dropped)) ? 1 : 0;
                float    conf = lib_type_conf;
                int64_t  pairs_in = n_merged + n_unmerged + n_orphan_r1 + n_orphan_r2 + n_dropped;
                fwrite(&ss8, 1, 1, bf); fwrite(&src8, 1, 1, bf); fwrite(&conf, 4, 1, bf);
                auto wstr = [&](const std::string& s) {
                    uint16_t n = (uint16_t)std::min(s.size(), (size_t)65535);
                    fwrite(&n, 2, 1, bf); if (n) fwrite(s.data(), 1, n, bf);
                };
                wstr(opts.adapter_5p_linker); wstr(opts.adapter1); wstr(opts.adapter2);
                int64_t cnts[6] = { pairs_in, n_merged, n_unmerged,
                                    n_orphan_r1, n_orphan_r2, n_dropped };
                fwrite(cnts, 8, 6, bf);
                fwrite(&bal8, 1, 1, bf);
            }
            // v2 extension: per-length-bin data
            int32_t n_len_bins = LenBins::N;
            fwrite(&n_len_bins, 4, 1, bf);
            for (int i = 0; i < LenBins::N; ++i) fwrite(&LenBins::EDGES[i],   4, 1, bf);
            for (int i = 0; i < LenBins::N; ++i) fwrite(&LenBins::EDGES[i+1], 4, 1, bf);
            fwrite(subst_ptr->bin_n_pairs, sizeof(subst_ptr->bin_n_pairs), 1, bf);
            fwrite(subst_ptr->bin_n_bases, sizeof(subst_ptr->bin_n_bases), 1, bf);
            fwrite(subst_ptr->fwd_len, sizeof(subst_ptr->fwd_len), 1, bf);
            fwrite(subst_ptr->rev_len, sizeof(subst_ptr->rev_len), 1, bf);
            fclose(bf);
            std::cerr << "Subst binary:     " << opts.subst_binary << "\n";
        }
    }

    double pct = n_pairs > 0 ? 100.0 * n_merged / n_pairs : 0.0;
    int64_t pairs_in = n_merged + n_unmerged + n_orphan_r1 + n_orphan_r2 + n_dropped;
    bool balanced = (n_orphan_r1 + n_orphan_r2 == n_orphan) && (pairs_in == n_pairs + n_dropped);
    std::cerr << "Pairs processed:  " << n_pairs    << "\n"
              << "Merged:           " << n_merged   << " (" << pct << "%)\n"
              << "Unmerged:         " << n_unmerged << "\n"
              << "Orphan:           " << n_orphan   << " (R1=" << n_orphan_r1
              << " R2=" << n_orphan_r2 << ")\n"
              << "Dropped:          " << n_dropped  << "\n"
              << "Balance:          " << (balanced ? "OK" : "MISMATCH")
              << " (pairs_in=" << pairs_in << ")\n";

    // Verification: the Myers adapter-fragment drop (highest FP risk). It is CLEAN if the
    // count is small and the length histogram peaks at adapter length (dimer-shaped ~15-24bp);
    // it is BIASING damage if it shaves the short real-insert tail. Length histogram of drops:
    {
        int64_t fd = g_frag_drop.load();
        std::cerr << "Adapter-fragment drops (Myers): " << fd;
        if (fd > 0) {
            std::cerr << "  len-hist[";
            for (int L = 12; L < 64; ++L) {
                int64_t c = g_frag_len_hist[L].load();
                if (c > 0) std::cerr << L << ":" << c << " ";
            }
            std::cerr << "]";
        }
        std::cerr << "\n";
        std::cerr << "Technical-drop reason: fragment=" << g_reason_frag.load() << "\n";
    }

    if (g_clip_dbg) {
        static const char* STEP[8] = {
            "wholeread@detected-adapter", "wholeread@universal", "linker-5p",
            "adapter5p-5p", "adapter3p-5p", "universal-5p", "polyG-3p", "already-short"};
        std::cerr << "Clip-death attribution (reads taken below min_length):\n";
        for (int m = 0; m < 2; ++m)
            for (int s = 0; s < 8; ++s) {
                int64_t v = g_clip_death[m][s].load();
                if (v) std::cerr << "  " << (m ? "R2" : "R1") << "  " << STEP[s] << " = " << v << "\n";
            }
    }

    // ---- comprehensive lossless merge-QC JSON ----
    if (!opts.json_out.empty()) {
        std::ofstream j(opts.json_out);
        if (!j) {
            std::cerr << "Error: cannot write --json " << opts.json_out << "\n";
        } else {
            j.precision(9);
            auto js = [](const std::string& s) {
                std::string o = "\"";
                for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
                o += '"'; return o;
            };
            const char* udg = det.is_udg ? "full" : (det.is_half_udg ? "half" : "none");
            j << "{\n";
            j << "  \"tool\": \"fqdup merge\",\n";
            j << "  \"library\": {\n";
            j << "    \"type\": \"" << (is_ss ? "ss" : "ds") << "\",\n";
            j << "    \"type_source\": " << js(lib_type_source) << ",\n";
            j << "    \"type_confidence\": " << lib_type_conf << ",\n";
            j << "    \"udg\": \"" << udg << "\",\n";
            j << "    \"skip_terminal\": " << det.skip_terminal << ",\n";
            j << "    \"damage_5prime\": " << det.damage_5p << ",\n";
            j << "    \"damage_3prime\": " << det.damage_3p << "\n";
            j << "  },\n";
            j << "  \"adapters\": {\n";
            j << "    \"adapter1\": " << js(opts.adapter1) << ",\n";
            j << "    \"adapter2\": " << js(opts.adapter2) << ",\n";
            j << "    \"linker_5prime\": " << js(opts.adapter_5p_linker) << ",\n";
            j << "    \"tech_constructs\": [";
            for (size_t i = 0; i < det.tech_info.size(); ++i) {
                const auto& t = det.tech_info[i];
                j << (i ? "," : "") << "\n      {\"seq\": " << js(t.seq)
                  << ", \"length\": " << t.seq.size()
                  << ", \"support\": " << t.support
                  << ", \"ic_body_bits\": " << t.ic_body << "}";
            }
            j << (det.tech_info.empty() ? "" : "\n    ") << "]\n";
            j << "  },\n";
            j << "  \"complexity_gate\": {\n";
            j << "    \"r1\": {\"entropy_lo\": " << det.complexity_entropy_lo_r1
              << ", \"dom_hi\": " << det.complexity_dom_hi_r1
              << ", \"ref_n\": " << det.complexity_ref_n_r1 << "},\n";
            j << "    \"r2\": {\"entropy_lo\": " << det.complexity_entropy_lo_r2
              << ", \"dom_hi\": " << det.complexity_dom_hi_r2
              << ", \"ref_n\": " << det.complexity_ref_n_r2 << "}\n";
            j << "  },\n";
            j << "  \"params\": {\"min_overlap\": " << opts.min_ov
              << ", \"max_mismatch_rate\": " << opts.max_mm_rate
              << ", \"min_length\": " << opts.min_length
              << ", \"max_length\": " << opts.max_length
              << ", \"poly_g_min_run\": " << opts.poly_g_min_run << "},\n";
            j << "  \"counts\": {\n";
            j << "    \"pairs_in\": " << pairs_in << ",\n";
            j << "    \"merged\": " << n_merged << ",\n";
            j << "    \"unmerged\": " << n_unmerged << ",\n";
            j << "    \"orphan_r1\": " << n_orphan_r1 << ",\n";
            j << "    \"orphan_r2\": " << n_orphan_r2 << ",\n";
            j << "    \"dropped\": " << n_dropped << ",\n";
            j << "    \"balanced\": " << (balanced ? "true" : "false") << "\n";
            j << "  },\n";
            j << "  \"drops\": {\n";
            j << "    \"technical_total\": " << g_frag_drop.load() << ",\n";
            j << "    \"by_reason\": {\"adapter_fragment\": " << g_reason_frag.load() << "},\n";
            j << "    \"length_histogram\": {";
            bool first = true;
            for (int L = 12; L < 64; ++L) {
                int64_t c = g_frag_len_hist[L].load();
                if (c > 0) { j << (first ? "" : ", ") << "\"" << L << "\": " << c; first = false; }
            }
            j << "}\n";
            j << "  }\n";
            j << "}\n";
            std::cerr << "Merge JSON:       " << opts.json_out << "\n";
        }
    }

    return 0;
}
