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
#include <cmath>
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
};

struct MergeBatch {
    uint64_t id = 0;
    std::vector<MergeRecord> records;
    int64_t n_merged   = 0;
    int64_t n_unmerged = 0;
    int64_t n_orphan   = 0;
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
    int   min_length    = 15;
    int   skip_terminal = 0;
    int   clip_5p       = 0;   // 0=disabled; hard-clip N bases from R1 5' end before merge
    int   poly_g_min_run = 0;   // 0=disabled; trim 3' poly-G runs >= this length
    float min_entropy   = 0.0f; // 0=disabled; Shannon entropy floor (bits, max=2.0)
    float max_n_rate    = 1.0f; // 1=disabled; max fraction of N bases
    std::string adapter1;
    std::string adapter2;
    std::vector<std::string> extra_adapters1;   // additional R1 adapters to try (from --adapter-fasta)
    std::string damage_out;     // path for paired damage profile JSON; empty=disabled
    std::string subst_out;      // path for overlap substitution matrix TSV; empty=disabled
    std::string subst_binary;   // path for binary .bsubst format; empty=disabled
};

// ============================================================================
// Pre-scan: auto-detect adapter sequences, library geometry, UDG status.
// Reads the first n_scan pairs; reader is positioned just after scan_buf on return.
// ============================================================================

struct DetectedParams {
    std::string adapter1;      // 24bp P7_RC suffix found on R1 past overlap
    std::string adapter2;      // 24bp P5_RC suffix found on RC(R2) past overlap
    bool   is_ss           = false;
    bool   is_udg          = false;
    bool   is_half_udg     = false;
    bool   poly_g_detected = false;  // NextSeq/NovaSeq dark cycle signal
    int    skip_terminal   = 0;
    float  damage_5p       = 0.f;
    float  damage_3p       = 0.f;
};

// Forward declaration (defined below)
static int find_adapter_in(const std::string& seq, const std::string& adapter,
                           int min_pos, int max_mm, int alen_req);

static DetectedParams detect_merge_params(
        std::vector<ReadPair>& scan_buf,
        int min_ov, float mm_loose = 0.15f) {

    // Per-position (0..29) mismatch tallies across all pairs with a found overlap
    std::vector<int64_t> pos_mm(30, 0), pos_total(30, 0);
    // 5' prefix agreement (ds geometry check): R1[0:8] vs rcR2[0:8]
    int64_t prefix_agree = 0, prefix_total = 0;
    // Adapter suffix k-mer frequency
    std::unordered_map<std::string, int> a1_freq, a2_freq;

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
        // adapter2 in RC(R2) at position r2s - 12 (just before the insert)
        if (r2s >= 12)
            a2_freq[rc2_seq.substr(r2s - 12, 12)]++;
    }

    DetectedParams d;
    if (prefix_total == 0) return d;

    // ---- library geometry ----
    float prefix_agree_rate = (float)prefix_agree / (float)prefix_total;
    d.is_ss = (prefix_agree_rate < 0.3f);

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

    float asym = (d.damage_3p > 0.001f) ? d.damage_5p / d.damage_3p : 10.f;
    if (d.damage_5p >= 0.02f && asym > 2.5f) {
        d.is_half_udg   = true;
        d.skip_terminal = 4;
    } else if (d.damage_5p >= 0.02f) {
        d.skip_terminal = 4;       // untreated: both ends damaged
    } else if (d.damage_5p < 0.005f && d.damage_3p < 0.005f) {
        d.is_udg = true;           // confidently UDG: very low at both ends
        d.skip_terminal = 0;
    } else if (asym > 2.5f && d.damage_5p >= 0.01f) {
        d.is_half_udg   = true;    // asymmetric even at low absolute rate → half-UDG
        d.skip_terminal = 4;
    }
    // else: low absolute damage, symmetric → ambiguous (low f_anc); skip_terminal=0

    // ---- poly-G tail detection (NextSeq/NovaSeq dark cycle: last 10 positions)----
    int64_t tail_g = 0, tail_n = 0;
    for (auto& pr : scan_buf) {
        int L1 = (int)pr.r1.seq.size();
        for (int i = std::max(0, L1 - 10); i < L1; ++i)
            { ++tail_n; if (pr.r1.seq[i] == 'G') ++tail_g; }
    }
    d.poly_g_detected = tail_n > 1000 && (float)tail_g / tail_n > 0.40f;

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
    // RC(R2) adapter: use detected if high confidence, else TruSeq default
    std::string a2 = best_kmer(a2_freq);
    if (a2.empty() || a2_freq[a2] < 100) a2 = TRUSEQ_RC2;
    d.adapter2 = a2;
    // adapter1 = RC(adapter2)
    d.adapter1 = revcomp(d.adapter2);  // revcomp from fastq_common.hpp

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

static bool passes_qc(const FastqRecord& rec, const MergeOpts& opts) {
    if ((int)rec.seq.size() < opts.min_length) return false;
    if (opts.max_n_rate < 1.0f) {
        int ns = (int)std::count(rec.seq.begin(), rec.seq.end(), 'N');
        if ((float)ns / (float)rec.seq.size() > opts.max_n_rate) return false;
    }
    if (opts.min_entropy > 0.f && shannon_entropy(rec.seq) < opts.min_entropy) return false;
    return true;
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

    PairBatch batch;
    while (in_q.pop(batch)) {
        MergeBatch out;
        out.id = batch.id;
        out.records.reserve(batch.pairs.size());

        for (auto& pr : batch.pairs) {
            FastqRecord& r1 = pr.r1;
            FastqRecord& r2 = pr.r2;

            if (opts.poly_g_min_run > 0) {
                trim_polybase(r1, 'G', opts.poly_g_min_run, opts.min_length);
                trim_polybase(r2, 'G', opts.poly_g_min_run, opts.min_length);
            }

            if (opts.clip_5p > 0 && (int)r1.seq.size() > opts.clip_5p) {
                r1.seq  = r1.seq.substr(opts.clip_5p);
                r1.qual = r1.qual.substr(opts.clip_5p);
            }

            int L1 = (int)r1.seq.size();
            int L2 = (int)r2.seq.size();

            if (L1 < opts.min_ov || L2 < opts.min_ov ||
                L1 > MAX_READ_LEN  || L2 > MAX_READ_LEN) {
                // Too short or too long: emit unmerged
                MergeRecord mr;
                mr.merged    = r1;
                mr.unmerged2 = r2;
                mr.is_merged = false;
                out.records.push_back(std::move(mr));
                ++out.n_unmerged;
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
                    if (passes_qc(mr.merged, opts)) {
                        ++out.n_merged;
                    } else {
                        mr.is_merged = false;
                        mr.merged    = r1;
                        mr.unmerged2 = r2;
                        ++out.n_unmerged;
                    }
                    out.records.push_back(std::move(mr));
                    continue;
                }
            }

            if (best_ov < opts.min_ov) {
                // No overlap found: trim adapters, apply QC, handle orphans
                FastqRecord t1 = r1, t2 = r2;
                if (!opts.adapter1.empty()) trim_adapter(t1, opts.adapter1, opts.min_length);
                if (!opts.adapter2.empty()) trim_adapter(t2, opts.adapter2, opts.min_length);
                bool ok1 = passes_qc(t1, opts);
                bool ok2 = passes_qc(t2, opts);
                MergeRecord mr;
                mr.is_merged = false;
                if (ok1 && ok2) {
                    mr.merged    = std::move(t1);
                    mr.unmerged2 = std::move(t2);
                    ++out.n_unmerged;
                } else if (ok1 || ok2) {
                    mr.is_orphan = true;
                    mr.merged    = ok1 ? std::move(t1) : std::move(t2);
                    ++out.n_orphan;
                }
                // both fail → drop (emit nothing; push empty record skipped by writer)
                if (ok1 || ok2) out.records.push_back(std::move(mr));
                continue;
            }

            // Merge via Bayesian consensus
            if (do_profile)
                taph::FrameSelector::update_sample_profile_paired(local_prof, r1.seq, r2.seq);
            if (do_subst)
                accum_overlap_subst(local_subst, r1.seq, rc2_seq, r2_offset, best_ov, opts.skip_terminal);
            MergeRecord mr;
            mr.is_merged = true;
            consensus_merge(r1, rc2_seq, rc2_qual, best_ov, r2_offset, mr.merged);
            // Trim adapter artifacts from merged output:
            // 5': adapter complement (CTCTTCCGATCT) from library-prep bleedthrough
            // 3': adapter sequence surviving in Phase-2 suffix or mismerges
            if (!opts.adapter2.empty()) trim_adapter_5p(mr.merged, opts.adapter2, opts.min_length);
            if (!opts.adapter1.empty()) trim_adapter(mr.merged, opts.adapter1, opts.min_length);

            if (!passes_qc(mr.merged, opts)) {
                mr.is_merged = false;
                mr.merged    = r1;
                mr.unmerged2 = r2;
                ++out.n_unmerged;
            } else {
                ++out.n_merged;
            }
            out.records.push_back(std::move(mr));
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
        "  --min-length N     Discard merged reads shorter than N bp (default: 15)\n"
        "  --clip-r1-5p N     Hard-clip N bases from R1 5' end before overlap (removes adapter stubs)\n"
        "  --min-entropy F    Discard low-complexity merged reads; Shannon entropy floor in bits\n"
        "                     (0=disabled; poly-G≈0, random≈2.0; default: 0)\n"
        "  --max-n-rate F     Discard merged reads with N fraction above F (default: 1.0=off)\n\n"
        "Adapter trimming (for unmerged pairs):\n"
        "  --adapter1 SEQ     R1 adapter sequence (Illumina P7 RC)\n"
        "  --adapter2 SEQ     R2 adapter sequence (Illumina P7)\n"
        "  --adapter-fasta F  FASTA with adapter pairs (odd=R1, even=R2); multiple pairs supported\n\n"
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
        else if (a == "--damage-out"       && i+1 < argc) opts.damage_out  = argv[++i];
        else if (a == "--subst-out"        && i+1 < argc) opts.subst_out   = argv[++i];
        else if (a == "--subst-binary"     && i+1 < argc) opts.subst_binary = argv[++i];
        else if ((a == "-p" || a == "--threads") && i+1 < argc) n_threads = std::stoi(argv[++i]);
        else if (a == "--min-overlap"      && i+1 < argc) opts.min_ov       = std::stoi(argv[++i]);
        else if (a == "--max-mm-rate"      && i+1 < argc) opts.max_mm_rate  = std::stof(argv[++i]);
        else if (a == "--min-length"       && i+1 < argc) opts.min_length   = std::stoi(argv[++i]);
        else if (a == "--adapter1"         && i+1 < argc) opts.adapter1       = argv[++i];
        else if (a == "--adapter2"         && i+1 < argc) opts.adapter2       = argv[++i];
        else if (a == "--adapter-fasta"    && i+1 < argc) adapter_fasta       = argv[++i];
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
        opts.skip_terminal = det.skip_terminal;
        if (opts.poly_g_min_run == 0 && det.poly_g_detected)
            opts.poly_g_min_run = 10;
    } else if (opts.adapter1.empty() && !opts.adapter2.empty()) {
        opts.adapter1 = revcomp(opts.adapter2);
    } else if (opts.adapter2.empty() && !opts.adapter1.empty()) {
        opts.adapter2 = revcomp(opts.adapter1);
    }

    std::string lib_type = det.is_ss ? "ss" : "ds";
    std::string udg_type = det.is_udg ? "UDG" : (det.is_half_udg ? "half-UDG" : "untreated");
    std::cerr << "fqdup merge: " << n_threads << " threads (workers=" << wrk_threads << ")\n"
              << "  library=" << lib_type << " damage=" << udg_type
              << " skip-terminal=" << opts.skip_terminal << "\n"
              << "  damage_5p=" << det.damage_5p << " damage_3p=" << det.damage_3p << "\n"
              << "  adapter1=" << (opts.adapter1.empty() ? "(none)" : opts.adapter1) << "\n"
              << "  adapter2=" << (opts.adapter2.empty() ? "(none)" : opts.adapter2) << "\n"
              << "  min-overlap=" << opts.min_ov
              << " max-mm-rate=" << opts.max_mm_rate
              << " min-length=" << opts.min_length;
    if (opts.clip_5p > 0)
        std::cerr << " clip-r1-5p=" << opts.clip_5p;
    if (opts.poly_g_min_run > 0)
        std::cerr << " poly-g=" << opts.poly_g_min_run;
    if (opts.min_entropy > 0.f)
        std::cerr << " min-entropy=" << opts.min_entropy;
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

    PairQueue    pair_q(2 * wrk_threads);
    MergeOutQueue out_q;

    int64_t n_pairs    = 0;
    int64_t n_merged   = 0;
    int64_t n_unmerged = 0;
    int64_t n_orphan   = 0;

    // Writer thread
    std::thread writer_thr([&] {
        uint64_t expected = 0;
        MergeBatch batch;
        while (out_q.pop_ordered(expected, batch)) {
            n_merged   += batch.n_merged;
            n_unmerged += batch.n_unmerged;
            n_orphan   += batch.n_orphan;
            for (auto& mr : batch.records) {
                ++n_pairs;
                if (mr.is_merged) {
                    if (!mr.merged.seq.empty())
                        merged_writer.write(mr.merged);
                } else if (mr.is_orphan) {
                    if (orphan_writer && !mr.merged.seq.empty())
                        orphan_writer->write(mr.merged);
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
            static const uint8_t MAGIC[8] = {'B','S','U','B','S','T',0x02,0x00};
            fwrite(MAGIC, 1, 8, bf);
            fwrite(&subst_ptr->n_pairs, 8, 1, bf);
            fwrite(&subst_ptr->n_bases, 8, 1, bf);
            int32_t npos = OverlapSubstCounts::MAX_POS;
            fwrite(&npos, 4, 1, bf);
            fwrite(subst_ptr->fwd, sizeof(subst_ptr->fwd), 1, bf);
            fwrite(subst_ptr->rev, sizeof(subst_ptr->rev), 1, bf);
            fwrite(subst_ptr->all, sizeof(subst_ptr->all), 1, bf);
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
    std::cerr << "Pairs processed:  " << n_pairs    << "\n"
              << "Merged:           " << n_merged   << " (" << pct << "%)\n"
              << "Unmerged:         " << n_unmerged << "\n"
              << "Orphan:           " << n_orphan   << "\n";

    return 0;
}
