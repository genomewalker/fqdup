// fqdup extend — self-referential k-mer graph extension for aDNA deduplication
//
// Replaces Tadpole in the pipeline.  Takes a FASTQ of merged reads (fastp output)
// and extends each read at both termini by walking unique paths in a de Bruijn graph
// built from the same input.  PCR copies of the same molecule extend identically →
// hash identically → deduplicated by derep_pairs.  Different molecules extend
// differently → kept separate.
//
// Ancient DNA damage masking: damaged terminal positions are excluded from k-mer
// construction so that deaminated copies build the same graph as undamaged copies.
// The walk can only infer across damaged terminals if other reads provide clean support
// for that transition; it never uses masked-terminal observations from the read being
// extended.
//
// Algorithm:
//   Pass 0 (optional): DART damage estimation (default: sample 500k reads)
//   Pass 1: Build canonical k-mer graph (only min(kmer, rc(kmer)) stored per position)
//   Pass 2: Extend each read at both termini and write extended FASTQ
//
// Memory design:
//   Canonical storage halves table entries vs oriented storage.
//   4096-shard ska::flat_hash_map for pass1 incremental insertion with per-shard locks.
//   After pass1, finalize() converts each shard to a sorted FlatEntry[] and frees the ska
//   maps — pass2 uses binary search at 16 bytes/entry with no hash-table overhead.
//   Fibonacci-mix shard routing gives uniform distribution.
//
// Short-read safe: reads with no clean interior k-mers are passed through unchanged.
//
// No stdin support: three-pass algorithm requires re-reading the file.

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <malloc.h>
#include <sys/stat.h>

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"

#ifdef HAVE_BBHASH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <BooPHF.h>
#pragma GCC diagnostic pop
#endif

namespace {

#ifdef HAVE_BBHASH
using ShardMPHF = boomphf::mphf<uint64_t, boomphf::SingleHashFunctor<uint64_t>>;
#endif

// ============================================================================
// 2-bit k-mer encoding (k ≤ 31, MSB-first: first base in highest bits)
// A=0, C=1, G=2, T=3
// ============================================================================

static auto make_enc4() {
    std::array<uint8_t, 256> t;
    t.fill(0xFF);
    t['A'] = t['a'] = 0;
    t['C'] = t['c'] = 1;
    t['G'] = t['g'] = 2;
    t['T'] = t['t'] = 3;
    return t;
}
static const auto kEnc4 = make_enc4();

static constexpr char kDec4[4] = {'A', 'C', 'G', 'T'};

inline uint8_t enc4(char c) { return kEnc4[static_cast<unsigned char>(c)]; }

// Complement in 2-bit space: A(0)↔T(3), C(1)↔G(2)
inline uint8_t comp4(uint8_t b) { return b ^ 3u; }

// Encode k consecutive bases starting at seq[pos].
// Returns false if any base is not in {A,C,G,T,a,c,g,t}.
static bool encode_kmer(const char* seq, int pos, int k, uint64_t& out) {
    uint64_t v = 0;
    for (int i = 0; i < k; ++i) {
        uint8_t b = enc4(seq[pos + i]);
        if (b == 0xFF) return false;
        v = (v << 2) | b;
    }
    out = v;
    return true;
}

// Shift k-mer left by 2 bits and append one base on the right.
static uint64_t advance_right(uint64_t kmer, uint8_t b, uint64_t mask) {
    return ((kmer << 2) | b) & mask;
}

// Reverse-complement of a k-mer.
static uint64_t kmer_rc(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; ++i) {
        rc = (rc << 2) | comp4(static_cast<uint8_t>(kmer & 3u));
        kmer >>= 2;
    }
    return rc;
}

// Count distinct 2-mers (dinucleotides) in a k-mer (MSB-first encoded).
// k-mers with < 3 distinct 2-mers are low-complexity and skipped.
static int count_distinct_dimers(uint64_t kmer, int k) {
    uint32_t seen = 0;
    for (int i = 0; i < k - 1; ++i) {
        int shift = 2 * (k - 2 - i);
        int dimer = static_cast<int>((kmer >> shift) & 0xFu);
        seen |= (1u << dimer);
    }
    return __builtin_popcount(seen);
}

// ============================================================================
// BoundedQueue — thread-safe MPSC/MPMC queue with capacity backpressure.
// push() blocks when full; pop() blocks when empty; close() unblocks all waiters.
// ============================================================================

template <class T>
class BoundedQueue {
public:
    explicit BoundedQueue(size_t cap) : cap_(cap), closed_(false) {}

    bool push(T item) {
        std::unique_lock<std::mutex> lk(mtx_);
        not_full_.wait(lk, [&] { return q_.size() < cap_ || closed_; });
        if (closed_) return false;
        q_.push_back(std::move(item));
        not_empty_.notify_one();
        return true;
    }

    bool pop(T& out) {
        std::unique_lock<std::mutex> lk(mtx_);
        not_empty_.wait(lk, [&] { return !q_.empty() || closed_; });
        if (q_.empty()) return false;
        out = std::move(q_.front());
        q_.pop_front();
        not_full_.notify_one();
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> lk(mtx_);
        closed_ = true;
        not_empty_.notify_all();
        not_full_.notify_all();
    }

private:
    size_t                  cap_;
    bool                    closed_;
    std::mutex              mtx_;
    std::condition_variable not_empty_, not_full_;
    std::deque<T>           q_;
};

// ============================================================================
// Canonical k-mer graph
//
// Only the canonical form (min(kmer, rc(kmer))) of each k-mer is stored.
// This halves table entries compared to oriented storage.
//
// Each entry stores edge counts for BOTH traversal orientations:
//   fwd[b]: times base b appeared after this k-mer when traversed in canonical direction.
//   rev[b]: times base b appeared after rc(canonical k-mer) when traversed in that direction.
//           Equivalently: times comp(b) appeared *before* the canonical k-mer.
//
// Edge counts saturate at uint8_t max (255).  At typical aDNA coverage (3-30×)
// counts stay well below 255; the rare saturation case does not affect walk
// correctness since we only branch on count ≥ min_count (default 2).
//
// Walk correctness:
//   At each walk step, let cur be the k-mer in the walking direction, canon = min(cur, rc(cur)).
//   is_fwd = (cur == canon).
//   Forward edges  = is_fwd ? entry.fwd : entry.rev
//   Predecessor edges (reciprocal check for next k-mer next_canon, next_is_fwd):
//                  = next_is_fwd ? next_entry.rev : next_entry.fwd
//   This gives the same number of table lookups (2 per step) as oriented storage,
//   but with a 2× smaller table for better cache utilization.
// ============================================================================

struct CanonKmerEntry {
    uint8_t fwd[4] = {};  // right-edge counts when traversing in canonical direction
    uint8_t rev[4] = {};  // right-edge counts when traversing in RC-of-canonical direction
};

// ============================================================================
// ShardedKmerTable
//
// 4096 shards.  Shard routing uses top 12 bits of a Fibonacci-mixed hash
// for uniform distribution (vs low-bit routing which has k-mer bias for
// short k where low bits = last few bases).
//
// Pass 1 design — KMC3-in-RAM (no hash map, no rehash peaks):
//   Each k-mer observation is packed into a single uint64_t:
//     bits 63..7 = canon (valid for k ≤ 28: 56 bits + 7 flag bits = 63 total)
//     bits 6..5  = lb (left base, 2 bits)
//     bits 4..3  = rb (right base, 2 bits)
//     bit  2     = is_fwd (canon == kmer_fwd)
//     bit  1     = lb_valid
//     bit  0     = rb_valid
//   Raw observations are appended to per-shard vectors under one lock per shard
//   per batch (same batch-sort approach as before, but appending not inserting).
//   finalize() sorts each shard's raw[] by the full uint64_t (= sort by canon
//   first, then flags) then reduces consecutive same-canon items into FlatEntry.
//   This eliminates the synchronized rehash that caused 179 GB peaks.
//
// Memory (DS4: 1.644B unique canonical k-mers from 198M reads, k=17):
//   Pass1 accumulation: 6.7B × 8 bytes = ~54 GB (demand-paged from reservation)
//   finalize() peak:    54 GB (raw) + 26 GB (building flat[]) ≈ 80 GB
//   After finalize:     26 GB (1.644B × 16 bytes FlatEntry) + malloc_trim release
//
// Pass 2 find() is lock-free (immutable sorted flat arrays, binary search).
// ============================================================================

struct ShardedKmerTable {
    static constexpr int N = 4096;

    // Fibonacci-mix hash: maps k-mer keys uniformly across shards.
    static uint32_t shard_of(uint64_t k) {
        k ^= k >> 30;
        k *= 0xbf58476d1ce4e5b9ULL;
        k ^= k >> 27;
        k *= 0x94d049bb133111ebULL;
        k ^= k >> 31;
        return static_cast<uint32_t>(k >> 52);  // top 12 bits → [0, 4095]
    }

    // Packed raw observation: (canon << 7) | flags.  Valid for k ≤ 28.
    // Sorting as uint64_t sorts by canon (high bits) first — correct for reduce.
    static constexpr int N_FLAG_BITS = 7;

    // Sorted flat entry for pass2 binary search (16 bytes, cache-line aligned run).
    struct FlatEntry {
        uint64_t       canon;
        CanonKmerEntry edges;
    };

    struct alignas(64) Shard {
        std::vector<uint64_t>  raw;     // pass1: packed observations, freed after finalize
        mutable std::mutex     mtx;
        std::vector<FlatEntry> flat;    // pass2 (binary search): sorted by canon
        std::vector<uint32_t>  prefix;  // pass2 (binary search): prefix index into flat[]
#ifdef HAVE_BBHASH
        std::unique_ptr<ShardMPHF>  mphf;        // pass2 (MPHF): built from flat[], then flat freed
        std::vector<CanonKmerEntry> mphf_edges;  // pass2 (MPHF): edges indexed by MPHF rank
        std::vector<uint32_t>       mphf_fp;     // pass2 (MPHF): upper-32-bit canon fingerprint
#endif
    };

    std::unique_ptr<Shard[]> shards;

    // Prefix-index parameters (set by init_lookup_layout before pass1).
    int      key_bits_    = 0;
    int      prefix_bits_ = 0;
    int      prefix_shift_= 0;
    uint32_t prefix_size_ = 0;

    ShardedKmerTable() : shards(std::make_unique<Shard[]>(N)) {}

    // Call once before pass1 with the chosen k.
    // Builds a 2^prefix_bits directory per shard, reducing binary search depth
    // from ~log2(flat.size()) to ~log2(flat.size() / 2^prefix_bits).
    void init_lookup_layout(int k) {
        key_bits_    = 2 * k;
        prefix_bits_ = std::min(12, key_bits_);
        prefix_shift_= key_bits_ - prefix_bits_;
        prefix_size_ = 1u << prefix_bits_;
    }

    // Pre-reserve raw[] buffers from file-size estimate.
    // Linux demand-paging: over-reservation is cheap — uncommitted pages have no RSS cost.
    void reserve_raw(const std::string& path, int k) {
        struct stat st{};
        if (stat(path.c_str(), &st) != 0) return;
        // gzip ratio ~4x; avg uncompressed FASTQ record ~150 bytes;
        // avg clean k-mers per read ~30; 2× safety margin (demand-paged = free if unused).
        const int64_t per_shard = static_cast<int64_t>(st.st_size) * 4LL * 30LL / 150LL
                                  * 2LL / N;
        for (int i = 0; i < N; ++i)
            shards[i].raw.reserve(static_cast<size_t>(std::max((int64_t)1024, per_shard)));
    }

    // Append a batch of packed observations (sorted by shard) to the raw[] buffers.
    // shard_start[sh] = first index in `items` for shard sh; shard_start[sh+1] = end.
    // One mutex acquisition per shard per batch.
    void batch_append(const uint64_t* items, const int* shard_start) {
        for (int sh = 0; sh < N; ++sh) {
            const int lo = shard_start[sh];
            const int hi = shard_start[sh + 1];
            if (lo == hi) continue;
            Shard& s = shards[sh];
            std::lock_guard<std::mutex> lk(s.mtx);
            s.raw.insert(s.raw.end(), items + lo, items + hi);
        }
    }

    // LSB radix sort for uint64_t using 8-bit digits.
    // n_passes: number of byte-passes — ceil((2k + N_FLAG_BITS) / 8).
    // After the call, sorted data is in `in`; `out` is a same-sized scratch buffer.
    // 3-5x faster than std::sort on large arrays where memory bandwidth dominates.
    static void radix_sort_u64(std::vector<uint64_t>& in, std::vector<uint64_t>& out,
                               int n_passes) {
        const size_t n = in.size();
        out.resize(n);
        uint32_t cnt[256];
        for (int pass = 0; pass < n_passes; ++pass) {
            const int shift = pass * 8;
            std::memset(cnt, 0, sizeof(cnt));
            for (size_t i = 0; i < n; ++i) ++cnt[(in[i] >> shift) & 0xFFu];
            uint32_t sum = 0;
            for (int b = 0; b < 256; ++b) { uint32_t c = cnt[b]; cnt[b] = sum; sum += c; }
            for (size_t i = 0; i < n; ++i) out[cnt[(in[i] >> shift) & 0xFFu]++] = in[i];
            // Swap contents: sorted result moves into `in` for the next pass.
            // std::swap on vectors is O(1) — swaps internal pointer/size/capacity only.
            std::swap(in, out);
        }
        // Invariant: after every pass the result lands in `in` (caller's first arg).
    }

    // Sort each shard's raw[] by packed uint64_t (= sort by canon), reduce
    // consecutive same-canon items into FlatEntry[], then free raw[].
    // Parallel over n_threads disjoint shard ranges.  Call once after all
    // batch_append() calls complete.
    void finalize(int n_threads) {
        const int nt      = std::max(1, n_threads);
        const int per     = (N + nt - 1) / nt;
        // Radix passes needed: ceil((2k + 7 flag bits) / 8).  Skips zero high bytes for
        // small k (e.g. k=17 → 6 passes instead of 8), saving ~25% sort work.
        const int n_passes = (key_bits_ + N_FLAG_BITS + 7) / 8;

        std::vector<std::thread> workers;
        workers.reserve(nt);
        for (int t = 0; t < nt; ++t) {
            const int lo = t * per;
            const int hi = std::min(lo + per, N);
            if (lo >= N) break;
            workers.emplace_back([this, lo, hi, n_passes]() {
                std::vector<uint64_t> tmp;  // per-thread scratch; reused across shards
                for (int i = lo; i < hi; ++i) {
                    Shard& s = shards[i];
                    if (s.raw.empty()) continue;

                    // Radix sort — result ends up in s.raw; tmp is scratch.
                    // Peak extra RAM per thread: sizeof(uint64_t) * largest_shard_size.
                    radix_sort_u64(s.raw, tmp, n_passes);

                    // Reduce: group consecutive same-canon observations into FlatEntry.
                    // Reserve raw/4: typical aDNA is 4-10× coverage so unique k-mers ≈
                    // raw/coverage.  push_back handles any under-estimate via doubling.
                    s.flat.reserve(s.raw.size() / 4);
                    uint64_t prev_canon = ~0ULL;
                    for (const uint64_t p : s.raw) {
                        const uint64_t canon  = p >> N_FLAG_BITS;
                        const uint8_t  flags  = static_cast<uint8_t>(p & 0x7F);
                        const bool     rb_v   = (flags >> 0) & 1;
                        const bool     lb_v   = (flags >> 1) & 1;
                        const bool     is_fwd = (flags >> 2) & 1;
                        const uint8_t  rb     = (flags >> 3) & 3;
                        const uint8_t  lb     = (flags >> 5) & 3;

                        if (canon != prev_canon) {
                            s.flat.push_back({canon, CanonKmerEntry{}});
                            prev_canon = canon;
                        }
                        CanonKmerEntry& e = s.flat.back().edges;
                        if (rb_v) {
                            uint8_t& c = is_fwd ? e.fwd[rb] : e.rev[rb];
                            if (c < 255) ++c;
                        }
                        if (lb_v) {
                            const uint8_t rev_edge = comp4(lb);
                            uint8_t& c = is_fwd ? e.rev[rev_edge] : e.fwd[rev_edge];
                            if (c < 255) ++c;
                        }
                    }
                    s.flat.shrink_to_fit();
                    build_prefix(s);
                    std::vector<uint64_t>().swap(s.raw);  // free raw pages
                }
            });
        }
        for (auto& w : workers) w.join();
        malloc_trim(0);  // return freed raw[] pages to OS
    }

    // Build per-shard prefix index after flat[] is sorted and finalized.
    // Each bucket covers one 1/prefix_size_ slice of the canon key space.
    void build_prefix(Shard& s) const {
        s.prefix.resize(prefix_size_ + 1);
        const uint32_t n = static_cast<uint32_t>(s.flat.size());
        uint32_t pos = 0;
        for (uint32_t bucket = 0; bucket < prefix_size_; ++bucket) {
            s.prefix[bucket] = pos;
            const uint64_t hi = static_cast<uint64_t>(bucket + 1) << prefix_shift_;
            while (pos < n && s.flat[pos].canon < hi) ++pos;
        }
        s.prefix[prefix_size_] = n;
    }

    // Mix all 64 bits of a canon k-mer key into a 32-bit fingerprint.
    // Necessary because k=17 k-mers are only 34 bits — canon >> 32 would give
    // only 2 useful bits (values 0-3), yielding a ~25% false positive rate.
    // This mix spreads all key bits uniformly into the lower 32 bits.
    static uint32_t fp32(uint64_t canon) {
        canon ^= canon >> 30;
        canon *= 0xbf58476d1ce4e5b9ULL;
        canon ^= canon >> 27;
        canon *= 0x94d049bb133111ebULL;
        canon ^= canon >> 31;
        return static_cast<uint32_t>(canon);
    }

    // Build per-shard BBHash MPHFs after finalize().  Replaces binary-search
    // lookup with O(1) MPHF + 32-bit fingerprint verification, reducing pass2
    // memory from ~26 GB (FlatEntry[]) to ~20.5 GB (MPHF + edges + fp).
    // Call once after finalize().  Clears flat[] and prefix[] on completion.
    void build_bbhash(int n_threads) {
#ifdef HAVE_BBHASH
        const int nt  = std::max(1, n_threads);
        const int per = (N + nt - 1) / nt;

        std::vector<std::thread> workers;
        workers.reserve(nt);
        for (int t = 0; t < nt; ++t) {
            const int lo = t * per;
            const int hi = std::min(lo + per, N);
            if (lo >= N) break;
            workers.emplace_back([this, lo, hi]() {
                for (int i = lo; i < hi; ++i) {
                    Shard& s = shards[i];
                    if (s.flat.empty()) continue;

                    const size_t n = s.flat.size();

                    // Extract keys for MPHF construction.
                    std::vector<uint64_t> keys;
                    keys.reserve(n);
                    for (const FlatEntry& e : s.flat)
                        keys.push_back(e.canon);

                    // Build MPHF (1 thread per shard; parallelism is at shard level).
                    auto range = boomphf::range(keys.begin(), keys.end());
                    s.mphf = std::make_unique<ShardMPHF>(n, range, 1, 2.0, false, false);

                    // Fill compact arrays indexed by MPHF rank.
                    s.mphf_edges.resize(n);
                    s.mphf_fp.resize(n);
                    for (const FlatEntry& e : s.flat) {
                        const uint64_t rank = s.mphf->lookup(e.canon);
                        s.mphf_edges[rank] = e.edges;
                        s.mphf_fp[rank]    = fp32(e.canon);
                    }

                    // Free binary-search structures; they're no longer needed.
                    std::vector<FlatEntry>().swap(s.flat);
                    std::vector<uint32_t>().swap(s.prefix);
                }
            });
        }
        for (auto& w : workers) w.join();
        malloc_trim(0);
#endif
    }

    // Lock-free; safe after finalize() — flat[] and prefix[] are immutable.
    // Also safe after build_bbhash() — MPHF arrays are immutable.
    const CanonKmerEntry* find(uint64_t canon) const {
        const Shard& s = shards[shard_of(canon)];

#ifdef HAVE_BBHASH
        if (s.mphf) {
            // O(1) MPHF path: lookup rank, verify 32-bit fingerprint.
            // fp catches both absent keys and any MPHF false associations.
            const uint64_t rank = s.mphf->lookup(canon);
            if (rank >= s.mphf_edges.size()) return nullptr;
            if (s.mphf_fp[rank] != fp32(canon)) return nullptr;
            return &s.mphf_edges[rank];
        }
#endif

        // Binary search fallback (default; BBHash path enabled with --bbhash).
        if (s.flat.empty()) return nullptr;

        const uint32_t bucket = static_cast<uint32_t>(canon >> prefix_shift_);
        const uint32_t lo = s.prefix[bucket];
        const uint32_t hi = s.prefix[bucket + 1];
        if (lo == hi) return nullptr;

        const FlatEntry* first = s.flat.data() + lo;
        const FlatEntry* last  = s.flat.data() + hi;
        const FlatEntry* it = std::lower_bound(first, last, canon,
            [](const FlatEntry& e, uint64_t k) { return e.canon < k; });
        return (it != last && it->canon == canon) ? &it->edges : nullptr;
    }

    bool empty() const {
        for (int i = 0; i < N; ++i) {
            const Shard& s = shards[i];
#ifdef HAVE_BBHASH
            if (s.mphf || !s.flat.empty()) return false;
#else
            if (!s.flat.empty()) return false;
#endif
        }
        return true;
    }

    size_t size() const {
        size_t n = 0;
        for (int i = 0; i < N; ++i) {
            const Shard& s = shards[i];
#ifdef HAVE_BBHASH
            n += s.mphf ? s.mphf_edges.size() : s.flat.size();
#else
            n += s.flat.size();
#endif
        }
        return n;
    }

    // Raw observation count (valid during pass1 accumulation, before finalize).
    size_t raw_size() const {
        size_t n = 0;
        for (int i = 0; i < N; ++i) n += shards[i].raw.size();
        return n;
    }
};

// ============================================================================
// Extension walk — canonical orientation tracking, rolling RC (O(1) per step).
// ============================================================================
//
// Walk right from start_kmer for up to max_steps steps.
// Stops at branch, dead end, or reciprocal non-uniqueness (unitig boundary).
// Returns the number of appended bases; written to out[].
//
// Canonical storage: we track both cur (the k-mer in walking direction) and
// its canonical form + orientation (is_fwd).  The forward-edge array and the
// predecessor-edge array swap based on orientation — see struct doc above.
//
// Damage-aware branch resolution: see inline comments.
// n_skip: retrace steps (bases inside the read between anchor k-mer and terminal).
// is_rc_walk: true for 5' extension via RC walk (use 5' damage params).
static int walk_right(const ShardedKmerTable& table,
                      uint64_t start_kmer,
                      int k, int max_steps, uint32_t min_count, char* out,
                      const DamageProfile* prof = nullptr,
                      bool is_rc_walk = false,
                      int n_skip = 0) {
    const uint64_t kmask     = (1ULL << (2 * k)) - 1;
    const int      msb_shift = 2 * (k - 1);

    // Initialise canonical orientation for the starting k-mer.
    uint64_t rc_cur    = kmer_rc(start_kmer, k);
    uint64_t canon_cur = (start_kmer < rc_cur) ? start_kmer : rc_cur;
    bool     is_fwd    = (start_kmer == canon_cur);
    uint64_t cur       = start_kmer;

    int n = 0;

    for (int step = 0; step < max_steps; ++step) {
        const CanonKmerEntry* entry = table.find(canon_cur);
        if (!entry) break;

        const uint8_t* edges = is_fwd ? entry->fwd : entry->rev;

        int     n_supported = 0;
        uint8_t best_b      = 0;
        for (uint8_t b = 0; b < 4; ++b) {
            if (edges[b] >= min_count) { ++n_supported; best_b = b; }
        }

        // Damage-aware branch resolution (retrace zone only, step < n_skip).
        bool damage_resolved = false;
        if (n_supported == 2 && prof && prof->enabled && step < n_skip) {
            uint8_t b0 = 0xFF, b1 = 0xFF;
            for (uint8_t b = 0; b < 4; ++b) {
                if (edges[b] >= min_count) {
                    if (b0 == 0xFF) b0 = b; else b1 = b;
                }
            }
            // Damage pairs: {A=0, G=2} (G→A deamination) or {C=1, T=3} (C→T deamination).
            if ((b0 == 0 && b1 == 2) || (b0 == 1 && b1 == 3)) {
                const uint8_t undmg_b = (b0 == 0) ? 2u : 1u;
                const uint8_t dmg_b   = (b0 == 0) ? 0u : 3u;
                const int c_undmg = edges[undmg_b];
                const int c_dmg   = edges[dmg_b];
                const int dist    = n_skip - 1 - step;
                const double d_max  = is_rc_walk ? prof->d_max_5prime  : prof->d_max_3prime;
                const double lambda = is_rc_walk ? prof->lambda_5prime : prof->lambda_3prime;
                const double p_dmg  = d_max * std::exp(-lambda * static_cast<double>(dist))
                                      + prof->background;
                const double tau    = std::min(0.25, p_dmg + 0.05);
                const double ratio  = static_cast<double>(c_dmg) / (c_undmg + c_dmg);
                if (ratio <= tau && c_undmg >= 2 * c_dmg) {
                    best_b        = undmg_b;
                    damage_resolved = true;
                }
            }
        }

        if (!damage_resolved && n_supported != 1) break;

        // Advance to next k-mer, maintaining canonical orientation.
        const uint64_t next_cur   = advance_right(cur, best_b, kmask);
        const uint64_t next_rc    = (rc_cur >> 2) |
                                    (static_cast<uint64_t>(comp4(best_b)) << msb_shift);
        const uint64_t next_canon = (next_cur < next_rc) ? next_cur : next_rc;
        const bool     next_is_fwd = (next_cur == next_canon);

        // Reciprocal uniqueness: predecessors of next k-mer must be unambiguous.
        // Predecessors of next_cur = reverse-direction edges of next_canon:
        //   if next_is_fwd: predecessors come from rev array (RC approaches)
        //   else:           predecessors come from fwd array (canonical approaches)
        if (!damage_resolved) {
            const CanonKmerEntry* next_entry = table.find(next_canon);
            if (!next_entry) break;
            const uint8_t* pred_edges = next_is_fwd ? next_entry->rev : next_entry->fwd;
            int n_pred = 0;
            for (uint8_t b = 0; b < 4; ++b)
                if (pred_edges[b] >= min_count) ++n_pred;
            if (n_pred != 1) break;
        }

        out[n++] = kDec4[best_b];
        cur       = next_cur;
        rc_cur    = next_rc;
        canon_cur = next_canon;
        is_fwd    = next_is_fwd;
    }
    return n;
}

// ============================================================================
// Damage mask helpers
// ============================================================================

static void compute_mask_lengths(const DamageProfile& prof, int L,
                                  int& mask5, int& mask3) {
    if (!prof.enabled) { mask5 = 0; mask3 = 0; return; }
    mask5 = 0;
    while (mask5 < L && mask5 < DamageProfile::MASK_POSITIONS && prof.mask_pos[mask5])
        ++mask5;
    mask3 = 0;
    while (mask3 < L - mask5 && mask3 < DamageProfile::MASK_POSITIONS && prof.mask_pos[mask3])
        ++mask3;
}

// ============================================================================
// ExtendEngine
// ============================================================================

struct ExtendConfig {
    int      k               = 17;   // k=17 fits short aDNA reads; k=25 fails for reads <37bp with masking
    uint32_t min_count       = 2;
    int      max_extend      = 100;  // Tadpole production uses el=er=100
    int      min_qual        = 20;
    bool     no_damage       = false;
    bool     use_bbhash      = false; // opt-in: build BBHash MPHF after finalize().
                                      // Saves ~6 GB RAM at DS4 scale but is slower than the
                                      // default prefix-indexed binary search for all tested
                                      // dataset sizes.  Enable with --bbhash when RAM is tight.
    int      mask_5_override = -1;   // >= 0: use directly, skip Pass 0
    int      mask_3_override = -1;
    double   mask_threshold  = 0.05;
    int64_t  damage_sample   = 500000;  // 0 = full scan; default samples 500k reads
    dart::SampleDamageProfile::LibraryType library_type =
        dart::SampleDamageProfile::LibraryType::UNKNOWN;
    int      n_threads       = 1;
};

class ExtendEngine {
public:
    ExtendEngine(const std::string& in_path, const std::string& out_path,
                 const ExtendConfig& cfg)
        : in_path_(in_path), out_path_(out_path), cfg_(cfg) {}

    int run();

private:
    std::string  in_path_, out_path_;
    ExtendConfig cfg_;
    ShardedKmerTable table_;
    DamageProfile profile_;
    // Asymmetric mask lengths; -1 = derive from profile_.
    int override_mask5_ = -1;
    int override_mask3_ = -1;
    // Walk scratch buffer size: max_extend + up to ~500bp read length for skip region.
    static constexpr int kWalkBufBase = 512;
    int walk_tmp_sz_ = 0;

    void get_mask_lengths(int L, int& m5, int& m3) const {
        if (override_mask5_ >= 0) {
            m5 = std::min(override_mask5_, L);
            m3 = std::min(override_mask3_, L - m5);
        } else {
            compute_mask_lengths(profile_, L, m5, m3);
        }
    }

    void pass0_estimate_damage();
    void pass1_build_graph();
    void pass2_extend_write();
};

void ExtendEngine::pass0_estimate_damage() {
    if (cfg_.no_damage) {
        log_info("extend: damage masking disabled (--no-damage)");
        return;
    }

    if (cfg_.mask_5_override >= 0 && cfg_.mask_3_override >= 0) {
        override_mask5_ = std::min(cfg_.mask_5_override, DamageProfile::MASK_POSITIONS);
        override_mask3_ = std::min(cfg_.mask_3_override, DamageProfile::MASK_POSITIONS);
        log_info("extend: manual mask override: 5'=" + std::to_string(override_mask5_) +
                 "  3'=" + std::to_string(override_mask3_) + " bp");
        return;
    }

    profile_ = estimate_damage(in_path_, cfg_.mask_threshold, cfg_.library_type,
                               cfg_.damage_sample);

    int m5 = 0, m3 = 0;
    get_mask_lengths(1000, m5, m3);
    log_info("extend: global mask lengths: 5'=" + std::to_string(m5) +
             "  3'=" + std::to_string(m3) + " bp");
}


void ExtendEngine::pass1_build_graph() {
    const int k = cfg_.k;
    if (2 * k + ShardedKmerTable::N_FLAG_BITS > 64)
        throw std::runtime_error("k=" + std::to_string(k) +
            " exceeds packed uint64_t capacity (max k=28 for sort-reduce build)");

    table_.init_lookup_layout(k);

    log_info("extend pass1: building canonical k-mer graph (k=" + std::to_string(k) +
             " min_count=" + std::to_string(cfg_.min_count) +
             " threads=" + std::to_string(cfg_.n_threads) + ")");

    // Pre-reserve raw[] per shard from file-size estimate.
    // Linux demand paging: uncommitted reserved pages cost no RSS.
    table_.reserve_raw(in_path_, k);

    const int      min_qual  = cfg_.min_qual;
    const int      nt        = cfg_.n_threads;
    const uint64_t kmask     = (1ULL << (2 * k)) - 1;
    const int      msb_shift = 2 * (k - 1);

    // Producer-consumer pipeline: reader fills batches, nt workers consume.
    const int BATCH_CAP = 4000;
    struct Batch {
        std::vector<FastqRecord> recs;
        int sz = 0;
        Batch() : recs(BATCH_CAP) {}
    };
    const int N_BUFS = nt + 4;
    std::vector<Batch> bufs(N_BUFS);

    std::vector<int>        free_slots(N_BUFS);
    std::iota(free_slots.begin(), free_slots.end(), 0);
    std::mutex              free_mtx;
    std::condition_variable free_cv;

    std::queue<int>         work_q;
    std::mutex              work_mtx;
    std::condition_variable work_cv;
    bool reader_done = false;

    std::atomic<int64_t> n_kmers{0};

    auto build_worker = [&]() {
        // Collect all k-mer observations as packed uint64_t (canon << 7 | flags),
        // counting-sort by shard, then batch_append() — one lock per shard per batch.
        // No hash map, no rehash.  Observations accumulate in per-shard raw[] vectors
        // and are sorted+reduced into FlatEntry[] by finalize() after all reads.
        std::vector<uint64_t> items;
        std::vector<uint64_t> sorted;
        items.reserve(BATCH_CAP * 42);
        sorted.reserve(BATCH_CAP * 42);
        std::array<int, ShardedKmerTable::N + 1> cnt{};

        for (;;) {
            int slot;
            {
                std::unique_lock<std::mutex> lk(work_mtx);
                work_cv.wait(lk, [&]{ return !work_q.empty() || reader_done; });
                if (work_q.empty()) return;
                slot = work_q.front(); work_q.pop();
            }
            Batch& b = bufs[slot];
            items.clear();
            int64_t local_kmers = 0;

            // ── Phase 1: collect all k-mers (no locking) ─────────────────
            for (int ri = 0; ri < b.sz; ++ri) {
                const FastqRecord& rec = b.recs[ri];
                const int L = static_cast<int>(rec.seq.size());
                if (L < k) continue;
                const char* seq  = rec.seq.data();
                const char* qual = rec.qual.data();

                int mask5, mask3;
                get_mask_lengths(L, mask5, mask3);
                const int first_kmer = mask5;
                const int last_kmer  = L - mask3 - k;
                if (first_kmer > last_kmer) continue;

                uint64_t kmer_fwd = 0, kmer_rev = 0;
                bool valid = false;
                int  last_invalid = -1;

                for (int p = first_kmer; p <= last_kmer; ++p) {
                    if (p == first_kmer) {
                        valid = true; kmer_fwd = 0;
                        for (int j = 0; j < k; ++j) {
                            uint8_t base = enc4(seq[p + j]);
                            int q = static_cast<int>(static_cast<unsigned char>(qual[p + j])) - 33;
                            if (base == 0xFF || q < min_qual) { valid = false; last_invalid = p + j; break; }
                            kmer_fwd = (kmer_fwd << 2) | base;
                        }
                        if (valid) kmer_rev = kmer_rc(kmer_fwd, k);
                    } else if (valid) {
                        int new_pos = p + k - 1;
                        uint8_t nb = enc4(seq[new_pos]);
                        int nq = static_cast<int>(static_cast<unsigned char>(qual[new_pos])) - 33;
                        if (nb == 0xFF || nq < min_qual) {
                            valid = false; last_invalid = new_pos;
                        } else {
                            kmer_fwd = ((kmer_fwd << 2) | nb) & kmask;
                            kmer_rev = (kmer_rev >> 2) |
                                       (static_cast<uint64_t>(comp4(nb)) << msb_shift);
                        }
                    } else if (last_invalid == p - 1) {
                        valid = true; kmer_fwd = 0;
                        for (int j = 0; j < k; ++j) {
                            uint8_t base = enc4(seq[p + j]);
                            int q = static_cast<int>(static_cast<unsigned char>(qual[p + j])) - 33;
                            if (base == 0xFF || q < min_qual) { valid = false; last_invalid = p + j; break; }
                            kmer_fwd = (kmer_fwd << 2) | base;
                        }
                        if (valid) kmer_rev = kmer_rc(kmer_fwd, k);
                    }
                    if (!valid) continue;
                    if (count_distinct_dimers(kmer_fwd, k) < 3) continue;
                    ++local_kmers;

                    uint8_t rb = 0; bool rb_valid = false;
                    const int right_pos = p + k;
                    if (right_pos < L - mask3) {
                        uint8_t rb_raw = enc4(seq[right_pos]);
                        if (rb_raw != 0xFF) {
                            int rq = static_cast<int>(static_cast<unsigned char>(qual[right_pos])) - 33;
                            if (rq >= min_qual) { rb = rb_raw; rb_valid = true; }
                        }
                    }
                    uint8_t lb = 0; bool lb_valid = false;
                    const int left_pos = p - 1;
                    if (left_pos >= mask5) {
                        uint8_t lb_raw = enc4(seq[left_pos]);
                        if (lb_raw != 0xFF) {
                            int lq = static_cast<int>(static_cast<unsigned char>(qual[left_pos])) - 33;
                            if (lq >= min_qual) { lb = lb_raw; lb_valid = true; }
                        }
                    }

                    const uint64_t canon  = (kmer_fwd < kmer_rev) ? kmer_fwd : kmer_rev;
                    const bool     is_fwd = (canon == kmer_fwd);
                    items.push_back((canon << ShardedKmerTable::N_FLAG_BITS) |
                                    (rb_valid ? 1ULL : 0ULL) |
                                    (lb_valid ? 2ULL : 0ULL) |
                                    (is_fwd   ? 4ULL : 0ULL) |
                                    (static_cast<uint64_t>(rb) << 3) |
                                    (static_cast<uint64_t>(lb) << 5));
                }
            }
            n_kmers.fetch_add(local_kmers, std::memory_order_relaxed);

            // ── Phase 2: counting sort by shard (O(N + S), no comparison) ─
            const int n = static_cast<int>(items.size());
            cnt.fill(0);
            for (uint64_t p : items)
                ++cnt[ShardedKmerTable::shard_of(p >> ShardedKmerTable::N_FLAG_BITS) + 1];
            for (int i = 1; i <= ShardedKmerTable::N; ++i) cnt[i] += cnt[i - 1];
            // Save start positions before scatter mutates cnt[].
            const std::array<int, ShardedKmerTable::N + 1> shard_start = cnt;
            sorted.resize(n);
            for (uint64_t p : items)
                sorted[cnt[ShardedKmerTable::shard_of(p >> ShardedKmerTable::N_FLAG_BITS)]++] = p;

            // ── Phase 3: append to per-shard raw[] (one lock per shard) ──
            table_.batch_append(sorted.data(), shard_start.data());

            {
                std::lock_guard<std::mutex> lk(free_mtx);
                free_slots.push_back(slot);
            }
            free_cv.notify_one();
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nt);
    for (int t = 0; t < nt; ++t) workers.emplace_back(build_worker);

    int64_t n_reads = 0;
    auto reader = make_fastq_reader(in_path_, static_cast<size_t>(cfg_.n_threads));
    for (;;) {
        int slot;
        {
            std::unique_lock<std::mutex> lk(free_mtx);
            free_cv.wait(lk, [&]{ return !free_slots.empty(); });
            slot = free_slots.back(); free_slots.pop_back();
        }
        Batch& b = bufs[slot];
        b.sz = 0;
        while (b.sz < BATCH_CAP && reader->read(b.recs[b.sz])) ++b.sz;
        n_reads += b.sz;
        if (b.sz == 0) {
            std::lock_guard<std::mutex> lk(free_mtx);
            free_slots.push_back(slot);
            break;
        }
        if (n_reads % 10000000 < BATCH_CAP)
            std::cerr << "\r[extend pass1] " << n_reads << " reads" << std::flush;
        {
            std::lock_guard<std::mutex> lk(work_mtx);
            work_q.push(slot);
        }
        work_cv.notify_one();
    }
    {
        std::lock_guard<std::mutex> lk(work_mtx);
        reader_done = true;
    }
    work_cv.notify_all();
    for (auto& w : workers) w.join();
    std::cerr << "\r";

    log_info("extend pass1: " + std::to_string(n_reads) + " reads"
             "  " + std::to_string(n_kmers.load()) + " k-mer positions accumulated"
             "  (" + std::to_string(table_.raw_size()) + " raw observations in shard buffers)");
}

void ExtendEngine::pass2_extend_write() {
    log_info("extend pass2: extending reads"
             " (max_extend=" + std::to_string(cfg_.max_extend) +
             " min_count=" + std::to_string(cfg_.min_count) +
             " threads=" + std::to_string(cfg_.n_threads) + ")");

    const int      k         = cfg_.k;
    const int      nt        = cfg_.n_threads;
    const uint64_t kmask     = (1ULL << (2 * k)) - 1;
    const int      msb_shift = 2 * (k - 1);

    // ── Chunk/task types ─────────────────────────────────────────────────────
    struct ExtOut {
        std::array<char, 256> left_ext{};
        std::array<char, 256> right_buf{};
        int left_len  = 0;
        int right_len = 0;
    };

    // CHUNK_READS: records per chunk; TASK_GRAIN: records per fine-grained task.
    // IN_FLIGHT: max chunks simultaneously in-flight across all queues.
    const int CHUNK_READS = 64000;
    const int TASK_GRAIN  = 2000;
    const int IN_FLIGHT   = std::max(4, nt / 2);

    struct Chunk {
        uint64_t                 batch_id = 0;
        uint32_t                 n        = 0;
        std::vector<FastqRecord> recs;
        std::vector<ExtOut>      outs;
        std::atomic<uint32_t>    tasks_remaining{0};
        Chunk() : recs(CHUNK_READS), outs(CHUNK_READS) {}
    };
    struct Task {
        std::shared_ptr<Chunk> chunk;
        uint32_t begin, end;
    };
    struct ChunkResult {
        uint64_t               batch_id;
        std::shared_ptr<Chunk> chunk;
    };

    // ── Queues ───────────────────────────────────────────────────────────────
    const int max_tasks = IN_FLIGHT * (CHUNK_READS / TASK_GRAIN + 1);
    BoundedQueue<std::shared_ptr<Chunk>> free_chunks(IN_FLIGHT + 2);
    BoundedQueue<Task>                   work_queue(max_tasks);
    BoundedQueue<ChunkResult>            done_queue(IN_FLIGHT + 2);

    for (int i = 0; i < IN_FLIGHT + 2; ++i)
        free_chunks.push(std::make_shared<Chunk>());

    // ── Stats (writer thread populates, main thread reads after join) ────────
    std::atomic<int64_t> N{0};
    std::atomic<int64_t> n_extended{0};
    std::atomic<int64_t> total_ext5{0};
    std::atomic<int64_t> total_ext3{0};

    // ── Reader thread ────────────────────────────────────────────────────────
    auto reader_fn = [&]() {
        uint64_t next_batch_id = 0;
        auto reader = make_fastq_reader(in_path_, static_cast<size_t>(cfg_.n_threads));
        std::shared_ptr<Chunk> chunk;

        for (;;) {
            if (!free_chunks.pop(chunk)) break;

            chunk->batch_id = next_batch_id++;
            chunk->n = 0;
            for (int i = 0; i < CHUNK_READS; ++i) {
                if (!reader->read(chunk->recs[i])) break;
                chunk->outs[i].left_len = 0;
                chunk->outs[i].right_len = 0;
                ++chunk->n;
            }
            if (chunk->n == 0) {
                free_chunks.push(std::move(chunk));
                break;
            }

            N.fetch_add(chunk->n, std::memory_order_relaxed);
            const uint32_t n_tasks = (chunk->n + TASK_GRAIN - 1) / TASK_GRAIN;
            chunk->tasks_remaining.store(n_tasks, std::memory_order_release);

            for (uint32_t t = 0; t < n_tasks; ++t) {
                uint32_t beg = t * static_cast<uint32_t>(TASK_GRAIN);
                uint32_t end = std::min(beg + static_cast<uint32_t>(TASK_GRAIN), chunk->n);
                work_queue.push(Task{chunk, beg, end});
            }
        }
        work_queue.close();
    };

    // ── Worker threads ───────────────────────────────────────────────────────
    std::atomic<int> workers_left{nt};

    auto worker_fn = [&]() {
        const int walk_tmp_sz = walk_tmp_sz_;
        std::vector<char> right_tmp(walk_tmp_sz);
        std::vector<char> left_tmp(walk_tmp_sz);
        std::array<char, 256> right_buf_local{};
        std::array<char, 256> left_ext_local{};

        Task task;
        while (work_queue.pop(task)) {
            Chunk& c = *task.chunk;

            for (uint32_t i = task.begin; i < task.end; ++i) {
                const FastqRecord& rec = c.recs[i];
                const int L = static_cast<int>(rec.seq.size());
                int right_len = 0, left_len = 0;

                if (L >= k && !table_.empty()) {
                    int mask5, mask3;
                    get_mask_lengths(L, mask5, mask3);

                    // ── 3' extension ──────────────────────────────────────
                    {
                        uint64_t kmer_fwd = 0;
                        bool valid = false;

                        for (int p = L - k; p >= 0; --p) {
                            if (p == L - k) {
                                valid = encode_kmer(rec.seq.data(), p, k, kmer_fwd);
                            } else if (valid) {
                                uint8_t new_l = enc4(rec.seq[p]);
                                if (new_l == 0xFF) { valid = false; }
                                else {
                                    kmer_fwd = (kmer_fwd >> 2) |
                                               (static_cast<uint64_t>(new_l) << msb_shift);
                                }
                            } else {
                                valid = encode_kmer(rec.seq.data(), p, k, kmer_fwd);
                            }
                            if (!valid) continue;

                            const int n_skip  = L - (p + k);
                            const int budget  = std::min(cfg_.max_extend + n_skip, walk_tmp_sz - 1);
                            const int n_walk  = walk_right(table_, kmer_fwd, k, budget,
                                                           cfg_.min_count, right_tmp.data(),
                                                           &profile_, false, n_skip);
                            right_len = std::max(0, n_walk - n_skip);
                            if (right_len > 0) {
                                std::memcpy(right_buf_local.data(), right_tmp.data() + n_skip,
                                            right_len);
                                break;
                            }
                        }
                    }

                    // ── 5' extension ──────────────────────────────────────
                    {
                        uint64_t kmer_fwd2 = 0, kmer_rv = 0;
                        bool valid = false;

                        for (int p = 0; p <= L - k; ++p) {
                            if (p == 0) {
                                valid = encode_kmer(rec.seq.data(), 0, k, kmer_fwd2);
                                if (valid) kmer_rv = kmer_rc(kmer_fwd2, k);
                            } else if (valid) {
                                uint8_t new_r = enc4(rec.seq[p + k - 1]);
                                if (new_r == 0xFF) { valid = false; }
                                else {
                                    kmer_fwd2 = ((kmer_fwd2 << 2) | new_r) & kmask;
                                    kmer_rv   = (kmer_rv >> 2) |
                                                (static_cast<uint64_t>(comp4(new_r)) << msb_shift);
                                }
                            } else {
                                valid = encode_kmer(rec.seq.data(), p, k, kmer_fwd2);
                                if (valid) kmer_rv = kmer_rc(kmer_fwd2, k);
                            }
                            if (!valid) continue;

                            const int n_skip  = p;
                            const int budget  = std::min(cfg_.max_extend + n_skip, walk_tmp_sz - 1);
                            const int n_walk  = walk_right(table_, kmer_rv, k, budget,
                                                           cfg_.min_count, left_tmp.data(),
                                                           &profile_, true, n_skip);
                            left_len = std::max(0, n_walk - n_skip);
                            if (left_len > 0) {
                                for (int ii = 0; ii < left_len; ++ii)
                                    left_ext_local[ii] =
                                        kDec4[comp4(enc4(left_tmp[n_walk - 1 - ii]))];
                                break;
                            }
                        }
                    }
                }

                c.outs[i].left_len  = left_len;
                c.outs[i].right_len = right_len;
                if (left_len  > 0) std::memcpy(c.outs[i].left_ext.data(), left_ext_local.data(), left_len);
                if (right_len > 0) std::memcpy(c.outs[i].right_buf.data(), right_buf_local.data(), right_len);
            }

            // Last task to finish signals chunk completion.
            if (task.chunk->tasks_remaining.fetch_sub(1, std::memory_order_acq_rel) == 1)
                done_queue.push(ChunkResult{task.chunk->batch_id, task.chunk});
        }

        if (workers_left.fetch_sub(1, std::memory_order_acq_rel) == 1)
            done_queue.close();
    };

    // ── Writer thread — ordered resequencer ─────────────────────────────────
    auto writer_fn = [&]() {
        const bool compress = (out_path_.size() > 3 &&
                               out_path_.substr(out_path_.size() - 3) == ".gz");
        FastqWriter fqw(out_path_, compress, cfg_.n_threads);

        int64_t local_ext = 0, local_ext5 = 0, local_ext3 = 0;
        uint64_t next_to_write = 0;
        int64_t  report_n      = 0;
        std::unordered_map<uint64_t, std::shared_ptr<Chunk>> pending;

        ChunkResult result;
        while (done_queue.pop(result)) {
            pending.emplace(result.batch_id, std::move(result.chunk));

            for (;;) {
                auto it = pending.find(next_to_write);
                if (it == pending.end()) break;

                Chunk& c = *it->second;
                for (uint32_t i = 0; i < c.n; ++i) {
                    const FastqRecord& rec = c.recs[i];
                    const ExtOut&      o   = c.outs[i];
                    const int L = static_cast<int>(rec.seq.size());

                    // Compute how many terminal positions are masked for this read.
                    // The original read's damaged bases at those positions vary between
                    // reads from the same molecule.  Replacing C/T (5') with \x01 and
                    // G/A (3') with \x02 makes the fingerprint identical for all copies,
                    // so derep_pairs can collapse them correctly.
                    int mask5 = 0, mask3 = 0;
                    if (profile_.enabled) {
                        while (mask5 < L && mask5 < DamageProfile::MASK_POSITIONS
                               && profile_.mask_pos[mask5]) ++mask5;
                        while (mask3 < L - mask5 && mask3 < DamageProfile::MASK_POSITIONS
                               && profile_.mask_pos[mask3]) ++mask3;
                    }

                    // Neutralize damage-masked positions in a sequence buffer.
                    // left_off: offset of the original read within the output buffer.
                    auto neutralize = [&](std::string& seq, int left_off) {
                        for (int p = 0; p < mask5; ++p) {
                            char& b = seq[left_off + p];
                            if (b == 'C' || b == 'T' || b == 'c' || b == 't') b = '\x01';
                        }
                        for (int p = 0; p < mask3; ++p) {
                            char& b = seq[left_off + L - mask3 + p];
                            if (b == 'G' || b == 'A' || b == 'g' || b == 'a') b = '\x02';
                        }
                    };

                    if (o.left_len > 0 || o.right_len > 0) {
                        ++local_ext;
                        local_ext5 += o.left_len;
                        local_ext3 += o.right_len;

                        FastqRecord out_rec;
                        out_rec.header = rec.header;
                        out_rec.plus   = rec.plus;
                        const int new_len = o.left_len + L + o.right_len;
                        out_rec.seq.resize(new_len);
                        out_rec.qual.resize(new_len);
                        std::memcpy(out_rec.seq.data(),                   o.left_ext.data(),  o.left_len);
                        std::memcpy(out_rec.seq.data() + o.left_len,      rec.seq.data(),     L);
                        std::memcpy(out_rec.seq.data() + o.left_len + L,  o.right_buf.data(), o.right_len);
                        std::memset(out_rec.qual.data(),                  '#',                o.left_len);
                        std::memcpy(out_rec.qual.data() + o.left_len,     rec.qual.data(),    L);
                        std::memset(out_rec.qual.data() + o.left_len + L, '#',                o.right_len);
                        neutralize(out_rec.seq, o.left_len);
                        fqw.write(out_rec);
                    } else if (mask5 > 0 || mask3 > 0) {
                        FastqRecord out_rec = rec;
                        neutralize(out_rec.seq, 0);
                        fqw.write(out_rec);
                    } else {
                        fqw.write(rec);
                    }
                }

                report_n += c.n;
                std::cerr << "\r[extend pass2] " << report_n << " reads processed" << std::flush;

                free_chunks.push(std::move(it->second));
                pending.erase(it);
                ++next_to_write;
            }
        }

        n_extended.store(local_ext,  std::memory_order_relaxed);
        total_ext5.store(local_ext5, std::memory_order_relaxed);
        total_ext3.store(local_ext3, std::memory_order_relaxed);
    };

    // ── Launch reader, workers, writer; wait ─────────────────────────────────
    std::thread reader_thread(reader_fn);

    std::vector<std::thread> worker_threads;
    worker_threads.reserve(nt);
    for (int t = 0; t < nt; ++t) worker_threads.emplace_back(worker_fn);

    std::thread writer_thread(writer_fn);

    reader_thread.join();
    for (auto& th : worker_threads) th.join();
    writer_thread.join();

    std::cerr << "\r";

    const int64_t Nv   = N.load();
    const int64_t next = n_extended.load();
    const int64_t e5   = total_ext5.load();
    const int64_t e3   = total_ext3.load();

    const double frac_ext = Nv   > 0 ? 100.0 * next / Nv   : 0.0;
    const double avg5     = next > 0 ? static_cast<double>(e5) / next : 0.0;
    const double avg3     = next > 0 ? static_cast<double>(e3) / next : 0.0;

    log_info("extend pass2: " + std::to_string(Nv) + " reads"
             "  extended=" + std::to_string(next)
             + " (" + std::to_string(static_cast<int>(frac_ext + 0.5)) + "%)"
             "  avg_5ext=" + std::to_string(avg5).substr(0, 4)
             + " avg_3ext=" + std::to_string(avg3).substr(0, 4) + " bp");
}

int ExtendEngine::run() {
    pass0_estimate_damage();
    pass1_build_graph();
    log_info("extend: finalizing k-mer table to sorted flat arrays...");
    table_.finalize(cfg_.n_threads);
    const size_t n_entries = table_.size();
    log_info("extend: finalized — " + std::to_string(n_entries) + " entries, memory released");
#ifdef HAVE_BBHASH
    if (cfg_.use_bbhash) {
        log_info("extend: building BBHash MPHF for O(1) pass2 lookup...");
        table_.build_bbhash(cfg_.n_threads);
        log_info("extend: BBHash ready — " + std::to_string(n_entries) + " entries");
    }
#endif
    walk_tmp_sz_ = cfg_.max_extend + kWalkBufBase;
    pass2_extend_write();
    return 0;
}

} // namespace

// ============================================================================
// CLI entry point
// ============================================================================

static void extend_usage(const char* prog) {
    std::cerr
        << "Usage: fqdup " << prog << " -i INPUT -o OUTPUT [options]\n"
        << "\nOptions:\n"
        << "  -i FILE                     Input FASTQ (.gz or plain); required\n"
        << "  -o FILE                     Output FASTQ (.gz or plain); required\n"
        << "  -k INT                      k-mer size, odd 5-31 (default: 17)\n"
        << "  --min-count INT             Min edge support to extend (default: 2)\n"
        << "  --max-extend INT            Max bases added per side (default: 100)\n"
        << "  --no-damage                 Skip damage estimation, no masking\n"
        << "  --bbhash                    Build BBHash MPHF for O(1) lookup; saves ~6 GB RAM at\n"
           << "                              DS4 scale but is slower for most datasets (default: off)\n"
        << "  --mask-5 INT                Override 5' mask length; must pair with --mask-3\n"
        << "  --mask-3 INT                Override 3' mask length; must pair with --mask-5\n"
        << "  --mask-threshold FLOAT      Excess damage threshold (default: 0.05)\n"
        << "  --damage-sample INT         Reads to sample for damage estimation (default: 500000; 0=all)\n"
        << "  --min-qual INT              Exclude bases below this Phred score (default: 20)\n"
        << "  --library-type auto|ds|ss   Library type for damage model (default: auto)\n"
        << "  --threads INT               Worker threads (default: 1)\n"
        << "\nNotes:\n"
        << "  - stdin is not supported (algorithm requires multiple passes over the file)\n"
        << "  - reads with no clean interior k-mers are written unchanged\n"
        << "  - extended bases get quality '#' (Phred 2); originals keep their quality\n"
        << "  - k must be odd (even k allows palindromic k-mers that corrupt orientation)\n";
}

int extend_main(int argc, char** argv) {
    std::string  in_path, out_path;
    ExtendConfig cfg;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            in_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            cfg.k = std::stoi(argv[++i]);
        } else if (arg == "--min-count" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 1) { std::cerr << "Error: --min-count must be >= 1, got " << v << "\n"; return 1; }
            cfg.min_count = static_cast<uint32_t>(v);
        } else if (arg == "--max-extend" && i + 1 < argc) {
            cfg.max_extend = std::stoi(argv[++i]);
        } else if (arg == "--no-damage") {
            cfg.no_damage = true;
        } else if (arg == "--bbhash") {
            cfg.use_bbhash = true;
        } else if (arg == "--mask-5" && i + 1 < argc) {
            cfg.mask_5_override = std::stoi(argv[++i]);
        } else if (arg == "--mask-3" && i + 1 < argc) {
            cfg.mask_3_override = std::stoi(argv[++i]);
        } else if (arg == "--mask-threshold" && i + 1 < argc) {
            cfg.mask_threshold = std::stod(argv[++i]);
        } else if (arg == "--damage-sample" && i + 1 < argc) {
            cfg.damage_sample = static_cast<int64_t>(std::stoll(argv[++i]));
        } else if (arg == "--min-qual" && i + 1 < argc) {
            cfg.min_qual = std::stoi(argv[++i]);
        } else if (arg == "--library-type" && i + 1 < argc) {
            std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                cfg.library_type = dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                cfg.library_type = dart::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt != "auto") {
                std::cerr << "Error: unknown --library-type: " << lt << "\n";
                return 1;
            }
        } else if (arg == "--threads" && i + 1 < argc) {
            cfg.n_threads = std::stoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            extend_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Error: unknown argument: " << arg << "\n";
            extend_usage(argv[0]);
            return 1;
        }
    }

    if (in_path.empty()) {
        std::cerr << "Error: -i INPUT required\n";
        extend_usage(argv[0]);
        return 1;
    }
    if (out_path.empty()) {
        std::cerr << "Error: -o OUTPUT required\n";
        extend_usage(argv[0]);
        return 1;
    }
    if (cfg.k < 5 || cfg.k > 31) {
        std::cerr << "Error: -k must be in [5, 31]\n";
        return 1;
    }
    if (cfg.k % 2 == 0) {
        std::cerr << "Error: -k must be odd (even k allows self-RC palindromes that "
                     "corrupt forward/RC edge orientation)\n";
        return 1;
    }
    if (cfg.max_extend < 1 || cfg.max_extend > 255) {
        std::cerr << "Error: --max-extend must be in [1, 255]\n";
        return 1;
    }
    if (cfg.n_threads < 1) {
        std::cerr << "Error: --threads must be >= 1\n";
        return 1;
    }
    if (cfg.mask_threshold <= 0.0 || cfg.mask_threshold >= 1.0) {
        std::cerr << "Error: --mask-threshold must be in (0, 1), got " << cfg.mask_threshold << "\n";
        return 1;
    }
    if ((cfg.mask_5_override >= 0) != (cfg.mask_3_override >= 0)) {
        std::cerr << "Error: --mask-5 and --mask-3 must be specified together\n";
        return 1;
    }

    init_logger("fqdup-extend.log");
    log_info("fqdup extend: " + in_path + " \xe2\x86\x92 " + out_path);

    try {
        ExtendEngine engine(in_path, out_path, cfg);
        int rc = engine.run();
        shutdown_logger();
        return rc;
    } catch (const std::exception& e) {
        log_info("Error: " + std::string(e.what()));
        std::cerr << "Error: " << e.what() << "\n";
        shutdown_logger();
        return 1;
    }
}
