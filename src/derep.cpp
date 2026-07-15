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
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include "taph/frame_selector_decl.hpp"

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

#include "fqdup/lr_score.hpp"
#include "fqdup/errcor_empirical.hpp"
#include "fqdup/bundle_key.hpp"

// fastq_common.hpp provides: FastqRecord, FastqReader, FastqReaderIgzip,
// FastqWriter, SequenceFingerprint, SequenceFingerprintHash, canonical_hash,
// revcomp, trim_id, shell_escape_path, GZBUF_SIZE — all in anonymous namespace.
#include "fqdup/cluster_format.hpp"
#include "fqdup/version.hpp"
#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"
#include "flat_hash_map.hpp"
#include "derep_detail/sharded_index.hpp"

#include "derep_detail/encoding.hpp"
#include "derep_detail/packed_ops.hpp"
#include "derep_detail/arena.hpp"
#include "derep_detail/damage_keys.hpp"
#include "derep_detail/banded_ed.hpp"
#include "derep_detail/syncmer_index.hpp"
#include "derep_detail/b3_deam_merge.hpp"
#include "fqdup/syncmer.hpp"

// All file-local types are in an anonymous namespace to give their member
// functions internal linkage, avoiding ODR violations with other TUs.
namespace cf = fqdup::clusterfmt;

namespace {

using namespace fqdup::derep_detail;

// ============================================================================
// DerepEngine — single-file two-pass deduplication
// ============================================================================

class Phase3Runner;  // re-parents phase3_error_correct's body (phase3_runner.hpp)

class DerepEngine {
    friend class Phase3Runner;
public:
    DerepEngine(bool use_revcomp, bool write_clusters,
                const DamageProfile& profile = DamageProfile{},
                const ErrCorParams& errcor = ErrCorParams{},
                bool /*errcor_mode_explicit*/ = false,
                bool bucket_cap_explicit = false,
                std::string fqcl_path = "",
                std::string input_fastq = "",
                std::string qc_json = "",
                std::string library_type_resolved = "auto",
                std::string prior_fqcl_path = "")
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          profile_(profile), errcor_(errcor),
          bucket_cap_explicit_(bucket_cap_explicit),
          fqcl_path_(std::move(fqcl_path)),
          fqcl_input_fastq_(std::move(input_fastq)),
          qc_json_(std::move(qc_json)),
          library_type_resolved_(std::move(library_type_resolved)),
          prior_fqcl_path_(std::move(prior_fqcl_path)),
          total_reads_(0), errcor_absorbed_(0), n_unique_clusters_(0) {}

    // Where the Phase-3 quality spill lives. $TMPDIR when set — on this cluster it
    // points at node-local scratch, and /tmp does not persist across nodes.
    // Otherwise sit next to the output, which is by definition a writable directory
    // the caller already sized for this run. The file is unlinked as soon as it is
    // created (arena.hpp), so nothing is left behind even if the run dies.
    static std::string spill_dir(const std::string& out_path) {
        if (const char* t = std::getenv("TMPDIR"))
            if (*t) return std::string(t);
        const auto slash = out_path.find_last_of('/');
        return slash == std::string::npos ? std::string(".") : out_path.substr(0, slash);
    }

    void process(const std::string& in_path,
                 const std::string& out_path,
                 const std::string& cluster_path) {
        using clk = std::chrono::steady_clock;
        auto fmt_secs = [](clk::duration d) {
            double s = std::chrono::duration<double>(d).count();
            char b[64];
            std::snprintf(b, sizeof(b), "%.1f s", s);
            return std::string(b);
        };

        log_info("=== Two-pass single-file deduplication ===");
        if (!prior_fqcl_path_.empty()) load_prior_fqcl(prior_fqcl_path_);
        log_info("Pass 1: Build lightweight index");
        // The "Decompression: ..." line moved into pass1(), where the reader actually
        // exists. Here it was derived from #ifdef HAVE_RAPIDGZIP -- whether rapidgzip was
        // COMPILED IN, not whether it was SELECTED -- so it printed "rapidgzip" on every
        // ISA-L run. make_fastq_reader() has defaulted to ISA-L since d1f732d.

        // Phase-3 qualities are spilled, not held: 21.46 GB of Phred bytes on the
        // 354 M-unique clay library existed to serve ~3.7 M single-byte q_at()
        // lookups. Opened only when Phase 3 will actually run, since that is the
        // only consumer (the emit pass re-reads qualities from the input).
        if (errcor_.enabled) qual_arena_.open_spill(spill_dir(out_path));

        auto t_pass1_begin = clk::now();
        pass1(in_path);
        if (errcor_.enabled) qual_arena_.seal();
        auto t_pass1_end = clk::now();
        log_info("Phase timer: Pass 1 = " + fmt_secs(t_pass1_end - t_pass1_begin));
        if (errcor_.enabled)
            log_info("Quality spill: " +
                     std::to_string(qual_arena_.bytes_spilled() / (1024ull * 1024ull)) +
                     " MB on disk, " +
                     std::to_string(qual_arena_.bytes_resident() / (1024ull * 1024ull)) +
                     " MB resident (offsets+lengths)");

        size_t index_mb = (index_.size() *
                           sizeof(std::pair<SequenceFingerprint, IndexEntry>)) /
                          (1024 * 1024);
        log_info("Index size: " + std::to_string(index_mb) + " MB for " +
                 std::to_string(total_reads_) + " reads");

        if (errcor_.enabled) {
            // Pass 1 has fixed the unique count, and Phase 3 has not allocated yet:
            // the only point where the run can refuse a doomed peak cheaply.
            preflight_memory();

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
            auto t_p3_begin = clk::now();
            phase3_error_correct();
            log_info("Phase timer: Phase 3 = " + fmt_secs(clk::now() - t_p3_begin));
        }

        // Free arena_ + qual_arena_ — only needed for phase3, not pass2.
        { SeqArena  empty; std::swap(arena_,      empty); }
        { QualArena empty; std::swap(qual_arena_, empty); }

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

        auto t_pass2_begin = clk::now();
        pass2(in_path, out_path, cluster_path);
        log_info("Phase timer: Pass 2 = " + fmt_secs(clk::now() - t_pass2_begin));

        print_stats();
    }

private:
    // Thread-safe canonical hash: all scratch is passed in, so worker threads
    // may call it concurrently during parallel Pass-1. Byte-identical to the
    // former member body (same LUTs, same masking, same canonical tiebreak).
    static XXH128_hash_t compute_hash_impl(const std::string& seq,
                                           const DamageProfile& profile,
                                           bool use_revcomp, bool& is_forward,
                                           std::vector<char>& s1, std::vector<char>& s2) {
        int L = static_cast<int>(seq.size());
        if (profile.enabled) {
            if (s1.size() < seq.size()) s1.resize(seq.size());
            if (s2.size() < seq.size()) s2.resize(seq.size());
            apply_damage_mask_inplace(seq, profile, s1.data());
            XXH128_hash_t h1 = XXH3_128bits(s1.data(), L);
            if (!use_revcomp) { is_forward = true; return h1; }
            for (int i = 0; i < L; ++i)
                s2[i] = static_cast<char>(
                    lut::kCompCase[static_cast<unsigned char>(seq[L - 1 - i])]);
            for (int i = 0; i < L; ++i) {
                char cu = static_cast<char>(lut::kUpper[static_cast<unsigned char>(s2[i])]);
                bool in_5zone = (i         < DamageProfile::MASK_POSITIONS) && profile.mask_pos[i];
                bool in_3zone = (L - 1 - i < DamageProfile::MASK_POSITIONS) && profile.mask_pos[L - 1 - i];
                if (profile.ss_mode) {
                    if ((in_5zone || in_3zone) && (cu == 'C' || cu == 'T')) s2[i] = '\x01';
                    else if ((in_5zone || in_3zone) && (cu == 'G' || cu == 'A')) s2[i] = '\x02';
                } else {
                    if (in_5zone && (cu == 'C' || cu == 'T')) s2[i] = '\x01';
                    else if (in_3zone && (cu == 'G' || cu == 'A')) s2[i] = '\x02';
                }
            }
            XXH128_hash_t h2 = XXH3_128bits(s2.data(), L);
            is_forward = (h1.high64 < h2.high64 ||
                          (h1.high64 == h2.high64 && h1.low64 <= h2.low64));
            return is_forward ? h1 : h2;
        }
        XXH128_hash_t h1 = XXH3_128bits(seq.data(), seq.size());
        if (!use_revcomp) { is_forward = true; return h1; }
        if (s1.size() < seq.size()) s1.resize(seq.size());
        for (int i = 0; i < L; ++i)
            s1[i] = static_cast<char>(
                lut::kCompCase[static_cast<unsigned char>(seq[L - 1 - i])]);
        XXH128_hash_t h2 = XXH3_128bits(s1.data(), L);
        is_forward = (h1.high64 < h2.high64 ||
                      (h1.high64 == h2.high64 && h1.low64 <= h2.low64));
        return is_forward ? h1 : h2;
    }

    XXH128_hash_t compute_hash(const std::string& seq, bool& is_forward) const {
        return compute_hash_impl(seq, profile_, use_revcomp_, is_forward,
                                 scratch1_, scratch2_);
    }

    // Count terminal damage markers: T at masked 5' positions + A at masked 3' positions.
    // Higher score = more terminal deamination signal = more likely to be authentic ancient DNA.
    // Used to select the most-damaged read as cluster representative, maximising the
    // damage signal in the output.
    // Per-read LLR damaged/non-damaged classifier. Delegates to split_model_ which
    // holds precomputed LOD tables (empirical per-bin or bulk exponential fallback).

    static uint8_t compute_damage_score(const std::string& seq,
                                        const DamageProfile& prof) {
        int score = 0;
        int L = static_cast<int>(seq.size());
        for (int p = 0; p < DamageProfile::MASK_POSITIONS && p < L; ++p) {
            if (!prof.mask_pos[p]) continue;
            if (lut::kUpper[static_cast<unsigned char>(seq[p])]     == 'T') ++score;
            if (lut::kUpper[static_cast<unsigned char>(seq[L-1-p])] == 'A') ++score;
        }
        return static_cast<uint8_t>(std::min(score, 255));
    }

    void load_prior_fqcl(const std::string& path) {
        cf::Reader rdr(path);
        prior_counts_.reserve(rdr.n_clusters());
        cf::ClusterRecord rec;
        for (std::uint64_t i = 0; i < rdr.n_clusters(); ++i) {
            rdr.read_cluster(i, rec);
            if (rec.n_members > 1)
                prior_counts_[rec.cluster_id] = rec.n_members;
        }
        log_info("Prior fqcl loaded: " + std::to_string(rdr.n_clusters()) +
                 " clusters, " + std::to_string(prior_counts_.size()) + " with count > 1");
    }

    // index_ slot = sherwood_v3_entry<pair<SequenceFingerprint,IndexEntry>> = 64 B.
    // At ska's default max_load_factor 0.5 a 355 M-unique library needs 2^30 slots
    // = 68.7 GB, only 33 % occupied — ~46 GB of empty table, the single largest term
    // in the derep peak (measured: 68.7 of 179 GB VmHWM). Robin-hood probing holds up
    // well past 0.5, so raise the ceiling and let 355 M entries live in 2^29 slots
    // (34.4 GB, load 0.66).
    //
    // Must be applied via rehash(), NOT reserve(): flat_hash_map.hpp:805 sizes
    // reserve() with min(0.5, _max_load_factor), so reserve() silently ignores any
    // load factor above 0.5. rehash() (…:626) and the growth check (…:827) both
    // honour it.
    //
    // Byte-identical: index_ is iterated at two sites, and neither lets slot order
    // reach the output. pass2 sorts the iterated result by record_index immediately.
    // Phase 3 (phase3_error_correct.inc, index_.for_each) writes only indexed slots
    // -- id_count[seq_id], id_fwd_count[seq_id], seq_entry[seq_id] -- so each entry
    // lands in its own cell and the visit order is irrelevant. Any layout/load-factor
    // change is therefore output-neutral; adding an order-DEPENDENT iteration (an
    // accumulation, a first-wins, a push_back) would break that and the gold md5s.
    static constexpr float kIndexMaxLoad = 0.75f;

    void pass1(const std::string& in_path) {
        index_.max_load_factor(kIndexMaxLoad);

        // Pre-reserve the hash map based on input file size. Without this,
        // ska::flat_hash_map regrows ~26 times to reach 100M entries, and each
        // regrow is a full re-hash of the existing data. On a 6.3 GB compressed
        // input this single change saves ~30 minutes of wall time.
        //
        // The arenas need no equivalent: they are block-chunked (arena.hpp), so
        // growth appends a block rather than reallocating, and an estimate this
        // rough would only mislead them. It is necessarily rough — the number of
        // *unique* reads depends on the duplication rate, which is not knowable
        // before Pass 1 runs. On the clay libraries it lands 1.7× low.
        //
        // Estimate: ~250 B per uncompressed FASTQ record (typical aDNA SE).
        // Compressed gz is ~5×, so compressed_bytes / 50 ≈ records.
        // Use 0.6× safety factor.
        try {
            std::error_code ec;
            auto fsize = std::filesystem::file_size(in_path, ec);
            if (!ec && fsize > 0) {
                bool gz = in_path.size() >= 3 &&
                          in_path.compare(in_path.size() - 3, 3, ".gz") == 0;
                size_t est_records = gz ? static_cast<size_t>(fsize / 50)
                                        : static_cast<size_t>(fsize / 250);
                est_records = static_cast<size_t>(est_records * 0.6);
                est_records = std::min(est_records, size_t{500'000'000});
                if (est_records >= 100'000) {
                    log_info("Pass 1: pre-reserving index capacity for ~" +
                             std::to_string(est_records) +
                             " entries (heuristic from input size " +
                             std::to_string(fsize / (1024 * 1024)) + " MB)");
                    index_.rehash(static_cast<size_t>(
                        std::ceil(est_records / static_cast<double>(kIndexMaxLoad))));
                }
            }
        } catch (...) {
            // reserve is best-effort; any failure just means we fall back to
            // the regrow path.
        }

        // Parallel Pass-1 offloads the pure per-read work (canonical hash, RC
        // build, damage score) to a worker pool; a single committer replays the
        // exact serial insert in record order -> byte-identical index_/arena_.
        // Disabled when a derep_pairs prior is loaded (that path needs raw seq
        // in the committer) or when only one thread is available.
        unsigned n_threads = errcor_.threads;
        if (n_threads == 0) n_threads = std::max(1u, std::thread::hardware_concurrency());
        bool want_par = n_threads > 1 && prior_counts_.empty();
        if (const char* e = std::getenv("FQDUP_PASS1_SERIAL"))
            if (e[0] == '1') want_par = false;

        if (want_par) pass1_parallel(in_path, n_threads);
        else          pass1_serial(in_path);

        if (n_below_min_length_ > 0)
            log_info("--min-length " + std::to_string(min_length_) + ": dropped " +
                     std::to_string(n_below_min_length_) + " reads shorter than " +
                     std::to_string(min_length_) + " bp");
        if (n_above_max_length_ > 0)
            log_info("--max-length " + std::to_string(max_length_) + ": dropped " +
                     std::to_string(n_above_max_length_) + " reads longer than " +
                     std::to_string(max_length_) + " bp");
    }

    void pass1_serial(const std::string& in_path) {
        auto reader = make_fastq_reader(in_path);
        log_info("Decompression: " + std::string(reader->backend_name()));
        FastqRecord rec;
        uint64_t record_idx = 0;

        while (reader->read(rec)) {
            // Length gate: advance the physical index (emit pass re-reads by position)
            // but skip indexing so dropped reads never become a cluster representative.
            const int Lrec = (int)rec.seq.size();
            if (Lrec < min_length_) { ++n_below_min_length_; record_idx++; continue; }
            if (max_length_ > 0 && Lrec > max_length_) { ++n_above_max_length_; record_idx++; continue; }
            bool is_forward = true;
            XXH128_hash_t h = compute_hash(rec.seq, is_forward);
            SequenceFingerprint fp(h, rec.seq.size());

            auto [it, inserted] = index_.emplace(fp, IndexEntry(record_idx));
            if (inserted) {
                // Seed count from derep_pairs prior when available.
                if (!prior_counts_.empty()) {
                    uint64_t h = XXH3_64bits(rec.seq.data(), rec.seq.size());
                    auto pit = prior_counts_.find(h);
                    if (pit != prior_counts_.end())
                        it->second.count = IndexEntry::narrow(pit->second);
                }
                if (profile_.enabled)
                    it->second.damage_score = compute_damage_score(rec.seq, profile_);
                if (errcor_.enabled) {
                    if (is_forward) {
                        it->second.seq_id = arena_.append(rec.seq);
                        if (!rec.qual.empty()) qual_arena_.append(rec.qual);
                        else                   qual_arena_.append_empty();
                    } else {
                        // Store in canonical orientation so Phase 3 comparisons are correct.
                        // Sequence: reverse-complement; quality: reverse only.
                        int L = static_cast<int>(rec.seq.size());
                        if (rc_buf_.size()  < static_cast<size_t>(L)) rc_buf_.resize(L);
                        if (rc_qbuf_.size() < static_cast<size_t>(L)) rc_qbuf_.resize(L);
                        for (int i = 0; i < L; ++i)
                            rc_buf_[i] = static_cast<char>(
                                lut::kCompUpper[static_cast<unsigned char>(rec.seq[L - 1 - i])]);
                        it->second.seq_id = arena_.append_chars(rc_buf_.data(), L);
                        if (!rec.qual.empty()) {
                            for (int i = 0; i < L; ++i) rc_qbuf_[i] = rec.qual[L - 1 - i];
                            qual_arena_.append_chars(rc_qbuf_.data(), L);
                        } else {
                            qual_arena_.append_empty();
                        }
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
                        it->second.record_index = IndexEntry::narrow(record_idx);
                        it->second.damage_score = score;
                    }
                }
            }
            if (is_forward) it->second.fwd_count++;

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

    // Byte-identical parallel variant of pass1_serial (see want_par guard).
    void pass1_parallel(const std::string& in_path, unsigned n_threads) {
        const size_t CHUNK = size_t{1} << 22;   // 4M records per chunk

        std::vector<char>          seqbuf, qualbuf, canon, canonq;
        std::vector<uint64_t>      off;          // n+1 prefix offsets (shared: seq==qual len)
        std::vector<uint16_t>      len;          // n
        std::vector<uint8_t>       hasq;         // n
        std::vector<XXH128_hash_t> hh;           // n
        std::vector<uint8_t>       fwd, dsc;     // n

        auto reader = make_fastq_reader(in_path);
        log_info("Decompression: " + std::string(reader->backend_name()));
        FastqRecord rec;
        uint64_t record_idx = 0;

        auto worker = [&](size_t lo, size_t hi) {
            std::vector<char> s1, s2;            // thread-local hash scratch
            for (size_t i = lo; i < hi; ++i) {
                int L = len[i];
                const char* sp = seqbuf.data() + off[i];
                std::string seq(sp, static_cast<size_t>(L));
                bool f = true;
                hh[i]  = compute_hash_impl(seq, profile_, use_revcomp_, f, s1, s2);
                fwd[i] = f ? 1 : 0;
                if (profile_.enabled) dsc[i] = compute_damage_score(seq, profile_);
                char* cs = canon.data() + off[i];
                if (f) std::memcpy(cs, sp, L);
                else   for (int k = 0; k < L; ++k)
                           cs[k] = static_cast<char>(
                               lut::kCompUpper[static_cast<unsigned char>(sp[L - 1 - k])]);
                if (hasq[i]) {
                    const char* qp = qualbuf.data() + off[i];
                    char* cq = canonq.data() + off[i];
                    if (f) std::memcpy(cq, qp, L);
                    else   for (int k = 0; k < L; ++k) cq[k] = qp[L - 1 - k];
                }
            }
        };

        while (true) {
            seqbuf.clear(); qualbuf.clear(); off.clear(); len.clear(); hasq.clear();
            off.push_back(0);
            size_t n = 0;
            while (n < CHUNK && reader->read(rec)) {
                uint16_t L = static_cast<uint16_t>(rec.seq.size());
                seqbuf.insert(seqbuf.end(), rec.seq.begin(), rec.seq.end());
                bool hq = !rec.qual.empty();
                if (hq) qualbuf.insert(qualbuf.end(), rec.qual.begin(), rec.qual.end());
                else    qualbuf.insert(qualbuf.end(), rec.seq.size(), '\0');
                len.push_back(L);
                hasq.push_back(hq ? 1 : 0);
                off.push_back(seqbuf.size());
                ++n;
            }
            if (n == 0) break;

            hh.resize(n); fwd.resize(n); dsc.assign(n, 0);
            canon.resize(seqbuf.size()); canonq.resize(qualbuf.size());

            unsigned nw = static_cast<unsigned>(std::min<size_t>(n_threads, n));
            std::vector<std::thread> ts; ts.reserve(nw);
            size_t per = (n + nw - 1) / nw;
            for (unsigned t = 0; t < nw; ++t) {
                size_t lo = static_cast<size_t>(t) * per, hi = std::min(n, lo + per);
                if (lo >= hi) break;
                ts.emplace_back(worker, lo, hi);
            }
            for (auto& th : ts) th.join();

            for (size_t i = 0; i < n; ++i) {
                int L = len[i];
                // Length gate (mirrors pass1_serial): advance physical index, skip indexing.
                // This merge loop runs single-threaded after join() -> counters need no atomic.
                if (L < min_length_) { ++n_below_min_length_; ++record_idx; continue; }
                if (max_length_ > 0 && L > max_length_) { ++n_above_max_length_; ++record_idx; continue; }
                SequenceFingerprint fp(hh[i], static_cast<size_t>(L));
                auto [it, inserted] = index_.emplace(fp, IndexEntry(record_idx));
                if (inserted) {
                    if (profile_.enabled) it->second.damage_score = dsc[i];
                    if (errcor_.enabled) {
                        it->second.seq_id = arena_.append_chars(canon.data() + off[i], L);
                        if (hasq[i]) qual_arena_.append_chars(canonq.data() + off[i], L);
                        else         qual_arena_.append_empty();
                    }
                } else {
                    it->second.count++;
                    if (profile_.enabled && dsc[i] > it->second.damage_score) {
                        it->second.record_index = IndexEntry::narrow(record_idx);
                        it->second.damage_score = dsc[i];
                    }
                }
                if (fwd[i]) it->second.fwd_count++;
                ++record_idx;
                ++total_reads_;
                if ((total_reads_ % 1000000) == 0) {
                    size_t uniq = index_.size();
                    double dup_pct = 100.0 * (1.0 - (double)uniq / total_reads_);
                    std::cerr << "\r[Pass 1] " << total_reads_ << " reads, "
                              << uniq << " unique, " << std::fixed << std::setprecision(1)
                              << dup_pct << "% dedup" << std::flush;
                }
            }
        }

        std::cerr << "\r";
        log_info("Pass 1 complete: " + std::to_string(total_reads_) +
                 " reads indexed (parallel x" + std::to_string(n_threads) + ")");
    }

    // Phase 3: parent-centric mismatch pattern detection.
    // Indexes sequences with count > min_parent_count as parents, then for each
    // potential child (count <= min_parent_count) finds all H=1 parent neighbours.
    // A child is absorbed unless the mismatch pattern is recurrent (SNP veto):
    //   sig_count_weighted >= snp_min_count AND
    //   sig_count_weighted / parent_count >= snp_threshold.
    void phase3_error_correct();  // out-of-line via Phase3Runner (phase3_runner.hpp)

    // Build per-cluster genealogy and emit .fqcl. Must run BEFORE arena_ is freed.
    void write_fqcl_(const std::vector<IndexEntry*>& seq_entry) {
        namespace cf = fqdup::clusterfmt;

        const uint32_t N = static_cast<uint32_t>(arena_.size());
        auto find_root = [&](uint32_t id) -> uint32_t {
            while (id < fqcl_parent_chain_.size() &&
                   fqcl_parent_chain_[id] != UINT32_MAX) {
                id = fqcl_parent_chain_[id];
            }
            return id;
        };

        // Group absorbed-edge mismatches by root cluster id, preserving order
        // (sorted by parent count desc → topological).
        std::unordered_map<uint32_t, std::vector<size_t>> per_root;
        per_root.reserve(fqcl_mismatches_.size() / 2 + 1);
        for (size_t i = 0; i < fqcl_mismatches_.size(); ++i) {
            uint32_t r = find_root(fqcl_mismatches_[i].child_id);
            per_root[r].push_back(i);
        }

        // Build cluster_id mapping over surviving roots: dense u64 in seq_id order.
        cf::WriterMetadata meta;
        meta.tool_version  = FQDUP_VERSION;
        meta.input_fastq   = fqcl_input_fastq_;
        // total_reads_ counts only indexed survivors (gates continue before the
        // increment); add length-filtered drops to report true physical input.
        meta.n_input_reads = total_reads_ + n_below_min_length_ + n_above_max_length_;
        // Bug 1 fix: report the resolved library type from the user flag
        // (or "auto" when truly auto-detected). Previously this read
        // profile_.enabled, which is false whenever deamination fit didn'''t converge (small
        // input, --no-damage-auto, manual dmax=0…), wiping the user's
        // explicit --library-type ds/ss request from the header.
        if (!library_type_resolved_.empty()) {
            meta.library_type = library_type_resolved_;
        } else {
            meta.library_type = profile_.ss_mode ? "ss"
                : (profile_.enabled ? "ds" : "auto");
        }
        meta.d_max_5       = profile_.d_max_5prime;
        meta.d_max_3       = profile_.d_max_3prime;
        meta.lambda_5      = profile_.lambda_5prime;
        meta.lambda_3      = profile_.lambda_3prime;
        // Count mask positions actually used.
        int mp5 = 0, mp3 = 0;
        for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p)
            if (profile_.mask_pos[p]) { mp5++; mp3++; }
        meta.mask_pos_5 = mp5;
        meta.mask_pos_3 = mp3;
        meta.snp_threshold = errcor_.snp_threshold;
        meta.snp_min_count = static_cast<int>(errcor_.snp_min_count);
        meta.bucket_cap    = static_cast<int>(errcor_.bucket_cap);
        meta.qc_json                 = qc_json_;
        meta.bucket_overflow_drops   = loss_bucket_overflow_drops_;
        meta.short_interior_skipped  = loss_short_interior_skipped_;
        meta.short_brute_evaluated   = loss_short_brute_evaluated_;
        meta.short_brute_found       = loss_short_brute_found_;
        meta.short_too_small_skipped = loss_short_too_small_skipped_;
        // v2: singletons are always written (as tinyblocks) — never silently dropped.

        cf::Writer writer(fqcl_path_, meta);
        writer.reserve_clusters(seq_entry.size());

        std::uint64_t cluster_id = 0;
        cf::ClusterRecord rec;
        for (uint32_t r = 0; r < N; ++r) {
            if (is_error_[r]) continue;             // absorbed root
            if (!seq_entry[r]) continue;            // no surviving entry pointer
            const IndexEntry* e = seq_entry[r];

            // v2: emit every surviving root, including pure singletons.
            // Pure singletons travel as tinyblocks (zero edges, zero damage).
            const uint16_t L      = arena_.length(r);
            const int      nbytes = (L + 3) / 4;

            rec.cluster_id     = cluster_id++;
            rec.flags          = cf::kFlagHasMemberIds;
            rec.n_members      = static_cast<uint32_t>(e->count);
            rec.n_after_damage = static_cast<uint32_t>(e->count);
            rec.parent_seq_len = L;
            rec.parent_seq.assign(arena_.data(r), arena_.data(r) + nbytes);
            rec.parent_qual.clear();
            rec.edges.clear();
            rec.member_ids.clear();
            // Damage track left zero in v0.1; per-cluster terminal counts not yet tracked.
            std::memset(rec.damage.term_5, 0, cf::kTermWindow);
            std::memset(rec.damage.term_3, 0, cf::kTermWindow);

            // Always record the root's seq_id as the first member (node 0 = root).
            rec.member_ids.push_back(std::to_string(r));

            auto it = per_root.find(r);
            if (it != per_root.end()) {
                // Mismatch positions are interior-relative; lift to full-read coords.
                auto [k5_off, k3_off] = damage_zone_bounds(L, profile_);
                (void)k3_off;
                const uint16_t k5u = static_cast<uint16_t>(k5_off);

                // Map seq_id → local node index. Node 0 = root.
                std::unordered_map<uint32_t, uint32_t> node_of;
                node_of.reserve(it->second.size() + 1);
                node_of[r] = 0;

                for (size_t midx : it->second) {
                    const ChildMismatch& cm = fqcl_mismatches_[midx];
                    auto pit = node_of.find(cm.parent_id);
                    if (pit == node_of.end()) continue;  // orphan (shouldn't happen)
                    uint32_t parent_node = pit->second;
                    uint32_t child_node  = static_cast<uint32_t>(node_of.size());
                    node_of[cm.child_id] = child_node;
                    // Record the absorbed child's seq_id so callers can enumerate
                    // every read that ended up under this root (singletons remain
                    // implicit: any seq_id never listed in any cluster's member_ids
                    // is a true unabsorbed singleton).
                    rec.member_ids.push_back(std::to_string(cm.child_id));

                    cf::Edge ed{};
                    ed.from_node = parent_node;
                    ed.to_node   = child_node;
                    ed.pos       = static_cast<uint16_t>(cm.mismatch_pos + k5u);
                    ed.from_base = cm.parent_base;
                    ed.to_base   = cm.alt_base;
                    ed.n_reads   = static_cast<uint32_t>(
                        cm.child_id < fqcl_parent_chain_.size() ? 1u : 1u);
                    ed.score     = cm.lr_score;
                    rec.edges.push_back(ed);

                    // H=2: emit a second edit hanging off the same child node.
                    if (cm.hamming == 2) {
                        uint32_t aux_node = static_cast<uint32_t>(node_of.size());
                        cf::Edge ed2{};
                        ed2.from_node = child_node;
                        ed2.to_node   = aux_node;
                        ed2.pos       = static_cast<uint16_t>(cm.mismatch_pos2 + k5u);
                        ed2.from_base = cm.parent_base2;
                        ed2.to_base   = cm.alt_base2;
                        ed2.n_reads   = 1u;
                        rec.edges.push_back(ed2);
                        // aux_node has no real seq_id — leave it out of node_of.
                    }
                }
            }

            writer.write_cluster(rec);
        }

        writer.close();
        log_info("Cluster-format file written: " + fqcl_path_ +
                 " (" + std::to_string(writer.n_clusters()) + " clusters; "
                 "singletons stored as tinyblocks)");

        // Free temporaries.
        std::vector<ChildMismatch>().swap(fqcl_mismatches_);
        std::vector<uint32_t>().swap(fqcl_parent_chain_);
    }

    void pass2(const std::string& in_path,
               const std::string& out_path,
               const std::string& cluster_path) {
        bool compress = (out_path.size() > 3 &&
                         out_path.substr(out_path.size() - 3) == ".gz");

        // Multi-thread compressed output via BGZF (htslib). Cap at 16: BGZF
        // throughput plateaus there and extra threads just steal from the read
        // loop. Without HAVE_BGZF or with threads<=1, FastqWriter falls back to
        // single-thread zlib.
        unsigned writer_threads = errcor_.threads;
        if (writer_threads == 0)
            writer_threads = std::max(1u, std::thread::hardware_concurrency());
        writer_threads = std::min(writer_threads, 16u);

        std::unique_ptr<FastqWriter> writer_dam, writer_und;
        if (!out_damaged_path_.empty()) {
            bool c = out_damaged_path_.size() > 3 &&
                     out_damaged_path_.substr(out_damaged_path_.size()-3) == ".gz";
            writer_dam = std::make_unique<FastqWriter>(
                out_damaged_path_, c, static_cast<int>(writer_threads));
        }
        if (!out_undamaged_path_.empty()) {
            bool c = out_undamaged_path_.size() > 3 &&
                     out_undamaged_path_.substr(out_undamaged_path_.size()-3) == ".gz";
            writer_und = std::make_unique<FastqWriter>(
                out_undamaged_path_, c, static_cast<int>(writer_threads));
        }
        bool do_split = writer_dam || writer_und;

        FastqWriter writer(out_path, compress, static_cast<int>(writer_threads));

        gzFile cluster_gz = nullptr;
#ifdef HAVE_BGZF
        BGZF*  cluster_bgz = nullptr;
#endif
        if (write_clusters_ && !cluster_path.empty()) {
#ifdef HAVE_BGZF
            if (writer_threads > 1) {
                cluster_bgz = bgzf_open(cluster_path.c_str(), "w");
                if (!cluster_bgz)
                    throw std::runtime_error("Cannot open cluster file: " + cluster_path);
                bgzf_mt(cluster_bgz, static_cast<int>(writer_threads), 0);
                static const char* hdr = "hash\tseq_len\tcount\n";
                if (bgzf_write(cluster_bgz, hdr, std::strlen(hdr)) < 0)
                    throw std::runtime_error("bgzf_write failed writing cluster header");
            } else
#endif
            {
                cluster_gz = gzopen(cluster_path.c_str(), "wb6");
                if (!cluster_gz)
                    throw std::runtime_error("Cannot open cluster file: " + cluster_path);
                gzbuffer(cluster_gz, GZBUF_SIZE);
                if (gzprintf(cluster_gz, "hash\tseq_len\tcount\n") < 0)
                    throw std::runtime_error("gzprintf failed while writing cluster header");
            }
        }

        struct WriteEntry {
            uint64_t            record_index;
            SequenceFingerprint fingerprint;
            uint64_t            count;
        };
        std::vector<WriteEntry> records_to_write;
        records_to_write.reserve(index_.size());
        index_.for_each([&](const SequenceFingerprint& fingerprint,
                            const IndexEntry& entry) {
            if (!is_error_.empty() && is_error_[entry.seq_id]) return;
            records_to_write.push_back({entry.record_index, fingerprint, entry.count});
        });
        std::sort(records_to_write.begin(), records_to_write.end(),
                  [](const WriteEntry& a, const WriteEntry& b) {
                      return a.record_index < b.record_index;
                  });
        const size_t n_to_write = records_to_write.size();
        n_unique_clusters_ = n_to_write;

        // Free index_ and is_error_ — no longer needed; records_to_write has everything.
        // index_ is the largest single structure in the run (32 GB measured on a
        // 355 M-unique library at kIndexMaxLoad), so this is the biggest drop before
        // the write loop begins.
        index_.release();
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
                if (do_split) {
                    float llr = split_model_.score(rec.seq);
                    if (llr >= split_threshold_) {
                        if (writer_dam) writer_dam->write(rec);
                    } else {
                        if (writer_und) writer_und->write(rec);
                    }
                }

                if (cluster_gz) {
                    if (gzprintf(cluster_gz, "%016lx%016lx\t%lu\t%lu\n",
                                 static_cast<unsigned long>(entry.fingerprint.hash_hi),
                                 static_cast<unsigned long>(entry.fingerprint.hash_lo),
                                 static_cast<unsigned long>(rec.seq.size()),
                                 static_cast<unsigned long>(entry.count)) < 0)
                        throw std::runtime_error(
                            "gzprintf failed while writing cluster record");
                }
#ifdef HAVE_BGZF
                else if (cluster_bgz) {
                    char buf[80];
                    int n = std::snprintf(buf, sizeof(buf), "%016lx%016lx\t%lu\t%lu\n",
                                          static_cast<unsigned long>(entry.fingerprint.hash_hi),
                                          static_cast<unsigned long>(entry.fingerprint.hash_lo),
                                          static_cast<unsigned long>(rec.seq.size()),
                                          static_cast<unsigned long>(entry.count));
                    if (n <= 0 || bgzf_write(cluster_bgz, buf, static_cast<size_t>(n)) < 0)
                        throw std::runtime_error("bgzf_write failed writing cluster record");
                }
#endif

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
#ifdef HAVE_BGZF
        if (cluster_bgz && bgzf_close(cluster_bgz) < 0)
            throw std::runtime_error("bgzf_close failed writing cluster file");
#endif
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
    std::string fqcl_path_;
    std::string fqcl_input_fastq_;
    std::string qc_json_;
    std::string library_type_resolved_;
    std::vector<ChildMismatch> fqcl_mismatches_;
    std::vector<uint32_t>      fqcl_parent_chain_;

    fqdup::derep_detail::ShardedIndex<SequenceFingerprint, IndexEntry,
                                      SequenceFingerprintHash> index_;

    // Phase 3 error correction state (populated only when errcor_.enabled)
    SeqArena          arena_;
    QualArena         qual_arena_;
    std::vector<bool> is_error_;

    // Reusable scratch buffers for allocation-free damage masking/hashing.
    // Grown lazily to max sequence length seen; never shrunk.
    mutable std::vector<char> scratch1_;
    mutable std::vector<char> scratch2_;
    mutable std::vector<char> rc_buf_;   // unmasked revcomp for canonical arena storage
    mutable std::vector<char> rc_qbuf_;  // reversed qual paired with rc_buf_

    uint64_t total_reads_;
    uint64_t errcor_absorbed_;
    size_t   n_unique_clusters_;  // saved before index_ is freed in pass2

    // Prior counts from a preceding derep_pairs --cluster-format run.
    // Maps XXH3_64(raw_seq) → n_members so Phase 1 seeds the correct count.
    std::string prior_fqcl_path_;
    ska::flat_hash_map<uint64_t, uint32_t> prior_counts_;

    // Optional split outputs: damaged / undamaged reads (LLR threshold = 0).
    std::string out_damaged_path_;
    std::string out_undamaged_path_;
    float split_threshold_ = 0.0f;
    DamageSplitModel split_model_;
    int min_length_ = 0;   // 0=disabled; CLI default 16. Reads shorter dropped before indexing.
    int max_length_ = 0;   // 0=disabled (no upper cap); reads longer dropped before indexing
    uint64_t n_below_min_length_ = 0;  // Pass-1 --min-length drops (never indexed/emitted)
    uint64_t n_above_max_length_ = 0;  // Pass-1 --max-length drops (never indexed/emitted)

    double max_mem_gb_ = 0.0;  // 0 = disabled

    // Preflight the run's peak once Pass 1 has fixed the unique count, and abort
    // there rather than letting the run discover the ceiling by being OOM-killed
    // several hours in. This is the only point where the count is known and Phase 3
    // has not allocated yet.
    //
    //   base = max( VmHWM already reached, VmRSS now + Phase-3 growth )
    //
    // Both halves are needed. VmHWM covers a peak the process has ALREADY survived
    // in Pass 1 (reader buffers, index growth) -- on a 6 M-read library that peak is
    // 2.27 GB while RSS at this point is only 1.14 GB, so an RSS-based estimate
    // would promise a run fits in less memory than it has already used. The RSS term
    // covers the opposite regime: at 355 M uniques Phase 3 is the peak and Pass 1
    // never came near it.
    //
    // Phase-3 growth counts only the allocations that are pure functions of N. The
    // rest -- the largest ilen's shard and postings, the mismatch records, the
    // bundle-occupancy maps, the fit accumulators -- is B1-derived and, with
    // bucket_cap=0, formally unbounded. kPeakFactor covers those plus allocator
    // slack: it is the ratio of measured VmHWM to this base on the reference runs
    // (355 M uniques: 83.68 / 69.4 = 1.21; 6 M: 2.27 / 2.27 = 1.00).
    //
    // ceiling: kPeakFactor is calibrated on two libraries. A library with an
    // unusually dominant single ilen, or --bucket-cap 0 on low-complexity data,
    // can exceed it -- this is an estimate, not a guarantee.
    // upgrade: log measured VmHWM / base from each run and widen to the observed max.
    void preflight_memory() const {
        if (max_mem_gb_ <= 0.0) return;
        constexpr double kPeakFactor = 1.25;
        constexpr double kGiB        = 1073741824.0;

        unsigned long rss_kb = 0, hwm_kb = 0;
        if (FILE* f = std::fopen("/proc/self/status", "r")) {
            char line[256];
            while (std::fgets(line, sizeof(line), f)) {
                std::sscanf(line, "VmRSS: %lu kB", &rss_kb);
                std::sscanf(line, "VmHWM: %lu kB", &hwm_kb);
            }
            std::fclose(f);
        }
        const uint64_t N = index_.size();
        // Phase-3 per-unique allocations, all sized N (phase3_error_correct.inc):
        //   id_count(4) + id_fwd_count(4) + seq_entry(8, IndexEntry*) + acc_count(4)
        //   + bundle_key_of(8) + rank_of(4)                                    = 32
        //   + rescue_bundle_key_of(8) when the rescue-indel path allocates it.
        // id_fwd_count and seq_entry were omitted, so this understated Phase 3 by
        // 12 B/unique -- ~7.4 GB on a 615 M-unique library. --max-mem then blessed a
        // run that OOMed anyway, which is the one thing this guard exists to prevent.
        const bool rescue = errcor_.rescue_indels && !errcor_.legacy_veto &&
                            errcor_.empirical;
        const uint64_t per_read = 32ull + (rescue ? 8ull : 0ull);
        const double rss_gb = static_cast<double>(rss_kb) * 1024.0 / kGiB;
        const double hwm_gb = static_cast<double>(hwm_kb) * 1024.0 / kGiB;
        const double p3_gb  = rss_gb + static_cast<double>(N * per_read) / kGiB;
        const double base_gb = std::max(hwm_gb, p3_gb);
        const double est_gb  = base_gb * kPeakFactor;

        char msg[256];
        std::snprintf(msg, sizeof(msg),
                      "--max-mem %.1f GB: estimated peak=%.2f GB "
                      "= max(VmHWM %.2f, resident %.2f + %llu B/read x %llu uniques "
                      "= %.2f) x %.2f",
                      max_mem_gb_, est_gb, hwm_gb, rss_gb,
                      static_cast<unsigned long long>(per_read),
                      static_cast<unsigned long long>(N), p3_gb, kPeakFactor);
        log_info(msg);
        if (est_gb > max_mem_gb_)
            throw std::runtime_error(
                "estimated peak " + std::to_string(est_gb).substr(0, 6) +
                " GB exceeds --max-mem " + std::to_string(max_mem_gb_).substr(0, 6) +
                " GB. Raise --max-mem, or lower the peak with --no-error-correct "
                "(skips Phase 3 entirely).");
    }

public:
    void set_max_mem(double gb) { max_mem_gb_ = gb; }
    void set_length_filter(int mn, int mx) { min_length_ = mn; max_length_ = mx; }
    void set_split_paths(std::string damaged, std::string undamaged,
                         float threshold = 0.0f) {
        out_damaged_path_   = std::move(damaged);
        out_undamaged_path_ = std::move(undamaged);
        split_threshold_    = threshold;
    }
    void set_split_model(DamageSplitModel model) {
        split_model_ = std::move(model);
    }
private:

    // Loss counters surfaced into .fqcl metadata (T1.3, T1.2, T3.1).
    uint64_t loss_bucket_overflow_drops_   = 0;
    uint64_t loss_short_interior_skipped_  = 0;
    uint64_t loss_short_brute_evaluated_   = 0;
    uint64_t loss_short_brute_found_       = 0;
    uint64_t loss_short_too_small_skipped_ = 0;
};

// Phase3Runner needs the complete DerepEngine type; include after the class,
// then define the thin wrapper. Body remains byte-identical (baseline af3ac4e1).
#include "derep_detail/phase3_runner.hpp"
void DerepEngine::phase3_error_correct() { Phase3Runner(*this).run(); }

}  // anonymous namespace

// ============================================================================
// Public entry point
// ============================================================================

static void print_usage(const char* prog, bool advanced = false) {
    std::cerr
        << "Usage: fqdup " << prog << " [OPTIONS]\n"
        << "\nSingle-file FASTQ deduplication with damage-aware hashing\n"
        << "and PCR error correction. Designed for sorted single-end FASTQ\n"
        << "(e.g. the non-extended output of 'fqdup derep_pairs').\n"
        << "\nInput MUST be sorted by read ID (use 'fqdup sort' first).\n"
        << "\nRequired:\n"
        << "  -i FILE              Input sorted FASTQ\n"
        << "  -o FILE              Output deduplicated FASTQ\n"
        << "\nOutput:\n"
        << "  -c FILE              Cluster statistics (gzipped TSV)\n"
        << "  --out-damaged FILE   Write LLR-classified damaged reads to FILE (.fq.gz)\n"
        << "  --out-undamaged FILE Write LLR-classified undamaged reads to FILE (.fq.gz)\n"
        << "  --split-threshold F  LLR threshold for split (default 0.0)\n"
        << "  --min-length N       Drop reads shorter than N bp before dedup (default: 16; 0 = off)\n"
        << "  --max-length N       Drop reads longer than N bp before dedup (0 = off; must be >= --min-length)\n"
        << "  --max-mem GB         Abort after Pass 1 if the estimated Phase-3 peak exceeds GB (0 = off)\n"
        << "  --split-model MODE   auto (default) | bulk | empirical\n"
        << "                         auto     — empirical per-bin if damage detected, else bulk\n"
        << "                         bulk     — bulk exponential only, no extra file pass\n"
        << "                         empirical — force per-bin scan (most accurate)\n"
        << "                       -o is optional when --out-damaged/--out-undamaged given\n"
        << "  --cluster-format F   Write .fqcl genealogy (requires error correction)\n"
        << "  --prior-fqcl F       Load cluster counts from a derep_pairs --cluster-format output;\n"
        << "                       seeds Phase 3 count weights for correct PCR-error scoring\n"
        << "\nCore:\n"
        << "  --no-revcomp         Disable reverse-complement matching\n"
        << "  --no-error-correct   Disable PCR error correction (Phase 3, ON by default)\n"
        << "  --errcor-rescue-indels  Enable syncmer-indexed indel rescue (ed<=2; default OFF)\n"
        << "  -h, --help           Show this help\n"
        << "  --help-advanced      Show advanced options (errcor/damage/PCR knobs)\n"
        << "\nDamage:\n"
        << "  --damage MODE        off | report | collapse  (default: report)\n"
        << "                         off      — no fit, no QC\n"
        << "                         report   — fit + QC, header only (no clustering use)\n"
        << "                         collapse — fit + QC + use damage when clustering\n"
        << "                                    WARNING: distorts per-position damage rates;\n"
        << "                                    do NOT use upstream of metaDMG/mapDamage.\n"
        << "  --library-type T     auto (default) | ds | ss | unknown\n"
        << "  --damage-dmax N[,N]  Manual d_max override (skips fit). \"0.21,0.13\" = 5',3'\n"
        << "  --damage-clip-pass M off | report (default) | refit\n"
        << "  --damage-scan N      Reads sampled for QC + adapter scan (default: 1000000; 0=all)\n"
        << "  --damage-deam-sample N  Reads sampled for deamination d_max/lambda fit (default: 5000000; 0=all)\n";

    if (!advanced) {
        std::cerr << "\nMemory: ~16 bytes per input read + 2-bit seq arena (~1 byte/4 bp) for error correction.\n";
        return;
    }

    std::cerr
        << "\n--- Advanced ---\n"
        << "\nError correction (Phase 3) tuning:\n"
        << "  --no-b1-damage-adjust        Widen the H=2 damage veto to ANY transition anywhere in the\n"
        << "                               read (legacy). Default vetoes only terminal-zone,\n"
        << "                               ss/ds-aware deaminations; this also blocks interior SNPs.\n"
        << "  --b1-kappa-min FLOAT         Interior CT/TV ratio threshold to enable b1_damage_adjust (default: 2.0)\n"
        << "  --b1-cs-min-total INT        Min child reads to apply cross-strand SNP veto (default: 3)\n"
        << "  --errcor-snp-threshold FLOAT SNP veto: sig/parent_count threshold (default: 0.20)\n"
        << "  --errcor-snp-min-count INT   SNP veto: min absolute sig_count (default: 1)\n"
        << "  --errcor-bucket-cap INT      Max pair-key bucket size (default: 0 = unlimited)\n"
        << "  --errcor-bucket-cap-lowcomplexity  Cap only low-complexity buckets; spare high-complexity (real abundance)\n"
        << "  --errcor-bucket-cap-complexity-frac FLOAT  Distinct-4mer frac below which an interior is low-complexity (default: 0.5)\n"
        << "  --errcor-snp-cutoff INT      Low-coverage cutoff (default: 10)\n"
        << "  --errcor-snp-factor FLOAT    Low-coverage SNP multiplier (default: 1.75)\n"
        << "  -p, --threads INT            Worker threads for Phase 3 + compressed I/O (default: 0 = HW, capped at 16 for writer)\n"
        << "  --errcor-empirical           Empirical posterior-odds rule (default; absorb iff S>0)\n"
        << "  --errcor-legacy-veto         Legacy SNP-veto rule\n"
        << "  --protect-transversions      Protect A↔T / C↔G (Channels H/G) from H=1/H=2 absorption\n"
        << "  --errcor-adj-len             Adjacent-length (L±1) indel probe\n"
        << "\nIndel rescue (T8) tuning:\n"
        << "  --rescue-min-hits INT        Min shared syncmers (default: auto: 3 cap=16, 2 cap=8)\n"
        << "  --rescue-topk INT            Per-child max parent candidates (default: 32)\n"
        << "  --rescue-bundle-hot INT      Skip bundles with occupancy >= N (default: 50)\n"
        << "  --rescue-hash-hot INT        Drop syncmer postings with list size >= N (default: 4096)\n"
        << "  --rescue-alpha-ins FLOAT     Per-insertion log-odds penalty (default: 2.0)\n"
        << "  --rescue-alpha-del FLOAT     Per-deletion log-odds penalty (default: 2.0)\n"
        << "  --rescue-mask-bonus FLOAT    Bonus when indel falls in damage mask (default: 1.0)\n"
        << "\nDamage expert overrides (rarely needed):\n"
        << "  --damage-dmax5 FLOAT         5'-only d_max (alternative to comma syntax)\n"
        << "  --damage-dmax3 FLOAT         3'-only d_max\n"
        << "  --damage-lambda FLOAT        Lambda (decay) for both ends\n"
        << "  --damage-lambda5 FLOAT       Lambda for 5' end\n"
        << "  --damage-lambda3 FLOAT       Lambda for 3' end\n"
        << "  --damage-bg FLOAT            Background deamination rate (default: 0.02)\n"
        << "  --mask-threshold FLOAT       Mask positions with excess damage > T (default: 0.05)\n"
        << "\nDeprecated aliases (still parsed):\n"
        << "  --collapse-damage            → --damage collapse\n"
        << "  --no-damage-qc               → --damage off\n"
        << "  --damage-qc-scan-reads N     → --damage-scan N\n"
        << "\nPCR (informational):\n"
        << "  --pcr-cycles INT             Number of PCR cycles (log only)\n"
        << "  --pcr-efficiency FLOAT       Efficiency per cycle, 0-1 (default: 1.0)\n"
        << "  --pcr-error-rate FLOAT       Error rate sub/base/doubling (default: 5.3e-7 = Q5)\n"
        << "\nMemory: ~16 bytes per input read + 2-bit seq arena (~1 byte/4 bp) for error correction.\n";
}

int derep_main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string in_path, out_path, cluster_path, fqcl_path, prior_fqcl_path;
    std::string out_damaged_path, out_undamaged_path;
    float split_threshold = 0.0f;
    int   min_length = 16;  // default 16 bp; 0=disabled (no lower cap)
    int   max_length = 0;   // 0=disabled (no upper cap)
    double max_mem_gb = 0.0;  // 0=disabled; preflight-abort ceiling for the Phase-3 peak
    enum class SplitModelMode { Auto, Bulk, Empirical } split_model_mode = SplitModelMode::Auto;
    bool use_revcomp = true;

    ErrCorParams errcor;
    errcor.enabled = true;
    bool bucket_cap_explicit  = false;   // default on: aDNA primary use case

    bool   damage_auto    = false;  // default off: damage-aware hashing distorts
                                    // downstream damage analysis (e.g. deamination fit) when
                                    // run on fqdup output. Use --damage-auto explicitly.
    bool    damage_qc_enabled    = true;
    int64_t damage_qc_scan_reads = 1'000'000;
    // deamination model fit sample cap. 0 = use all reads (was the silent default,
    // which scans the entire input before Pass 1 even starts). 5M is plenty
    // for stable d_max/lambda estimates and saves a full input pass on large
    // libraries.
    int64_t damage_deam_max_reads = 5'000'000;
    DamageClipPolicy damage_clip_policy = DamageClipPolicy::ReportOnly;
    taph::SampleDamageProfile::LibraryType forced_library_type =
        taph::SampleDamageProfile::LibraryType::UNKNOWN;  // auto-detect by default
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
            print_usage(argv[0], /*advanced=*/false);
            return 0;
        } else if (arg == "--help-advanced") {
            print_usage(argv[0], /*advanced=*/true);
            return 0;
        } else if (arg == "-i" && i + 1 < argc) {
            in_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg == "--out-damaged" && i + 1 < argc) {
            out_damaged_path = argv[++i];
        } else if (arg == "--out-undamaged" && i + 1 < argc) {
            out_undamaged_path = argv[++i];
        } else if (arg == "--split-threshold" && i + 1 < argc) {
            split_threshold = std::stof(argv[++i]);
        } else if (arg == "--min-length" && i + 1 < argc) {
            min_length = std::stoi(argv[++i]);
        } else if (arg == "--max-length" && i + 1 < argc) {
            max_length = std::stoi(argv[++i]);
        } else if (arg == "--max-mem" && i + 1 < argc) {
            max_mem_gb = std::stod(argv[++i]);
        } else if (arg == "--split-model" && i + 1 < argc) {
            std::string m(argv[++i]);
            if      (m == "auto")     split_model_mode = SplitModelMode::Auto;
            else if (m == "bulk")     split_model_mode = SplitModelMode::Bulk;
            else if (m == "empirical") split_model_mode = SplitModelMode::Empirical;
            else { std::cerr << "Error: --split-model must be auto|bulk|empirical\n"; return 1; }
        } else if (arg == "-c" && i + 1 < argc) {
            cluster_path = argv[++i];
        } else if (arg == "--prior-fqcl" && i + 1 < argc) {
            prior_fqcl_path = argv[++i];
        } else if (arg == "--cluster-format" && i + 1 < argc) {
            fqcl_path = argv[++i];
        } else if (arg == "--no-revcomp") {
            use_revcomp = false;
        } else if (arg == "--error-correct") {
            errcor.enabled = true;
        } else if (arg == "--no-error-correct") {
            errcor.enabled = false;
        } else if (arg == "--no-b1-damage-adjust") {
            errcor.b1_damage_adjust = false;
        } else if (arg == "--b1-kappa-min" && i + 1 < argc) {
            errcor.b1_kappa_min = std::stod(argv[++i]);
            if (errcor.b1_kappa_min < 0.0) {
                std::cerr << "Error: --b1-kappa-min must be >= 0, got " << errcor.b1_kappa_min << "\n";
                return 1;
            }
        } else if (arg == "--b1-cs-min-total" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 1) {
                std::cerr << "Error: --b1-cs-min-total must be >= 1, got " << v << "\n";
                return 1;
            }
            errcor.b1_cs_min_total = static_cast<uint32_t>(v);
        } else if (arg == "--errcor-snp-threshold" && i + 1 < argc) {
            errcor.snp_threshold = std::stod(argv[++i]);
        } else if (arg == "--errcor-snp-min-count" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 0) { std::cerr << "Error: --errcor-snp-min-count must be >= 0, got " << v << "\n"; return 1; }
            errcor.snp_min_count = static_cast<uint32_t>(v);
        } else if (arg == "--errcor-bucket-cap" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 0) { std::cerr << "Error: --errcor-bucket-cap must be >= 0 (0 = unlimited), got " << v << "\n"; return 1; }
            errcor.bucket_cap = static_cast<uint32_t>(v);
            bucket_cap_explicit = true;
        } else if (arg == "--errcor-bucket-cap-lowcomplexity") {
            errcor.bucket_cap_lowcomplexity_only = true;
        } else if (arg == "--errcor-bucket-cap-complexity-frac" && i + 1 < argc) {
            errcor.bucket_cap_complexity_frac = std::stod(argv[++i]);
        } else if (arg == "--errcor-snp-cutoff" && i + 1 < argc) {
            long v = std::stol(argv[++i]);
            if (v < 1) { std::cerr << "Error: --errcor-snp-cutoff must be >= 1, got " << v << "\n"; return 1; }
            errcor.snp_low_cov_cutoff = static_cast<uint32_t>(v);
        } else if (arg == "--errcor-snp-factor" && i + 1 < argc) {
            errcor.snp_low_cov_factor = std::stod(argv[++i]);
        } else if (arg == "--errcor-empirical") {
            errcor.empirical   = true;
            errcor.legacy_veto = false;
        } else if (arg == "--errcor-legacy-veto") {
            errcor.legacy_veto = true;
        } else if (arg == "--protect-transversions") {
            errcor.protect_transversions = true;
        } else if (arg == "--errcor-adj-len") {
            errcor.adj_len_probe = true;
        } else if (arg == "--errcor-rescue-indels") {
            errcor.rescue_indels = true;
        } else if (arg == "--rescue-min-hits" && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 1) { std::cerr << "Error: --rescue-min-hits must be >= 1, got " << v << "\n"; return 1; }
            errcor.rescue_min_hits = v;
        } else if (arg == "--rescue-topk" && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 1) { std::cerr << "Error: --rescue-topk must be >= 1, got " << v << "\n"; return 1; }
            errcor.rescue_topk = static_cast<uint32_t>(v);
        } else if (arg == "--rescue-bundle-hot" && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 2) { std::cerr << "Error: --rescue-bundle-hot must be >= 2, got " << v << "\n"; return 1; }
            errcor.rescue_bundle_hot = static_cast<uint32_t>(v);
        } else if (arg == "--rescue-hash-hot" && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 1) { std::cerr << "Error: --rescue-hash-hot must be >= 1, got " << v << "\n"; return 1; }
            errcor.rescue_hash_hot = static_cast<uint32_t>(v);
        } else if (arg == "--rescue-alpha-ins" && i + 1 < argc) {
            errcor.rescue_alpha_ins = std::stod(argv[++i]);
        } else if (arg == "--rescue-alpha-del" && i + 1 < argc) {
            errcor.rescue_alpha_del = std::stod(argv[++i]);
        } else if (arg == "--rescue-mask-bonus" && i + 1 < argc) {
            errcor.rescue_mask_bonus = std::stod(argv[++i]);
        } else if (arg == "--b3-disable") {
            errcor.b3_enabled = false;
        } else if (arg == "--b3-min-n-elig" && i + 1 < argc) {
            errcor.b3_min_n_elig = std::stoi(argv[++i]);
        } else if (arg == "--b3-min-mass" && i + 1 < argc) {
            errcor.b3_min_mass = static_cast<float>(std::stod(argv[++i]));
        } else if (arg == "--b3-deam-threshold" && i + 1 < argc) {
            errcor.b3_deam_threshold = static_cast<float>(std::stod(argv[++i]));
        } else if (arg == "--b3-max-hamming" && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 3 || v > fqdup::b3::kMaxDiffs) {
                std::cerr << "Error: --b3-max-hamming must be 3.." << fqdup::b3::kMaxDiffs
                          << ", got " << v << "\n";
                return 1;
            }
            errcor.b3_max_hamming = v;
        } else if (arg == "--b3-count-ratio" && i + 1 < argc) {
            errcor.b3_count_ratio = std::stod(argv[++i]);
        } else if ((arg == "-p" || arg == "--threads") && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 0) { std::cerr << "Error: " << arg << " must be >= 0 (0 = auto), got " << v << "\n"; return 1; }
            errcor.threads = static_cast<unsigned>(v);
        } else if (arg == "--library-type" && i + 1 < argc) {
            std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                forced_library_type = taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                forced_library_type = taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt == "auto")
                forced_library_type = taph::SampleDamageProfile::LibraryType::UNKNOWN;
            else { std::cerr << "Error: Unknown --library-type: " << lt << " (use auto, ds, ss)\n"; return 1; }
        } else if (arg == "--damage" && i + 1 < argc) {
            std::string m = argv[++i];
            if      (m == "off")      { damage_qc_enabled = false; damage_auto = false; }
            else if (m == "report")   { damage_qc_enabled = true;  damage_auto = false; }
            else if (m == "collapse") { damage_qc_enabled = true;  damage_auto = true;  }
            else { std::cerr << "Error: --damage must be off|report|collapse, got " << m << "\n"; return 1; }
        } else if (arg == "--collapse-damage") {
            // Deprecated alias for --damage collapse
            damage_qc_enabled = true; damage_auto = true;
        } else if (arg == "--damage-dmax" && i + 1 < argc) {
            std::string v = argv[++i];
            auto comma = v.find(',');
            if (comma != std::string::npos) {
                damage_dmax5 = std::stod(v.substr(0, comma));
                damage_dmax3 = std::stod(v.substr(comma + 1));
            } else {
                damage_dmax5 = damage_dmax3 = std::stod(v);
            }
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
        } else if (arg == "--no-damage-qc") {
            damage_qc_enabled = false;
        } else if (arg == "--damage-deam-sample" && i + 1 < argc) {
            long long v = std::stoll(argv[++i]);
            if (v < 0) {
                std::cerr << "Error: --damage-deam-sample must be >= 0 (0 = scan all), got "
                          << v << "\n";
                return 1;
            }
            damage_deam_max_reads = static_cast<int64_t>(v);
        } else if ((arg == "--damage-scan" || arg == "--damage-qc-scan-reads") && i + 1 < argc) {
            long long v = std::stoll(argv[++i]);
            if (v < 0) {
                std::cerr << "Error: --damage-qc-scan-reads must be >= 0 (0 = scan all), got "
                          << v << "\n";
                return 1;
            }
            damage_qc_scan_reads = static_cast<int64_t>(v);
        } else if (arg == "--damage-clip-pass" && i + 1 < argc) {
            std::string v(argv[++i]);
            if      (v == "off")    damage_clip_policy = DamageClipPolicy::Off;
            else if (v == "report") damage_clip_policy = DamageClipPolicy::ReportOnly;
            else if (v == "refit")  damage_clip_policy = DamageClipPolicy::Refit;
            else {
                std::cerr << "Error: --damage-clip-pass must be off|report|refit, got " << v << "\n";
                return 1;
            }
        } else if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Error: Unknown argument: " << arg << "\n\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (out_path.empty() && out_damaged_path.empty() && out_undamaged_path.empty())
        out_path = "/dev/null";  // allow omitting -o when only split outputs wanted
    if (in_path.empty() || out_path.empty()) {
        std::cerr << "Error: Missing required arguments (-i and -o)\n\n";
        print_usage(argv[0]);
        return 1;
    }
    if (min_length < 0) {
        std::cerr << "Error: --min-length must be >= 0 (0 = disabled)\n";
        return 1;
    }
    if (max_length != 0 && max_length < min_length) {
        std::cerr << "Error: --max-length (" << max_length << ") must be 0 (disabled) or >= --min-length ("
                  << min_length << ")\n";
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
    log_info(std::string("Damage mode: ") +
             (damage_auto ? "collapse" :
              (damage_qc_enabled ? "report" : "off")) +
             (damage_dmax5 >= 0.0 ? " (manual d_max override)" : ""));
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
    DamageQcReport qc_report;
    LengthStratifiedDamageProfile lsd_data;

    try {
        // Always fit deamination + run QC so the .fqcl header reports full damage
        // stats (d_max, BIC, library evidence, adapter/hexamer QC) for
        // downstream consumers. Only gate the *use* of damage in clustering
        // on --collapse-damage via profile.enabled below.
        if (damage_auto || damage_qc_enabled) {
            DamageEstimateOptions opts;
            opts.mask_threshold     = mask_threshold;
            opts.forced_lib         = forced_library_type;
            opts.qc_enabled         = damage_qc_enabled;
            opts.adapter_scan_reads = damage_qc_scan_reads;
            opts.clip_policy        = damage_qc_enabled ? damage_clip_policy
                                                        : DamageClipPolicy::Off;
            opts.skip_deam_fit      = false;
            opts.max_reads          = damage_deam_max_reads;
            // Forward the derep thread count to the taph deamination+QC pass;
            // it defaulted to 1 (single-threaded), making Pass 0 ~2x slower than
            // the standalone `profile -p N` path doing identical work. Match
            // derep's 0->hardware_concurrency convention (see lines 688-689).
            opts.threads            = errcor.threads ? errcor.threads
                                                     : std::max(1u, std::thread::hardware_concurrency());
            if (damage_deam_max_reads > 0)
                log_info("Damage deamination fit: sampling first " +
                         std::to_string(damage_deam_max_reads) +
                         " reads (--damage-deam-sample 0 to scan all)");
            auto t_taph_begin = std::chrono::steady_clock::now();
            DamageEstimate est = estimate_damage_with_qc(in_path, opts);
            {
                double s = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t_taph_begin).count();
                char b[64];
                std::snprintf(b, sizeof(b), "%.1f s", s);
                log_info("Phase timer: taph (deamination+QC) = " + std::string(b));
            }
            profile = est.profile;
            profile.pcr_error_rate = pcr_total_rate;
            // profile.enabled controls whether derep CONSUMES damage for
            // clustering decisions. Default off — preserve damage signal.
            if (!damage_auto) profile.enabled = false;
            qc_report = est.qc;

            // Build per-bin empirical split model when split output is requested.
            // auto  = empirical if damage detected (d_max>0.01), else bulk fallback.
            // bulk  = skip extra pass entirely.
            // empirical = always run the extra pass.
            const bool want_split = !out_damaged_path.empty() || !out_undamaged_path.empty();
            const bool run_empirical =
                want_split &&
                split_model_mode != SplitModelMode::Bulk &&
                (split_model_mode == SplitModelMode::Empirical ||
                 // max(d5,d3): 3'-only DS damage was gated out of the empirical model (see split.cpp).
                 std::max(est.profile.d_max_5prime, est.profile.d_max_3prime) > 0.01);
            if (run_empirical) {
                const std::vector<uint64_t>* hist_ptr =
                    est.lsd_hist.empty() ? nullptr : &est.lsd_hist;
                auto t_lsd = std::chrono::steady_clock::now();
                lsd_data = estimate_damage_split_model(
                    in_path, est.profile, hist_ptr);
                double s = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t_lsd).count();
                char tb[64];
                std::snprintf(tb, sizeof(tb), "%.1f s", s);
                log_info("Phase timer: length-stratified LLR model = " +
                         std::string(tb) + " (" +
                         std::to_string(lsd_data.bins.size()) + " bins)");
            } else if (want_split) {
                log_info("Split model: bulk exponential (--split-model bulk or no damage detected)");
            }
        }
        if (damage_dmax5 >= 0.0) {
            profile.d_max_5prime  = damage_dmax5;
            profile.d_max_3prime  = (damage_dmax3 >= 0.0) ? damage_dmax3 : damage_dmax5;
            profile.lambda_5prime = damage_lambda5;
            profile.lambda_3prime = damage_lambda3;
            profile.background    = damage_bg;
            profile.mask_threshold = mask_threshold;
            profile.pcr_error_rate = pcr_total_rate;
            profile.enabled       = (damage_dmax5 > 0.0);
            profile.ss_mode       = (forced_library_type ==
                                     taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
            profile.populate_mask_from_model();
            profile.print_info(/*typical_read_length=*/0);
        }
        if (damage_auto && errcor.enabled && !profile.enabled) {
            log_warn("WARNING: --error-correct is ON but no ancient DNA damage detected. "
                     "On modern DNA this may absorb genuine low-frequency variants. "
                     "Use --no-error-correct to disable.");
        }

        if (!fqcl_path.empty() && !errcor.enabled) {
            log_warn("--cluster-format requires PCR error correction; ignoring (run without --no-error-correct).");
            fqcl_path.clear();
        }

        // Build QC JSON object string (opaque payload for cluster_format).
        std::string qc_json;
        if (qc_report.enabled) {
            std::ostringstream q;
            q.precision(6);
            const char* policy_str =
                damage_clip_policy == DamageClipPolicy::Off    ? "off" :
                damage_clip_policy == DamageClipPolicy::Refit  ? "refit" : "report";
            q << "{\"enabled\":true"
              << ",\"clip_policy\":\"" << policy_str << "\""
              << ",\"profile_clipped\":" << (qc_report.profile_clipped ? "true" : "false")
              << ",\"adapter_reads_scanned\":" << qc_report.adapter_reads_scanned
              << ",\"adapter_stubs_5prime\":[";
            for (size_t i = 0; i < qc_report.adapter.stubs5.size(); ++i) {
                if (i) q << ",";
                q << "\"" << qc_report.adapter.stubs5[i] << "\"";
            }
            q << "],\"adapter_stubs_3prime\":[";
            for (size_t i = 0; i < qc_report.adapter.stubs3.size(); ++i) {
                if (i) q << ",";
                q << "\"" << qc_report.adapter.stubs3[i] << "\"";
            }
            q << "],\"hexamer_entropy_5prime\":" << qc_report.hex_stats.entropy_terminal
              << ",\"hexamer_terminal_interior_jsd\":" << qc_report.hex_stats.jsd
              << ",\"d5_hexamer_corrected\":" << qc_report.preservation.d5_hexamer_corrected
              << ",\"flags\":[";
            for (size_t i = 0; i < qc_report.flag_names.size(); ++i) {
                if (i) q << ",";
                q << "\"" << qc_report.flag_names[i] << "\"";
            }
            q << "]}";
            qc_json = q.str();
        }

        // Bug 1: resolve library_type from the user-provided CLI flag.
        // Order: explicit --library-type ds/ss > auto-detected by the deamination fit
        // (when damage_auto produced an evaluable result) > "auto".
        std::string library_type_str = "auto";
        if (forced_library_type ==
                taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED) {
            library_type_str = "ds";
        } else if (forced_library_type ==
                taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED) {
            library_type_str = "ss";
        } else if (profile.enabled) {
            library_type_str = profile.ss_mode ? "ss" : "ds";
        }

        if (!prior_fqcl_path.empty())
            log_info("Prior fqcl: " + prior_fqcl_path);
        DerepEngine engine(use_revcomp, !cluster_path.empty(), profile, errcor,
                           false, bucket_cap_explicit, fqcl_path, in_path,
                           std::move(qc_json), library_type_str, prior_fqcl_path);
        float effective_split_threshold = split_threshold;
        if (!out_damaged_path.empty() || !out_undamaged_path.empty()) {
            if (profile.d_max_5prime > 0.01f) {
                LsdClassifyParams cls = make_lsd_classify_params(profile);
                if (cls.d_anc > 0.01) {
                    double pi = std::clamp(
                        static_cast<double>(profile.d_max_5prime) / cls.d_anc,
                        1e-6, 1.0 - 1e-6);
                    double log_prior_odds = std::log(pi / (1.0 - pi));
                    effective_split_threshold = static_cast<float>(split_threshold - log_prior_odds);
                    char b[80];
                    std::snprintf(b, sizeof(b), "%.4f  (pi=%.4f  log_prior_odds=%.4f)",
                                  effective_split_threshold, pi, log_prior_odds);
                    log_info("Posterior-adjusted split threshold: " + std::string(b));
                }
            }
            engine.set_split_paths(out_damaged_path, out_undamaged_path, effective_split_threshold);
        }
        engine.set_split_model(DamageSplitModel::build(lsd_data, profile));
        engine.set_length_filter(min_length, max_length);
        engine.set_max_mem(max_mem_gb);
        engine.process(in_path, out_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
