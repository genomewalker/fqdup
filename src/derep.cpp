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
#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"
#include "flat_hash_map.hpp"

#include "derep_detail/encoding.hpp"
#include "derep_detail/packed_ops.hpp"
#include "derep_detail/arena.hpp"
#include "derep_detail/damage_keys.hpp"

// All file-local types are in an anonymous namespace to give their member
// functions internal linkage, avoiding ODR violations with other TUs.
namespace {

using namespace fqdup::derep_detail;

// ============================================================================
// DerepEngine — single-file two-pass deduplication
// ============================================================================

class DerepEngine {
public:
    DerepEngine(bool use_revcomp, bool write_clusters,
                const DamageProfile& profile = DamageProfile{},
                const ErrCorParams& errcor = ErrCorParams{},
                bool /*errcor_mode_explicit*/ = false,
                bool bucket_cap_explicit = false,
                std::string fqcl_path = "",
                std::string input_fastq = "")
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          profile_(profile), errcor_(errcor),
          bucket_cap_explicit_(bucket_cap_explicit),
          fqcl_path_(std::move(fqcl_path)),
          fqcl_input_fastq_(std::move(input_fastq)),
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
                        if (!rec.qual.empty()) qual_arena_.append(rec.qual);
                        else                   qual_arena_.append_empty();
                    } else {
                        // Store in canonical orientation so Phase 3 comparisons are correct.
                        // Sequence: reverse-complement; quality: reverse only.
                        int L = static_cast<int>(rec.seq.size());
                        if (rc_buf_.size()  < static_cast<size_t>(L)) rc_buf_.resize(L);
                        if (rc_qbuf_.size() < static_cast<size_t>(L)) rc_qbuf_.resize(L);
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
        if (errcor_.adj_len_probe && errcor_.legacy_veto)
            log_warn("--errcor-adj-len: requires the empirical posterior-odds rule (legacy SNP veto has no indel notion); probe disabled.");
        if (errcor_.adj_len_probe && !errcor_.legacy_veto && !fqcl_path_.empty())
            log_warn("--errcor-adj-len + --cluster-format: indel-edge absorptions are NOT recorded in the .fqcl genealogy "
                     "(ChildMismatch encodes substitutions only). Dedup output is correct, but the genealogy will under-report "
                     "edges. Disable --errcor-adj-len if exact edge accounting is required.");
        log_info(errcor_.legacy_veto
                 ? "Phase 3 decision rule: legacy SNP veto"
                 : "Phase 3 decision rule: empirical posterior odds (S > 0)");
        if (arena_.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("Phase 3: arena size exceeds uint32_t range");
        const uint32_t N = static_cast<uint32_t>(arena_.size());
        if (N == 0) return;

        // Trailing padding for safe extraction in extract_packed_part
        arena_.packed.push_back(0);

        // T5.3: precompute bundle_key + bundle occupancy for LR prior. The
        // bundle_key (start_kmer ⊕ end_kmer) is our reference-free locus
        // proxy — clusters sharing a key landed at the same capture site, so
        // a high-occupancy bundle is evidence that mismatches in this group
        // could be PCR siblings rather than independent capture events.
        std::vector<uint64_t> bundle_key_of(N, 0);
        // For the LR prior we only need bundle occupancy (a count). The full
        // member list is only needed by the T5.6 adj-len probe, which uses a
        // separate length-agnostic map. Keep this one count-only — saves a
        // vector header per bundle (~24 B × ~N/avg_occ on capture data).
        ska::flat_hash_map<uint64_t, uint32_t> bundle_occ_map;
        bundle_occ_map.reserve(N);
        {
            std::vector<uint8_t> dec;
            std::vector<char>    asc;
            for (uint32_t id = 0; id < N; ++id) {
                if (!arena_.is_eligible(id)) continue;
                int L = arena_.length(id);
                if ((int)dec.size() < L) { dec.resize(L); asc.resize(L); }
                arena_.decode_range(id, 0, L, dec.data());
                for (int i = 0; i < L; ++i) asc[i] = static_cast<char>(dec[i]);
                uint64_t k = fqdup::bundlekey::from_decoded(
                    asc.data(), L, fqdup::bundlekey::kDefaultEndK);
                bundle_key_of[id] = k;
                ++bundle_occ_map[k];
            }
        }
        auto bundle_occ_of = [&](uint32_t id) -> uint64_t {
            auto it = bundle_occ_map.find(bundle_key_of[id]);
            return it == bundle_occ_map.end() ? 1u : it->second;
        };

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
        // Below this, even brute-force pairwise comparison is unsafe (random match
        // probability for a single mismatch at length 7 is non-negligible).
        constexpr int kBruteforceMinLen = 8;
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

        // Short-interior fallback: for ilen in [kBruteforceMinLen, kMinInteriorLen),
        // pigeonhole splits degenerate. Index parents in same-length buckets and
        // brute-force scan in B1 — slower per child but bounded by bucket size.
        ska::flat_hash_map<int, std::vector<uint32_t>> short_parents;

        uint64_t n_parents = 0;
        for (uint32_t id = 0; id < N; ++id) {
            if (!arena_.is_eligible(id)) continue;
            int L = arena_.length(id);
            const auto& lay = get_layout(L);
            if (lay.ilen < kMinInteriorLen) {
                if (lay.ilen >= kBruteforceMinLen)
                    short_parents[lay.ilen].push_back(id);
                continue;
            }

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

        uint64_t n_short_parents = 0;
        for (const auto& [ilen, v] : short_parents) n_short_parents += v.size();
        stats.short_brute_parents = n_short_parents;
        if (n_short_parents > 0)
            log_info("Phase 3: short-interior fallback indexes " +
                     std::to_string(n_short_parents) + " parents across " +
                     std::to_string(short_parents.size()) + " length buckets [" +
                     std::to_string(kBruteforceMinLen) + "," +
                     std::to_string(kMinInteriorLen - 1) + "]");

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
                        if (errcor_.bucket_cap == 0 || cnt < errcor_.bucket_cap) {
                            pi.ids.push_back(ev[i].id); cnt++;
                        } else {
                            stats.bucket_overflow_drops++;
                        }
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

        // ── Phase B1: collect ChildMismatch records (parallel) ──────────────
        // Pre-warm layout_cache before workers start so it becomes read-only.
        for (uint32_t id = 0; id < N; ++id) {
            if (!arena_.is_eligible(id)) continue;
            (void)get_layout(arena_.length(id));
        }

        unsigned n_threads = errcor_.threads;
        if (n_threads == 0) n_threads = std::max(1u, std::thread::hardware_concurrency());

        struct B1LocalStats {
            uint64_t children_scanned       = 0;
            uint64_t total_candidates       = 0;
            uint64_t children_found         = 0;
            uint64_t short_interior_skipped = 0;
            uint64_t short_brute_evaluated  = 0;
            uint64_t short_brute_found      = 0;
            uint64_t short_too_small_skipped = 0;
            double   decode_hash_child_ms   = 0;
            double   query_ms               = 0;
            double   check_ms               = 0;
        };

        std::vector<std::vector<ChildMismatch>> per_thread_mm(n_threads);
        std::vector<B1LocalStats>              per_thread_stats(n_threads);

        constexpr uint32_t kChunkSize = 256;
        const uint32_t n_chunks = (N + kChunkSize - 1) / kChunkSize;
        std::atomic<uint32_t> next_chunk{0};

        auto b1_worker = [&](unsigned tid) {
            auto& local_mm = per_thread_mm[tid];
            auto& ls       = per_thread_stats[tid];
            local_mm.reserve(std::min(static_cast<size_t>(N) * 2 / n_threads + 1024,
                                      static_cast<size_t>(1 << 22)));

            std::vector<uint8_t>  scratch;
            scratch.reserve(65535);
            std::vector<uint32_t> cand_storage;
            if (errcor_.bucket_cap > 0)
                cand_storage.resize(12u * errcor_.bucket_cap);
            else
                cand_storage.reserve(256);

            while (true) {
                uint32_t chunk = next_chunk.fetch_add(1, std::memory_order_relaxed);
                if (chunk >= n_chunks) break;
                uint32_t cid_lo = chunk * kChunkSize;
                uint32_t cid_hi = std::min(cid_lo + kChunkSize, N);

                for (uint32_t cid = cid_lo; cid < cid_hi; ++cid) {
                    if (is_error_[cid]) continue;
                    if (!arena_.is_eligible(cid)) continue;

                    int L = arena_.length(cid);
                    const auto& lay = get_layout(L);
                    if (lay.ilen < kMinInteriorLen) {
                        ls.short_interior_skipped++;
                        if (lay.ilen < kBruteforceMinLen) {
                            ls.short_too_small_skipped++;
                            continue;
                        }
                        // Brute-force same-length parent scan.
                        auto sp_it = short_parents.find(lay.ilen);
                        if (sp_it == short_parents.end()) continue;
                        const auto& parents_vec = sp_it->second;

                        int nf = (lay.ilen + 3) / 4;
                        size_t need = static_cast<size_t>(3 * nf);
                        if (scratch.size() < need) scratch.resize(need);
                        uint8_t* ci_full_s  = scratch.data();
                        uint8_t* crc_full_s = ci_full_s  + nf;
                        uint8_t* pi_buf_s   = crc_full_s + nf;

                        extract_packed_part(arena_.data(cid), lay.k5, lay.ilen, ci_full_s);
                        compute_interior_rc(ci_full_s, lay.ilen, crc_full_s);

                        ls.short_brute_evaluated++;
                        for (uint32_t pid : parents_vec) {
                            if (is_error_[pid]) continue;
                            if (pid == cid) continue;
                            if (id_count[pid] < id_count[cid]) continue;
                            if (id_count[pid] == id_count[cid] && pid >= cid) continue;
                            extract_packed_part(arena_.data(pid), lay.k5, lay.ilen, pi_buf_s);

                            MismatchInfo mm = packed_find_mismatch(
                                pi_buf_s, ci_full_s, 0, lay.ilen, profile_.enabled);
                            if (mm.found) {
                                local_mm.push_back({pid, cid, mm.position,
                                                    mm.base_b, mm.base_a,
                                                    0, 0, 0, 1, {}});
                                ls.short_brute_found++;
                                ls.children_found++;
                                continue;
                            }
                            mm = packed_find_mismatch(
                                pi_buf_s, crc_full_s, 0, lay.ilen, profile_.enabled);
                            if (mm.found) {
                                uint16_t cpos = static_cast<uint16_t>(lay.ilen - 1 - mm.position);
                                local_mm.push_back({pid, cid, cpos,
                                    static_cast<uint8_t>(mm.base_b ^ 0x3u),
                                    static_cast<uint8_t>(mm.base_a ^ 0x3u),
                                    0, 0, 0, 1, {}});
                                ls.short_brute_found++;
                                ls.children_found++;
                            }
                        }
                        continue;
                    }

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
            ls.decode_hash_child_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

            ls.children_scanned++;

            // Collect unique parent candidates — query both canonical and RC hash keys.
            // With bucket_cap=0 (default), the storage grows; with bucket_cap>0 it is
            // pre-sized and any overflow is counted (paired with bucket_overflow_drops).
            uint32_t n_cands = 0;
            const bool unbounded = (errcor_.bucket_cap == 0);
            auto collect = [&](uint32_t pid) {
                if (unbounded) {
                    if (n_cands >= cand_storage.size()) cand_storage.resize(n_cands + 1);
                    cand_storage[n_cands++] = pid;
                } else if (n_cands < static_cast<uint32_t>(cand_storage.size())) {
                    cand_storage[n_cands++] = pid;
                }
            };
            uint32_t* cand_buf = cand_storage.data();

            auto tq0 = clk::now();
            auto& sh = shard_it->second;
            for (int t = 0; t < 6; ++t) {
                sh.pi[t].query(pair_key(h[kPairA[t]],  h[kPairB[t]],  t, lay.ilen), collect);
                sh.pi[t].query(pair_key(rh[kPairA[t]], rh[kPairB[t]], t, lay.ilen), collect);
            }
            auto tq1 = clk::now();
            ls.query_ms += std::chrono::duration<double, std::milli>(tq1 - tq0).count();

            // Refresh: collect() may have reallocated cand_storage in unbounded mode.
            cand_buf = cand_storage.data();
            std::sort(cand_buf, cand_buf + n_cands);
            n_cands = static_cast<uint32_t>(
                std::unique(cand_buf, cand_buf + n_cands) - cand_buf);
            ls.total_candidates += n_cands;

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
                    local_mm.push_back({pid, cid, mm.position, mm.base_b, mm.base_a,
                                          0, 0, 0, 1, {}});
                    ls.children_found++;
                    continue;
                }
                // Try RC comparison (H=1)
                mm = packed_find_mismatch(pi_buf, crc_full, 0, lay.ilen, profile_.enabled);
                if (mm.found) {
                    uint16_t cpos = static_cast<uint16_t>(lay.ilen - 1 - mm.position);
                    local_mm.push_back({pid, cid, cpos,
                                          static_cast<uint8_t>(mm.base_b ^ 0x3u),
                                          static_cast<uint8_t>(mm.base_a ^ 0x3u),
                                          0, 0, 0, 1, {}});
                    ls.children_found++;
                    continue;
                }
                // H=2 path: only for reads with count ≤ max_h2_count.
                // Both mismatches must be non-damage transversions (A↔T or C↔G, xr==3).
                if (id_count[cid] <= errcor_.max_h2_count) {
                    MismatchInfo2 mm2 = packed_find_mismatches2(pi_buf, ci_full, 0, lay.ilen);
                    if (mm2.count == 2 &&
                        !is_damage_sub_packed(mm2.base_a[0], mm2.base_b[0], false) &&
                        !is_damage_sub_packed(mm2.base_a[1], mm2.base_b[1], false)) {
                        local_mm.push_back({pid, cid,
                                              mm2.pos[0], mm2.base_b[0], mm2.base_a[0],
                                              mm2.pos[1], mm2.base_b[1], mm2.base_a[1],
                                              2, {}});
                        ls.children_found++;
                        continue;
                    }
                    mm2 = packed_find_mismatches2(pi_buf, crc_full, 0, lay.ilen);
                    if (mm2.count == 2 &&
                        !is_damage_sub_packed(mm2.base_a[0], mm2.base_b[0], false) &&
                        !is_damage_sub_packed(mm2.base_a[1], mm2.base_b[1], false)) {
                        // Convert RC positions/bases to child canonical frame
                        uint16_t p0 = static_cast<uint16_t>(lay.ilen - 1 - mm2.pos[0]);
                        uint16_t p1 = static_cast<uint16_t>(lay.ilen - 1 - mm2.pos[1]);
                        local_mm.push_back({pid, cid,
                                              p0, static_cast<uint8_t>(mm2.base_b[0] ^ 0x3u),
                                                  static_cast<uint8_t>(mm2.base_a[0] ^ 0x3u),
                                              p1, static_cast<uint8_t>(mm2.base_b[1] ^ 0x3u),
                                                  static_cast<uint8_t>(mm2.base_a[1] ^ 0x3u),
                                              2, {}});
                        ls.children_found++;
                    }
                }
            }
            ls.check_ms += std::chrono::duration<double, std::milli>(clk::now() - tc0).count();
                }  // end for cid in chunk
            }  // end while chunk dispatcher
        };  // end b1_worker

        if (n_threads == 1) {
            b1_worker(0);
        } else {
            std::vector<std::thread> ts;
            ts.reserve(n_threads - 1);
            for (unsigned t = 1; t < n_threads; ++t) ts.emplace_back(b1_worker, t);
            b1_worker(0);
            for (auto& th : ts) th.join();
        }

        // Merge per-thread state. B2 sorts mismatches deterministically, but
        // we concat in tid order so the same -j N gives byte-identical output.
        size_t total_mm = 0;
        for (const auto& v : per_thread_mm) total_mm += v.size();
        std::vector<ChildMismatch> mismatches;
        mismatches.reserve(total_mm);
        for (auto& v : per_thread_mm) {
            mismatches.insert(mismatches.end(), v.begin(), v.end());
            std::vector<ChildMismatch>().swap(v);
        }
        for (const auto& ls : per_thread_stats) {
            stats.children_scanned       += ls.children_scanned;
            stats.total_candidates       += ls.total_candidates;
            stats.children_found         += ls.children_found;
            stats.decode_hash_child_ms   += ls.decode_hash_child_ms;
            stats.query_ms               += ls.query_ms;
            stats.check_ms               += ls.check_ms;
            short_interior_skipped       += ls.short_interior_skipped;
            stats.short_brute_evaluated  += ls.short_brute_evaluated;
            stats.short_brute_found      += ls.short_brute_found;
            stats.short_too_small_skipped += ls.short_too_small_skipped;
        }

        log_info("Phase 3 B1: " + std::to_string(n_threads) + " thread" +
                 (n_threads == 1 ? "" : "s") + ", found " +
                 std::to_string(stats.children_found) +
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

        // ── T5.8: empirical posterior-odds model ────────────────────────────
        // Build EdgeCandidate features for every B1 edge and fit the model
        // (per-bin P_real from cross-bundle recurrence; per-occ-bin
        // log(π_pcr/π_real) from singleton-vs-recurring split). The fitted
        // model replaces lr_threshold + kP_real_uniform + kBundlePrior_alpha.
        constexpr int kDamageEdgeMargin_T58 = 5;
        auto build_edge = [&](const ChildMismatch& cm) -> fqdup::errcor_emp::EdgeCandidate {
            const int L  = static_cast<int>(arena_.length(cm.child_id));
            const int k5 = get_layout(L).k5;
            const int ilen = get_layout(L).ilen;
            fqdup::errcor_emp::EdgeCandidate e{};
            e.child_id   = cm.child_id;
            e.parent_id  = cm.parent_id;
            e.bundle_key = bundle_key_of[cm.parent_id];
            e.bundle_occ = bundle_occ_of(cm.parent_id);
            e.child_len  = static_cast<uint16_t>(L);
            e.n_mm       = cm.hamming;
            auto fill = [&](int i, uint16_t mm_pos, uint8_t alt, uint8_t pb) {
                uint16_t pos_full = static_cast<uint16_t>(mm_pos + k5);
                e.mm[i].pos        = pos_full;
                e.mm[i].parent_2b  = pb;
                e.mm[i].alt_2b     = alt;
                e.mm[i].qual       = qual_arena_.q_at(cm.child_id, pos_full);
                bool dmg_xor       = ((alt ^ pb) == 2u);
                bool in_zone       = (mm_pos < kDamageEdgeMargin_T58 ||
                                      mm_pos >= static_cast<uint16_t>(ilen - kDamageEdgeMargin_T58));
                e.mm[i].damage_chan = (dmg_xor && in_zone) ? 1 : 0;
                e.mm[i].p_damage    = (dmg_xor && profile_.enabled)
                                      ? profile_.p_damage_at(pos_full, L) : 0.0;
            };
            fill(0, cm.mismatch_pos, cm.alt_base, cm.parent_base);
            if (cm.hamming == 2)
                fill(1, cm.mismatch_pos2, cm.alt_base2, cm.parent_base2);
            return e;
        };

        fqdup::errcor_emp::ErrCorEmpiricalModel emp_model;
        if (errcor_.empirical && !errcor_.legacy_veto && !mismatches.empty()) {
            std::vector<fqdup::errcor_emp::EdgeCandidate> all_edges;
            all_edges.reserve(mismatches.size());
            std::unordered_set<uint64_t> distinct_bundles;
            distinct_bundles.reserve(mismatches.size());
            for (const auto& cm : mismatches) {
                all_edges.emplace_back(build_edge(cm));
                distinct_bundles.insert(bundle_key_of[cm.parent_id]);
            }
            emp_model.fit(all_edges, distinct_bundles.size());
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.4g", emp_model.p_real_global);
            log_info("Phase 3 empirical model: " + std::to_string(all_edges.size()) +
                     " edges, " + std::to_string(distinct_bundles.size()) +
                     " distinct bundles, p_real_global=" + buf);
            std::snprintf(buf, sizeof(buf), "%.3f", emp_model.log_pi_ratio_global);
            log_info("  log_pi_ratio_global=" + std::string(buf) + " nats; per-occ:");
            for (int o = 0; o < fqdup::errcor_emp::kNumOccBins; ++o) {
                std::snprintf(buf, sizeof(buf), "%.3f", emp_model.log_pi_ratio_occ[o]);
                log_info("    occ_bin[" + std::to_string(o) + "] = " + buf);
            }
        }

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
                    // T5.8: empirical posterior-odds for H=2.
                bool   absorb_h2;
                if (errcor_.empirical && !errcor_.legacy_veto) {
                    auto e = build_edge(cm);
                    e.bundle_key = bundle_key_of[eff_pid];
                    e.bundle_occ = bundle_occ_of(eff_pid);
                    double S = emp_model.score(e);
                    absorb_h2 = (S > 0.0);
                    if (snp_veto_h2) { stats.lr_sum_protected += S; ++stats.lr_n_protected; }
                    else             { stats.lr_sum_absorbed  += S; ++stats.lr_n_absorbed;  }
                    if (absorb_h2) ++stats.edge_lr_absorbed;
                    else           ++stats.edge_lr_protected;
                    int ob = fqdup::errcor_emp::occ_bin(e.bundle_occ);
                    if (absorb_h2) ++stats.edge_lr_absorbed_by_occ[ob];
                    else           ++stats.edge_lr_protected_by_occ[ob];
                } else {
                    absorb_h2 = !snp_veto_h2;
                    if (snp_veto_h2) ++stats.edge_legacy_veto;
                    else             ++stats.edge_legacy_absorb;
                }

                if (!absorb_h2) {
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
                            if (!fqcl_path_.empty()) fqcl_mismatches_.push_back(cm);
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
                if (damage_bypass) ++stats.edge_damage_bypass;

                // T5.8: empirical posterior-odds for H=1.
                bool absorb_h1;
                if (errcor_.empirical && !errcor_.legacy_veto) {
                    auto e = build_edge(cm);
                    e.bundle_key = bundle_key_of[eff_pid];
                    e.bundle_occ = bundle_occ_of(eff_pid);
                    double S = emp_model.score(e);
                    absorb_h1 = (S > 0.0);
                    if (snp_veto) { stats.lr_sum_protected += S; ++stats.lr_n_protected; }
                    else          { stats.lr_sum_absorbed  += S; ++stats.lr_n_absorbed;  }
                    if (absorb_h1) ++stats.edge_lr_absorbed;
                    else           ++stats.edge_lr_protected;
                    int ob = fqdup::errcor_emp::occ_bin(e.bundle_occ);
                    if (absorb_h1) ++stats.edge_lr_absorbed_by_occ[ob];
                    else           ++stats.edge_lr_protected_by_occ[ob];
                } else {
                    absorb_h1 = !snp_veto;
                    if (snp_veto) ++stats.edge_legacy_veto;
                    else          ++stats.edge_legacy_absorb;
                }

                if (!absorb_h1) {
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
                        if (!fqcl_path_.empty()) fqcl_mismatches_.push_back(cm);
                    }
                }
            }
        }

        // ── T5.6: Adjacent-length (L±1) probe ──────────────────────────────
        // Catches PCR siblings differing by a single insertion/deletion at any
        // position (damage-truncated reads, polymerase slippage). Restricted to
        // bundle members so candidate-set is locus-bounded — never N² over all
        // sequences. Per-pair cost: O(L) banded check; absorbs only on exact
        // 1-indel + 0-substitution match. LR rule applies as for H=1.
        if (errcor_.adj_len_probe && !errcor_.legacy_veto) {
            // Build a length-agnostic bundle map so L vs L±1 reads share keys.
            // Only constructed when the probe is enabled — zero overhead off-path.
            ska::flat_hash_map<uint64_t, std::vector<uint32_t>> adj_bundles;
            adj_bundles.reserve(N);
            {
                std::vector<uint8_t> dec;
                std::vector<char>    asc;
                for (uint32_t id = 0; id < N; ++id) {
                    if (!arena_.is_eligible(id)) continue;
                    int L = arena_.length(id);
                    if ((int)dec.size() < L) { dec.resize(L); asc.resize(L); }
                    arena_.decode_range(id, 0, L, dec.data());
                    for (int i = 0; i < L; ++i) asc[i] = static_cast<char>(dec[i]);
                    uint64_t k = fqdup::bundlekey::from_decoded_no_len(
                        asc.data(), L, fqdup::bundlekey::kDefaultEndK);
                    adj_bundles[k].push_back(id);
                }
            }

            // Per-bundle parent shortlist cap. On capture hotspots a bundle can
            // hold hundreds of clusters; pairing all-vs-all is O(M²). We restrict
            // candidate parents to the top-K members by id_count (the only ones
            // that can plausibly absorb others), then test every potential child
            // of adjacent length against those parents. Work is O(K·M_adj).
            constexpr size_t kAdjLenParentCap = 16;

            std::vector<uint8_t> abuf, bbuf;
            for (auto& kv : adj_bundles) {
                auto& members = kv.second;
                if (members.size() < 2) continue;
                std::vector<std::pair<int, uint32_t>> by_len;
                by_len.reserve(members.size());
                for (uint32_t mid : members)
                    by_len.emplace_back(static_cast<int>(arena_.length(mid)), mid);
                std::sort(by_len.begin(), by_len.end());
                const size_t M = by_len.size();

                // Build parent shortlist: top-K members by id_count.
                std::vector<uint32_t> parents;
                if (M <= kAdjLenParentCap) {
                    parents.reserve(M);
                    for (auto& p : by_len) parents.push_back(p.second);
                } else {
                    parents.reserve(M);
                    for (auto& p : by_len) parents.push_back(p.second);
                    std::partial_sort(
                        parents.begin(),
                        parents.begin() + kAdjLenParentCap,
                        parents.end(),
                        [&](uint32_t a, uint32_t b) {
                            if (id_count[a] != id_count[b]) return id_count[a] > id_count[b];
                            return a < b;  // stable tie-break for determinism
                        });
                    parents.resize(kAdjLenParentCap);
                }
                ska::flat_hash_set<uint32_t> parent_set(parents.begin(), parents.end());

                for (size_t i = 0; i < M; ++i) {
                    int Li = by_len[i].first;
                    for (size_t j = i + 1; j < M; ++j) {
                        int Lj = by_len[j].first;
                        if (Lj - Li > 1) break;
                        if (Lj - Li != 1) continue;  // only adjacent (skip same-length, handled by B2)
                        // At least one side must be in the parent shortlist —
                        // otherwise neither could plausibly be the absorber.
                        if (!parent_set.count(by_len[i].second) &&
                            !parent_set.count(by_len[j].second)) continue;
                        uint32_t a_id = by_len[i].second;
                        uint32_t b_id = by_len[j].second;  // longer
                        ++stats.adj_len_evaluated;

                        if ((int)abuf.size() < Li) abuf.resize(Li);
                        if ((int)bbuf.size() < Lj) bbuf.resize(Lj);
                        arena_.decode_range(a_id, 0, Li, abuf.data());
                        arena_.decode_range(b_id, 0, Lj, bbuf.data());

                        // Test: is `a` (shorter) exactly a 1-deletion of `b` (longer)?
                        int d = 0;
                        while (d < Li && bbuf[d] == abuf[d]) ++d;
                        bool ok = true;
                        for (int p = d; p < Li; ++p)
                            if (bbuf[p + 1] != abuf[p]) { ok = false; break; }
                        if (!ok) continue;
                        ++stats.adj_len_matched;

                        // Higher count wins parent role (matches B2 convention).
                        uint32_t parent_id, child_id;
                        if (id_count[a_id] >= id_count[b_id]) {
                            parent_id = a_id; child_id = b_id;
                        } else {
                            parent_id = b_id; child_id = a_id;
                        }
                        if (is_error_[child_id]) continue;
                        uint32_t eff_pid = parent_id;
                        if (is_error_[eff_pid]) eff_pid = find_root_chain(eff_pid);

                        // T5.8: posterior odds for adj-len indel edge. Treated
                        // as a single-mismatch edge at the indel position.
                        if (errcor_.empirical && !errcor_.legacy_veto) {
                            int L_child = static_cast<int>(arena_.length(child_id));
                            uint16_t pos16 = static_cast<uint16_t>(std::min(d, L_child - 1));
                            fqdup::errcor_emp::EdgeCandidate e{};
                            e.child_id   = child_id;
                            e.parent_id  = parent_id;
                            e.bundle_key = bundle_key_of[eff_pid];
                            e.bundle_occ = members.size();
                            e.child_len  = static_cast<uint16_t>(L_child);
                            e.n_mm       = 1;
                            e.mm[0].pos        = pos16;
                            e.mm[0].parent_2b  = 0;
                            e.mm[0].alt_2b     = 1;  // arbitrary non-equal pair: flagged as transversion
                            e.mm[0].qual       = qual_arena_.q_at(child_id, pos16);
                            e.mm[0].damage_chan = 0;
                            e.mm[0].p_damage    = 0.0;
                            double S = emp_model.score(e);
                            if (S <= 0.0) {
                                ++stats.adj_len_protected;
                                continue;
                            }
                        }

                        is_error_[child_id] = true;
                        parent_chain[child_id] = eff_pid;
                        acc_count[eff_pid] += acc_count[child_id];
                        acc_count[child_id] = 0;
                        if (seq_entry[child_id] && seq_entry[eff_pid] &&
                            seq_entry[child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                            seq_entry[eff_pid]->record_index = seq_entry[child_id]->record_index;
                            seq_entry[eff_pid]->damage_score = seq_entry[child_id]->damage_score;
                        }
                        ++stats.absorbed;
                        ++errcor_absorbed_;
                        ++stats.adj_len_absorbed;
                        // No fqcl record: ChildMismatch encodes substitutions only.
                        // The .fqcl genealogy will lack the indel edge but the
                        // dedup output is correct (is_error_ drops the read).
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
        // Persist for write_fqcl_ → metadata.
        loss_bucket_overflow_drops_   = stats.bucket_overflow_drops;
        loss_short_brute_evaluated_   = stats.short_brute_evaluated;
        loss_short_brute_found_       = stats.short_brute_found;
        loss_short_too_small_skipped_ = stats.short_too_small_skipped;
        // Only the truly-unrecoverable (ilen < kBruteforceMinLen) reads remain a "loss".
        loss_short_interior_skipped_ = stats.short_too_small_skipped;
        if (short_interior_skipped > 0) {
            double pct_short = 100.0 * static_cast<double>(short_interior_skipped) /
                               static_cast<double>(N);
            log_info("Phase 3: " + std::to_string(short_interior_skipped) +
                     " reads have interior < " + std::to_string(kMinInteriorLen) +
                     " bp (" + std::to_string(pct_short).substr(0, 4) + "%); " +
                     std::to_string(stats.short_brute_evaluated) + " ran brute-force fallback (" +
                     std::to_string(stats.short_brute_found) + " edges), " +
                     std::to_string(stats.short_too_small_skipped) +
                     " unrecoverable (ilen < " + std::to_string(kBruteforceMinLen) + ")");
            if (pct_short > 10.0)
                log_warn("Over 10% of reads have interiors too short for the pigeonhole index "
                         "after damage masking. Consider reducing --mask-threshold or "
                         "adjusting damage zone parameters.");
        }
        log_info("Phase 3 complete: absorbed " + std::to_string(errcor_absorbed_) +
                 " sequences (H=1 directed-ascent + H=2 direct-detection, snp_threshold=" +
                 std::to_string(errcor_.snp_threshold) + ")");

        if (!fqcl_path_.empty()) {
            fqcl_parent_chain_ = std::move(parent_chain);
            write_fqcl_(seq_entry);
        }
    }

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
        meta.tool_version  = "0.1";
        meta.input_fastq   = fqcl_input_fastq_;
        meta.n_input_reads = total_reads_;
        meta.library_type  = profile_.ss_mode ? "ss" : (profile_.enabled ? "ds" : "auto");
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
        meta.bucket_overflow_drops   = loss_bucket_overflow_drops_;
        meta.short_interior_skipped  = loss_short_interior_skipped_;
        meta.short_brute_evaluated   = loss_short_brute_evaluated_;
        meta.short_brute_found       = loss_short_brute_found_;
        meta.short_too_small_skipped = loss_short_too_small_skipped_;
        // v2: singletons are always written (as tinyblocks) — never silently dropped.

        cf::Writer writer(fqcl_path_, meta);

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
    std::string fqcl_path_;
    std::string fqcl_input_fastq_;
    std::vector<ChildMismatch> fqcl_mismatches_;
    std::vector<uint32_t>      fqcl_parent_chain_;

    ska::flat_hash_map<SequenceFingerprint, IndexEntry, SequenceFingerprintHash> index_;

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

    // Loss counters surfaced into .fqcl metadata (T1.3, T1.2, T3.1).
    uint64_t loss_bucket_overflow_drops_   = 0;
    uint64_t loss_short_interior_skipped_  = 0;
    uint64_t loss_short_brute_evaluated_   = 0;
    uint64_t loss_short_brute_found_       = 0;
    uint64_t loss_short_too_small_skipped_ = 0;
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
        << "  --errcor-bucket-cap INT      Max pair-key bucket size (default: 0 = unlimited).\n"
        << "                               Set >0 only as an OOM safety valve; any drops are\n"
        << "                               counted and warned (Bucket overflow drops).\n"
        << "  --errcor-snp-cutoff INT      Parent count below which SNP threshold is tightened (default: 10)\n"
        << "  --errcor-snp-factor FLOAT    SNP threshold multiplier at low coverage (default: 1.75)\n"
        << "  --errcor-threads INT         Phase 3 B1 worker threads (default: 1, 0 = hardware concurrency)\n"
        << "  --errcor-empirical           Use empirical posterior-odds rule (default; absorb iff S>0).\n"
        << "  --errcor-legacy-veto         Use legacy SNP-veto rule instead of the empirical model.\n"
        << "  --errcor-adj-len             Enable adjacent-length (L±1) indel probe (uses the empirical rule).\n"
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

    std::string in_path, out_path, cluster_path, fqcl_path;
    bool use_revcomp = true;

    ErrCorParams errcor;
    errcor.enabled = true;
    bool bucket_cap_explicit  = false;   // default on: aDNA primary use case

    bool   damage_auto    = false;  // default off: damage-aware hashing distorts
                                    // downstream damage analysis (e.g. DART) when
                                    // run on fqdup output. Use --damage-auto explicitly.
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
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-i" && i + 1 < argc) {
            in_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg == "-c" && i + 1 < argc) {
            cluster_path = argv[++i];
        } else if (arg == "--cluster-format" && i + 1 < argc) {
            fqcl_path = argv[++i];
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
            if (v < 0) { std::cerr << "Error: --errcor-bucket-cap must be >= 0 (0 = unlimited), got " << v << "\n"; return 1; }
            errcor.bucket_cap = static_cast<uint32_t>(v);
            bucket_cap_explicit = true;
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
        } else if (arg == "--errcor-adj-len") {
            errcor.adj_len_probe = true;
        } else if (arg == "--errcor-threads" && i + 1 < argc) {
            int v = std::stoi(argv[++i]);
            if (v < 0) { std::cerr << "Error: --errcor-threads must be >= 0 (0 = auto), got " << v << "\n"; return 1; }
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
        DerepEngine engine(use_revcomp, !cluster_path.empty(), profile, errcor,
                           false, bucket_cap_explicit, fqcl_path, in_path);
        engine.process(in_path, out_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
