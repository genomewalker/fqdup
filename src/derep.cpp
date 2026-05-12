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

class DerepEngine {
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
        log_info("Decompression: " + std::string(
#ifdef HAVE_RAPIDGZIP
            "rapidgzip (parallel multi-threaded)"
#elif defined(HAVE_ISAL)
            "ISA-L (hardware-accelerated)"
#else
            "zlib"
#endif
        ));

        auto t_pass1_begin = clk::now();
        pass1(in_path);
        auto t_pass1_end = clk::now();
        log_info("Phase timer: Pass 1 = " + fmt_secs(t_pass1_end - t_pass1_begin));

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

    void pass1(const std::string& in_path) {
        // Pre-reserve hash map + arenas based on input file size. Without this,
        // ska::flat_hash_map (max_load_factor=0.5) regrows ~26 times to reach
        // 100M entries, and SeqArena::packed doubles ~30 times — each regrow is
        // a full re-hash + memcpy of the existing data. On a 6.3 GB compressed
        // input this single change saves ~30 minutes of wall time.
        //
        // Estimate: ~250 B per uncompressed FASTQ record (typical aDNA SE).
        // Compressed gz is ~5×, so compressed_bytes / 50 ≈ records.
        // Use 0.6× safety factor (we'd rather one regrow at the end than
        // over-reserve memory by 2×).
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
                    log_info("Pass 1: pre-reserving capacity for ~" +
                             std::to_string(est_records) +
                             " entries (heuristic from input size " +
                             std::to_string(fsize / (1024 * 1024)) + " MB)");
                    index_.reserve(est_records);
                    if (errcor_.enabled) {
                        arena_.offsets.reserve(est_records);
                        arena_.lengths.reserve(est_records);
                        arena_.eligible.reserve(est_records);
                        // ~100 bp avg sequence → 25 bytes packed
                        arena_.packed.reserve(est_records * 25);
                        qual_arena_.offsets.reserve(est_records);
                        qual_arena_.lengths.reserve(est_records);
                        qual_arena_.bytes.reserve(est_records * 100);
                    }
                }
            }
        } catch (...) {
            // reserve is best-effort; any failure just means we fall back to
            // the regrow path.
        }

        auto reader = make_fastq_reader(in_path);
        FastqRecord rec;
        uint64_t record_idx = 0;

        while (reader->read(rec)) {
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
                        it->second.count = pit->second;
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

        // T8 (rescue): length-agnostic parallel bundle key + occupancy. Built
        // only when --errcor-rescue-indels is set so the standard path pays
        // nothing. Indel pairs (L vs L±1, L vs L±2) collapse to the same key
        // here, so the within-bundle gate in the rescue driver actually fires
        // for them. Occupancy is computed against this key so bundle_hot sees
        // the true PCR cloud size that spans neighbouring lengths.
        std::vector<uint64_t> rescue_bundle_key_of;
        ska::flat_hash_map<uint64_t, uint32_t> rescue_bundle_occ_map;
        if (errcor_.rescue_indels && !errcor_.legacy_veto && errcor_.empirical) {
            rescue_bundle_key_of.assign(N, 0);
            rescue_bundle_occ_map.reserve(N);
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
                rescue_bundle_key_of[id] = k;
                ++rescue_bundle_occ_map[k];
            }
        }
        auto rescue_bundle_occ_of = [&](uint32_t id) -> uint64_t {
            if (rescue_bundle_key_of.empty()) return 1u;
            auto it = rescue_bundle_occ_map.find(rescue_bundle_key_of[id]);
            return it == rescue_bundle_occ_map.end() ? 1u : it->second;
        };

        is_error_.assign(N, false);

        std::vector<uint64_t> id_count(N, 0);
        std::vector<uint32_t> id_fwd_count(N, 0);
        // Reverse map: seq_id → IndexEntry* for representative propagation during absorption.
        std::vector<IndexEntry*> seq_entry(N, nullptr);
        for (auto& [fp, entry] : index_) {
            id_count[entry.seq_id] = entry.count;
            id_fwd_count[entry.seq_id] = entry.fwd_count;
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
                                pi_buf_s, ci_full_s, 0, lay.ilen, errcor_.protect_transversions);
                            if (mm.found) {
                                local_mm.push_back({pid, cid, mm.position,
                                                    mm.base_b, mm.base_a,
                                                    0, 0, 0, 1, {}});
                                ls.short_brute_found++;
                                ls.children_found++;
                                continue;
                            }
                            mm = packed_find_mismatch(
                                pi_buf_s, crc_full_s, 0, lay.ilen, errcor_.protect_transversions);
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
                MismatchInfo mm = packed_find_mismatch(pi_buf, ci_full, 0, lay.ilen, errcor_.protect_transversions);
                if (mm.found) {
                    local_mm.push_back({pid, cid, mm.position, mm.base_b, mm.base_a,
                                          0, 0, 0, 1, {}});
                    ls.children_found++;
                    continue;
                }
                // Try RC comparison (H=1)
                mm = packed_find_mismatch(pi_buf, crc_full, 0, lay.ilen, errcor_.protect_transversions);
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
                // With b1_damage_adjust=true (default): terminal damage mismatches
                // (C→T at 5', G→A at 3' for DS; both channels at both ends for SS)
                // are not counted toward H_adj. h_adj=2 → record both; h_adj=1 →
                // record only the non-damage mismatch as H=1; h_adj=0 → skip (both
                // terminal damage, handled by Phase 2 tag-folding).
                // With b1_damage_adjust=false: reverts to old transversion-only rule.
                if (id_count[cid] <= errcor_.max_h2_count) {
                    auto admit_h2 = [&](MismatchInfo2 mm2, bool rc) -> bool {
                        if (mm2.count != 2) return false;
                        // Transversion guard: A↔T / C↔G (xr==3) protected regardless
                        // of b1_damage_adjust when --protect-transversions is active.
                        if (errcor_.protect_transversions) {
                            if ((mm2.base_a[0] ^ mm2.base_b[0]) == 3u) return false;
                            if ((mm2.base_a[1] ^ mm2.base_b[1]) == 3u) return false;
                        }
                        bool dmg0, dmg1;
                        if (errcor_.b1_damage_adjust) {
                            dmg0 = is_terminal_damage_mismatch(
                                       mm2.base_a[0], mm2.base_b[0],
                                       mm2.pos[0], lay.ilen, profile_.ss_mode);
                            dmg1 = is_terminal_damage_mismatch(
                                       mm2.base_a[1], mm2.base_b[1],
                                       mm2.pos[1], lay.ilen, profile_.ss_mode);
                        } else {
                            dmg0 = is_damage_sub_packed(mm2.base_a[0], mm2.base_b[0], errcor_.protect_transversions);
                            dmg1 = is_damage_sub_packed(mm2.base_a[1], mm2.base_b[1], errcor_.protect_transversions);
                            if (dmg0 || dmg1) return false; // old: reject any damage type
                        }
                        int h_adj = 2 - (int)dmg0 - (int)dmg1;
                        if (h_adj == 0) return false; // both terminal damage → skip (Phase 2)
                        auto cvt = [&](int i, uint16_t p, uint8_t ab, uint8_t pb) {
                            if (rc) {
                                p  = static_cast<uint16_t>(lay.ilen - 1 - p);
                                ab = static_cast<uint8_t>(ab ^ 0x3u);
                                pb = static_cast<uint8_t>(pb ^ 0x3u);
                            }
                            return std::make_tuple(p, ab, pb);
                        };
                        if (h_adj == 1) {
                            int ri = dmg0 ? 1 : 0;
                            auto [p, ab, pb] = cvt(ri, mm2.pos[ri], mm2.base_b[ri], mm2.base_a[ri]);
                            local_mm.push_back({pid, cid, p, ab, pb, 0, 0, 0, 1, {}});
                        } else {
                            auto [p0, ab0, pb0] = cvt(0, mm2.pos[0], mm2.base_b[0], mm2.base_a[0]);
                            auto [p1, ab1, pb1] = cvt(1, mm2.pos[1], mm2.base_b[1], mm2.base_a[1]);
                            local_mm.push_back({pid, cid, p0, ab0, pb0, p1, ab1, pb1, 2, {}});
                        }
                        ls.children_found++;
                        return true;
                    };
                    MismatchInfo2 mm2 = packed_find_mismatches2(pi_buf, ci_full, 0, lay.ilen);
                    if (admit_h2(mm2, false)) continue;
                    mm2 = packed_find_mismatches2(pi_buf, crc_full, 0, lay.ilen);
                    if (admit_h2(mm2, true)) { /* admitted */ }
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
                const uint32_t dmg_xr = alt ^ pb;
                const bool is_ct_e  = (dmg_xr == 2u) && ((alt & pb) == 1u);
                const bool is_ga_e  = (dmg_xr == 2u) && ((alt & pb) == 0u);
                const bool at_5_e   = mm_pos < kDamageEdgeMargin_T58;
                const bool at_3_e   = mm_pos >= static_cast<uint16_t>(ilen - kDamageEdgeMargin_T58);
                e.mm[i].damage_chan = (profile_.ss_mode
                    ? (is_ct_e && (at_5_e || at_3_e))
                    : ((is_ct_e && at_5_e) || (is_ga_e && at_3_e))) ? 1 : 0;
                e.mm[i].p_damage    = (dmg_xr && profile_.enabled)
                                      ? profile_.p_damage_at(pos_full, L) : 0.0;
            };
            fill(0, cm.mismatch_pos, cm.alt_base, cm.parent_base);
            if (cm.hamming == 2)
                fill(1, cm.mismatch_pos2, cm.alt_base2, cm.parent_base2);
            return e;
        };

        // Cross-bundle recurrence map for interior-transition candidates.
        // Key: mismatch_pos * 4 + alt_base. Value: number of distinct bundle_keys
        // that have an interior-transition H=1 mismatch at (pos, alt).
        // Built once before B2; used in LRT path to veto absorption when ≥2
        // independent fragments (bundles) share the same interior mismatch.
        std::unordered_map<uint32_t, uint32_t> interior_trans_bundle_count;
        uint64_t n_interior_ct = 0;
        uint64_t n_interior_tv = 0;
        if (errcor_.b1_damage_adjust) {
            std::unordered_map<uint32_t, std::unordered_set<uint64_t>> key_to_bundles;
            for (const auto& cm : mismatches) {
                if (cm.hamming != 1) continue;
                const int pid_len = get_layout(static_cast<int>(arena_.length(cm.parent_id))).ilen;
                if (is_terminal_damage_mismatch(cm.parent_base, cm.alt_base,
                                                static_cast<int>(cm.mismatch_pos),
                                                pid_len, profile_.ss_mode))
                    continue;
                const uint32_t xr = cm.alt_base ^ cm.parent_base;
                if (xr == 2u) {
                    ++n_interior_ct;
                    uint32_t key = static_cast<uint32_t>(cm.mismatch_pos) * 4u + cm.alt_base;
                    key_to_bundles[key].insert(bundle_key_of[cm.parent_id]);
                } else {
                    ++n_interior_tv;
                }
            }
            interior_trans_bundle_count.reserve(key_to_bundles.size());
            for (const auto& [k, bset] : key_to_bundles)
                interior_trans_bundle_count[k] = static_cast<uint32_t>(bset.size());
        }
        // κ = interior C↔T rate / interior transversion rate.
        // κ≈1: no excess above PCR error floor → absorption would destroy real variants.
        // κ≥b1_kappa_min: genuine interior deamination signal → absorption is safe.
        double kappa = (n_interior_tv > 0)
            ? static_cast<double>(n_interior_ct) / static_cast<double>(n_interior_tv)
            : (n_interior_ct > 0 ? std::numeric_limits<double>::infinity() : 0.0);
        bool b1_active = errcor_.b1_damage_adjust && (kappa >= errcor_.b1_kappa_min);
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.3f", kappa);
            std::string verdict = !errcor_.b1_damage_adjust ? "" :
                (b1_active ? " → active" : " → suppressed (κ < " +
                    std::to_string(errcor_.b1_kappa_min).substr(0, 4) + ")");
            log_info("Phase 3: interior κ=" + std::string(buf) +
                     " (CT=" + std::to_string(n_interior_ct) +
                     " TV=" + std::to_string(n_interior_tv) + ")" + verdict);
        }

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
        auto joint_dispersion_adj = [](
            const std::vector<std::pair<uint32_t,uint64_t>>& pac,
            uint64_t veto_count,
            double snp_threshold) -> double
        {
            if (pac.empty() || veto_count < 4) return 0.0;
            uint64_t total_w = 0, max_w = 0;
            for (const auto& [k, w] : pac) { total_w += w; max_w = std::max(max_w, w); }
            if (total_w == 0) return 0.0;
            double max_frac = static_cast<double>(max_w) / static_cast<double>(total_w);
            if (max_frac >= snp_threshold) return 0.0;
            if (max_frac > 0.10)
                return -1.5 * (max_frac - 0.10) / (snp_threshold - 0.10);
            if (pac.size() >= 3) return 0.3;
            return 0.0;
        };

        size_t i = 0;
        uint64_t n_joint_adj_pos = 0;
        uint64_t n_joint_adj_neg = 0;
        uint64_t n_parents_processed = 0;
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
            // Key: mismatch_pos * 4 + alt_base. Sort+collapse a small vector
            // (typical fanout ≤ 32) — beats unordered_map allocations + cache misses.
            std::vector<std::pair<uint32_t, uint64_t>> pos_alt_counts;
            pos_alt_counts.reserve(2 * (parent_end - parent_start));
            for (size_t j = parent_start; j < parent_end; ++j) {
                uint64_t cc = id_count[mismatches[j].child_id];
                pos_alt_counts.emplace_back(
                    static_cast<uint32_t>(mismatches[j].mismatch_pos) * 4u
                        + mismatches[j].alt_base,
                    cc);
                if (mismatches[j].hamming == 2) {
                    pos_alt_counts.emplace_back(
                        static_cast<uint32_t>(mismatches[j].mismatch_pos2) * 4u
                            + mismatches[j].alt_base2,
                        cc);
                }
            }
            std::sort(pos_alt_counts.begin(), pos_alt_counts.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            // Collapse runs of equal key in-place.
            size_t pac_w = 0;
            for (size_t r = 0; r < pos_alt_counts.size(); ++r) {
                if (pac_w > 0 && pos_alt_counts[pac_w - 1].first == pos_alt_counts[r].first)
                    pos_alt_counts[pac_w - 1].second += pos_alt_counts[r].second;
                else
                    pos_alt_counts[pac_w++] = pos_alt_counts[r];
            }
            pos_alt_counts.resize(pac_w);
            double joint_adj = joint_dispersion_adj(pos_alt_counts, veto_count, errcor_.snp_threshold);
            ++n_parents_processed;
            if (joint_adj > 0.0)  ++n_joint_adj_pos;
            else if (joint_adj < 0.0) ++n_joint_adj_neg;
            auto find_pac = [&](uint32_t k) -> uint64_t {
                auto it = std::lower_bound(
                    pos_alt_counts.begin(), pos_alt_counts.end(), k,
                    [](const std::pair<uint32_t, uint64_t>& p, uint32_t v) {
                        return p.first < v;
                    });
                return (it != pos_alt_counts.end() && it->first == k) ? it->second : 0ull;
            };

            // For each child, decide: SNP-protected or absorb
            const int pid_ilen = get_layout(static_cast<int>(arena_.length(pid))).ilen;
            for (size_t j = parent_start; j < parent_end; ++j) {
                const ChildMismatch& cm = mismatches[j];
                uint32_t key = static_cast<uint32_t>(cm.mismatch_pos) * 4u + cm.alt_base;
                uint64_t sig = find_pac(key);

                // H=2 path: both mismatches already verified non-damage in B1.
                // SNP veto: protect if either mismatch position has population support.
                if (cm.hamming == 2) {
                    uint32_t key2 = static_cast<uint32_t>(cm.mismatch_pos2) * 4u + cm.alt_base2;
                    uint64_t sig2 = find_pac(key2);

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
                float  fqcl_lr_h2 = std::numeric_limits<float>::quiet_NaN();
                if (errcor_.empirical && !errcor_.legacy_veto) {
                    auto e = build_edge(cm);
                    e.bundle_key = bundle_key_of[eff_pid];
                    e.bundle_occ = bundle_occ_of(eff_pid);
                    double S = emp_model.score(e);
                    // Count-asymmetry prior: PCR error families have
                    // parent_count >> child_count (the error is rare per
                    // amplification, so the descendant lineage is small).
                    // Real biology has parent_count ~ child_count (two
                    // independent molecules drawn from the same coverage).
                    // log((p+1)/(c+1)) shifts the log-odds toward PCR when
                    // the parent dominates. Scale by n_mm: each independent
                    // PCR error compounds the parent-vs-child likelihood
                    // ratio (probability of N errors going the same direction
                    // is the product), so the H=2 hypothesis needs 2× the
                    // count-asymmetry support that H=1 needs to overcome the
                    // 2× per-mismatch penalty.
                    S += static_cast<double>(e.n_mm) *
                         std::log(static_cast<double>(parent_count + 1) /
                                  static_cast<double>(id_count[cm.child_id] + 1));
                    S += joint_adj;
                    absorb_h2 = (S > 0.0);
                    fqcl_lr_h2 = static_cast<float>(S);
                    // H=3 shadow: track near-boundary LR scores before veto logic.
                    constexpr double kH3ShadowNat = 1.0;
                    if      (S >  0.0 && S <  kH3ShadowNat) ++stats.h3_shadow_absorbed_near_boundary;
                    else if (S <= 0.0 && S > -kH3ShadowNat) ++stats.h3_shadow_protected_near_boundary;
                    // Count-based co-occurrence clamp. The ratio test in
                    // snp_veto_h2 fails at hyper-cluster scale (sig must clear
                    // ~snp_threshold * inflated parent_count). When ≥K
                    // independent reads share the same (pos, alt), that is
                    // direct population evidence of a real SNP — block
                    // absorption regardless of S. K = max(snp_min_count, 2):
                    // sig==1 is the child alone; sig≥2 is genuine co-occurrence.
                    // H=2 mismatches were pre-filtered against damage in B1, so
                    // no damage_bypass exemption applies here.
                    uint32_t hard_veto_min_h2 = std::max(errcor_.snp_min_count, 2u);
                    // Singleton high-Q clamp. H=2 mismatches were pre-filtered
                    // against damage in B1, so any high-Q mismatch is variant
                    // evidence regardless of co-occurrence count.
                    bool hi_q_singleton_h2 =
                        (e.mm[0].qual >= errcor_.singleton_qual_min) ||
                        (e.mm[1].qual >= errcor_.singleton_qual_min);
                    if (absorb_h2 && (sig >= hard_veto_min_h2 ||
                                      sig2 >= hard_veto_min_h2 ||
                                      hi_q_singleton_h2))
                        absorb_h2 = false;
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
                            if (!fqcl_path_.empty()) {
                                auto cm_s = cm; cm_s.lr_score = fqcl_lr_h2;
                                fqcl_mismatches_.push_back(cm_s);
                            }
                        }
                    }
                    continue;  // skip the H=1 block below
                }

                // Damage-aware bypass: residual deamination at terminal positions.
                // DS libraries: C→T at 5' end, G→A at 3' end.
                // SS libraries: C→T at both ends (G→A is not a damage pattern in SS).
                // Distinguish C↔T (xr=2, AND=1) from G↔A (xr=2, AND=0).
                const uint32_t dmg_xr_h1 = cm.alt_base ^ cm.parent_base;
                const bool is_ct_h1 = (dmg_xr_h1 == 2u) && ((cm.alt_base & cm.parent_base) == 1u);
                const bool is_ga_h1 = (dmg_xr_h1 == 2u) && ((cm.alt_base & cm.parent_base) == 0u);
                const bool at_5_h1  = cm.mismatch_pos < kDamageEdgeMargin;
                const bool at_3_h1  = cm.mismatch_pos >=
                                          static_cast<uint16_t>(pid_ilen - kDamageEdgeMargin);
                bool damage_bypass = profile_.enabled && (
                    profile_.ss_mode
                    ? (is_ct_h1 && (at_5_h1 || at_3_h1))
                    : ((is_ct_h1 && at_5_h1) || (is_ga_h1 && at_3_h1)));

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

                // Interior transition: C↔T or G↔A at a non-terminal position.
                // Admitted by b1_damage_adjust but requires calibration-free LRT
                // (not the empirical model, which was fit on transversion candidates).
                const bool is_interior_trans_h1 = !damage_bypass && (dmg_xr_h1 == 2u) && b1_active;

                // T5.8: empirical posterior-odds for H=1.
                bool absorb_h1;
                float fqcl_lr_h1 = std::numeric_limits<float>::quiet_NaN();
                if (errcor_.empirical && !errcor_.legacy_veto) {
                    auto e = build_edge(cm);
                    e.bundle_key = bundle_key_of[eff_pid];
                    e.bundle_occ = bundle_occ_of(eff_pid);

                    if (is_interior_trans_h1) {
                        // === Calibration-free LRT path for interior transitions ===
                        uint64_t n_child_count = id_count[cm.child_id];
                        double log_lambda = interior_transition_lrt(
                            parent_count, n_child_count, errcor_.lrt_f0, errcor_.lrt_f1);
                        bool lrt_absorb = (log_lambda > std::log(errcor_.lrt_T));
                        bool cross_strand_veto = false;
                        bool hard_absorb = (parent_count >
                                            static_cast<uint64_t>(errcor_.lrt_hard_ratio * n_child_count))
                                        && (e.mm[0].qual < errcor_.singleton_qual_min);
                        // Wilson CI SNP veto: require sig ≥ 2 (same as empirical hard_veto_min).
                        double wlb = wilson_lower95(sig, veto_count);
                        bool wilson_veto = (sig >= 2 && wlb >= errcor_.snp_cp_lb_threshold);
                        // Cross-bundle recurrence veto: if ≥2 independent bundles show
                        // this same (pos, alt) interior mismatch it is biological, not PCR.
                        auto cb_it = interior_trans_bundle_count.find(key);
                        bool cross_bundle_veto = (cb_it != interior_trans_bundle_count.end()
                                                  && cb_it->second >= 2);
                        // Cross-strand veto: ancient deamination is strand-specific; a
                        // deamination artifact child has reads from one orientation only.
                        // A true SNP child has reads from both fwd and RC orientations.
                        // Require ≥b1_cs_min_total total child reads to avoid false vetoes
                        // where orientation balance is random at very low coverage.
                        {
                            uint64_t cf = id_fwd_count[cm.child_id];
                            uint64_t ct = id_count[cm.child_id];
                            uint64_t cr = (ct >= cf) ? ct - cf : 0;
                            cross_strand_veto = (cf > 0 && cr > 0 &&
                                                 ct >= errcor_.b1_cs_min_total);
                        }
                        bool hi_q = (e.mm[0].qual >= errcor_.singleton_qual_min);
                        absorb_h1 = (hard_absorb || lrt_absorb)
                                    && !wilson_veto && !cross_bundle_veto
                                    && !cross_strand_veto && !hi_q;
                        ++stats.edge_lr_absorbed;  // reuse counter for logging
                    } else {
                    double S = emp_model.score(e);
                    // See H=2 site for rationale on the count-asymmetry term.
                    S += std::log(static_cast<double>(parent_count + 1) /
                                  static_cast<double>(id_count[cm.child_id] + 1));
                    S += joint_adj;
                    absorb_h1 = (S > 0.0);
                    fqcl_lr_h1 = static_cast<float>(S);
                    // Count-based co-occurrence clamp. snp_veto's ratio test
                    // fails at hyper-cluster scale (sig must clear ~snp_threshold
                    // * inflated parent_count). When ≥K independent reads share
                    // the same (pos, alt), that is direct population evidence of
                    // a real SNP — block absorption regardless of S.
                    // K = max(snp_min_count, 2): sig==1 is the child alone;
                    // sig≥2 is genuine co-occurrence. damage_bypass mismatches
                    // at terminal positions remain exempt (residual deamination).
                    uint32_t hard_veto_min = std::max(errcor_.snp_min_count, 2u);
                    // Singleton high-Q clamp. A single high-Q non-damage
                    // mismatch is direct sequence-divergence evidence; the
                    // empirical S>0 path otherwise allows absorption of
                    // singleton variant alleles at low-coverage loci.
                    // damage_bypass remains exempt (residual deamination).
                    bool hi_q_singleton =
                        (e.mm[0].qual >= errcor_.singleton_qual_min);
                    if (absorb_h1 && !damage_bypass &&
                        (sig >= hard_veto_min || hi_q_singleton))
                        absorb_h1 = false;
                    if (snp_veto) { stats.lr_sum_protected += S; ++stats.lr_n_protected; }
                    else          { stats.lr_sum_absorbed  += S; ++stats.lr_n_absorbed;  }
                    if (absorb_h1) ++stats.edge_lr_absorbed;
                    else           ++stats.edge_lr_protected;
                    int ob = fqdup::errcor_emp::occ_bin(e.bundle_occ);
                    if (absorb_h1) ++stats.edge_lr_absorbed_by_occ[ob];
                    else           ++stats.edge_lr_protected_by_occ[ob];
                    }  // end empirical path (non-interior-transition)
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
                        if (!fqcl_path_.empty()) {
                            auto cm_s = cm; cm_s.lr_score = fqcl_lr_h1;
                            fqcl_mismatches_.push_back(cm_s);
                        }
                    }
                }
            }
        }

        // ── Phase B3: damage-aware H>2 merge ───────────────────────────────
        // For heavily damaged ancient DNA (d_max ≥ b3_min_dmax), PCR copies of
        // the same molecule can accumulate 3–5 deamination differences in the
        // near-terminal interior zone (just beyond the mask). Phase 3 pigeonhole
        // only covers H≤2 — this pass fills that gap.
        // Algorithm: K_deam hash buckets (normalize damage-probable positions)
        // → same-length pairs within bucket → all-deam-consistent diff check
        // → count LRT (same threshold as B2).
        if (errcor_.b3_enabled && profile_.enabled && !errcor_.legacy_veto) {
            // Profile-derived gate: count interior positions where P(damage) >
            // b3_deam_threshold using the median read length from the arena.
            // Skip if n_elig < b3_min_n_elig (B3 cannot produce H>=3 absorptions)
            // or sum(P) < b3_min_mass (events too rare to matter).
            int rep_L = 0;
            {
                std::vector<uint16_t> ls(arena_.lengths.begin(), arena_.lengths.end());
                if (!ls.empty()) {
                    std::nth_element(ls.begin(), ls.begin() + ls.size() / 2, ls.end());
                    rep_L = static_cast<int>(ls[ls.size() / 2]);
                }
            }
            int    b3_n_elig = 0;
            double b3_mass   = 0.0;
            if (rep_L > 0) {
                auto [k5, k3] = damage_zone_bounds(rep_L, profile_);
                int interior_end = rep_L - k3;
                for (int p = k5; p < interior_end; ++p) {
                    double q = profile_.p_damage_at(p, rep_L);
                    if (q > static_cast<double>(errcor_.b3_deam_threshold)) {
                        ++b3_n_elig;
                        b3_mass += q;
                    }
                }
            }
            double dmax_avg = 0.5 * (profile_.d_max_5prime + profile_.d_max_3prime);
            if (b3_n_elig < errcor_.b3_min_n_elig ||
                b3_mass   < static_cast<double>(errcor_.b3_min_mass)) {
                log_info("Phase B3: skipped (n_elig=" + std::to_string(b3_n_elig) +
                         " mass=" + std::to_string(b3_mass) +
                         " d_max_avg=" + std::to_string(dmax_avg) + ")");
            } else {
                log_info("Phase B3: running (n_elig=" + std::to_string(b3_n_elig) +
                         " mass=" + std::to_string(b3_mass) +
                         " d_max_avg=" + std::to_string(dmax_avg) + ")");
                // Build K_deam buckets.
                ska::flat_hash_map<uint64_t, std::vector<uint32_t>> b3_buckets;
                b3_buckets.reserve(N / 2);
                {
                    std::vector<uint8_t> scratch;
                    for (uint32_t id = 0; id < N; ++id) {
                        if (is_error_[id]) continue;
                        if (!arena_.is_eligible(id)) continue;
                        int L = arena_.length(id);
                        const auto& lay = get_layout(L);
                        if (lay.ilen < kMinInteriorLen) continue;
                        scratch.resize(static_cast<size_t>((lay.ilen + 3) / 4));
                        uint64_t kd = fqdup::b3::kdamage_hash(
                            arena_.data(id), lay.k5, lay.ilen, L,
                            profile_, errcor_.b3_deam_threshold, scratch.data());
                        // Gate by locus: combine bundle_key (same fragment endpoints)
                        // with kdamage_hash so only reads from the same capture site
                        // are ever compared.
                        uint64_t b3k = XXH3_64bits_withSeed(&kd, sizeof(kd), bundle_key_of[id]);
                        b3_buckets[b3k].push_back(id);
                    }
                }
                // Pairwise check within each bucket.
                std::vector<uint8_t> pi_buf, ci_buf;
                for (auto& [kd, bucket] : b3_buckets) {
                    if (bucket.size() < 2) continue;
                    std::sort(bucket.begin(), bucket.end(),
                              [&](uint32_t a, uint32_t b_) {
                                  return acc_count[a] > acc_count[b_];
                              });
                    for (size_t pidx = 0; pidx < bucket.size(); ++pidx) {
                        uint32_t pid = bucket[pidx];
                        if (is_error_[pid]) continue;
                        uint32_t eff_pid = find_root_chain(pid);
                        int Lp = arena_.length(eff_pid);
                        const auto& layp = get_layout(Lp);
                        if (layp.ilen < kMinInteriorLen) continue;
                        int nbp = (layp.ilen + 3) / 4;
                        pi_buf.resize(static_cast<size_t>(nbp));
                        extract_packed_part(arena_.data(eff_pid), layp.k5, layp.ilen,
                                            pi_buf.data());
                        uint64_t cp = acc_count[eff_pid];
                        for (size_t cidx = pidx + 1; cidx < bucket.size(); ++cidx) {
                            uint32_t cid = bucket[cidx];
                            if (is_error_[cid]) continue;
                            if (arena_.length(cid) != Lp) continue;
                            uint64_t cc = acc_count[cid];
                            if (cc == 0) continue;
                            if (cp > 0 && static_cast<double>(cc) /
                                static_cast<double>(cp) > errcor_.b3_count_ratio)
                                continue;
                            ++stats.b3_candidates;
                            ci_buf.resize(static_cast<size_t>(nbp));
                            extract_packed_part(arena_.data(cid), layp.k5, layp.ilen,
                                                ci_buf.data());
                            fqdup::b3::DiffList diffs;
                            fqdup::b3::packed_diff_list(pi_buf.data(), ci_buf.data(),
                                                         layp.ilen, diffs);
                            if (diffs.count < 3 || diffs.count > errcor_.b3_max_hamming)
                                continue;
                            if (!fqdup::b3::all_diffs_deam_consistent(
                                    diffs, layp.k5, Lp, profile_,
                                    errcor_.b3_deam_threshold))
                                continue;
                            double lrt = fqdup::b3::b3_count_lrt(
                                cp, cc, errcor_.lrt_f0, errcor_.lrt_f1);
                            if (lrt <= std::log(errcor_.lrt_T)) {
                                ++stats.b3_protected;
                                continue;
                            }
                            if (!is_error_[cid]) {
                                is_error_[cid] = true;
                                parent_chain[cid] = eff_pid;
                                acc_count[eff_pid] += cc;
                                acc_count[cid] = 0;
                                if (seq_entry[cid] && seq_entry[eff_pid] &&
                                    seq_entry[cid]->damage_score >
                                    seq_entry[eff_pid]->damage_score) {
                                    seq_entry[eff_pid]->record_index =
                                        seq_entry[cid]->record_index;
                                    seq_entry[eff_pid]->damage_score =
                                        seq_entry[cid]->damage_score;
                                }
                                ++stats.b3_absorbed;
                                ++errcor_absorbed_;
                            }
                        }
                    }
                }
                log_info("Phase B3: candidates=" + std::to_string(stats.b3_candidates) +
                         " absorbed=" + std::to_string(stats.b3_absorbed) +
                         " protected=" + std::to_string(stats.b3_protected));
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
                            // Count-asymmetry term (see B2 H=2 site).
                            S += std::log(static_cast<double>(id_count[parent_id] + 1) /
                                          static_cast<double>(id_count[child_id]  + 1));
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

        // ── T8.7: Phase3RescueIndels ──────────────────────────────────────────
        // Opt-in indel-rescue path. T8 rescue is strictly within-bundle
        // (the bundle filter is a hard correctness constraint, not a
        // scoring preference), so candidate generation iterates bundles
        // directly. For each bundle with ≥2 eligible members:
        //   - direct path (default): pairwise banded-ed over same-bundle
        //     parents with |Δlen|≤1
        //   - syncmer path (only above ~32k directed pairs): bundle-local
        //     ephemeral SyncmerIndex
        // Decisions go to per-thread vectors and are merged in (child_id,
        // score desc, parent_id asc) order, applying at most one
        // absorption per child.
        if (errcor_.rescue_indels && !errcor_.legacy_veto && errcor_.empirical) {
            using fqdup::derep_detail::SyncmerIndex;
            using fqdup::derep_detail::build_syncmer_index_from_pids;
            using fqdup::derep_detail::EditScript;
            using fqdup::derep_detail::banded_edit_distance_le2;
            using fqdup::errcor_emp::EdgeCandidate;
            using fqdup::errcor_emp::IndelScoreParams;

            IndelScoreParams ip;
            ip.alpha_ins  = errcor_.rescue_alpha_ins;
            ip.alpha_del  = errcor_.rescue_alpha_del;
            ip.mask_bonus = errcor_.rescue_mask_bonus;

            auto child_bundle_occ_of = rescue_bundle_occ_of;

            // T8 Step 1 — Pre-group eligible ids by rescue_bundle_key. Bundle
            // membership is the hard correctness gate; iterating bundles
            // directly makes the SyncmerIndex unnecessary for the common case
            // (mean bundle size ~3 → direct banded-ed beats build+query+sort
            // overhead). The pre-grouping pass also (a) accounts the bundle-
            // hot skip into stats.rescue_skip_bundle_hot (previously silent;
            // GPT-5.5 audit), and (b) computes pairs_est used to drive LPT
            // scheduling and the syncmer-path threshold.
            ska::flat_hash_map<uint64_t, std::vector<uint32_t>> bundle_members;
            bundle_members.reserve(N / 4 + 16);
            for (uint32_t id = 0; id < N; ++id) {
                if (!arena_.is_eligible(id)) continue;
                if (is_error_[id]) continue;
                uint32_t occ = static_cast<uint32_t>(child_bundle_occ_of(id));
                if (occ < 2u) continue;
                if (occ >= errcor_.rescue_bundle_hot) {
                    // Fix silent recall loss: bundle-hot skip was uncounted.
                    ++stats.rescue_skip_bundle_hot;
                    continue;
                }
                bundle_members[rescue_bundle_key_of[id]].push_back(id);
            }

            struct BundleTask {
                uint64_t                key;
                std::vector<uint32_t>*  members;     // borrowed
                uint32_t                M;
                uint64_t                pairs_est;   // est. directed pairs |Δlen|≤2
            };
            std::vector<BundleTask> tasks;
            tasks.reserve(bundle_members.size());
            for (auto& kv : bundle_members) {
                if (kv.second.size() < 2) continue;
                std::sort(kv.second.begin(), kv.second.end());
                const auto& bm = kv.second;
                const uint32_t M = static_cast<uint32_t>(bm.size());
                // Length histogram. M is bounded by rescue_bundle_hot (≤50
                // by default) so a flat hashmap is overkill; sort the lens
                // and run-length-encode.
                std::vector<int> lens; lens.reserve(M);
                for (uint32_t i = 0; i < M; ++i) lens.push_back(arena_.length(bm[i]));
                std::sort(lens.begin(), lens.end());
                uint64_t pairs_est = 0;
                size_t i = 0;
                while (i < lens.size()) {
                    int Lv = lens[i];
                    size_t j = i + 1;
                    while (j < lens.size() && lens[j] == Lv) ++j;
                    uint64_t cnt = j - i;
                    // Window neighbours within ±2 (inclusive) of Lv.
                    uint64_t neigh = 0;
                    for (size_t k = 0; k < lens.size(); ++k)
                        if (std::abs(lens[k] - Lv) <= 2) ++neigh;
                    // Exclude self from neighbour count.
                    pairs_est += cnt * (neigh - 1);
                    i = j;
                }
                tasks.push_back(BundleTask{kv.first, &kv.second, M, pairs_est});
            }
            // LPT scheduling: biggest-pair-count bundles start first so tail
            // latency is minimised. Tie-break by key for determinism.
            std::sort(tasks.begin(), tasks.end(),
                      [](const BundleTask& a, const BundleTask& b) {
                          if (a.pairs_est != b.pairs_est) return a.pairs_est > b.pairs_est;
                          return a.key < b.key;
                      });

            struct RescueDecision {
                uint32_t child_id;
                uint32_t parent_id;
                double   score;
                uint8_t  ed;
            };

            unsigned n_threads = errcor_.threads;
            if (n_threads == 0) n_threads = std::max(1u, std::thread::hardware_concurrency());

            struct alignas(64) T8Local {
                std::vector<RescueDecision> dec;
                uint64_t children_examined = 0;
                uint64_t index_queries     = 0;
                uint64_t topk_truncated    = 0;
                uint64_t pairs_banded      = 0;
                uint64_t banded_reject     = 0;
                uint64_t protected_        = 0;
                uint64_t banded_ed[3]      = {0, 0, 0};
                char _pad[8] = {};
            };
            std::vector<T8Local> per_thread(n_threads);
            std::vector<SyncmerIndex::QueryScratch> per_thread_qs(n_threads);

            // Above this directed-pair count per bundle, switch to a bundle-
            // local SyncmerIndex. Below, direct banded-ed is cheaper than
            // build+query overhead. Long-tail bundles trigger the syncmer
            // path; mean-bundle case never does.
            constexpr uint64_t kSyncmerPairThreshold = 32000;

            // T8 Step 2 — atomic dynamic scheduler (replaces static round-robin).
            // Combined with LPT ordering this keeps tail-thread idle time small
            // even when bundle work is heavily skewed.
            std::atomic<size_t> next_bundle{0};

            auto worker = [&](unsigned tid) {
                T8Local& T = per_thread[tid];
                SyncmerIndex::QueryScratch& qs = per_thread_qs[tid];
                // T8 Step 3 / Step 7 — per-thread reusable scratch.
                thread_local std::vector<uint8_t>  decoded;
                thread_local std::vector<uint32_t> off;          // size = M+1
                thread_local std::vector<int>      L;            // per-member length
                thread_local std::vector<uint32_t> order_by_len; // local indices sorted by L
                thread_local std::vector<uint64_t> ch_h;
                thread_local std::vector<uint16_t> ch_p;
                thread_local ska::flat_hash_map<uint32_t, uint32_t> pid_to_local;
                // Step 7 — reusable bundle-local SyncmerIndex (zero allocs after warmup).
                thread_local SyncmerIndex bidx;
                thread_local std::vector<std::pair<uint64_t, uint32_t>> pairs_scratch;

                while (true) {
                    size_t bi = next_bundle.fetch_add(1, std::memory_order_relaxed);
                    if (bi >= tasks.size()) break;
                    const BundleTask& task = tasks[bi];
                    const auto& bm = *task.members;
                    const size_t M = task.M;
                    if (M < 2) continue;

                    // Step 3 — cache lengths and compute decode offsets in one pass.
                    L.assign(M, 0);
                    off.assign(M + 1, 0);
                    for (size_t i = 0; i < M; ++i) {
                        L[i] = arena_.length(bm[i]);
                        off[i + 1] = off[i] + L[i];
                    }
                    if (decoded.size() < off.back()) decoded.resize(off.back());
                    for (size_t i = 0; i < M; ++i)
                        arena_.decode_range(bm[i], 0, L[i], decoded.data() + off[i]);

                    // Local indices sorted by (length asc, parent_id asc) — used
                    // for the 2-pointer / binary-search Δlen window in the
                    // direct path. Determinism via the parent_id tiebreak.
                    order_by_len.resize(M);
                    for (uint32_t i = 0; i < M; ++i) order_by_len[i] = i;
                    std::sort(order_by_len.begin(), order_by_len.end(),
                              [&](uint32_t a, uint32_t b) {
                                  if (L[a] != L[b]) return L[a] < L[b];
                                  return bm[a] < bm[b];
                              });

                    // T8 Step 1 / Step 3 — pairs_est drives the path switch.
                    // Old code recomputed an O(M²) probe here; deleted entirely.
                    const bool over_threshold = (task.pairs_est > kSyncmerPairThreshold);

                    if (over_threshold) {
                        // Step 7 — build into reusable thread_local index.
                        build_syncmer_index_from_decoded(
                            bidx,
                            decoded.data(), off.data(), L.data(),
                            bm.data(), M,
                            pairs_scratch,
                            errcor_.rescue_hash_hot);
                        pid_to_local.clear();
                        if (pid_to_local.bucket_count() < M * 2)
                            pid_to_local.reserve(M * 2);
                        for (size_t i = 0; i < M; ++i)
                            pid_to_local[bm[i]] = static_cast<uint32_t>(i);
                    }

                    // T8 Step 5 — refactored try_pair. Returns through the
                    // best_* refs; the per-child best is pushed to T.dec ONCE
                    // after the parent loop (Step 4 / GPT-5.5 P1.6).
                    auto try_pair = [&](uint32_t child_id, uint32_t parent_id,
                                        int Lc, int Lp,
                                        const uint8_t* child_dec_p,
                                        const uint8_t* parent_dec_p,
                                        int k5_c, int k3_c,
                                        double& best_score, uint32_t& best_parent,
                                        uint8_t& best_ed) {
                        // |Δlen|≤2 already filtered by callers; defensive.
                        if (std::abs(Lp - Lc) > 2) return;

                        EditScript es{};
                        int ed = banded_edit_distance_le2(
                            parent_dec_p, Lp, child_dec_p, Lc, &es);
                        ++T.pairs_banded;
                        if (ed < 0) { ++T.banded_reject; return; }
                        if (ed >= 0 && ed < 3) ++T.banded_ed[ed];

                        EdgeCandidate e{};
                        e.child_id   = child_id;
                        e.parent_id  = parent_id;
                        e.bundle_key = bundle_key_of[parent_id];
                        e.bundle_occ = bundle_occ_of(parent_id);
                        e.child_len  = static_cast<uint16_t>(Lc);
                        e.n_mm       = 0;
                        e.n_ins      = static_cast<uint8_t>(es.n_ins);
                        e.n_del      = static_cast<uint8_t>(es.n_del);
                        e.indel_in_mask = 0;
                        for (int s = 0; s < 2; ++s) {
                            if (es.kind[s] != 1 && es.kind[s] != 2) continue;
                            int p = es.pos[s];
                            if (p < k5_c || p >= Lc - k3_c) {
                                e.indel_in_mask = 1; break;
                            }
                        }
                        for (int s = 0; s < 2 && e.n_mm < 2; ++s) {
                            if (es.kind[s] != 0) continue;
                            int p = es.pos[s];
                            auto& mm = e.mm[e.n_mm];
                            mm.pos        = static_cast<uint16_t>(p);
                            mm.parent_2b  = es.ref_base[s];
                            mm.alt_2b     = es.alt_base[s];
                            mm.qual       = qual_arena_.q_at(child_id, p);
                            mm.damage_chan = 0;
                            mm.p_damage    = 0.0;
                            ++e.n_mm;
                        }

                        double S = emp_model.score_with_indel(e, ip);
                        if (!(S > 0.0)) { ++T.protected_; return; }
                        // Track per-child best (score desc, parent_id asc).
                        if (S > best_score ||
                            (S == best_score && parent_id < best_parent)) {
                            best_score  = S;
                            best_parent = parent_id;
                            best_ed     = static_cast<uint8_t>(ed);
                        }
                    };

                    // id_count gate hoisted out of try_pair (cheap; pre-prunes
                    // banded-ed work).
                    auto id_count_admits = [&](uint32_t child_id, uint32_t parent_id) {
                        if (parent_id == child_id) return false;
                        if (id_count[parent_id] < id_count[child_id]) return false;
                        if (id_count[parent_id] == id_count[child_id] &&
                            parent_id > child_id) return false;
                        return true;
                    };

                    for (size_t ci_local = 0; ci_local < M; ++ci_local) {
                        const uint32_t child_id = bm[ci_local];
                        ++T.children_examined;
                        const int Lc = L[ci_local];
                        const uint8_t* child_dec_p = decoded.data() + off[ci_local];
                        auto [k5_c, k3_c] = damage_zone_bounds(Lc, profile_);

                        double   best_score  = 0.0;
                        uint32_t best_parent = UINT32_MAX;
                        uint8_t  best_ed     = 0;

                        if (!over_threshold) {
                            // T8 Step 4 — direct path. |Δlen|≤2 (was ≤1; the
                            // banded_ed cap supports 2 — silent recall miss
                            // flagged by GPT-5.5). Iterate via order_by_len
                            // window; M ≤ rescue_bundle_hot (≤50 default), so
                            // std::lower_bound is fine.
                            auto lo = std::lower_bound(
                                order_by_len.begin(), order_by_len.end(), Lc - 2,
                                [&](uint32_t a, int v) { return L[a] < v; });
                            auto hi = std::upper_bound(
                                order_by_len.begin(), order_by_len.end(), Lc + 2,
                                [&](int v, uint32_t a) { return v < L[a]; });
                            for (auto it = lo; it != hi; ++it) {
                                const uint32_t pi_local = *it;
                                if (pi_local == ci_local) continue;
                                const uint32_t parent_id = bm[pi_local];
                                if (!id_count_admits(child_id, parent_id)) continue;
                                try_pair(child_id, parent_id,
                                         Lc, L[pi_local],
                                         child_dec_p,
                                         decoded.data() + off[pi_local],
                                         k5_c, k3_c,
                                         best_score, best_parent, best_ed);
                            }
                        } else {
                            fqdup::SyncmerParams sp = fqdup::syncmer_params_for_length(Lc);
                            uint16_t min_hits = errcor_.rescue_min_hits > 0
                                ? static_cast<uint16_t>(errcor_.rescue_min_hits)
                                : (sp.cap >= 16 ? uint16_t{3} : uint16_t{2});
                            fqdup::compute_sketch(child_dec_p, Lc, ch_h, ch_p);
                            if (ch_h.empty()) continue;
                            ++T.index_queries;
                            // T8 Step 6 — filtered query. Pre-tally rejection
                            // by length AND id_count avoids same-bundle
                            // wrong-length parents stealing topk slots from
                            // valid same-length ones (silent recall loss).
                            auto accept = [&](uint32_t pid) -> bool {
                                auto pit = pid_to_local.find(pid);
                                if (pit == pid_to_local.end()) return false;
                                int Lp = L[pit->second];
                                if (std::abs(Lp - Lc) > 2) return false;
                                return id_count_admits(child_id, pid);
                            };
                            const auto& hits = bidx.query_filtered(
                                qs, ch_h.data(), (int)ch_h.size(),
                                min_hits, errcor_.rescue_topk, accept);
                            if (hits.size() == errcor_.rescue_topk) ++T.topk_truncated;
                            for (auto& h : hits) {
                                auto pit = pid_to_local.find(h.parent_id);
                                if (pit == pid_to_local.end()) continue;
                                const uint32_t pi_local = pit->second;
                                try_pair(child_id, h.parent_id,
                                         Lc, L[pi_local],
                                         child_dec_p,
                                         decoded.data() + off[pi_local],
                                         k5_c, k3_c,
                                         best_score, best_parent, best_ed);
                            }
                        }

                        // Push at most ONE decision per child — the global
                        // sort + first-wins below reduces to identity.
                        if (best_parent != UINT32_MAX) {
                            T.dec.push_back(RescueDecision{
                                child_id, best_parent, best_score, best_ed});
                        }
                    }
                }
            };

            if (n_threads == 1) {
                worker(0);
            } else {
                std::vector<std::thread> ts;
                ts.reserve(n_threads - 1);
                for (unsigned t = 1; t < n_threads; ++t) ts.emplace_back(worker, t);
                worker(0);
                for (auto& t : ts) t.join();
            }

            // Fold per-thread sub-counters and decisions into globals.
            std::vector<RescueDecision> decisions;
            for (auto& T : per_thread) {
                stats.rescue_children_examined += T.children_examined;
                stats.rescue_index_queries     += T.index_queries;
                stats.rescue_topk_truncated    += T.topk_truncated;
                stats.rescue_pairs_banded      += T.pairs_banded;
                stats.rescue_banded_reject     += T.banded_reject;
                stats.rescue_protected         += T.protected_;
                for (int e = 0; e < 3; ++e)
                    stats.rescue_banded_ed[e]  += T.banded_ed[e];
                decisions.insert(decisions.end(),
                                 T.dec.begin(), T.dec.end());
            }

            // Deterministic merge: at most one absorption per child, picking
            // the highest-scoring parent (ties broken by parent_id asc).
            std::sort(decisions.begin(), decisions.end(),
                      [](const RescueDecision& a, const RescueDecision& b) {
                          if (a.child_id != b.child_id) return a.child_id < b.child_id;
                          if (a.score    != b.score)    return a.score    > b.score;
                          return a.parent_id < b.parent_id;
                      });
            uint32_t prev_child = UINT32_MAX;
            for (auto& d : decisions) {
                if (d.child_id == prev_child) continue;
                prev_child = d.child_id;
                if (is_error_[d.child_id]) continue;
                uint32_t eff_pid = d.parent_id;
                if (is_error_[eff_pid]) eff_pid = find_root_chain(eff_pid);
                if (eff_pid == d.child_id) continue;
                is_error_[d.child_id] = true;
                parent_chain[d.child_id] = eff_pid;
                acc_count[eff_pid]      += acc_count[d.child_id];
                acc_count[d.child_id]    = 0;
                if (seq_entry[d.child_id] && seq_entry[eff_pid] &&
                    seq_entry[d.child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                    seq_entry[eff_pid]->record_index = seq_entry[d.child_id]->record_index;
                    seq_entry[eff_pid]->damage_score = seq_entry[d.child_id]->damage_score;
                }
                ++stats.absorbed;
                ++stats.rescue_absorbed;
                int ob = static_cast<int>(fqdup::errcor_emp::occ_bin(rescue_bundle_occ_of(eff_pid)));
                if (ob >= 0 && ob < 6) ++stats.rescue_absorbed_by_occ[ob];
                ++errcor_absorbed_;
                // No fqcl record (indel edge — same constraint as adj-len).
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
        log_info("  Joint adj fired: neg=" + std::to_string(n_joint_adj_neg) +
                 " pos=" + std::to_string(n_joint_adj_pos) +
                 " (of " + std::to_string(n_parents_processed) + " parents)");

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
        meta.tool_version  = FQDUP_VERSION;
        meta.input_fastq   = fqcl_input_fastq_;
        meta.n_input_reads = total_reads_;
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

    // Prior counts from a preceding derep_pairs --cluster-format run.
    // Maps XXH3_64(raw_seq) → n_members so Phase 1 seeds the correct count.
    std::string prior_fqcl_path_;
    ska::flat_hash_map<uint64_t, uint32_t> prior_counts_;

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
        << "  --b1-kappa-min FLOAT         Interior CT/TV ratio threshold to enable b1_damage_adjust (default: 2.0)\n"
        << "  --b1-cs-min-total INT        Min child reads to apply cross-strand SNP veto (default: 3)\n"
        << "  --errcor-snp-threshold FLOAT SNP veto: sig/parent_count threshold (default: 0.20)\n"
        << "  --errcor-snp-min-count INT   SNP veto: min absolute sig_count (default: 1)\n"
        << "  --errcor-bucket-cap INT      Max pair-key bucket size (default: 0 = unlimited)\n"
        << "  --errcor-snp-cutoff INT      Low-coverage cutoff (default: 10)\n"
        << "  --errcor-snp-factor FLOAT    Low-coverage SNP multiplier (default: 1.75)\n"
        << "  -t, --threads INT            Worker threads for Phase 3 + compressed I/O (default: 0 = HW, capped at 16 for writer)\n"
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
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
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
        engine.process(in_path, out_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
