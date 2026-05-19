// derep_pairs.cpp
// Paired-end FASTQ deduplication on pre-sorted files.
// No damage modeling or PCR error correction — for those, run 'fqdup derep'
// on the output afterwards.
//
// REQUIRES: Both input files sorted by read ID (use 'fqdup sort' first).
//
// Strategy:
//   Pass 1: Stream both files in lockstep, build hash → position index
//   Pass 2: Stream again, write representative pairs
//
// Representative selection: longest non-extended mate per cluster.
// Memory: ~24 bytes per input read pair.

#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>
#include <zlib.h>

#ifdef __linux__
#include <malloc.h>
#endif

#ifdef __GLIBC__
extern "C" {
    int mallctl(const char *name, void *oldp, size_t *oldlenp,
                void *newp, size_t newlen) __attribute__((weak));
}
#endif

#include "fqdup/cluster_format.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"
#include "fqdup/version.hpp"
#include "derep_detail/damage_keys.hpp"
#include "flat_hash_map.hpp"
#include <xxhash.h>

namespace cf = fqdup::clusterfmt;

static inline std::string_view trim_id_view(std::string_view header) {
    size_t pos = (!header.empty() && header[0] == '@') ? 1 : 0;
    size_t sp = header.find_first_of(" \t", pos);
    return (sp == std::string_view::npos) ? header.substr(pos) : header.substr(pos, sp - pos);
}

namespace {

// ============================================================================
// Index entry — no damage/error state needed here
// ============================================================================

struct IndexEntry {
    uint64_t record_index;  // sequential record number of the representative
    uint64_t count;         // total reads in this cluster
    uint32_t non_len;       // non-extended seq length of the representative

    IndexEntry() : record_index(0), count(1), non_len(0) {}
    IndexEntry(uint64_t idx, uint32_t len)
        : record_index(idx), count(1), non_len(len) {}
};

// ============================================================================
// DerepPairsEngine — two-pass paired deduplication
// ============================================================================

static std::vector<std::uint8_t> pack_2bit(const std::string& seq) {
    static const auto enc = []() {
        std::array<std::uint8_t, 256> t{};
        t['A'] = t['a'] = 0; t['C'] = t['c'] = 1;
        t['G'] = t['g'] = 2; t['T'] = t['t'] = 3;
        return t;
    }();
    std::uint32_t L = static_cast<std::uint32_t>(seq.size());
    std::vector<std::uint8_t> out((L + 3) / 4, 0u);
    for (std::uint32_t i = 0; i < L; ++i)
        out[i >> 2] |= enc[static_cast<unsigned char>(seq[i])] << (6 - 2 * (i & 3));
    return out;
}

class DerepPairsEngine {
public:
    DerepPairsEngine(bool use_revcomp, bool write_clusters,
                     std::string fqcl_path, bool allow_id_mismatch = false,
                     size_t threads = 0)
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          fqcl_path_(std::move(fqcl_path)),
          allow_id_mismatch_(allow_id_mismatch), total_reads_(0) {
        if (threads == 0)
            threads = std::max(1u, std::thread::hardware_concurrency());
        // Two readers run concurrently per pass; split budget evenly, cap at 8.
        decomp_threads_ = std::min(static_cast<size_t>(8), std::max(static_cast<size_t>(1), threads / 2));
        // BGZF output: independent of read threads, cap at 16 (plateaus there).
        write_threads_ = static_cast<int>(std::min(threads, static_cast<size_t>(16)));
    }

    void process(const std::string& ext_path,
                 const std::string& non_path,
                 const std::string& out_ext_path,
                 const std::string& out_non_path,
                 const std::string& cluster_path) {

        log_info("=== Two-pass paired deduplication ===");
        log_info("Pass 1: Build lightweight index");
        log_info("Decompression threads per reader: " + std::to_string(decomp_threads_));
        log_info("Decompression: " + std::string(
#ifdef HAVE_RAPIDGZIP
            "rapidgzip (parallel multi-threaded)"
#elif defined(HAVE_ISAL)
            "ISA-L (hardware-accelerated)"
#else
            "zlib"
#endif
        ));

        pass1(ext_path, non_path);

        size_t index_mb = (index_.size() *
                           sizeof(std::pair<SequenceFingerprint, IndexEntry>)) /
                          (1024 * 1024);
        log_info("Index size: " + std::to_string(index_mb) + " MB for " +
                 std::to_string(total_reads_) + " reads");

        log_info("Pass 2: Write unique records");

        pass2(ext_path, non_path, out_ext_path, out_non_path, cluster_path, fqcl_path_);

        print_stats();
    }

private:
    bool read_paired_checked(FastqReaderBase& ext_reader, FastqReaderBase& non_reader,
                             FastqRecord& ext_rec, FastqRecord& non_rec,
                             uint64_t record_idx, const char* phase) {
        bool ext_ok = ext_reader.read(ext_rec);
        bool non_ok = non_reader.read(non_rec);

        if (ext_ok != non_ok) {
            throw std::runtime_error(
                "Paired FASTQ files have different lengths in " +
                std::string(phase) + " at record " +
                std::to_string(record_idx + 1));
        }
        if (!ext_ok) return false;

        std::string_view ext_id = trim_id_view(ext_rec.header);
        std::string_view non_id = trim_id_view(non_rec.header);
        if (ext_id != non_id) {
            if (allow_id_mismatch_) {
                log_warn("Paired FASTQ ID mismatch in " + std::string(phase) +
                         " at record " + std::to_string(record_idx + 1) +
                         ": ext=" + std::string(ext_id) + ", non=" + std::string(non_id));
            } else {
                throw std::runtime_error(
                    "Paired FASTQ ID mismatch in " + std::string(phase) +
                    " at record " + std::to_string(record_idx + 1) +
                    ": ext=" + std::string(ext_id) + ", non=" + std::string(non_id) +
                    " (ensure both files are sorted identically, or use --allow-id-mismatch)");
            }
        }
        return true;
    }

    void pass1(const std::string& ext_path, const std::string& non_path) {
        auto ext_reader = make_fastq_reader(ext_path, decomp_threads_);
        auto non_reader = make_fastq_reader(non_path, decomp_threads_);
        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;

        while (read_paired_checked(*ext_reader, *non_reader, ext_rec, non_rec,
                                   record_idx, "pass1")) {
            if (ext_rec.seq.size() > rc_scratch_.size())
                rc_scratch_.resize(ext_rec.seq.size());
            XXH128_hash_t h = fqdup::derep_detail::canonical_hash_noalloc(
                ext_rec.seq, use_revcomp_, rc_scratch_.data());
            SequenceFingerprint fp(h, ext_rec.seq.size());

            auto it = index_.find(fp);
            if (it == index_.end()) {
                index_.emplace(fp, IndexEntry(record_idx,
                    static_cast<uint32_t>(non_rec.seq.size())));
            } else {
                it->second.count++;
                if (non_rec.seq.size() > it->second.non_len) {
                    it->second.record_index = record_idx;
                    it->second.non_len = static_cast<uint32_t>(non_rec.seq.size());
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

    void pass2(const std::string& ext_path, const std::string& non_path,
               const std::string& out_ext_path, const std::string& out_non_path,
               const std::string& cluster_path, const std::string& fqcl_path) {
        bool compress_ext = (out_ext_path.size() > 3 &&
                             out_ext_path.substr(out_ext_path.size() - 3) == ".gz");
        bool compress_non = (out_non_path.size() > 3 &&
                             out_non_path.substr(out_non_path.size() - 3) == ".gz");

        FastqWriter ext_writer(out_ext_path, compress_ext, write_threads_);
        FastqWriter non_writer(out_non_path, compress_non, write_threads_);

        gzFile cluster_gz = nullptr;
        if (write_clusters_ && !cluster_path.empty()) {
            cluster_gz = gzopen(cluster_path.c_str(), "wb6");
            if (!cluster_gz)
                throw std::runtime_error("Cannot open cluster file: " + cluster_path);
            gzbuffer(cluster_gz, GZBUF_SIZE);
            if (gzprintf(cluster_gz, "hash\text_len\tpair_count\tnon_len\n") < 0)
                throw std::runtime_error("gzprintf failed while writing cluster header");
        }

        // fqcl writer — optional, opened only when --cluster-format was given.
        std::unique_ptr<cf::Writer> fqcl_writer;
        if (!fqcl_path.empty()) {
            cf::WriterMetadata meta;
            meta.tool         = "fqdup-pairs";
            meta.tool_version = FQDUP_VERSION;
            meta.input_fastq  = non_path;
            meta.n_input_reads = total_reads_;
            fqcl_writer = std::make_unique<cf::Writer>(fqcl_path, std::move(meta));
            fqcl_writer->reserve_clusters(index_.size());
        }
        uint64_t fqcl_cluster_id = 0;

        ska::flat_hash_map<uint64_t, SequenceFingerprint> records_to_write;
        records_to_write.reserve(index_.size());
        for (const auto& [fingerprint, entry] : index_)
            records_to_write[entry.record_index] = fingerprint;
        const size_t n_to_write = records_to_write.size();

        auto ext_reader = make_fastq_reader(ext_path, decomp_threads_);
        auto non_reader = make_fastq_reader(non_path, decomp_threads_);
        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;
        size_t written = 0;

        while (read_paired_checked(*ext_reader, *non_reader, ext_rec, non_rec,
                                   record_idx, "pass2")) {
            auto it = records_to_write.find(record_idx);
            if (it != records_to_write.end()) {
                ext_writer.write(ext_rec);
                non_writer.write(non_rec);

                if (fqcl_writer || cluster_gz) {
                    const SequenceFingerprint& fp = it->second;
                    auto ei = index_.find(fp);
                    if (ei == index_.end())
                        throw std::runtime_error(
                            "Internal error: missing fingerprint for cluster output");
                    const IndexEntry& entry = ei->second;

                    if (fqcl_writer) {
                        cf::ClusterRecord rec;
                        rec.cluster_id    = XXH3_64bits(non_rec.seq.data(), non_rec.seq.size());
                        rec.n_members     = static_cast<uint32_t>(entry.count);
                        rec.n_after_damage = rec.n_members;
                        rec.parent_seq_len = static_cast<uint32_t>(non_rec.seq.size());
                        rec.parent_seq    = pack_2bit(non_rec.seq);
                        rec.flags         = 0;
                        (void)fqcl_cluster_id++;
                        fqcl_writer->write_cluster(rec);
                    }
                    if (cluster_gz &&
                        gzprintf(cluster_gz, "%016lx%016lx\t%lu\t%lu\t%u\n",
                                 static_cast<unsigned long>(fp.hash_hi),
                                 static_cast<unsigned long>(fp.hash_lo),
                                 static_cast<unsigned long>(ext_rec.seq.size()),
                                 static_cast<unsigned long>(entry.count),
                                 entry.non_len) < 0)
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

        if (fqcl_writer) fqcl_writer->close();
        if (cluster_gz && gzclose(cluster_gz) != Z_OK)
            throw std::runtime_error("gzclose failed writing cluster file");
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
    std::string fqcl_path_;
    bool allow_id_mismatch_;
    size_t decomp_threads_;
    int    write_threads_;
    std::vector<char> rc_scratch_;

    ska::flat_hash_map<SequenceFingerprint, IndexEntry, SequenceFingerprintHash> index_;
    uint64_t total_reads_;
};

}  // anonymous namespace

// ============================================================================
// Public entry point
// ============================================================================

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: fqdup " << prog << " [OPTIONS]\n"
        << "\nPaired-end FASTQ deduplication on pre-sorted files.\n"
        << "No damage modeling or error correction — use 'fqdup derep' for those.\n"
        << "\nBoth input files MUST be sorted by read ID (use 'fqdup sort' first).\n"
        << "\nRequired:\n"
        << "  -n FILE      Sorted non-extended FASTQ\n"
        << "  -e FILE      Sorted extended FASTQ\n"
        << "  -o-non FILE  Output non-extended FASTQ\n"
        << "  -o-ext FILE  Output extended FASTQ\n"
        << "\nOptional:\n"
        << "  -c FILE             Write cluster statistics to gzipped TSV\n"
        << "  --cluster-format F  Write cluster genealogy to .fqcl (use with derep --prior-fqcl)\n"
        << "  --no-revcomp        Disable reverse-complement matching (default: enabled)\n"
        << "  --allow-id-mismatch Warn instead of failing on read ID mismatches\n"
        << "  -t N, --threads N   Decompression threads (default: auto, capped at 8/reader)\n"
        << "  -h, --help          Show this help\n"
        << "\nMemory: ~24 bytes per input read pair.\n";
}

int derep_pairs_main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string nonext_path, ext_path, out_ext_path, out_non_path, cluster_path, fqcl_path;
    bool use_revcomp = true;
    bool allow_id_mismatch = false;
    size_t threads = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-n" && i + 1 < argc) {
            nonext_path = argv[++i];
        } else if (arg == "-e" && i + 1 < argc) {
            ext_path = argv[++i];
        } else if (arg == "-o-ext" && i + 1 < argc) {
            out_ext_path = argv[++i];
        } else if (arg == "-o-non" && i + 1 < argc) {
            out_non_path = argv[++i];
        } else if (arg == "-c" && i + 1 < argc) {
            cluster_path = argv[++i];
        } else if (arg == "--cluster-format" && i + 1 < argc) {
            fqcl_path = argv[++i];
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            threads = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--no-revcomp") {
            use_revcomp = false;
        } else if (arg == "--allow-id-mismatch") {
            allow_id_mismatch = true;
        }
    }

    if (nonext_path.empty() || ext_path.empty() ||
        out_ext_path.empty() || out_non_path.empty()) {
        std::cerr << "Error: Missing required arguments\n\n";
        print_usage(argv[0]);
        return 1;
    }

    init_logger("fqdup-derep-pairs.log");
    log_info("=== fqdup derep_pairs: Paired-end two-pass deduplication ===");
    log_info("Extended input (sorted): " + ext_path);
    log_info("Non-extended input (sorted): " + nonext_path);
    log_info("Extended output: " + out_ext_path);
    log_info("Non-extended output: " + out_non_path);
    if (!cluster_path.empty())
        log_info("Cluster output: " + cluster_path);
    log_info("Reverse-complement: " + std::string(use_revcomp ? "enabled" : "disabled"));
    if (allow_id_mismatch)
        log_info("ID mismatch handling: warn (--allow-id-mismatch)");

    try {
        if (!fqcl_path.empty())
            log_info("Cluster fqcl output: " + fqcl_path);
        DerepPairsEngine engine(use_revcomp, !cluster_path.empty(), fqcl_path, allow_id_mismatch, threads);
        engine.process(ext_path, nonext_path, out_ext_path, out_non_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
