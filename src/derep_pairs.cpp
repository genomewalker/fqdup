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

#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"
#include "flat_hash_map.hpp"

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

class DerepPairsEngine {
public:
    DerepPairsEngine(bool use_revcomp, bool write_clusters,
                     bool use_pigz, bool use_isal)
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          use_pigz_(use_pigz), use_isal_(use_isal), total_reads_(0) {}

    void process(const std::string& ext_path,
                 const std::string& non_path,
                 const std::string& out_ext_path,
                 const std::string& out_non_path,
                 const std::string& cluster_path) {

        log_info("=== Two-pass paired deduplication ===");
        log_info("Pass 1: Build lightweight index");
        log_info("Decompression: " + std::string(
            use_isal_ ? "ISA-L (hardware-accelerated)" :
            use_pigz_ ? "pigz (parallel)" : "zlib (standard)"));

#ifdef HAVE_ISAL
        if (use_isal_) {
            pass1_isal(ext_path, non_path);
        } else
#endif
        {
            pass1_standard(ext_path, non_path);
        }

        size_t index_mb = (index_.size() *
                           sizeof(std::pair<SequenceFingerprint, IndexEntry>)) /
                          (1024 * 1024);
        log_info("Index size: " + std::to_string(index_mb) + " MB for " +
                 std::to_string(total_reads_) + " reads");

        log_info("Pass 2: Write unique records");

#ifdef HAVE_ISAL
        if (use_isal_) {
            pass2_isal(ext_path, non_path, out_ext_path, out_non_path, cluster_path);
        } else
#endif
        {
            pass2_standard(ext_path, non_path, out_ext_path, out_non_path, cluster_path);
        }

        print_stats();
    }

private:
    template <typename TExtReader, typename TNonReader>
    bool read_paired_checked(TExtReader& ext_reader, TNonReader& non_reader,
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

        const std::string ext_id = trim_id(ext_rec.header);
        const std::string non_id = trim_id(non_rec.header);
        if (ext_id != non_id) {
            log_warn("Paired FASTQ ID mismatch in " + std::string(phase) +
                     " at record " + std::to_string(record_idx + 1) +
                     ": ext=" + ext_id + ", non=" + non_id);
        }
        return true;
    }

    void pass1_standard(const std::string& ext_path, const std::string& non_path) {
        FastqReader ext_reader(ext_path, use_pigz_);
        FastqReader non_reader(non_path, use_pigz_);
        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec,
                                   record_idx, "pass1")) {
            uint64_t h = canonical_hash(ext_rec.seq, use_revcomp_);
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

#ifdef HAVE_ISAL
    void pass1_isal(const std::string& ext_path, const std::string& non_path) {
        FastqReaderIgzip ext_reader(ext_path);
        FastqReaderIgzip non_reader(non_path);
        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec,
                                   record_idx, "pass1")) {
            uint64_t h = canonical_hash(ext_rec.seq, use_revcomp_);
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
#endif

    void write_pass2(FastqReader& ext_reader, FastqReader& non_reader,
                     const std::string& out_ext_path, const std::string& out_non_path,
                     const std::string& cluster_path) {
        bool compress_ext = (out_ext_path.size() > 3 &&
                             out_ext_path.substr(out_ext_path.size() - 3) == ".gz");
        bool compress_non = (out_non_path.size() > 3 &&
                             out_non_path.substr(out_non_path.size() - 3) == ".gz");

        FastqWriter ext_writer(out_ext_path, compress_ext);
        FastqWriter non_writer(out_non_path, compress_non);

        gzFile cluster_gz = nullptr;
        if (write_clusters_ && !cluster_path.empty()) {
            cluster_gz = gzopen(cluster_path.c_str(), "wb6");
            if (cluster_gz) {
                gzbuffer(cluster_gz, GZBUF_SIZE);
                if (gzprintf(cluster_gz,
                             "hash\text_len\text_count\tnon_len\tnon_count\n") < 0)
                    throw std::runtime_error(
                        "gzprintf failed while writing cluster header");
            }
        }

        ska::flat_hash_map<uint64_t, SequenceFingerprint> records_to_write;
        records_to_write.reserve(index_.size());
        for (const auto& [fingerprint, entry] : index_)
            records_to_write[entry.record_index] = fingerprint;
        const size_t n_to_write = records_to_write.size();

        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;
        size_t written = 0;

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec,
                                   record_idx, "pass2")) {
            auto it = records_to_write.find(record_idx);
            if (it != records_to_write.end()) {
                ext_writer.write(ext_rec);
                non_writer.write(non_rec);

                if (cluster_gz) {
                    const SequenceFingerprint& fp = it->second;
                    auto ei = index_.find(fp);
                    if (ei == index_.end())
                        throw std::runtime_error(
                            "Internal error: missing fingerprint for cluster output");
                    const IndexEntry& entry = ei->second;
                    if (gzprintf(cluster_gz, "%016lx\t%lu\t%lu\t%u\t%lu\n",
                                 static_cast<unsigned long>(fp.hash),
                                 static_cast<unsigned long>(ext_rec.seq.size()),
                                 static_cast<unsigned long>(entry.count),
                                 entry.non_len,
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

    void pass2_standard(const std::string& ext_path, const std::string& non_path,
                        const std::string& out_ext_path, const std::string& out_non_path,
                        const std::string& cluster_path) {
        FastqReader ext_reader(ext_path, use_pigz_);
        FastqReader non_reader(non_path, use_pigz_);
        write_pass2(ext_reader, non_reader, out_ext_path, out_non_path, cluster_path);
    }

#ifdef HAVE_ISAL
    void pass2_isal(const std::string& ext_path, const std::string& non_path,
                    const std::string& out_ext_path, const std::string& out_non_path,
                    const std::string& cluster_path) {
        FastqReaderIgzip ext_reader(ext_path);
        FastqReaderIgzip non_reader(non_path);

        bool compress_ext = (out_ext_path.size() > 3 &&
                             out_ext_path.substr(out_ext_path.size() - 3) == ".gz");
        bool compress_non = (out_non_path.size() > 3 &&
                             out_non_path.substr(out_non_path.size() - 3) == ".gz");

        FastqWriter ext_writer(out_ext_path, compress_ext);
        FastqWriter non_writer(out_non_path, compress_non);

        gzFile cluster_gz = nullptr;
        if (write_clusters_ && !cluster_path.empty()) {
            cluster_gz = gzopen(cluster_path.c_str(), "wb6");
            if (cluster_gz) {
                gzbuffer(cluster_gz, GZBUF_SIZE);
                if (gzprintf(cluster_gz,
                             "hash\text_len\text_count\tnon_len\tnon_count\n") < 0)
                    throw std::runtime_error(
                        "gzprintf failed while writing cluster header");
            }
        }

        ska::flat_hash_map<uint64_t, SequenceFingerprint> records_to_write;
        records_to_write.reserve(index_.size());
        for (const auto& [fingerprint, entry] : index_)
            records_to_write[entry.record_index] = fingerprint;
        const size_t n_to_write = records_to_write.size();

        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;
        size_t written = 0;

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec,
                                   record_idx, "pass2")) {
            auto it = records_to_write.find(record_idx);
            if (it != records_to_write.end()) {
                ext_writer.write(ext_rec);
                non_writer.write(non_rec);

                if (cluster_gz) {
                    const SequenceFingerprint& fp = it->second;
                    auto ei = index_.find(fp);
                    if (ei == index_.end())
                        throw std::runtime_error(
                            "Internal error: missing fingerprint for cluster output");
                    const IndexEntry& entry = ei->second;
                    if (gzprintf(cluster_gz, "%016lx\t%lu\t%lu\t%u\t%lu\n",
                                 static_cast<unsigned long>(fp.hash),
                                 static_cast<unsigned long>(ext_rec.seq.size()),
                                 static_cast<unsigned long>(entry.count),
                                 entry.non_len,
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
#endif

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
    bool use_pigz_;
    bool use_isal_;

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
        << "  -c FILE      Write cluster statistics to gzipped TSV\n"
        << "  --no-revcomp Disable reverse-complement matching (default: enabled)\n"
        << "  --pigz       Use pigz for parallel decompression\n"
        << "  --isal       Use ISA-L for hardware-accelerated decompression\n"
        << "  -h, --help   Show this help\n"
        << "\nMemory: ~24 bytes per input read pair.\n";
}

int derep_pairs_main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string nonext_path, ext_path, out_ext_path, out_non_path, cluster_path;
    bool use_revcomp = true;
    bool use_pigz    = false;
    bool use_isal    = false;

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
        } else if (arg == "--no-revcomp") {
            use_revcomp = false;
        } else if (arg == "--pigz") {
            use_pigz = true;
        } else if (arg == "--isal") {
            use_isal = true;
        }
    }

    if (nonext_path.empty() || ext_path.empty() ||
        out_ext_path.empty() || out_non_path.empty()) {
        std::cerr << "Error: Missing required arguments\n\n";
        print_usage(argv[0]);
        return 1;
    }

#ifndef HAVE_ISAL
    if (use_isal) {
        std::cerr << "Warning: ISA-L requested but not available.\n"
                  << "Falling back to " << (use_pigz ? "pigz" : "zlib") << ".\n";
        use_isal = false;
    }
#endif

    init_logger("fqdup-derep-pairs.log");
    log_info("=== fqdup derep_pairs: Paired-end two-pass deduplication ===");
    log_info("Extended input (sorted): " + ext_path);
    log_info("Non-extended input (sorted): " + nonext_path);
    log_info("Extended output: " + out_ext_path);
    log_info("Non-extended output: " + out_non_path);
    if (!cluster_path.empty())
        log_info("Cluster output: " + cluster_path);
    log_info("Reverse-complement: " + std::string(use_revcomp ? "enabled" : "disabled"));

    try {
        DerepPairsEngine engine(use_revcomp, !cluster_path.empty(), use_pigz, use_isal);
        engine.process(ext_path, nonext_path, out_ext_path, out_non_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
