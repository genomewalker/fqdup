// derep.cpp
// Ultra-fast, MEMORY-EFFICIENT FASTQ deduplication for PRE-SORTED paired files
// Two-pass algorithm: stores only hashes + positions, not full records
//
// REQUIRES: Both input files MUST be sorted by read ID (use fastq_sort first!)
//
// Strategy:
// Pass 1: Stream files, build hash → position index (16 bytes per record!)
// Pass 2: Stream again, write unique records using index
//
// Memory: ~16 bytes per INPUT record (not per unique cluster!)
// For 400M reads: ~6.4 GB (vs 78 GB for storing full records)
//
// Key optimizations:
// - ISA-L: Hardware-accelerated decompression (4-6× faster)
// - Pigz: Parallel decompression
// - Minimal memory: Only hash + position stored

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <zlib.h>

#ifdef HAVE_ISAL
#include <isa-l/igzip_lib.h>
#endif

#include "fqdup/logger.hpp"
#include "flat_hash_map.hpp"
#include <xxhash.h>

#ifdef __linux__
#include <malloc.h>
#endif

#ifdef __GLIBC__
extern "C" {
    int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen) __attribute__((weak));
}
#endif

// All file-local types are in an anonymous namespace to give their member
// functions internal linkage and avoid ODR violations with sort.cpp (which
// defines identically-named classes with different layouts).
namespace {

// ============================================================================
// Configuration
// ============================================================================

constexpr size_t GZBUF_SIZE = 4 * 1024 * 1024;  // 4MB

// ============================================================================
// Minimal index entry - only 16 bytes per record!
// ============================================================================

struct IndexEntry {
    uint64_t record_index;  // Sequential record number
    uint32_t non_len;       // Non-extended sequence length
    uint32_t count;         // Number of duplicates

    IndexEntry() : record_index(0), non_len(0), count(1) {}
    IndexEntry(uint64_t idx, uint32_t len) : record_index(idx), non_len(len), count(1) {}
};

struct SequenceFingerprint {
    uint64_t hash;
    uint32_t seq_len;

    SequenceFingerprint() : hash(0), seq_len(0) {}
    SequenceFingerprint(uint64_t h, size_t len) : hash(h) {
        if (len > std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("Sequence length exceeds supported fingerprint range");
        }
        seq_len = static_cast<uint32_t>(len);
    }

    bool operator==(const SequenceFingerprint& other) const {
        return hash == other.hash && seq_len == other.seq_len;
    }
};

struct SequenceFingerprintHash {
    size_t operator()(const SequenceFingerprint& fp) const {
        uint64_t mixed = fp.hash ^ (static_cast<uint64_t>(fp.seq_len) * 0x9e3779b97f4a7c15ULL);
        return static_cast<size_t>(mixed ^ (mixed >> 32));
    }
};

// ============================================================================
// FASTQ Record
// ============================================================================

struct FastqRecord {
    std::string header;
    std::string seq;
    std::string plus;
    std::string qual;

    void clear() {
        header.clear();
        seq.clear();
        plus.clear();
        qual.clear();
    }
};

// ============================================================================
// Utilities
// ============================================================================

static inline std::string trim_id(const std::string& header) {
    size_t pos = 0;
    if (!header.empty() && (header[0] == '@' || header[0] == '>'))
        pos = 1;
    size_t sp = header.find_first_of(" \t", pos);
    if (sp == std::string::npos)
        return header.substr(pos);
    return header.substr(pos, sp - pos);
}

static inline std::string revcomp(const std::string& s) {
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        unsigned char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': case 'a': out[i] = (c == 'A') ? 'T' : 't'; break;
            case 'C': case 'c': out[i] = (c == 'C') ? 'G' : 'g'; break;
            case 'G': case 'g': out[i] = (c == 'G') ? 'C' : 'c'; break;
            case 'T': case 't': out[i] = (c == 'T') ? 'A' : 'a'; break;
            default: out[i] = 'N'; break;
        }
    }
    return out;
}

static inline uint64_t canonical_hash(const std::string& seq, bool use_revcomp) {
    uint64_t h1 = XXH3_64bits(seq.data(), seq.size());
    if (!use_revcomp)
        return h1;
    std::string rc = revcomp(seq);
    uint64_t h2 = XXH3_64bits(rc.data(), rc.size());
    return std::min(h1, h2);
}

static std::string shell_escape_path(const std::string& path) {
    std::string escaped;
    escaped.reserve(path.size() + 2);
    escaped.push_back('\'');
    for (char c : path) {
        if (c == '\'') {
            escaped += "'\\''";
        } else {
            escaped.push_back(c);
        }
    }
    escaped.push_back('\'');
    return escaped;
}

// ============================================================================
// FASTQ Reader with ISA-L support
// ============================================================================

#ifdef HAVE_ISAL
class FastqReaderIgzip {
public:
    FastqReaderIgzip(const std::string& path)
        : path_(path), eof_(false), decomp_buffer_pos_(0), decomp_buffer_used_(0),
          record_count_(0) {

        fp_ = fopen(path.c_str(), "rb");
        if (!fp_) throw std::runtime_error("Cannot open: " + path);

        setvbuf(fp_, nullptr, _IOFBF, GZBUF_SIZE);

        input_buffer_ = new uint8_t[GZBUF_SIZE];
        decomp_buffer_ = new uint8_t[GZBUF_SIZE];
        state_ = new inflate_state();

        memset(state_, 0, sizeof(*state_));
        isal_inflate_init(state_);
        state_->crc_flag = ISAL_GZIP;
    }

    ~FastqReaderIgzip() {
        if (fp_) fclose(fp_);
        delete[] input_buffer_;
        delete[] decomp_buffer_;
        delete state_;
    }

    bool read(FastqRecord& rec) {
        if (!getline_igzip(rec.header)) return false;
        if (!getline_igzip(rec.seq)) return false;
        if (!getline_igzip(rec.plus)) return false;
        if (!getline_igzip(rec.qual)) return false;
        record_count_++;
        return true;
    }

    uint64_t record_count() const { return record_count_; }

private:
    bool getline_igzip(std::string& line) {
        line.clear();
        line.reserve(256);

        while (true) {
            if (decomp_buffer_pos_ >= decomp_buffer_used_) {
                if (!refill_decomp_buffer()) {
                    return !line.empty();
                }
            }

            while (decomp_buffer_pos_ < decomp_buffer_used_) {
                char c = decomp_buffer_[decomp_buffer_pos_++];
                if (c == '\n') {
                    return true;
                }
                line.push_back(c);
            }
        }
    }

    bool refill_decomp_buffer() {
        if (eof_ && state_->avail_in == 0) {
            return false;
        }

        if (state_->avail_in == 0 && !eof_) {
            size_t bytes_read = fread(input_buffer_, 1, GZBUF_SIZE, fp_);
            if (bytes_read == 0) {
                eof_ = true;
                if (state_->avail_in == 0) {
                    return false;
                }
            } else {
                state_->next_in = input_buffer_;
                state_->avail_in = bytes_read;
            }
        }

        state_->next_out = decomp_buffer_;
        state_->avail_out = GZBUF_SIZE;

        int ret = isal_inflate(state_);
        if (ret < 0) {
            throw std::runtime_error("ISA-L decompression error: " + std::to_string(ret));
        }

        decomp_buffer_used_ = GZBUF_SIZE - state_->avail_out;
        decomp_buffer_pos_ = 0;

        return decomp_buffer_used_ > 0;
    }

    std::string path_;
    bool eof_;
    FILE* fp_;
    inflate_state* state_;

    uint8_t* input_buffer_;
    uint8_t* decomp_buffer_;
    size_t decomp_buffer_pos_;
    size_t decomp_buffer_used_;
    uint64_t record_count_;
};
#endif

// Standard reader
class FastqReader {
public:
    FastqReader(const std::string& path, bool use_pigz = false)
        : path_(path), use_pigz_(use_pigz), pipe_(nullptr), gzfp_(nullptr), record_count_(0) {

        if (use_pigz_) {
            std::string cmd = "pigz -dc -p 4 -- " + shell_escape_path(path);
            pipe_ = popen(cmd.c_str(), "r");
            if (!pipe_) {
                throw std::runtime_error("Cannot open pigz pipe: " + path);
            }
            setvbuf(pipe_, nullptr, _IOFBF, GZBUF_SIZE);
        } else {
            gzfp_ = gzopen(path.c_str(), "rb");
            if (!gzfp_) {
                throw std::runtime_error("Cannot open file: " + path);
            }
            gzbuffer(gzfp_, GZBUF_SIZE);
        }
    }

    ~FastqReader() {
        if (gzfp_) gzclose(gzfp_);
        if (pipe_) {
            int rc = pclose(pipe_);
            if (rc != 0) {
                std::cerr << "Fatal: pclose failed for pigz reader pipe (" << path_
                          << "), status " << rc << "\n";
                std::terminate();
            }
        }
    }

    bool read(FastqRecord& rec) {
        if (!getline_gz(rec.header)) return false;
        if (!getline_gz(rec.seq)) return false;
        if (!getline_gz(rec.plus)) return false;
        if (!getline_gz(rec.qual)) return false;
        record_count_++;
        return true;
    }

    uint64_t record_count() const { return record_count_; }

private:
    bool getline_gz(std::string& line) {
        line.clear();
        line.reserve(256);
        char buffer[8192];

        if (pipe_) {
            while (true) {
                if (fgets(buffer, sizeof(buffer), pipe_) == nullptr) {
                    return !line.empty();
                }
                size_t len = strlen(buffer);
                if (len > 0 && buffer[len - 1] == '\n') {
                    line.append(buffer, len - 1);
                    return true;
                }
                line.append(buffer, len);
            }
        } else {
            while (true) {
                if (gzgets(gzfp_, buffer, sizeof(buffer)) == nullptr) {
                    return !line.empty();
                }
                size_t len = strlen(buffer);
                if (len > 0 && buffer[len - 1] == '\n') {
                    line.append(buffer, len - 1);
                    return true;
                }
                line.append(buffer, len);
            }
        }
    }

    std::string path_;
    bool use_pigz_;
    FILE* pipe_;
    gzFile gzfp_;
    uint64_t record_count_;
};

// ============================================================================
// FASTQ Writer
// ============================================================================

class FastqWriter {
public:
    FastqWriter(const std::string& path, bool compress)
        : path_(path), compress_(compress), gzfp_(nullptr) {

        if (compress_) {
            // Try to use pigz for parallel compression first (popen -> write)
            // Match reader pattern: use 4 threads for pigz here as well
            std::string pigz_cmd = "pigz -c -p 4 > " + shell_escape_path(path) + " 2>/dev/null";
            pigz_pipe_ = popen(pigz_cmd.c_str(), "w");
            if (pigz_pipe_) {
                // Set a large buffer for pipe writes
                setvbuf(pigz_pipe_, nullptr, _IOFBF, GZBUF_SIZE);
            } else {
                // Fallback to gzFile
                gzfp_ = gzopen(path.c_str(), "wb6");
                if (!gzfp_) {
                    throw std::runtime_error("Cannot open output: " + path);
                }
                gzbuffer(gzfp_, GZBUF_SIZE);
            }
        } else {
            plain_out_.open(path);
            if (!plain_out_.good()) {
                throw std::runtime_error("Cannot open output: " + path);
            }
        }
    }

    ~FastqWriter() {
        if (pigz_pipe_) {
            int rc = pclose(pigz_pipe_);
            if (rc != 0) {
                std::cerr << "Fatal: pclose failed for pigz writer pipe (" << path_
                          << "), status " << rc << "\n";
                std::terminate();
            }
        }
        else if (gzfp_) gzclose(gzfp_);
    }

    void write(const FastqRecord& rec) {
        if (compress_) {
            if (pigz_pipe_) {
                fprintf(pigz_pipe_, "%s\n%s\n%s\n%s\n",
                        rec.header.c_str(), rec.seq.c_str(),
                        rec.plus.c_str(), rec.qual.c_str());
            } else {
                if (gzprintf(gzfp_, "%s\n%s\n%s\n%s\n",
                             rec.header.c_str(), rec.seq.c_str(),
                             rec.plus.c_str(), rec.qual.c_str()) < 0) {
                    throw std::runtime_error("gzprintf failed while writing compressed FASTQ record");
                }
            }
        } else {
            plain_out_ << rec.header << '\n'
                      << rec.seq << '\n'
                      << rec.plus << '\n'
                      << rec.qual << '\n';
        }
    }

private:
    std::string path_;
    bool compress_;
    gzFile gzfp_;
    FILE* pigz_pipe_ = nullptr;
    std::ofstream plain_out_;
};

// ============================================================================
// Deduplication Engine - Two-pass with minimal memory
// ============================================================================

class DerepEngine {
public:
    DerepEngine(bool use_revcomp, bool write_clusters, bool use_pigz, bool use_isal)
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          use_pigz_(use_pigz), use_isal_(use_isal), total_reads_(0) {}

    void process(const std::string& ext_path,
                 const std::string& non_path,
                 const std::string& out_ext_path,
                 const std::string& out_non_path,
                 const std::string& cluster_path) {

        log_info("=== Two-pass deduplication (memory-efficient) ===");
        log_info("Pass 1: Build lightweight index");
        log_info("Decompression: " + std::string(
            use_isal_ ? "ISA-L (hardware-accelerated)" :
            use_pigz_ ? "pigz (parallel)" : "zlib (standard)"));

        // Pass 1: Build index (hash → record index + metadata)
#ifdef HAVE_ISAL
        if (use_isal_) {
            pass1_isal(ext_path, non_path);
        } else
#endif
        {
            pass1_standard(ext_path, non_path);
        }

        size_t index_mb = (index_.size() * sizeof(std::pair<SequenceFingerprint, IndexEntry>)) / (1024 * 1024);
        log_info("Index size: " + std::to_string(index_mb) + " MB for " +
                std::to_string(total_reads_) + " reads");

        log_info("Pass 2: Write unique records");

        // Pass 2: Write unique records
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
    bool read_paired_checked(TExtReader& ext_reader,
                             TNonReader& non_reader,
                             FastqRecord& ext_rec,
                             FastqRecord& non_rec,
                             uint64_t record_idx,
                             const char* phase) {
        bool ext_ok = ext_reader.read(ext_rec);
        bool non_ok = non_reader.read(non_rec);

        if (ext_ok != non_ok) {
            throw std::runtime_error("Paired FASTQ files have different lengths in " +
                                     std::string(phase) + " at record " + std::to_string(record_idx + 1));
        }
        if (!ext_ok) {
            return false;
        }

        const std::string ext_id = trim_id(ext_rec.header);
        const std::string non_id = trim_id(non_rec.header);
        if (ext_id != non_id) {
            log_warn("Paired FASTQ ID mismatch in " + std::string(phase) + " at record " +
                     std::to_string(record_idx + 1) + ": ext=" + ext_id + ", non=" + non_id);
        }
        return true;
    }

    void pass1_standard(const std::string& ext_path, const std::string& non_path) {
        FastqReader ext_reader(ext_path, use_pigz_);
        FastqReader non_reader(non_path, use_pigz_);

        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec, record_idx, "pass1")) {
            uint64_t h = canonical_hash(ext_rec.seq, use_revcomp_);
            SequenceFingerprint fp(h, ext_rec.seq.size());

            auto it = index_.find(fp);
            if (it == index_.end()) {
                // New cluster - store minimal info
                index_.emplace(fp, IndexEntry(record_idx, non_rec.seq.size()));
            } else {
                // Existing cluster - update count and check if this is longer
                it->second.count++;
                if (non_rec.seq.size() > it->second.non_len) {
                    it->second.record_index = record_idx;
                    it->second.non_len = non_rec.seq.size();
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

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec, record_idx, "pass1")) {
            uint64_t h = canonical_hash(ext_rec.seq, use_revcomp_);
            SequenceFingerprint fp(h, ext_rec.seq.size());

            auto it = index_.find(fp);
            if (it == index_.end()) {
                index_.emplace(fp, IndexEntry(record_idx, non_rec.seq.size()));
            } else {
                it->second.count++;
                if (non_rec.seq.size() > it->second.non_len) {
                    it->second.record_index = record_idx;
                    it->second.non_len = non_rec.seq.size();
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

    void pass2_standard(const std::string& ext_path,
                       const std::string& non_path,
                       const std::string& out_ext_path,
                       const std::string& out_non_path,
                       const std::string& cluster_path) {

        FastqReader ext_reader(ext_path, use_pigz_);
        FastqReader non_reader(non_path, use_pigz_);

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
                if (gzprintf(cluster_gz, "hash\text_len\text_count\tnon_len\tnon_count\n") < 0) {
                    throw std::runtime_error("gzprintf failed while writing cluster header");
                }
            }
        }

        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;
        size_t written = 0;

        // Create set of record indices to write
        ska::flat_hash_map<uint64_t, SequenceFingerprint> records_to_write;  // record_idx → fingerprint
        for (const auto& [fingerprint, entry] : index_) {
            records_to_write[entry.record_index] = fingerprint;
        }

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec, record_idx, "pass2")) {
            auto it = records_to_write.find(record_idx);
            if (it != records_to_write.end()) {
                // This is a representative record - write it
                ext_writer.write(ext_rec);
                non_writer.write(non_rec);

                if (cluster_gz) {
                    const SequenceFingerprint& fingerprint = it->second;
                    auto entry_it = index_.find(fingerprint);
                    if (entry_it == index_.end()) {
                        throw std::runtime_error("Internal error: missing fingerprint for cluster output");
                    }
                    const IndexEntry& entry = entry_it->second;
                    if (gzprintf(cluster_gz, "%016lx\t%lu\t%u\t%lu\t%u\n",
                                 fingerprint.hash, ext_rec.seq.size(), entry.count,
                                 non_rec.seq.size(), entry.count) < 0) {
                        throw std::runtime_error("gzprintf failed while writing cluster record");
                    }
                }

                written++;
                if ((written % 100000) == 0) {
                    std::cerr << "\r[Pass 2] " << written << " / " << index_.size()
                             << " unique records written" << std::flush;
                }
            }

            record_idx++;
        }

        if (cluster_gz) gzclose(cluster_gz);

        std::cerr << "\r";
        log_info("Pass 2 complete: " + std::to_string(written) + " unique reads written");
    }

#ifdef HAVE_ISAL
    void pass2_isal(const std::string& ext_path,
                   const std::string& non_path,
                   const std::string& out_ext_path,
                   const std::string& out_non_path,
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
                if (gzprintf(cluster_gz, "hash\text_len\text_count\tnon_len\tnon_count\n") < 0) {
                    throw std::runtime_error("gzprintf failed while writing cluster header");
                }
            }
        }

        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;
        size_t written = 0;

        ska::flat_hash_map<uint64_t, SequenceFingerprint> records_to_write;
        for (const auto& [fingerprint, entry] : index_) {
            records_to_write[entry.record_index] = fingerprint;
        }

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec, record_idx, "pass2")) {
            auto it = records_to_write.find(record_idx);
            if (it != records_to_write.end()) {
                ext_writer.write(ext_rec);
                non_writer.write(non_rec);

                if (cluster_gz) {
                    const SequenceFingerprint& fingerprint = it->second;
                    auto entry_it = index_.find(fingerprint);
                    if (entry_it == index_.end()) {
                        throw std::runtime_error("Internal error: missing fingerprint for cluster output");
                    }
                    const IndexEntry& entry = entry_it->second;
                    if (gzprintf(cluster_gz, "%016lx\t%lu\t%u\t%lu\t%u\n",
                                 fingerprint.hash, ext_rec.seq.size(), entry.count,
                                 non_rec.seq.size(), entry.count) < 0) {
                        throw std::runtime_error("gzprintf failed while writing cluster record");
                    }
                }

                written++;
                if ((written % 100000) == 0) {
                    std::cerr << "\r[Pass 2] " << written << " / " << index_.size()
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
        log_info("Unique clusters: " + std::to_string(index_.size()));

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

    // Minimal index: (hash + ext_len) → (record_index, non_len, count)
    // Secondary fingerprint (ext_len) helps detect hash collisions.
    ska::flat_hash_map<SequenceFingerprint, IndexEntry, SequenceFingerprintHash> index_;

    uint64_t total_reads_;
};

// ============================================================================
// Main
// ============================================================================

}  // namespace

static void print_usage(const char* prog) {
    std::cerr << "Usage: fqdup " << prog << " [OPTIONS]\n"
              << "\n**MEMORY-EFFICIENT VERSION**: Uses two-pass algorithm\n"
              << "  Pass 1: Build lightweight index (~16 bytes per read)\n"
              << "  Pass 2: Write unique records\n"
              << "\n**IMPORTANT**: Both input files MUST be sorted by read ID first!\n"
              << "\nRequired:\n"
              << "  -n FILE     Sorted non-extended FASTQ (.gz required)\n"
              << "  -e FILE     Sorted extended FASTQ (.gz required)\n"
              << "  -o-ext FILE Output extended FASTQ (.gz for compression)\n"
              << "  -o-non FILE Output non-extended FASTQ (.gz for compression)\n"
              << "\nOptional:\n"
              << "  -c FILE     Write cluster info to gzipped TSV\n"
              << "  --no-revcomp Disable reverse-complement matching (default: enabled)\n"
              << "  --pigz      Use pigz for parallel decompression\n"
              << "  --isal      Use ISA-L for hardware-accelerated decompression (FASTEST!)\n"
              << "  -h, --help  Show this help\n"
              << "\nMemory usage:\n"
              << "  ~16 bytes per input read (NOT per unique cluster!)\n"
              << "  400M reads = ~6.4 GB RAM\n"
              << "\nPerformance:\n"
              << "  --isal: 4-6× faster decompression\n"
              << "  --pigz: 2-3× faster decompression\n";
}

int derep_main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string nonext_path, ext_path, out_ext_path, out_non_path, cluster_path;
    bool use_revcomp = true;
    bool use_pigz = false;
    bool use_isal = false;

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
        std::cerr << "Warning: ISA-L requested but not available.\n";
        std::cerr << "Falling back to " << (use_pigz ? "pigz" : "zlib") << ".\n";
        use_isal = false;
    }
#endif

    init_logger("fqdup-derep.log");
    log_info("=== fqdup derep: Memory-efficient two-pass deduplication ===");
    log_info("Extended input (sorted): " + ext_path);
    log_info("Non-extended input (sorted): " + nonext_path);
    log_info("Extended output: " + out_ext_path);
    log_info("Non-extended output: " + out_non_path);
    if (!cluster_path.empty()) {
        log_info("Cluster output: " + cluster_path);
    }
    log_info("Reverse-complement: " + std::string(use_revcomp ? "enabled" : "disabled"));

    try {
        DerepEngine engine(use_revcomp, !cluster_path.empty(), use_pigz, use_isal);
        engine.process(ext_path, nonext_path, out_ext_path, out_non_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
