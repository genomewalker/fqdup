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
#include <cmath>
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
// Ancient DNA Damage Model
//
// DART-style per-position deamination model:
//   P_CT(p) = d_max_5 * exp(-lambda_5 * p) + background   [C→T at 5' terminus]
//   P_GA(p) = d_max_3 * exp(-lambda_3 * p) + background   [G→A at 3' terminus]
// where p is distance (0-indexed) from the relevant end.
//
// Two PCR-duplicate reads damaged independently have expected Hamming distance:
//   E[k] = Σ_p [ 2*P_CT(p)*(1-P_CT(p)) + 2*P_GA(L-1-p)*(1-P_GA(L-1-p)) ]
//          + 2 * L * pcr_error_rate
//
// We handle damage by masking positions whose damage probability exceeds
// mask_threshold before hashing, so damaged and undamaged copies of the same
// molecule map to the same fingerprint.
// ============================================================================

struct DamageProfile {
    double d_max_5prime   = 0.0;   // Max C→T rate at 5' end
    double d_max_3prime   = 0.0;   // Max G→A rate at 3' end
    double lambda_5prime  = 0.5;   // Exponential decay, 5' end
    double lambda_3prime  = 0.5;   // Exponential decay, 3' end
    double background     = 0.02;  // Background deamination across all positions
    double mask_threshold = 0.05;  // Mask when P(deamination) > this value
    double pcr_error_rate = 0.0;   // n_cycles * phi (per bp per cycle)
    bool   enabled        = false;

    double p_ct(int p) const {
        return d_max_5prime * std::exp(-lambda_5prime * p) + background;
    }
    double p_ga(int p) const {
        return d_max_3prime * std::exp(-lambda_3prime * p) + background;
    }

    // Expected Hamming distance between two independently-damaged duplicate
    // reads of length L (Poisson approximation, see header comment above).
    double expected_mismatches(int L) const {
        double e = 0.0;
        for (int p = 0; p < L; ++p) {
            double pct = p_ct(p);
            double pga = p_ga(L - 1 - p);
            e += 2.0 * pct * (1.0 - pct);
            e += 2.0 * pga * (1.0 - pga);
        }
        e += 2.0 * L * pcr_error_rate;
        return e;
    }

    // 99th-percentile mismatch tolerance under a Poisson approximation.
    int mismatch_tolerance(int L) const {
        double lam = expected_mismatches(L);
        if (lam < 1e-9) return 0;
        double cumP = 0.0, pois = std::exp(-lam);
        int k = 0;
        while (cumP + pois < 0.99 && k < 100) {
            cumP += pois;
            ++k;
            pois *= lam / k;
        }
        return k;
    }

    void print_info(int typical_read_length) const {
        log_info("--- Damage-Aware Deduplication ---");
        log_info("  5'-end d_max:  " + std::to_string(d_max_5prime));
        log_info("  3'-end d_max:  " + std::to_string(d_max_3prime));
        log_info("  5'-end lambda: " + std::to_string(lambda_5prime));
        log_info("  3'-end lambda: " + std::to_string(lambda_3prime));
        log_info("  Background:    " + std::to_string(background));
        log_info("  Mask threshold:" + std::to_string(mask_threshold));
        if (pcr_error_rate > 0.0) {
            log_info("  PCR error rate:" + std::to_string(pcr_error_rate) + " per bp");
        }
        if (typical_read_length > 0) {
            double e   = expected_mismatches(typical_read_length);
            int    tol = mismatch_tolerance(typical_read_length);
            log_info("  Expected mismatches (L=" + std::to_string(typical_read_length) +
                     "): " + std::to_string(e) +
                     ", 99th-pct tolerance: " + std::to_string(tol));
        }
    }
};

// Mask damage-prone positions before hashing.
// C and T at 5'-end positions with P_CT > threshold → sentinel 0x01.
// G and A at 3'-end positions with P_GA > threshold → sentinel 0x02.
// swap_ends=true is used for the reverse-complement strand: the 5' end of the
// RC corresponds to the 3' end of the original molecule, so the 5' and 3'
// damage rates are exchanged.
static std::string apply_damage_mask(const std::string& seq,
                                     const DamageProfile& prof,
                                     bool swap_ends) {
    std::string m = seq;
    int L = static_cast<int>(seq.size());
    for (int i = 0; i < L; ++i) {
        char cu = static_cast<char>(
            std::toupper(static_cast<unsigned char>(m[i])));

        double rate5 = swap_ends ? prof.p_ga(i)         : prof.p_ct(i);
        double rate3 = swap_ends ? prof.p_ct(L - 1 - i) : prof.p_ga(L - 1 - i);

        if (rate5 > prof.mask_threshold && (cu == 'C' || cu == 'T')) {
            m[i] = '\x01';
        } else if (rate3 > prof.mask_threshold && (cu == 'G' || cu == 'A')) {
            m[i] = '\x02';
        }
    }
    return m;
}

// Forward declaration — defined in the Utilities section below.
static inline std::string revcomp(const std::string& s);

// Damage-aware canonical hash: mask terminal damage positions so that reads
// from the same molecule that differ only by deamination hash identically.
// For the reverse-complement strand the 5'/3' damage profiles are swapped
// (the 5' end of the RC corresponds to the 3' overhang of the original).
static uint64_t damage_canonical_hash(const std::string& seq,
                                      const DamageProfile& prof,
                                      bool use_revcomp) {
    std::string masked = apply_damage_mask(seq, prof, /*swap_ends=*/false);
    uint64_t h1 = XXH3_64bits(masked.data(), masked.size());
    if (!use_revcomp) return h1;

    std::string rc        = revcomp(seq);
    std::string masked_rc = apply_damage_mask(rc, prof, /*swap_ends=*/true);
    uint64_t h2 = XXH3_64bits(masked_rc.data(), masked_rc.size());
    return std::min(h1, h2);
}

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
// Damage estimation (DART-inspired, needs FastqReader defined above)
// ============================================================================

// Sample n_reads records uniformly from a (possibly gzipped) FASTQ file,
// accumulate per-position T/(C+T) and A/(A+G) frequencies, estimate a
// background level from positions 20-30, then fit an exponential decay
//   excess(p) = d_max * exp(-lambda * p)
// via ordinary least squares on the log-linearised data.
//
// Uses systematic stride sampling (every stride-th record) so that samples
// are spread throughout the file rather than drawn from the beginning only.
// This matters when the input is sorted — position-sorted reads would make
// the first N records an unrepresentative subsample of one chromosome.
static DamageProfile estimate_damage(const std::string& path,
                                     int    n_reads,
                                     double mask_threshold,
                                     int    stride = 1000) {
    FastqReader reader(path, /*use_pigz=*/false);
    FastqRecord rec;

    constexpr int MAX_POS = 35;
    double T_5[MAX_POS]  = {};  // T  count at position p from 5' end
    double CT_5[MAX_POS] = {};  // C+T count at position p from 5' end
    double A_3[MAX_POS]  = {};  // A  count at position p from 3' end
    double AG_3[MAX_POS] = {};  // A+G count at position p from 3' end

    int      reads_scanned = 0;
    int      typical_len   = 0;
    uint64_t record_pos    = 0;

    // Read up to n_reads * stride total records so that n_reads samples
    // are spread across the first (n_reads * stride) records of the file.
    // For a 400 M-read file with defaults (n_reads=100k, stride=1000) this
    // scans 100 M records — roughly the whole file for typical library sizes.
    const uint64_t max_to_scan = static_cast<uint64_t>(n_reads) * stride;

    while (record_pos < max_to_scan && reader.read(rec)) {
        if (record_pos % stride == 0) {
            int L = static_cast<int>(rec.seq.size());
            if (L > typical_len) typical_len = L;

            for (int p = 0; p < std::min(L, MAX_POS); ++p) {
                char c5 = static_cast<char>(
                    std::toupper(static_cast<unsigned char>(rec.seq[p])));
                if      (c5 == 'T') { T_5[p]++;  CT_5[p]++; }
                else if (c5 == 'C') {             CT_5[p]++; }

                char c3 = static_cast<char>(
                    std::toupper(static_cast<unsigned char>(rec.seq[L - 1 - p])));
                if      (c3 == 'A') { A_3[p]++;  AG_3[p]++; }
                else if (c3 == 'G') {             AG_3[p]++; }
            }
            reads_scanned++;
        }
        record_pos++;
    }

    // Per-position frequencies (-1 when count is too low)
    double freq_T[MAX_POS] = {};
    double freq_A[MAX_POS] = {};
    for (int p = 0; p < MAX_POS; ++p) {
        freq_T[p] = (CT_5[p]  >= 10.0) ? T_5[p]  / CT_5[p]  : -1.0;
        freq_A[p] = (AG_3[p] >= 10.0) ? A_3[p] / AG_3[p] : -1.0;
    }

    // Background: average of positions 20-30 (minimal damage zone)
    double bg_sum = 0.0;
    int    bg_n   = 0;
    for (int p = 20; p < std::min(30, MAX_POS); ++p) {
        if (freq_T[p] >= 0.0) { bg_sum += freq_T[p]; bg_n++; }
        if (freq_A[p] >= 0.0) { bg_sum += freq_A[p]; bg_n++; }
    }
    double background = (bg_n > 0) ? bg_sum / bg_n : 0.02;
    background = std::max(0.005, std::min(0.15, background));

    // OLS fit: excess(p) ≈ d_max * exp(-lambda * p)
    // → log(excess / d_max) = -lambda * p   (simple regression, slope = -lambda)
    auto fit_exp = [&](const double* freq,
                       double& d_max_out, double& lambda_out) {
        d_max_out  = (freq[0] >= 0.0) ? std::max(0.0, freq[0] - background) : 0.0;
        lambda_out = 0.5;
        if (d_max_out < 0.01) return;

        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int n = 0;
        for (int p = 0; p < std::min(20, MAX_POS); ++p) {
            if (freq[p] < 0.0) continue;
            double excess = freq[p] - background;
            if (excess < 0.005) continue;
            double y = std::log(excess / d_max_out);
            sx += p; sy += y; sxx += p * p; sxy += p * y;
            n++;
        }
        if (n < 2) return;
        double denom = n * sxx - sx * sx;
        if (std::abs(denom) < 1e-10) return;
        double slope = (n * sxy - sx * sy) / denom;
        lambda_out = std::max(0.05, std::min(10.0, -slope));
    };

    double d_max_5 = 0.0, lambda_5 = 0.5;
    double d_max_3 = 0.0, lambda_3 = 0.5;
    fit_exp(freq_T, d_max_5, lambda_5);
    fit_exp(freq_A, d_max_3, lambda_3);

    DamageProfile profile;
    profile.d_max_5prime   = d_max_5;
    profile.d_max_3prime   = d_max_3;
    profile.lambda_5prime  = lambda_5;
    profile.lambda_3prime  = lambda_3;
    profile.background     = background;
    profile.mask_threshold = mask_threshold;
    profile.enabled        = (d_max_5 > 0.02 || d_max_3 > 0.02);

    log_info("Pass 0: damage estimation — sampled " +
             std::to_string(reads_scanned) + " reads (every " +
             std::to_string(stride) + "th) from " +
             std::to_string(record_pos) + " scanned in " + path);
    log_info("  5'-end: d_max=" + std::to_string(d_max_5) +
             " lambda=" + std::to_string(lambda_5));
    log_info("  3'-end: d_max=" + std::to_string(d_max_3) +
             " lambda=" + std::to_string(lambda_3));
    log_info("  Background: " + std::to_string(background));
    if (profile.enabled && typical_len > 0) {
        profile.print_info(typical_len);
    } else {
        log_info("  Damage below threshold — standard exact hashing will be used");
    }
    return profile;
}

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
    DerepEngine(bool use_revcomp, bool write_clusters, bool use_pigz, bool use_isal,
                const DamageProfile& profile = DamageProfile{})
        : use_revcomp_(use_revcomp), write_clusters_(write_clusters),
          use_pigz_(use_pigz), use_isal_(use_isal), profile_(profile), total_reads_(0) {}

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

    uint64_t compute_hash(const std::string& seq) const {
        if (profile_.enabled)
            return damage_canonical_hash(seq, profile_, use_revcomp_);
        return canonical_hash(seq, use_revcomp_);
    }

    void pass1_standard(const std::string& ext_path, const std::string& non_path) {
        FastqReader ext_reader(ext_path, use_pigz_);
        FastqReader non_reader(non_path, use_pigz_);

        FastqRecord ext_rec, non_rec;
        uint64_t record_idx = 0;

        while (read_paired_checked(ext_reader, non_reader, ext_rec, non_rec, record_idx, "pass1")) {
            uint64_t h = compute_hash(ext_rec.seq);
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
            uint64_t h = compute_hash(ext_rec.seq);
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
    DamageProfile profile_;

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
              << "\nAncient DNA damage-aware deduplication:\n"
              << "  --damage-auto        Estimate damage parameters from input (Pass 0)\n"
              << "  --damage-stride INT  Sampling stride for Pass 0 (default: 1000)\n"
              << "                       Samples every N-th read to cover sorted files uniformly.\n"
              << "                       Total records scanned = 100000 * stride.\n"
              << "  --damage-dmax  FLOAT Set d_max for both 5' and 3' ends manually\n"
              << "  --damage-dmax5 FLOAT Set d_max for 5' end only\n"
              << "  --damage-dmax3 FLOAT Set d_max for 3' end only\n"
              << "  --damage-lambda FLOAT  Set lambda (decay) for both ends\n"
              << "  --damage-lambda5 FLOAT Set lambda for 5' end only\n"
              << "  --damage-lambda3 FLOAT Set lambda for 3' end only\n"
              << "  --damage-bg FLOAT    Background deamination rate (default: 0.02)\n"
              << "  --mask-threshold FLOAT  Mask positions with P(deamination) > T (default: 0.05)\n"
              << "  --pcr-cycles INT     Number of PCR cycles (for error modeling)\n"
              << "  --pcr-efficiency FLOAT  PCR efficiency per cycle, 0-1 (default: 1.0)\n"
              << "  --pcr-error-rate FLOAT  Error rate in sub/base/doubling (default: 5.3e-7 = Q5/HiFi)\n"
              << "                          Typical values: Q5=5.3e-7, Phusion=3.9e-6, Taq=1.5e-4\n"
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
    bool use_pigz    = false;
    bool use_isal    = false;

    // Damage model options
    bool   damage_auto     = false;
    double damage_dmax5    = -1.0;  // < 0 means "not set"
    double damage_dmax3    = -1.0;
    double damage_lambda5  = 0.5;
    double damage_lambda3  = 0.5;
    double damage_bg       = 0.02;
    double mask_threshold  = 0.05;
    int    damage_stride   = 1000;  // sample every N-th read during damage estimation
    int    pcr_cycles      = 0;
    double pcr_efficiency  = 1.0;     // PCR amplification efficiency per cycle (0–1)
    // Error rate in sub/base/doubling (Potapov & Ong 2017):
    //   Q5/HiFi:  5.3e-7   Phusion/Pfu: 3.9e-6   KOD: 1.2e-5   Taq: 1.5e-4
    // Thermocycling damage adds ~1.4e-6/base/cycle (97% C→T, any polymerase).
    // Default: Q5-class HiFi polymerase — appropriate for well-controlled libraries.
    double pcr_phi         = 5.3e-7;  // sub/base/doubling

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
        } else if (arg == "--damage-auto") {
            damage_auto = true;
        } else if (arg == "--damage-stride" && i + 1 < argc) {
            damage_stride = std::stoi(argv[++i]);
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

    // Build damage profile
    // Effective doublings: D = n * log2(1+E)  (Potapov & Ong 2017, Eq. 6 framework)
    // PCR error contribution per base (for two independent reads in a duplicate):
    //   pcr_error_rate = 2 * phi * D  (factor 2: both reads accumulate errors)
    // We store the per-base-per-read rate; expected_mismatches() multiplies by L and 2.
    double D_eff = pcr_cycles * std::log2(1.0 + pcr_efficiency);
    double pcr_total_rate = pcr_phi * D_eff;  // per base, single read

    DamageProfile profile;
    profile.mask_threshold = mask_threshold;
    profile.pcr_error_rate = pcr_total_rate;

    try {
        if (damage_auto) {
            // Pass 0: estimate from the non-extended file (typically shorter reads,
            // more representative of capture/library prep damage)
            profile = estimate_damage(nonext_path, 100000, mask_threshold, damage_stride);
            profile.pcr_error_rate = pcr_total_rate;
        } else if (damage_dmax5 >= 0.0) {
            // Manual specification
            profile.d_max_5prime  = damage_dmax5;
            profile.d_max_3prime  = (damage_dmax3 >= 0.0) ? damage_dmax3 : damage_dmax5;
            profile.lambda_5prime = damage_lambda5;
            profile.lambda_3prime = damage_lambda3;
            profile.background    = damage_bg;
            profile.pcr_error_rate = pcr_total_rate;
            profile.enabled       = (damage_dmax5 > 0.0);
            profile.print_info(/*typical_read_length=*/0);
        }

        DerepEngine engine(use_revcomp, !cluster_path.empty(), use_pigz, use_isal, profile);
        engine.process(ext_path, nonext_path, out_ext_path, out_non_path, cluster_path);
        log_info("=== Deduplication complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
