#pragma once
// Shared FASTQ I/O primitives for derep_pairs.cpp and derep.cpp.
//
// FastqRecord, FastqReaderBase, and make_fastq_reader() are declared in
// fastq_types.hpp (global scope, ODR-safe, includeable from any TU).
//
// Reader/writer implementations with inline method bodies are in an anonymous
// namespace below (internal linkage per TU), avoiding ODR violations.

#include "fqdup/fastq_types.hpp"
#include "fqdup/fastq_writer.hpp"  // shared FastqWriter (bgzf / isa-l / zlib)

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <zlib.h>

#ifdef HAVE_ISAL
#include <isa-l/igzip_lib.h>
#endif

#ifdef HAVE_BGZF
#include <htslib/bgzf.h>
#endif

#include <xxhash.h>

namespace {

constexpr size_t GZBUF_SIZE = 4 * 1024 * 1024;  // 4 MB

// ============================================================================
// Sequence fingerprint — composite key for the deduplication index
//
// Uses XXH3_128bits (two 64-bit words) to make hash collisions negligible.
// At 100M reads the birthday-paradox collision probability is ~3×10^-24 with
// a 128-bit hash vs ~3×10^-10 with 64-bit. seq_len is retained as an extra
// discriminator and to avoid recomputing length on lookup.
// ============================================================================

struct SequenceFingerprint {
    uint64_t hash_lo;   // XXH3_128bits low 64 bits
    uint64_t hash_hi;   // XXH3_128bits high 64 bits
    uint32_t seq_len;

    SequenceFingerprint() : hash_lo(0), hash_hi(0), seq_len(0) {}
    SequenceFingerprint(XXH128_hash_t h128, size_t len)
        : hash_lo(h128.low64), hash_hi(h128.high64) {
        if (len > std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("Sequence length exceeds supported fingerprint range");
        seq_len = static_cast<uint32_t>(len);
    }

    bool operator==(const SequenceFingerprint& other) const {
        return hash_lo == other.hash_lo && hash_hi == other.hash_hi &&
               seq_len == other.seq_len;
    }
};

struct SequenceFingerprintHash {
    size_t operator()(const SequenceFingerprint& fp) const {
        // Combine the two halves with a mixing step
        uint64_t mixed = fp.hash_lo ^ (fp.hash_hi * 0x9e3779b97f4a7c15ULL)
                                    ^ (static_cast<uint64_t>(fp.seq_len) * 0x517cc1b727220a95ULL);
        return static_cast<size_t>(mixed ^ (mixed >> 32));
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

static inline XXH128_hash_t canonical_hash(const std::string& seq, bool use_revcomp) {
    XXH128_hash_t h1 = XXH3_128bits(seq.data(), seq.size());
    if (!use_revcomp)
        return h1;
    std::string rc = revcomp(seq);
    XXH128_hash_t h2 = XXH3_128bits(rc.data(), rc.size());
    // Canonical = lexicographically smaller of the two 128-bit hashes
    if (h1.high64 < h2.high64 || (h1.high64 == h2.high64 && h1.low64 <= h2.low64))
        return h1;
    return h2;
}

// ============================================================================
// FASTQ Reader — ISA-L hardware-accelerated decompression (optional)
// ============================================================================

#ifdef HAVE_ISAL
class FastqReaderIgzip : public FastqReaderBase {
public:
    FastqReaderIgzip(const std::string& path)
        : path_(path), eof_(false), decomp_buffer_pos_(0), decomp_buffer_used_(0),
          record_count_(0) {

        fp_ = fopen(path.c_str(), "rb");
        if (!fp_) throw std::runtime_error("Cannot open: " + path);

        setvbuf(fp_, nullptr, _IOFBF, GZBUF_SIZE);

        input_buffer_  = new uint8_t[GZBUF_SIZE];
        decomp_buffer_ = new uint8_t[GZBUF_SIZE];
        state_         = new inflate_state();

        memset(state_, 0, sizeof(*state_));
        isal_inflate_init(state_);
        state_->crc_flag = ISAL_GZIP;
    }

    ~FastqReaderIgzip() override {
        if (fp_) fclose(fp_);
        delete[] input_buffer_;
        delete[] decomp_buffer_;
        delete state_;
    }

    bool read(FastqRecord& rec) override {
        if (!getline_igzip(rec.header)) return false;
        if (!getline_igzip(rec.seq))
            throw std::runtime_error("Truncated FASTQ: missing sequence line after header '" +
                                     rec.header + "' (record " +
                                     std::to_string(record_count_ + 1) + ")");
        if (!getline_igzip(rec.plus))
            throw std::runtime_error("Truncated FASTQ: missing '+' line after sequence in record " +
                                     std::to_string(record_count_ + 1));
        if (!getline_igzip(rec.qual))
            throw std::runtime_error("Truncated FASTQ: missing quality line in record " +
                                     std::to_string(record_count_ + 1));
        record_count_++;
        return true;
    }

    uint64_t record_count() const override { return record_count_; }

private:
    bool getline_igzip(std::string& line) {
        line.clear();
        line.reserve(256);

        while (true) {
            if (decomp_buffer_pos_ >= decomp_buffer_used_) {
                if (!refill_decomp_buffer())
                    return !line.empty();
            }
            while (decomp_buffer_pos_ < decomp_buffer_used_) {
                char c = static_cast<char>(decomp_buffer_[decomp_buffer_pos_++]);
                if (c == '\n') return true;
                line.push_back(c);
            }
        }
    }

    bool refill_decomp_buffer() {
        if (eof_ && state_->avail_in == 0) return false;

        if (state_->avail_in == 0 && !eof_) {
            size_t bytes_read = fread(input_buffer_, 1, GZBUF_SIZE, fp_);
            if (bytes_read == 0) {
                eof_ = true;
                if (state_->avail_in == 0) return false;
            } else {
                state_->next_in  = input_buffer_;
                state_->avail_in = static_cast<uint32_t>(bytes_read);
            }
        }

        state_->next_out  = decomp_buffer_;
        state_->avail_out = GZBUF_SIZE;

        int ret = isal_inflate(state_);
        if (ret < 0)
            throw std::runtime_error("ISA-L decompression error: " + std::to_string(ret));

        decomp_buffer_used_ = GZBUF_SIZE - state_->avail_out;
        decomp_buffer_pos_  = 0;
        return decomp_buffer_used_ > 0;
    }

    std::string    path_;
    bool           eof_;
    FILE*          fp_;
    inflate_state* state_;
    uint8_t*       input_buffer_;
    uint8_t*       decomp_buffer_;
    size_t         decomp_buffer_pos_;
    size_t         decomp_buffer_used_;
    uint64_t       record_count_;
};
#endif  // HAVE_ISAL

// ============================================================================
// FASTQ Reader — standard (zlib)
// ============================================================================

class FastqReader : public FastqReaderBase {
public:
    explicit FastqReader(const std::string& path)
        : path_(path), gzfp_(nullptr), record_count_(0) {

        if (path == "-" || path == "/dev/stdin")
            gzfp_ = gzdopen(fileno(stdin), "rb");
        else
            gzfp_ = gzopen(path.c_str(), "rb");
        if (!gzfp_)
            throw std::runtime_error("Cannot open file: " + path);
        gzbuffer(gzfp_, GZBUF_SIZE);
    }

    ~FastqReader() override {
        if (gzfp_) gzclose(gzfp_);
    }

    bool read(FastqRecord& rec) override {
        if (!getline_gz(rec.header)) return false;
        if (!getline_gz(rec.seq))
            throw std::runtime_error("Truncated FASTQ: missing sequence line after header '" +
                                     rec.header + "' (record " +
                                     std::to_string(record_count_ + 1) + ")");
        if (!getline_gz(rec.plus))
            throw std::runtime_error("Truncated FASTQ: missing '+' line after sequence in record " +
                                     std::to_string(record_count_ + 1));
        if (!getline_gz(rec.qual))
            throw std::runtime_error("Truncated FASTQ: missing quality line in record " +
                                     std::to_string(record_count_ + 1));
        record_count_++;
        return true;
    }

    uint64_t record_count() const override { return record_count_; }

private:
    bool getline_gz(std::string& line) {
        line.clear();
        line.reserve(256);
        char buffer[8192];

        while (true) {
            if (gzgets(gzfp_, buffer, sizeof(buffer)) == nullptr)
                return !line.empty();
            size_t len = strlen(buffer);
            if (len > 0 && buffer[len - 1] == '\n') {
                line.append(buffer, len - 1);
                return true;
            }
            line.append(buffer, len);
        }
    }

    std::string path_;
    gzFile      gzfp_;
    uint64_t    record_count_;
};

// FastqWriter now lives in fqdup/fastq_writer.hpp (included above) so sort and
// every other subcommand share one implementation.

}  // anonymous namespace
