// fastq_io_backend.cpp
// Implements make_fastq_reader() — the only translation unit that includes
// rapidgzip headers, keeping its non-inline symbols out of all other TUs.

#ifdef HAVE_RAPIDGZIP
#include <rapidgzip/rapidgzip.hpp>
#include <filereader/Standard.hpp>
#endif

#include "fqdup/fastq_common.hpp"
#include <memory>
#include <sys/stat.h>
#include <stdexcept>

namespace {

#ifdef HAVE_RAPIDGZIP
class FastqReaderRapidgzip : public FastqReaderBase {
public:
    explicit FastqReaderRapidgzip(const std::string& path, size_t threads = 0)
        : path_(path), buffer_pos_(0), buffer_used_(0), eof_(false), record_count_(0) {
        buffer_.resize(GZBUF_SIZE);
        reader_ = std::make_unique<rapidgzip::ParallelGzipReader<>>(
            std::make_unique<rapidgzip::StandardFileReader>(path),
            threads,    // 0 = auto-detect from hardware_concurrency
            GZBUF_SIZE
        );
        struct stat st{};
        compressed_size_ = (::stat(path.c_str(), &st) == 0)
                               ? static_cast<uint64_t>(st.st_size) : 0;
    }

    bool read(FastqRecord& rec) override {
        if (!readline(rec.header)) return false;
        if (!readline(rec.seq))
            throw std::runtime_error("Truncated FASTQ: missing sequence line after header '" +
                                     rec.header + "' (record " +
                                     std::to_string(record_count_ + 1) + ")");
        if (!readline(rec.plus))
            throw std::runtime_error("Truncated FASTQ: missing '+' line after sequence in record " +
                                     std::to_string(record_count_ + 1));
        if (!readline(rec.qual))
            throw std::runtime_error("Truncated FASTQ: missing quality line in record " +
                                     std::to_string(record_count_ + 1));
        record_count_++;
        return true;
    }

    uint64_t record_count() const override { return record_count_; }

private:
    bool readline(std::string& line) {
        line.clear();
        while (true) {
            for (size_t i = buffer_pos_; i < buffer_used_; ++i) {
                if (buffer_[i] == '\n') {
                    line.append(buffer_.data() + buffer_pos_, i - buffer_pos_);
                    buffer_pos_ = i + 1;
                    return true;
                }
            }
            if (buffer_pos_ < buffer_used_)
                line.append(buffer_.data() + buffer_pos_, buffer_used_ - buffer_pos_);
            buffer_pos_ = 0;
            buffer_used_ = 0;
            if (eof_) return !line.empty();
            const size_t n = reader_->read(buffer_.data(), buffer_.size());
            if (n == 0) {
                eof_ = true;
                // Detect truncated gzip: if rapidgzip stops decoding before the
                // last compressed byte, the underlying stream lacked a proper
                // gzip footer (CRC32+ISIZE) — i.e. it was truncated mid-stream.
                // tellCompressed() returns bits; compare with file size in bits.
                if (compressed_size_ > 0) {
                    const uint64_t bits_consumed = reader_->tellCompressed();
                    const uint64_t bits_total    = compressed_size_ * 8ULL;
                    // Allow up to 7 bits of padding within the last byte.
                    if (bits_consumed + 7 < bits_total)
                        throw std::runtime_error(
                            "Truncated gzip input '" + path_ +
                            "': decoded " + std::to_string(bits_consumed / 8) +
                            "/" + std::to_string(compressed_size_) +
                            " bytes (missing gzip footer)");
                }
                return !line.empty();
            }
            buffer_used_ = n;
        }
    }

    std::unique_ptr<rapidgzip::ParallelGzipReader<>> reader_;
    std::string       path_;
    std::vector<char> buffer_;
    size_t            buffer_pos_;
    size_t            buffer_used_;
    bool              eof_;
    uint64_t          record_count_;
    uint64_t          compressed_size_ = 0;
};
#endif  // HAVE_RAPIDGZIP

}  // anonymous namespace

static bool is_gzip(const std::string& path) {
    if (path == "/dev/stdin" || path == "-") return false;  // stdin: FastqReader uses gzdopen(fileno(stdin))
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    unsigned char magic[2] = {0, 0};
    bool gz = (fread(magic, 1, 2, f) == 2) &&
              (magic[0] == 0x1f && magic[1] == 0x8b);
    fclose(f);
    return gz;
}

std::unique_ptr<FastqReaderBase> make_fastq_reader(const std::string& path,
                                                    size_t threads) {
    if (is_gzip(path)) {
#ifdef HAVE_RAPIDGZIP
        return std::make_unique<FastqReaderRapidgzip>(path, threads);
#elif defined(HAVE_ISAL)
        return std::make_unique<FastqReaderIgzip>(path);
#endif
    }
    // Plain text (or non-gz): zlib gzopen handles both transparently
    return std::make_unique<FastqReader>(path);
}
