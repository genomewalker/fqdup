// fastq_io_backend.cpp
// Implements make_fastq_reader() — the only translation unit that includes
// rapidgzip headers, keeping its non-inline symbols out of all other TUs.

#ifdef HAVE_RAPIDGZIP
#include <rapidgzip/rapidgzip.hpp>
#include <filereader/Standard.hpp>
#endif

#include "fqdup/fastq_common.hpp"
#include <memory>

namespace {

#ifdef HAVE_RAPIDGZIP
class FastqReaderRapidgzip : public FastqReaderBase {
public:
    explicit FastqReaderRapidgzip(const std::string& path, size_t threads = 0)
        : buffer_pos_(0), buffer_used_(0), eof_(false), record_count_(0) {
        buffer_.resize(GZBUF_SIZE);
        reader_ = std::make_unique<rapidgzip::ParallelGzipReader<>>(
            std::make_unique<rapidgzip::StandardFileReader>(path),
            threads,    // 0 = auto-detect from hardware_concurrency
            GZBUF_SIZE
        );
    }

    bool read(FastqRecord& rec) override {
        if (!readline(rec.header)) return false;
        if (!readline(rec.seq))    return false;
        if (!readline(rec.plus))   return false;
        if (!readline(rec.qual))   return false;
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
            if (n == 0) { eof_ = true; return !line.empty(); }
            buffer_used_ = n;
        }
    }

    std::unique_ptr<rapidgzip::ParallelGzipReader<>> reader_;
    std::vector<char> buffer_;
    size_t            buffer_pos_;
    size_t            buffer_used_;
    bool              eof_;
    uint64_t          record_count_;
};
#endif  // HAVE_RAPIDGZIP

}  // anonymous namespace

static bool is_gzip(const std::string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    unsigned char magic[2] = {0, 0};
    bool gz = (fread(magic, 1, 2, f) == 2) &&
              (magic[0] == 0x1f && magic[1] == 0x8b);
    fclose(f);
    return gz;
}

std::unique_ptr<FastqReaderBase> make_fastq_reader(const std::string& path) {
    if (is_gzip(path)) {
#ifdef HAVE_RAPIDGZIP
        return std::make_unique<FastqReaderRapidgzip>(path);
#elif defined(HAVE_ISAL)
        return std::make_unique<FastqReaderIgzip>(path);
#endif
    }
    // Plain text (or non-gz): zlib gzopen handles both transparently
    return std::make_unique<FastqReader>(path);
}
