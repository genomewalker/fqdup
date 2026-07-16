// fastq_io_backend.cpp
// Implements make_fastq_reader(): ISA-L (igzip) for gz, zlib for plain/stdin.
//
// rapidgzip was removed 2026-07-16. Its parallel-decode path raced on large gz over NFS and
// SILENTLY truncated a 30 GB sort input, costing 115,146,482 reads with no error anywhere;
// `gzip -t` cannot see it (a race-truncated member is a complete, valid gzip). It was also
// MEASURED 9.7x SLOWER than ISA-L on the exact KapK profile (job 20189414: 2427.85 s vs
// 250.98 s, cpu 132% vs 1046% — the "parallel" decoder just stalled on the filer). Both worse,
// no upside, so retire it rather than fix it. If decode ever becomes the bottleneck, parallelise
// by gzip MEMBER (the merged files are ~1.32M independent 64 KB members) — race-free by
// construction, unlike rapidgzip's shared-state inflate.

#include "fqdup/fastq_common.hpp"
#include <algorithm>
#include <cstdlib>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <stdexcept>

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
    (void)threads;  // ISA-L decode is single-threaded; the compute -p never drove decode
    if (is_gzip(path)) {
        // Single-threaded ISA-L: rapidgzip's parallel-decode race cannot exist by construction,
        // and it was measured 9.7x faster anyway (see file header). igzip decodes ~98 MB/s here,
        // I/O-bound on the filer, so parallel decode buys nothing on this workload.
#if defined(HAVE_ISAL)
        return std::make_unique<FastqReaderIgzip>(path);
#else
        return std::make_unique<FastqReader>(path);  // zlib fallback when ISA-L is unavailable
#endif
    }
    // Plain text (or non-gz): zlib gzopen handles both transparently
    return std::make_unique<FastqReader>(path);
}

namespace {
class ChainedFastqReader : public FastqReaderBase {
public:
    ChainedFastqReader(std::vector<std::string> paths, size_t threads)
        : paths_(std::move(paths)), idx_(0), threads_(threads), total_(0) {
        if (!paths_.empty())
            cur_ = make_fastq_reader(paths_[0], threads_);
    }

    bool read(FastqRecord& rec) override {
        while (cur_) {
            if (cur_->read(rec)) { ++total_; return true; }
            if (++idx_ < paths_.size())
                cur_ = make_fastq_reader(paths_[idx_], threads_);
            else
                cur_.reset();
        }
        return false;
    }

    uint64_t record_count() const override { return total_; }

    // Every chained file goes through make_fastq_reader, so they all share a backend.
    const char* backend_name() const override {
        return cur_ ? cur_->backend_name() : "(no input)";
    }

private:
    std::vector<std::string> paths_;
    size_t idx_, threads_;
    uint64_t total_;
    std::unique_ptr<FastqReaderBase> cur_;
};
}  // namespace

std::unique_ptr<FastqReaderBase> make_chained_fastq_reader(
    const std::vector<std::string>& paths, size_t threads) {
    if (paths.size() == 1) return make_fastq_reader(paths[0], threads);
    return std::make_unique<ChainedFastqReader>(paths, threads);
}
