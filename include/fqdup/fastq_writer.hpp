#pragma once
// Single shared FASTQ writer for every subcommand (sort, derep, split, merge…).
//
// Backend precedence on the compressed path:
//   bgzf (htslib, multi-threaded) when n_threads>1,
//   else isa-l igzip (fast single-thread deflate) when built with HAVE_ISAL,
//   else zlib.
// isa-l manages the full gzip wrapper (header + CRC32/ISIZE trailer) via
// IGZIP_GZIP — never hand-roll the header/trailer (doing so once dropped the
// trailer and produced "missing trailer" output).
//
// Header-only, all methods implicitly inline; the class has external linkage and
// an identical definition in every TU, so it is ODR-safe to include anywhere.

#include "fqdup/fastq_types.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

#ifdef HAVE_ISAL
#include <isa-l/igzip_lib.h>
#endif

#ifdef HAVE_BGZF
#include <htslib/bgzf.h>
#endif

// Should output at `path` be gzip-compressed? Compression is chosen from the NAME, so this is the
// one place that decision lives -- callers must not re-derive it with an ad-hoc `endswith(".gz")`,
// because that silently writes PLAIN fastq behind a .gz name when the path carries an atomic-rename
// temp suffix (the classic write-to-`out.fq.gz.tmp`-then-rename pattern). Strip one known temp
// suffix, then honour the real extension.
inline bool output_wants_gzip(std::string path) {
    static const char* const kTempSuffixes[] = {
        ".tmp", ".part", ".partial", ".inprogress", ".building", ".incomplete"};
    for (const char* suf : kTempSuffixes) {
        const size_t n = std::strlen(suf);
        if (path.size() > n && path.compare(path.size() - n, n, suf) == 0) {
            path.resize(path.size() - n);
            break;  // strip at most one temp suffix
        }
    }
    return path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

class FastqWriter {
public:
    FastqWriter(const std::string& path, bool compress, int n_threads = 1)
        : path_(path), compress_(compress), gzfp_(nullptr)
#ifdef HAVE_BGZF
        , bgzfp_(nullptr)
#endif
    {
        if (compress_) {
#ifdef HAVE_BGZF
            if (n_threads > 1) {
                bgzfp_ = bgzf_open(path.c_str(), "w");
                if (!bgzfp_)
                    throw std::runtime_error("Cannot open output: " + path);
                bgzf_mt(bgzfp_, n_threads, 0);
                return;
            }
#endif
#ifdef HAVE_ISAL
            fp_ = fopen(path.c_str(), "wb");
            if (!fp_) throw std::runtime_error("Cannot open output: " + path);
            out_buf_.resize(OUT_BUF_SIZE);
            level_buf_.resize(ISAL_DEF_LVL1_DEFAULT);
            isal_deflate_init(&stream_);
            stream_.level          = 1;
            stream_.level_buf      = level_buf_.data();
            stream_.level_buf_size = static_cast<uint32_t>(level_buf_.size());
            stream_.next_out       = out_buf_.data();
            stream_.avail_out      = static_cast<uint32_t>(out_buf_.size());
            stream_.gzip_flag      = IGZIP_GZIP;  // isa-l writes header + trailer
#else
            gzfp_ = gzopen(path.c_str(), "wb6");
            if (!gzfp_)
                throw std::runtime_error("Cannot open output: " + path);
            gzbuffer(gzfp_, OUT_BUF_SIZE);
#endif
        } else {
            plain_out_.open(path);
            if (!plain_out_.good())
                throw std::runtime_error("Cannot open output: " + path);
        }
    }

    ~FastqWriter() {
        try { close(); }
        catch (const std::exception& e) {
            std::cerr << "Fatal: " << e.what() << "\n";
            std::exit(1);
        }
        catch (...) {
            std::cerr << "Fatal: unknown error closing output\n";
            std::exit(1);
        }
    }

    void close() {
        if (closed_) return;
        closed_ = true;
#ifdef HAVE_BGZF
        if (bgzfp_) {
            if (bgzf_close(bgzfp_) < 0)
                throw std::runtime_error("bgzf_close failed writing " + path_);
            bgzfp_ = nullptr;
            return;
        }
#endif
#ifdef HAVE_ISAL
        if (fp_) {
            stream_.end_of_stream = 1;
            stream_.next_in  = nullptr;
            stream_.avail_in = 0;
            // Finalize per isa-l contract: call isal_deflate until ZSTATE_END.
            // Only then has the final block + gzip trailer (CRC32+ISIZE) been
            // emitted. Terminating early drops the trailer -> "missing trailer".
            do {
                if (stream_.avail_out == 0) flush_isal_buf();
                if (isal_deflate(&stream_) != COMP_OK)
                    throw std::runtime_error("isal_deflate failed finalizing " + path_);
            } while (stream_.internal_state.state != ZSTATE_END);
            flush_isal_buf();
            if (fclose(fp_) != 0)
                throw std::runtime_error("fclose failed writing " + path_);
            fp_ = nullptr;
            return;
        }
#endif
        if (gzfp_) {
            if (gzclose(gzfp_) != Z_OK)
                throw std::runtime_error("gzclose failed writing " + path_);
            gzfp_ = nullptr;
        }
        if (plain_out_.is_open()) {
            plain_out_.flush();
            plain_out_.close();
            if (plain_out_.fail())
                throw std::runtime_error("close failed writing " + path_);
        }
    }

    void write(const FastqRecord& rec) {
        write_fields(rec.header.data(), rec.header.size(),
                     rec.seq.data(),    rec.seq.size(),
                     rec.plus.data(),   rec.plus.size(),
                     rec.qual.data(),   rec.qual.size());
    }

    // Backend-agnostic record write from raw field pointers+lengths. Lets
    // callers with arena/zero-copy records (e.g. sort's FastqRecordArena) use
    // this writer without a dependency on their record type.
    void write_fields(const char* h, size_t hl, const char* s, size_t sl,
                      const char* p, size_t pl, const char* q, size_t ql) {
        // One call per record, not eight. Both backends buffer internally and take a lock
        // (or a branch) per call, so a 188M-read merge was paying 1.5 billion of them to
        // move ~4 short fields. Assemble the record once and hand it over whole: the byte
        // stream is unchanged, only its call granularity is.
        rec_buf_.clear();
        rec_buf_.append(h, hl); rec_buf_.push_back('\n');
        rec_buf_.append(s, sl); rec_buf_.push_back('\n');
        rec_buf_.append(p, pl); rec_buf_.push_back('\n');
        rec_buf_.append(q, ql); rec_buf_.push_back('\n');

#ifdef HAVE_BGZF
        if (bgzfp_) {
            if (bgzf_write(bgzfp_, rec_buf_.data(), rec_buf_.size()) < 0)
                throw std::runtime_error("bgzf_write failed writing " + path_);
            return;
        }
#endif
#ifdef HAVE_ISAL
        if (fp_) {
            isal_write(rec_buf_.data(), rec_buf_.size());
            return;
        }
#endif
        if (compress_) {
            if (gzprintf(gzfp_, "%.*s\n%.*s\n%.*s\n%.*s\n",
                         (int)hl, h, (int)sl, s, (int)pl, p, (int)ql, q) < 0)
                throw std::runtime_error("gzprintf failed writing FASTQ record");
        } else {
            // sort's temp chunks land here (written uncompressed by design).
            plain_out_.write(rec_buf_.data(), static_cast<std::streamsize>(rec_buf_.size()));
            if (!plain_out_.good())
                throw std::runtime_error("write failed writing " + path_);
        }
    }

private:
    // Scratch for the one-call-per-record assembly above. A member, not a local, so its
    // capacity survives across records. Not thread-safe -- but a FastqWriter never was:
    // two threads sharing one would already interleave record bytes.
    std::string rec_buf_;

#ifdef HAVE_ISAL
    void isal_write(const char* data, size_t len) {
        stream_.next_in  = const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(data));
        stream_.avail_in = static_cast<uint32_t>(len);
        while (stream_.avail_in > 0) {
            if (stream_.avail_out == 0) flush_isal_buf();
            isal_deflate(&stream_);
        }
    }

    void flush_isal_buf() {
        size_t n = out_buf_.size() - stream_.avail_out;
        if (n == 0) return;
        if (fwrite(out_buf_.data(), 1, n, fp_) != n)
            throw std::runtime_error("fwrite failed writing " + path_);
        stream_.next_out  = out_buf_.data();
        stream_.avail_out = static_cast<uint32_t>(out_buf_.size());
    }
#endif

    static constexpr size_t OUT_BUF_SIZE = 4 * 1024 * 1024;  // 4 MB
    std::string   path_;
    bool          compress_;
    bool          closed_ = false;
    gzFile        gzfp_;
#ifdef HAVE_ISAL
    FILE*                fp_ = nullptr;
    isal_zstream         stream_{};
    std::vector<uint8_t> out_buf_;
    std::vector<uint8_t> level_buf_;
#endif
#ifdef HAVE_BGZF
    BGZF*         bgzfp_;
#endif
    std::ofstream plain_out_;
};
