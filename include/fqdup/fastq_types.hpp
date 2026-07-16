#pragma once
// fastq_types.hpp — global-scope FASTQ types and factory declaration.
//
// Only FastqRecord, FastqReaderBase, and make_fastq_reader() live here.
// These are at global scope so they can safely cross TU boundaries without
// ODR violations. This header has no anonymous namespace — it is safe to
// include from any translation unit, including sort.cpp.
//
// The reader/writer implementations with inline method bodies are in
// fastq_common.hpp (anonymous namespace, one definition per TU).

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// ============================================================================
// FASTQ Record — global scope, safe to cross TU boundaries
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
// Abstract reader — virtual dispatch across TU boundaries
// ============================================================================

class FastqReaderBase {
public:
    virtual ~FastqReaderBase() = default;
    virtual bool read(FastqRecord& rec) = 0;
    virtual uint64_t record_count() const = 0;
    // The reader names itself. Call sites used to log the backend from #ifdef
    // HAVE_RAPIDGZIP, i.e. "was it compiled in" -- not "was it selected". derep
    // therefore printed "Decompression: rapidgzip" on runs that decoded with
    // ISA-L, which is precisely the line an operator reads to confirm the
    // parallel-decode race is not in play.
    virtual const char* backend_name() const = 0;
};

// Factory: ISA-L (igzip) for gz when available, else zlib; zlib for plain/stdin.
// `threads` is ignored for decode (ISA-L is single-threaded) and kept only for signature
// compatibility. rapidgzip was removed 2026-07-16 (raced on NFS + 9.7x slower than ISA-L;
// see fastq_io_backend.cpp). Implemented in src/fastq_io_backend.cpp.
std::unique_ptr<FastqReaderBase> make_fastq_reader(const std::string& path,
                                                    size_t threads = 0);

// Chains multiple files into a single sequential reader.
// If paths.size() == 1 this is identical to make_fastq_reader(paths[0], threads).
std::unique_ptr<FastqReaderBase> make_chained_fastq_reader(
    const std::vector<std::string>& paths, size_t threads = 0);
