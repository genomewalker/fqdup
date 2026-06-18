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
};

// Factory: returns the best available backend (rapidgzip > ISA-L > zlib).
// threads: decompression thread budget for rapidgzip (0 = auto from hardware_concurrency).
// Callers that already have worker pools should pass their own thread count to avoid
// oversubscription. Implemented in src/fastq_io_backend.cpp.
std::unique_ptr<FastqReaderBase> make_fastq_reader(const std::string& path,
                                                    size_t threads = 0);

// Chains multiple files into a single sequential reader.
// If paths.size() == 1 this is identical to make_fastq_reader(paths[0], threads).
std::unique_ptr<FastqReaderBase> make_chained_fastq_reader(
    const std::vector<std::string>& paths, size_t threads = 0);
