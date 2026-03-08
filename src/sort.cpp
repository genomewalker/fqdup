// sort.cpp
// Optimized external FASTQ sorter with parallel chunk processing
//
// Strategy:
// 1. Single-threaded chunk creation (write to disk immediately)
// 2. Parallel chunk sorting (read from disk, sort, write back)
// 3. K-way merge
//
// Memory: Predictable - only 1 chunk in RAM at a time during creation
// Speed: Parallel sorting of chunks for maximum throughput

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>
#include <zlib.h>

#include "fqdup/fastq_types.hpp"  // FastqRecord, FastqReaderBase, make_fastq_reader
#include "fqdup/logger.hpp"

// For memory release to OS
#ifdef __linux__
#include <malloc.h>  // For malloc_trim (glibc)
#endif

// jemalloc memory release
#ifdef __GLIBC__
extern "C" {
    // jemalloc specific - declare mallctl if using jemalloc
    int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen) __attribute__((weak));
}
#endif

// SortReaderBridge: wraps make_fastq_reader() (rapidgzip/ISA-L/zlib) and bridges
// its FastqRecord interface to the arena-based read used by Phase 1 chunk creation.
// Lives outside the anonymous namespace — uses global FastqRecord from fastq_types.hpp.
class SortReaderBridge {
    std::unique_ptr<FastqReaderBase> inner_;
    FastqRecord tmp_;
public:
    explicit SortReaderBridge(const std::string& path)
        : inner_(make_fastq_reader(path)) {}

    // Arena-based read for Phase 1: zero heap allocation per record.
    // StringArena and FastqRecordArena are defined below in the anonymous namespace;
    // this method is defined after them (see bottom of anonymous namespace section).
    template<typename Arena, typename ArenRec>
    bool read(ArenRec& rec, Arena& arena) {
        if (!inner_->read(tmp_)) return false;
        auto store = [&](const std::string& s) { return arena.store(s.data(), s.size()); };
        rec.header     = store(tmp_.header); rec.header_len = tmp_.header.size();
        rec.seq        = store(tmp_.seq);    rec.seq_len    = tmp_.seq.size();
        rec.plus       = store(tmp_.plus);   rec.plus_len   = tmp_.plus.size();
        rec.qual       = store(tmp_.qual);   rec.qual_len   = tmp_.qual.size();
        return true;
    }

    // Direct read for merge phase (no arena needed).
    bool read(FastqRecord& rec) { return inner_->read(rec); }
};

// All remaining file-local types are in an anonymous namespace to give their
// member functions internal linkage and avoid ODR violations.
namespace {

constexpr size_t GZBUF_SIZE = 4 * 1024 * 1024;  // 4MB for better I/O performance

// Parse memory size string (e.g., "4G", "2048M")
static size_t parse_memory_size(const std::string& str) {
    if (str.empty()) return 2048;

    size_t value = 0;
    size_t i = 0;
    while (i < str.size() && std::isdigit(str[i])) {
        value = value * 10 + (str[i] - '0');
        i++;
    }

    if (i < str.size()) {
        char unit = std::toupper(str[i]);
        switch (unit) {
            case 'K': {
                size_t kb = value;
                if (kb < 1024)
                    throw std::runtime_error(
                        "Memory value too small: '" + str + "' (minimum 1M / 1024K)");
                return kb / 1024;
            }
            case 'M': return value;
            case 'G': return value * 1024;
            case 'T': return value * 1024 * 1024;
            default: throw std::runtime_error("Unknown memory unit '" +
                         std::string(1, unit) + "' in '" + str +
                         "' (use K, M, G, or T)");
        }
    }
    return value;
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

// String arena - single large allocation for all string data in a chunk
// Eliminates millions of individual malloc calls and fragmentation
class StringArena {
public:
    explicit StringArena(size_t capacity) : capacity_(capacity), used_(0) {
        data_ = new char[capacity_];
    }

    ~StringArena() {
        delete[] data_;
    }

    // Allocate string from arena, return pointer into arena
    const char* store(const std::string& str) {
        return store(str.data(), str.size());
    }

    const char* store(const char* data, size_t len) {
        if (used_ + len > capacity_) {
            throw std::runtime_error("StringArena overflow - increase chunk size");
        }
        char* ptr = data_ + used_;
        std::memcpy(ptr, data, len);
        used_ += len;
        return ptr;
    }

    void reset() {
        used_ = 0;
    }

    size_t bytes_used() const { return used_; }
    size_t capacity() const { return capacity_; }

private:
    char* data_;
    size_t capacity_;
    size_t used_;
};

// FastqRecord is provided by fastq_types.hpp (global scope).

// FastqRecordArena using raw pointers + lengths (zero allocations for chunk creation!)
struct FastqRecordArena {
    const char* header = nullptr;
    const char* seq = nullptr;
    const char* plus = nullptr;
    const char* qual = nullptr;

    size_t header_len = 0;
    size_t seq_len = 0;
    size_t plus_len = 0;
    size_t qual_len = 0;

    void clear() {
        header = seq = plus = qual = nullptr;
        header_len = seq_len = plus_len = qual_len = 0;
    }

    size_t memory_size() const {
        return header_len + seq_len + plus_len + qual_len;
    }

    // Get string_view for each field
    std::string_view get_header() const { return std::string_view(header, header_len); }
    std::string_view get_seq() const { return std::string_view(seq, seq_len); }
    std::string_view get_plus() const { return std::string_view(plus, plus_len); }
    std::string_view get_qual() const { return std::string_view(qual, qual_len); }
};

// Chunk buffer: vector of arena records + string arena
struct ChunkBuffer {
    std::vector<FastqRecordArena> records;
    StringArena arena;

    explicit ChunkBuffer(size_t chunk_size_bytes)
        : arena(chunk_size_bytes * 3) {  // 3x for safety (strings are bigger than raw size)
        records.reserve(chunk_size_bytes / 300);  // Estimate ~300 bytes per record
    }

    void reset() {
        records.clear();
        arena.reset();
    }
};

// Extract read ID from header (remove @ and everything after first space/tab)
static std::string_view trim_id(std::string_view header) {
    size_t pos = (!header.empty() && header[0] == '@') ? 1 : 0;
    size_t sp = header.find_first_of(" \t", pos);
    return (sp == std::string_view::npos) ? header.substr(pos) : header.substr(pos, sp - pos);
}

// Parse numeric suffix for natural ordering (e.g., "read_123" -> 123)
static bool parse_numeric_suffix(const std::string& id, std::string& prefix, long long& num) {
    size_t i = id.size();

    // Find last run of digits
    while (i > 0 && std::isdigit(id[i - 1])) {
        i--;
    }

    if (i == id.size()) {
        // No digits at end
        prefix = id;
        num = 0;
        return false;
    }

    prefix = id.substr(0, i);
    num = std::stoll(id.substr(i));
    return true;
}

// Comparators now work with string_view (FastqRecord stores pointers)
struct NaturalOrderComparator {
    bool operator()(const FastqRecord& a, const FastqRecord& b) const {
        std::string_view id_a = trim_id(a.header);
        std::string_view id_b = trim_id(b.header);

        // Find numeric suffix
        size_t i = id_a.size();
        while (i > 0 && std::isdigit(id_a[i-1])) i--;
        size_t j = id_b.size();
        while (j > 0 && std::isdigit(id_b[j-1])) j--;

        // Compare prefixes
        std::string_view prefix_a = id_a.substr(0, i);
        std::string_view prefix_b = id_b.substr(0, j);

        if (prefix_a != prefix_b) return prefix_a < prefix_b;
        if (i == id_a.size() && j == id_b.size()) return false;
        if (i == id_a.size()) return true;
        if (j == id_b.size()) return false;

        // Parse and compare numbers
        long long num_a = std::stoll(std::string(id_a.substr(i)));
        long long num_b = std::stoll(std::string(id_b.substr(j)));
        return num_a < num_b;
    }
};

struct LexicographicComparator {
    bool operator()(const FastqRecord& a, const FastqRecord& b) const {
        return trim_id(a.header) < trim_id(b.header);
    }
};

// Legacy zlib reader — used by the dead-code ParallelChunkSorter only.
// Live paths (create_chunks, merge_chunks) now use SortReaderBridge above.
class FastqReader {
public:
    FastqReader(const std::string& path) : path_(path), use_stdin_(false), pipe_(nullptr) {
        // Check if reading from stdin
        if (path == "/dev/stdin" || path == "-") {
            use_stdin_ = true;
            gzfp_ = gzdopen(fileno(stdin), "rb");
            if (!gzfp_) throw std::runtime_error("Cannot open stdin");
        } else {
            gzfp_ = gzopen(path.c_str(), "rb");
            if (!gzfp_) throw std::runtime_error("Cannot open: " + path);
        }
        gzbuffer(gzfp_, GZBUF_SIZE);
    }

    // Constructor for reading from FILE* pipe
    FastqReader(FILE* pipe, bool is_pipe) : path_("pipe"), use_stdin_(false), pipe_(pipe), gzfp_(nullptr) {
        if (!pipe_) throw std::runtime_error("Invalid pipe");
    }

    ~FastqReader() {
        if (gzfp_) gzclose(gzfp_);
        // Note: Don't close pipe_ here, caller will handle it with pclose()
    }

    // Old interface (for merge phase that needs std::string)
    bool read(std::string& header, std::string& seq, std::string& plus, std::string& qual) {
        if (!getline_gz(header)) return false;
        if (!getline_gz(seq))
            throw std::runtime_error("Truncated FASTQ: missing sequence line after header '" + header + "'");
        if (!getline_gz(plus))
            throw std::runtime_error("Truncated FASTQ: missing '+' line after sequence");
        if (!getline_gz(qual))
            throw std::runtime_error("Truncated FASTQ: missing quality line");
        return true;
    }

    bool read(FastqRecord& rec) {
        if (!getline_gz(rec.header)) return false;
        if (!getline_gz(rec.seq))
            throw std::runtime_error("Truncated FASTQ: missing sequence line after header '" + rec.header + "'");
        if (!getline_gz(rec.plus))
            throw std::runtime_error("Truncated FASTQ: missing '+' line after sequence");
        if (!getline_gz(rec.qual))
            throw std::runtime_error("Truncated FASTQ: missing quality line");
        return true;
    }

private:
    bool getline_gz(std::string& line) {
        line.clear();
        line.reserve(256);
        char buffer[8192];

        if (pipe_) {
            // Read from uncompressed pipe (from pigz)
            while (true) {
                if (fgets(buffer, sizeof(buffer), pipe_) == nullptr) return !line.empty();
                size_t len = strlen(buffer);
                if (len > 0 && buffer[len - 1] == '\n') {
                    line.append(buffer, len - 1);
                    return true;
                }
                line.append(buffer, len);
            }
        } else {
            // Read from gzFile
            while (true) {
                if (gzgets(gzfp_, buffer, sizeof(buffer)) == nullptr) return !line.empty();
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
    bool use_stdin_;
    FILE* pipe_;
    gzFile gzfp_;
};

// FastqReaderIgzip removed from sort.cpp — ISA-L and rapidgzip decompression
// are now accessed through SortReaderBridge → make_fastq_reader().


class FastqWriter {
public:
    FastqWriter(const std::string& path, bool compress) : compress_(compress) {
        if (compress_) {
            gzfp_ = gzopen(path.c_str(), "wb6");
            if (!gzfp_) throw std::runtime_error("Cannot open: " + path);
            gzbuffer(gzfp_, GZBUF_SIZE);
        } else {
            out_.open(path);
            if (!out_.good()) throw std::runtime_error("Cannot open: " + path);
        }
    }
    ~FastqWriter() { if (gzfp_) gzclose(gzfp_); }

    // Write FastqRecord (std::string - for merge phase)
    void write(const FastqRecord& rec) {
        if (compress_) {
            if (gzprintf(gzfp_, "%s\n%s\n%s\n%s\n",
                         rec.header.c_str(), rec.seq.c_str(),
                         rec.plus.c_str(), rec.qual.c_str()) < 0) {
                throw std::runtime_error("gzprintf failed while writing compressed FASTQ record");
            }
        } else {
            out_ << rec.header << '\n' << rec.seq << '\n'
                 << rec.plus << '\n' << rec.qual << '\n';
        }
    }

    // Write FastqRecordArena (raw pointers - for chunk writing)
    void write(const FastqRecordArena& rec) {
        if (compress_) {
            if (gzprintf(gzfp_, "%.*s\n%.*s\n%.*s\n%.*s\n",
                         (int)rec.header_len, rec.header,
                         (int)rec.seq_len, rec.seq,
                         (int)rec.plus_len, rec.plus,
                         (int)rec.qual_len, rec.qual) < 0) {
                throw std::runtime_error("gzprintf failed while writing compressed FASTQ record");
            }
        } else {
            out_.write(rec.header, rec.header_len) << '\n';
            out_.write(rec.seq, rec.seq_len) << '\n';
            out_.write(rec.plus, rec.plus_len) << '\n';
            out_.write(rec.qual, rec.qual_len) << '\n';
        }
    }

private:
    bool compress_;
    gzFile gzfp_ = nullptr;
    std::ofstream out_;
};

// Parallel chunk sorter - sorts multiple chunks concurrently
class ParallelChunkSorter {
public:
    ParallelChunkSorter(int num_threads, bool natural_order, bool fast_mode)
        : num_threads_(num_threads), natural_order_(natural_order), fast_mode_(fast_mode), done_(false), active_workers_(0) {

        for (int i = 0; i < num_threads_; ++i) {
            workers_.emplace_back(&ParallelChunkSorter::worker_thread, this);
        }
    }

    ~ParallelChunkSorter() {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            done_ = true;
        }
        cv_.notify_all();

        for (auto& w : workers_) {
            if (w.joinable()) w.join();
        }
    }

    // Submit unsorted chunk file for sorting
    void submit(const std::string& unsorted_path, const std::string& sorted_path) {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            queue_.push({unsorted_path, sorted_path});
            chunks_queued_++;
        }
        cv_.notify_one();
    }

    // Wait for all chunks to be sorted
    void wait() {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_done_.wait(lock, [this]{ return queue_.empty() && active_workers_ == 0; });
    }

    size_t get_chunks_sorted() const { return chunks_sorted_; }

private:
    struct ChunkTask {
        std::string unsorted_path;
        std::string sorted_path;
    };

    void worker_thread() {
        while (true) {
            ChunkTask task;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [this]{ return !queue_.empty() || done_; });

                if (queue_.empty() && done_) break;
                if (!queue_.empty()) {
                    task = queue_.front();
                    queue_.pop();
                    active_workers_++;
                }
            }

            if (!task.unsorted_path.empty()) {
                sort_chunk(task.unsorted_path, task.sorted_path);

                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    active_workers_--;
                    chunks_sorted_++;

                    // Delete unsorted file
                    std::remove(task.unsorted_path.c_str());

                    if (queue_.empty() && active_workers_ == 0) {
                        cv_done_.notify_all();
                    }
                }
            }
        }
    }

    void sort_chunk(const std::string& input, const std::string& output) {
        // Read chunk - use standard FastqReader (handles both .gz and uncompressed)
        // Note: chunks are often uncompressed in fast mode, so using zlib is fine
        FastqReader reader(input);
        std::vector<FastqRecord> records;
        records.reserve(50000);  // Pre-allocate for typical chunk
        FastqRecord rec;
        while (reader.read(rec)) {
            records.push_back(std::move(rec));
        }

        const uint32_t N = static_cast<uint32_t>(records.size());
        std::vector<uint32_t> indices(N);
        for (uint32_t i = 0; i < N; ++i) indices[i] = i;

        if (natural_order_) {
            struct NatKey { std::string_view prefix; uint64_t suffix; bool has_suffix; };
            std::vector<NatKey> keys;
            keys.reserve(N);
            for (const auto& r : records) {
                std::string_view id = trim_id(r.header);
                size_t i = id.size();
                while (i > 0 && static_cast<unsigned>(id[i-1] - '0') <= 9u) --i;
                uint64_t suf = 0;
                for (size_t p = i; p < id.size(); ++p)
                    suf = suf * 10u + static_cast<unsigned>(id[p] - '0');
                keys.push_back({id.substr(0, i), suf, i < id.size()});
            }
            std::sort(indices.begin(), indices.end(), [&keys](uint32_t a, uint32_t b) {
                const NatKey& ka = keys[a]; const NatKey& kb = keys[b];
                if (ka.prefix != kb.prefix) return ka.prefix < kb.prefix;
                if (!ka.has_suffix && !kb.has_suffix) return false;
                if (!ka.has_suffix) return true;
                if (!kb.has_suffix) return false;
                return ka.suffix < kb.suffix;
            });
        } else {
            std::vector<std::string_view> keys;
            keys.reserve(N);
            for (const auto& r : records)
                keys.push_back(trim_id(r.header));
            std::sort(indices.begin(), indices.end(), [&keys](uint32_t a, uint32_t b) {
                return keys[a] < keys[b];
            });
        }

        // Apply permutation to records (write in sorted order)
        // In fast mode, write uncompressed for speed
        bool compress = !fast_mode_;
        FastqWriter writer(output, compress);
        for (uint32_t idx : indices) {
            writer.write(records[idx]);
        }
    }

    int num_threads_;
    bool natural_order_;
    bool fast_mode_;
    std::vector<std::thread> workers_;
    std::queue<ChunkTask> queue_;
    std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable cv_done_;
    std::atomic<bool> done_;
    std::atomic<int> active_workers_;
    std::atomic<size_t> chunks_queued_{0};
    std::atomic<size_t> chunks_sorted_{0};
};

// Main sorter with parallel chunk processing
class OptimizedSorter {
public:
    OptimizedSorter(size_t max_memory_mb, const std::string& temp_dir, int threads, bool natural_order, bool fast_mode)
        : max_memory_mb_(max_memory_mb), temp_dir_(temp_dir), threads_(threads),
          natural_order_(natural_order), fast_mode_(fast_mode), chunks_created_(0) {

        // Memory model: C++ overhead is significant
        // - std::string overhead: 32 bytes per string × 4 = 128 bytes per record
        // - std::vector capacity overgrowth: ~1.5-2x
        // - malloc fragmentation: ~1.2x
        // - Total multiplier: ~6-8x raw data
        size_t target_chunks = threads_;
        chunk_size_bytes_ = (max_memory_mb_ * 1024ULL * 1024ULL) / (target_chunks * 10);

        // Use 1GB chunks to reduce memory pressure
        const size_t MIN_CHUNK_MB = 1024;
        if (chunk_size_bytes_ < MIN_CHUNK_MB * 1024ULL * 1024ULL) {
            chunk_size_bytes_ = MIN_CHUNK_MB * 1024ULL * 1024ULL;
        }
    }

    ~OptimizedSorter() {
        for (const auto& f : sorted_files_) std::remove(f.c_str());
    }

    void sort(const std::string& input, const std::string& output) {
        // Estimate input size and adjust chunk size accordingly
        adjust_chunk_size_for_input(input);

        // sort.cpp cannot use rapidgzip (ODR constraint — fastq_io_backend.cpp owns
        // all rapidgzip symbols). Reader selection is ISA-L or zlib only.
        log_info("Decompression: " + std::string(
#ifdef HAVE_ISAL
            "ISA-L (hardware-accelerated)"
#else
            "zlib"
#endif
        ));
        log_info("Phase 1: Creating sorted chunks (" + std::to_string(threads_) + " writer threads)");
        create_chunks(input);

        log_info("Phase 2: Merging " + std::to_string(sorted_files_.size()) + " sorted chunks");
        merge_chunks(output);
    }

private:
    void adjust_chunk_size_for_input(const std::string& input) {
        // Skip for stdin (can't determine size)
        if (input == "/dev/stdin" || input == "-") {
            return;
        }

        // Get input file size
        std::ifstream file(input, std::ios::binary | std::ios::ate);
        if (!file.good()) {
            return;
        }
        size_t compressed_size_bytes = file.tellg();
        file.close();

        double compressed_size_gb = compressed_size_bytes / (1024.0 * 1024.0 * 1024.0);

        // Estimate uncompressed size (gzip typically 3-4× compression for FASTQ)
        double compression_ratio = 3.5;
        double estimated_uncompressed_gb = compressed_size_gb * compression_ratio;

        // Calculate optimal number of chunks based on data size
        // Use fewer, larger chunks for better parallelism in Phase 2
        size_t target_chunks = threads_;

        if (estimated_uncompressed_gb < 10) {
            target_chunks = threads_ / 2;  // Very small: 8 chunks for 16 threads
            if (target_chunks < 4) target_chunks = 4;  // Minimum 4 chunks
        } else if (estimated_uncompressed_gb > 100) {
            target_chunks = threads_ * 2;  // Large: 32 chunks for 16 threads
        } else {
            target_chunks = threads_;  // Medium: 16 chunks for 16 threads
        }

        // Calculate chunk size from estimated data size
        size_t estimated_chunk_mb = (estimated_uncompressed_gb * 1024) / target_chunks;

        // Apply constraints
        const size_t MIN_CHUNK_MB = 2048;   // Minimum 2GB chunks
        const size_t MAX_CHUNK_MB = max_memory_mb_ / 20;  // Safety limit

        estimated_chunk_mb = std::max(MIN_CHUNK_MB,
                             std::min(MAX_CHUNK_MB, estimated_chunk_mb));

        // Round to nearest GB for cleaner numbers
        estimated_chunk_mb = ((estimated_chunk_mb + 512) / 1024) * 1024;

        chunk_size_bytes_ = estimated_chunk_mb * 1024ULL * 1024ULL;

        log_info("Chunk size: " + std::to_string(estimated_chunk_mb) + " MB (" +
                std::to_string(estimated_chunk_mb / 1024.0) + " GB)");
    }

    void create_chunks(const std::string& input) {
        auto phase1_start = std::chrono::steady_clock::now();

        // SortReaderBridge dispatches to the best available backend via make_fastq_reader():
        // rapidgzip (parallel, multi-threaded) > ISA-L (SIMD) > zlib.
        SortReaderBridge reader(input);

        // Now with string arenas, we can use maximum concurrency!
        // Each buffer has its own arena - zero malloc fragmentation
        int num_writers = threads_ / 2;  // Use half threads for writing
        if (num_writers < 2) num_writers = 2;
        const size_t NUM_CHUNK_BUFFERS = num_writers + 2;

        const size_t max_memory_bytes = max_memory_mb_ * 1024ULL * 1024ULL;
        const size_t max_chunk_bytes = max_memory_bytes / (NUM_CHUNK_BUFFERS * 3ULL);
        if (max_chunk_bytes == 0) {
            throw std::runtime_error("Insufficient --max-memory for configured chunk buffers");
        }
        if (chunk_size_bytes_ > max_chunk_bytes) {
            chunk_size_bytes_ = max_chunk_bytes;
            log_info("Adjusted chunk size to " + std::to_string(chunk_size_bytes_ / (1024ULL * 1024ULL)) +
                     " MB to stay within --max-memory preallocation budget");
        }

        // Pre-allocate all chunk buffers with arenas
        std::vector<ChunkBuffer> chunk_pool;
        chunk_pool.reserve(NUM_CHUNK_BUFFERS);
        for (size_t i = 0; i < NUM_CHUNK_BUFFERS; ++i) {
            chunk_pool.emplace_back(chunk_size_bytes_);
        }

        // Track which buffers are available for reuse
        std::queue<size_t> available_buffers;
        for (size_t i = 0; i < NUM_CHUNK_BUFFERS; ++i) {
            available_buffers.push(i);
        }

        std::queue<size_t> write_queue;  // Queue of buffer indices, not actual data
        std::mutex queue_mutex;
        std::condition_variable queue_cv;
        std::condition_variable reader_cv;  // For backpressure
        std::atomic<bool> done_reading{false};
        std::atomic<int> active_writers{0};
        std::atomic<bool> writer_failed{false};
        std::exception_ptr writer_exception;
        std::mutex writer_exception_mutex;

        size_t chunk_bytes = 0;
        size_t total_reads = 0;
        size_t current_buffer_idx = 0;  // Will be set from available_buffers

        // Time spent reading vs writing
        size_t read_time_ms = 0;
        auto last_read_start = std::chrono::steady_clock::now();

        auto writer_func = [&]() {
            while (true) {
                size_t buffer_idx = 0;
                bool has_work = false;
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    queue_cv.wait(lock, [&]{ return !write_queue.empty() || done_reading.load(); });

                    if (write_queue.empty() && done_reading.load()) break;

                    if (!write_queue.empty()) {
                        buffer_idx = write_queue.front();
                        write_queue.pop();
                        active_writers++;
                        has_work = true;
                    }
                }

                if (has_work) {
                    try {
                        // Sort and write the chunk from the buffer pool
                        write_sorted_chunk(std::move(chunk_pool[buffer_idx]));
                    } catch (...) {
                        active_writers--;
                        {
                            std::lock_guard<std::mutex> lock(writer_exception_mutex);
                            if (!writer_exception) {
                                writer_exception = std::current_exception();
                            }
                        }
                        writer_failed = true;
                        done_reading = true;
                        queue_cv.notify_all();
                        reader_cv.notify_all();
                        break;
                    }
                    active_writers--;

                    // Return buffer to available pool
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex);
                        chunk_pool[buffer_idx].reset();  // Reset for reuse (arena + vector)
                        available_buffers.push(buffer_idx);
                    }

                    // Notify reader that a buffer is available
                    reader_cv.notify_one();
                }
            }
        };

        // Start writer threads
        std::vector<std::thread> writer_threads;
        for (int i = 0; i < num_writers; ++i) {
            writer_threads.emplace_back(writer_func);
        }

        // Get initial buffer from pool
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            reader_cv.wait(lock, [&]{ return !available_buffers.empty(); });
            current_buffer_idx = available_buffers.front();
            available_buffers.pop();
        }

        FastqRecordArena rec;
        while (!writer_failed.load() && reader.read(rec, chunk_pool[current_buffer_idx].arena)) {
            auto read_end = std::chrono::steady_clock::now();
            read_time_ms += std::chrono::duration_cast<std::chrono::milliseconds>(read_end - last_read_start).count();

            chunk_bytes += rec.memory_size();
            chunk_pool[current_buffer_idx].records.push_back(rec);
            total_reads++;

            last_read_start = std::chrono::steady_clock::now();

            // Queue chunk for writing when size limit reached
            if (chunk_bytes >= chunk_size_bytes_) {
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);

                    // Submit current buffer for writing
                    write_queue.push(current_buffer_idx);

                    // Wait for an available buffer (backpressure)
                    reader_cv.wait(lock, [&]{ return !available_buffers.empty() || writer_failed.load(); });
                    if (writer_failed.load()) {
                        break;
                    }

                    // Get next buffer from pool
                    current_buffer_idx = available_buffers.front();
                    available_buffers.pop();
                }
                queue_cv.notify_one();

                chunk_bytes = 0;
            }

            if (writer_failed.load()) {
                break;
            }

            if ((total_reads % 1000000) == 0) {
                double avg_write_speed = 0;
                {
                    std::lock_guard<std::mutex> lock(stats_mutex_);
                    if (total_write_time_ms_ > 0) {
                        avg_write_speed = (total_bytes_written_ / 1024.0 / 1024.0) / (total_write_time_ms_ / 1000.0);
                    }
                }

                size_t pending_writes = 0;
                {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    pending_writes = write_queue.size();
                }

                int num_active = active_writers.load();

                std::cerr << "\r[Phase 1] " << total_reads << " reads, "
                         << chunks_created_ << " chunks written, "
                         << pending_writes << " queued + " << num_active << " writing, "
                         << "avg: " << std::fixed << std::setprecision(0) << avg_write_speed << " MB/s" << std::flush;
            }
        }

        // Queue final chunk if it has data
        if (!writer_failed.load()) {
            if (!chunk_pool[current_buffer_idx].records.empty()) {
                std::unique_lock<std::mutex> lock(queue_mutex);
                write_queue.push(current_buffer_idx);
                queue_cv.notify_one();
            } else {
                // Return unused buffer
                std::unique_lock<std::mutex> lock(queue_mutex);
                available_buffers.push(current_buffer_idx);
            }
        }

        // Wait for all writers to finish
        done_reading = true;
        queue_cv.notify_all();  // Wake all writers
        for (auto& t : writer_threads) {
            t.join();
        }

        if (writer_exception) {
            std::rethrow_exception(writer_exception);
        }

        std::cerr << "\r" << std::string(100, ' ') << "\r" << std::flush;

        auto phase1_end = std::chrono::steady_clock::now();
        auto phase1_total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(phase1_end - phase1_start).count();

        double avg_write_speed = 0;
        size_t write_time_ms = 0;
        {
            std::lock_guard<std::mutex> lock(stats_mutex_);
            write_time_ms = total_write_time_ms_;
            if (total_write_time_ms_ > 0) {
                avg_write_speed = (total_bytes_written_ / 1024.0 / 1024.0) / (total_write_time_ms_ / 1000.0);
            }
        }

        double read_pct = (read_time_ms * 100.0) / phase1_total_ms;
        double write_pct = (write_time_ms * 100.0) / phase1_total_ms;

        log_info("Phase 1 complete: " + std::to_string(total_reads) + " reads, " +
                std::to_string(chunks_created_) + " chunks");

    }

    void write_sorted_chunk(ChunkBuffer&& chunk) {
        size_t chunk_id = chunks_created_++;
        // Temp files are always uncompressed — they are written, read once, and deleted.
        // Gzip here burns CPU with zero user benefit; final output compression is
        // determined by the output filename extension in merge_chunks.
        std::string sorted_file = temp_dir_ + "/chunk_" + std::to_string(chunk_id) + ".sorted.fq";

        {
            std::lock_guard<std::mutex> lock(sorted_mutex_);
            sorted_files_.push_back(sorted_file);
        }

        const uint32_t N = static_cast<uint32_t>(chunk.records.size());

        // uint32_t indices: 4 bytes vs 8 — fits 2× more entries per cache line
        std::vector<uint32_t> indices(N);
        for (uint32_t i = 0; i < N; ++i) indices[i] = i;

        if (natural_order_) {
            // Precompute (prefix_view, numeric_suffix) once — zero allocs, O(1) per comparison.
            // Keys point into chunk.arena which lives for the duration of this function.
            struct NatKey { std::string_view prefix; uint64_t suffix; bool has_suffix; };
            std::vector<NatKey> keys;
            keys.reserve(N);
            for (const auto& r : chunk.records) {
                std::string_view id = trim_id(r.get_header());
                size_t i = id.size();
                while (i > 0 && static_cast<unsigned>(id[i-1] - '0') <= 9u) --i;
                uint64_t suf = 0;
                for (size_t p = i; p < id.size(); ++p)
                    suf = suf * 10u + static_cast<unsigned>(id[p] - '0');
                keys.push_back({id.substr(0, i), suf, i < id.size()});
            }
            std::sort(indices.begin(), indices.end(), [&keys](uint32_t a, uint32_t b) {
                const NatKey& ka = keys[a]; const NatKey& kb = keys[b];
                if (ka.prefix != kb.prefix) return ka.prefix < kb.prefix;
                if (!ka.has_suffix && !kb.has_suffix) return false;
                if (!ka.has_suffix) return true;
                if (!kb.has_suffix) return false;
                return ka.suffix < kb.suffix;
            });
        } else {
            // string_view points into chunk.arena — zero allocs, comparison = fast memcmp
            std::vector<std::string_view> keys;
            keys.reserve(N);
            for (const auto& r : chunk.records)
                keys.push_back(trim_id(r.get_header()));
            std::sort(indices.begin(), indices.end(), [&keys](uint32_t a, uint32_t b) {
                return keys[a] < keys[b];
            });
        }

        // Write sorted chunk
        auto start = std::chrono::steady_clock::now();

        FastqWriter writer(sorted_file, false);
        for (uint32_t idx : indices) {
            writer.write(chunk.records[idx]);
        }

        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // Update cumulative stats
        {
            std::lock_guard<std::mutex> lock(stats_mutex_);
            total_write_time_ms_ += duration;
            total_bytes_written_ += chunk.records.size() * 300;  // Estimate
        }
    }


    void merge_chunks(const std::string& output) {
        if (sorted_files_.size() == 1) {
            if (std::rename(sorted_files_[0].c_str(), output.c_str()) != 0) {
                throw std::runtime_error("rename failed from " + sorted_files_[0] + " to " + output +
                                         ": " + std::strerror(errno));
            }
            sorted_files_.clear();
            return;
        }

        // Temp files are always uncompressed .fq — make_fastq_reader detects plain text.
        std::vector<std::unique_ptr<FastqReaderBase>> readers;
        for (const auto& f : sorted_files_) {
            readers.push_back(make_fastq_reader(f));
        }

        // Use lambda comparator to avoid static member in local struct
        auto compare_records = [this](const FastqRecord& a, const FastqRecord& b) -> bool {
            if (natural_order_) {
                return NaturalOrderComparator()(a, b);
            } else {
                return trim_id(a.header) < trim_id(b.header);
            }
        };

        struct MergeEntry {
            FastqRecord record;
            size_t idx;
        };

        auto compare_entries = [&compare_records](const MergeEntry& a, const MergeEntry& b) -> bool {
            return compare_records(b.record, a.record);  // Inverted for min-heap
        };

        std::priority_queue<MergeEntry, std::vector<MergeEntry>, decltype(compare_entries)> pq(compare_entries);

        for (size_t i = 0; i < readers.size(); ++i) {
            MergeEntry e;
            e.idx = i;
            if (readers[i]->read(e.record)) pq.push(std::move(e));
        }

        bool compress = (output.size() > 3 && output.substr(output.size() - 3) == ".gz");

        // For compressed output, try to use pigz for parallel compression
        FILE* pigz_pipe = nullptr;
        if (compress) {
            // Try to use pigz with multiple threads for faster compression
            std::string pigz_cmd = "pigz -c -p " + std::to_string(threads_) + " > " + shell_escape_path(output) + " 2>/dev/null";
            pigz_pipe = popen(pigz_cmd.c_str(), "w");
        }

        // Fall back to FastqWriter if pigz not available or uncompressed output
        std::unique_ptr<FastqWriter> writer;
        if (!pigz_pipe) {
            writer = std::make_unique<FastqWriter>(output, compress);
        }

        size_t merged = 0;
        while (!pq.empty()) {
            MergeEntry e = pq.top();
            pq.pop();

            if (pigz_pipe) {
                // Write directly to pigz pipe
                fprintf(pigz_pipe, "%s\n%s\n%s\n%s\n",
                       e.record.header.c_str(), e.record.seq.c_str(),
                       e.record.plus.c_str(), e.record.qual.c_str());
            } else {
                writer->write(e.record);
            }
            merged++;

            if (readers[e.idx]->read(e.record)) pq.push(std::move(e));

            if ((merged % 1000000) == 0) {
                std::cerr << "\r[Phase 2] Merged " << merged << " reads" << std::flush;
            }
        }

        if (pigz_pipe) {
            int rc = pclose(pigz_pipe);
            if (rc != 0) {
                throw std::runtime_error("pclose failed for pigz writer pipe with status " + std::to_string(rc));
            }
        }

        std::cerr << "\r" << std::string(100, ' ') << "\r" << std::flush;
        log_info("Phase 2 complete: " + std::to_string(merged) + " reads");
    }

    size_t max_memory_mb_;
    size_t chunk_size_bytes_;
    std::string temp_dir_;
    int threads_;
    bool natural_order_;
    bool fast_mode_;
    std::vector<std::string> sorted_files_;
    std::atomic<size_t> chunks_created_{0};

    // Stats for progress display
    std::mutex stats_mutex_;
    std::mutex sorted_mutex_;
    size_t total_write_time_ms_{0};
    size_t total_bytes_written_{0};
};

}  // namespace

static void print_usage(const char* prog) {
    std::cerr << "Usage: fqdup " << prog << " -i INPUT -o OUTPUT --max-memory SIZE [-p THREADS] [-t TMPDIR] [-N] [--fast]\n"
              << "\nOptimized FASTQ sorter with parallel chunk processing\n"
              << "  -i FILE            Input FASTQ (.gz supported)\n"
              << "  -o FILE            Output FASTQ (.gz for compression)\n"
              << "  --max-memory SIZE  Memory limit (e.g., 64G, 32G)\n"
              << "  -p THREADS         Number of sorting threads (default: auto)\n"
              << "  -t DIR             Temp directory (default: .)\n"
              << "  -N                 Natural order sorting (numeric suffixes)\n"
              << "  --fast             Fast mode: uncompressed intermediates (3× faster, more disk)\n"
              << "\nAlgorithm:\n"
              << "  Phase 1: Create sorted chunks (parallel writer threads)\n"
              << "  Phase 2: K-way merge (streaming)\n"
#ifdef HAVE_ISAL
              << "\nDecompression: ISA-L igzip (hardware-accelerated)\n"
#else
              << "\nDecompression: zlib\n"
#endif
              ;
}

int sort_main(int argc, char** argv) {
    // Configure jemalloc for aggressive memory release
    #ifdef __GLIBC__
    if (mallctl) {
        // Set dirty page decay to 1 second (default is 10 seconds)
        ssize_t decay_ms = 1000;
        mallctl("arenas.dirty_decay_ms", nullptr, nullptr, &decay_ms, sizeof(decay_ms));
        mallctl("arenas.muzzy_decay_ms", nullptr, nullptr, &decay_ms, sizeof(decay_ms));
    }
    #endif

    std::string input, output, temp_dir = ".";
    size_t max_memory_mb = 0;
    int threads = 0;
    bool natural_order = false;
    bool fast_mode = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "-i" && i + 1 < argc) input = argv[++i];
        else if (arg == "-o" && i + 1 < argc) output = argv[++i];
        else if (arg == "--max-memory" && i + 1 < argc) max_memory_mb = parse_memory_size(argv[++i]);
        else if (arg == "-p" && i + 1 < argc) threads = std::atoi(argv[++i]);
        else if (arg == "-t" && i + 1 < argc) temp_dir = argv[++i];
        else if (arg == "-N") natural_order = true;
        else if (arg == "--fast") fast_mode = true;
        else if (arg == "-h" || arg == "--help") { print_usage(argv[0]); return 0; }
    }

    if (input.empty() || output.empty() || max_memory_mb == 0) {
        std::cerr << "Error: Missing required arguments\n";
        print_usage(argv[0]);
        return 1;
    }

    if (threads <= 0) {
        threads = std::thread::hardware_concurrency();
        if (threads < 1) threads = 4;
    }

    std::string log_file = temp_dir + "/fqdup-sort.log";
    init_logger(log_file.c_str());

    log_info("Input: " + input);
    log_info("Output: " + output);
    log_info("Memory limit: " + std::to_string(max_memory_mb) + " MB (" +
            std::to_string(max_memory_mb / 1024.0) + " GB)");
    log_info("Threads: " + std::to_string(threads));
    log_info("Temp dir: " + temp_dir);
    log_info("Fast mode: " + std::string(fast_mode ? "enabled (uncompressed temps)" : "disabled"));
    log_info("Natural order: " + std::string(natural_order ? "enabled" : "disabled"));

    try {
        OptimizedSorter sorter(max_memory_mb, temp_dir, threads, natural_order, fast_mode);
        sorter.sort(input, output);
        log_info("=== Complete ===");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
