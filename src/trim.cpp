// fqdup trim — detect and remove 5'/3' adapter stubs from collapsed FASTQ
//
// Pass 1 (pre-scan, --scan-reads): hexamer frequency analysis → detect stubs
// Pass 2 (all reads):              parallel clip pipeline → ordered BGZF output
//
// Pipeline: producer (reader) → bounded clip queue → N clip workers
//                                                  → ordered output queue → writer thread

#include "fqdup/fastq_common.hpp"
#include "taph/frame_selector_decl.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/sample_damage_profile.hpp"

#include <array>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

static constexpr int LSD_L_MIN = 30;
static constexpr int BATCH_SZ  = 4096;

// ---- batch -----------------------------------------------------------------
struct TrimBatch {
    uint64_t              id        = 0;
    std::vector<FastqRecord> records;
    int64_t               n_clipped5 = 0;
    int64_t               n_clipped3 = 0;
    int64_t               n_dropped  = 0;
};

// ---- bounded clip queue (producer → workers) --------------------------------
struct ClipQueue {
    std::mutex              mtx;
    std::condition_variable cv_ne, cv_nf;
    std::deque<TrimBatch>   q;
    bool                    done      = false;
    int                     max_depth;

    explicit ClipQueue(int d) : max_depth(d) {}

    void push(TrimBatch&& b) {
        std::unique_lock lk(mtx);
        cv_nf.wait(lk, [&]{ return (int)q.size() < max_depth || done; });
        q.push_back(std::move(b));
        cv_ne.notify_one();
    }

    bool pop(TrimBatch& b) {
        std::unique_lock lk(mtx);
        cv_ne.wait(lk, [&]{ return !q.empty() || done; });
        if (q.empty()) return false;
        b = std::move(q.front());
        q.pop_front();
        cv_nf.notify_one();
        return true;
    }

    void set_done() {
        std::unique_lock lk(mtx);
        done = true;
        cv_ne.notify_all();
        cv_nf.notify_all();
    }
};

// ---- ordered output queue (workers → writer thread) -------------------------
struct OutQueue {
    std::mutex                      mtx;
    std::condition_variable         cv;
    std::map<uint64_t, TrimBatch>   pending;
    bool                            done = false;

    void push(TrimBatch&& b) {
        std::unique_lock lk(mtx);
        pending.emplace(b.id, std::move(b));
        cv.notify_one();
    }

    // Pop batch with exactly `expected_id`; returns false when done+empty.
    bool pop_ordered(uint64_t expected_id, TrimBatch& out) {
        std::unique_lock lk(mtx);
        cv.wait(lk, [&]{ return pending.count(expected_id) || done; });
        auto it = pending.find(expected_id);
        if (it == pending.end()) return false;
        out = std::move(it->second);
        pending.erase(it);
        return true;
    }

    void set_done() {
        std::unique_lock lk(mtx);
        done = true;
        cv.notify_all();
    }
};

// ---- clip worker ------------------------------------------------------------
static void clip_worker(ClipQueue& in_q, OutQueue& out_q,
                        const std::vector<std::string>& stubs5,
                        const std::vector<std::string>& stubs3,
                        int min_length) {
    TrimBatch batch;
    while (in_q.pop(batch)) {
        for (auto& rec : batch.records) {
            if (!stubs5.empty() && (int)rec.seq.size() >= 6) {
                for (const auto& s : stubs5) {
                    if (rec.seq.compare(0, 6, s) == 0) {
                        rec.seq.erase(0, 6);
                        rec.qual.erase(0, 6);
                        ++batch.n_clipped5;
                        break;
                    }
                }
            }
            if (!stubs3.empty()) {
                bool trimmed;
                do {
                    trimmed = false;
                    int L = (int)rec.seq.size();
                    if (L < 12) break;
                    for (const auto& s : stubs3) {
                        if (rec.seq.compare(L - 6, 6, s) == 0) {
                            rec.seq.erase(L - 6, 6);
                            rec.qual.erase(L - 6, 6);
                            trimmed = true;
                            ++batch.n_clipped3;
                            break;
                        }
                    }
                } while (trimmed);
            }
            if ((int)rec.seq.size() < min_length) {
                ++batch.n_dropped;
                rec.seq.clear();  // sentinel: writer skips empty-seq records
            }
        }
        out_q.push(std::move(batch));
    }
}

// ---- usage -----------------------------------------------------------------
static void usage() {
    std::cerr <<
        "Usage: fqdup trim -i FILE -o FILE [options]\n\n"
        "Detect and remove 5'/3' adapter stubs from collapsed FASTQ.\n"
        "Uses hexamer frequency analysis on the first --scan-reads reads\n"
        "to detect stubs, then trims all reads via a parallel clip pipeline.\n\n"
        "Required:\n"
        "  -i FILE          Input FASTQ (.gz or plain)\n"
        "  -o FILE          Output FASTQ (.gz)\n\n"
        "Options:\n"
        "  -p N             Threads (default: all cores)\n"
        "  --scan-reads N   Reads sampled for stub detection (default: 1000000; 0=all)\n"
        "  --min-length N   Discard trimmed reads shorter than N bp (default: 15)\n"
        "  --stub5 SEQ      Force 5' stub (skip detection)\n"
        "  --stub3 SEQ      Force 3' stub (skip detection)\n"
        "  -h, --help       Show this help\n";
}

// ---- main ------------------------------------------------------------------
int trim_main(int argc, char** argv) {
    std::string in_path, out_path;
    int64_t scan_reads = 1'000'000;
    int     min_length = 15;
    int     n_threads  = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    std::string forced_stub5, forced_stub3;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      ((a == "-i" || a == "--input")   && i+1 < argc) { in_path  = argv[++i]; }
        else if ((a == "-o" || a == "--output")  && i+1 < argc) { out_path = argv[++i]; }
        else if ((a == "-p" || a == "--threads") && i+1 < argc) { n_threads = std::stoi(argv[++i]); }
        else if (a == "--scan-reads"             && i+1 < argc) { scan_reads = std::stoll(argv[++i]); }
        else if (a == "--min-length"             && i+1 < argc) { min_length = std::stoi(argv[++i]); }
        else if (a == "--stub5"                  && i+1 < argc) { forced_stub5 = argv[++i]; }
        else if (a == "--stub3"                  && i+1 < argc) { forced_stub3 = argv[++i]; }
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else { std::cerr << "Error: unknown option '" << a << "'\n"; usage(); return 1; }
    }

    if (in_path.empty() || out_path.empty()) {
        std::cerr << "Error: -i and -o are required\n";
        usage();
        return 1;
    }

    // ---- stub detection or manual override ---------------------------------
    std::vector<std::string> stubs5, stubs3;

    // Single reader opened before detection and kept alive for the trim pass.
    // scan_buf holds records buffered during hexamer analysis (single-pass —
    // no second file open needed). Empty in the forced-stubs branch.
    std::unique_ptr<FastqReaderBase> rdr;
    std::vector<FastqRecord> scan_buf;

    if (!forced_stub5.empty() || !forced_stub3.empty()) {
        if (!forced_stub5.empty()) stubs5.push_back(forced_stub5);
        if (!forced_stub3.empty()) stubs3.push_back(forced_stub3);
        std::cerr << "Stubs (manual):";
        if (!stubs5.empty()) std::cerr << " 5'=" << stubs5[0];
        if (!stubs3.empty()) std::cerr << " 3'=" << stubs3[0];
        std::cerr << "\n";
        rdr = make_fastq_reader(in_path, static_cast<size_t>(n_threads));
    } else {
        std::cerr << "Scanning " << (scan_reads ? std::to_string(scan_reads) : "all")
                  << " reads for adapter stubs...\n";
        taph::SampleDamageProfile scan_profile;
        std::array<uint32_t, 4096> hex3_terminal{};
        uint64_t n_hex3 = 0;
        int64_t  n      = 0;

        if (scan_reads > 0)
            scan_buf.reserve(static_cast<size_t>(scan_reads));
        rdr = make_fastq_reader(in_path, static_cast<size_t>(n_threads));
        FastqRecord rec;
        while (rdr->read(rec) && (scan_reads == 0 || n < scan_reads)) {
            int L = static_cast<int>(rec.seq.size());
            if (L >= LSD_L_MIN)
                taph::FrameSelector::update_sample_profile(scan_profile, rec.seq);
            if (L >= 12) {
                int code = taph::encode_hex_at(rec.seq, L - 6);
                if (code >= 0) { ++hex3_terminal[code]; ++n_hex3; }
            }
            ++n;
            scan_buf.push_back(std::move(rec));
        }
        // rdr stays open — remaining records consumed in producer below.

        auto detected = taph::detect_adapter_stubs(scan_profile, hex3_terminal.data(), n_hex3);
        stubs5 = detected.stubs5;
        stubs3 = detected.stubs3;

        if (stubs5.empty() && stubs3.empty()) {
            std::cerr << "No adapter stubs detected — output will be identical to input.\n";
        } else {
            std::cerr << "Detected stubs:";
            for (const auto& s : stubs5) std::cerr << " 5'=" << s;
            for (const auto& s : stubs3) std::cerr << " 3'=" << s;
            std::cerr << "\n";
        }
    }

    // ---- parallel trim pass ------------------------------------------------
    int clip_threads  = std::max(1, n_threads - 1);  // leave 1 for writer
    int write_threads = std::min(n_threads, 16);

    ClipQueue clip_q(2 * clip_threads);
    OutQueue  out_q;

    // Writer thread: drains out_q in batch-id order, writes to BGZF output.
    int64_t n_in = 0, n_out = 0, n_clipped5 = 0, n_clipped3 = 0, n_dropped = 0;
    uint64_t total_batches = 0;  // set by producer after pushing all batches

    bool compress = out_path.size() >= 3 &&
                    out_path.compare(out_path.size() - 3, 3, ".gz") == 0;
    FastqWriter writer(out_path, compress, write_threads);

    std::thread writer_thread([&] {
        uint64_t expected = 0;
        TrimBatch batch;
        // Loop until we've consumed all batches (total_batches set by producer).
        // total_batches is written before out_q.set_done(), so the memory order
        // is guaranteed via the mutex inside set_done/pop_ordered.
        while (true) {
            if (!out_q.pop_ordered(expected, batch)) break;
            for (auto& rec : batch.records) {
                ++n_in;
                if (rec.seq.empty()) continue;  // dropped
                writer.write(rec);
                ++n_out;
            }
            n_clipped5 += batch.n_clipped5;
            n_clipped3 += batch.n_clipped3;
            n_dropped  += batch.n_dropped;
            ++expected;
        }
    });

    // Worker threads: clip batches in parallel.
    std::vector<std::thread> workers;
    workers.reserve(clip_threads);
    for (int t = 0; t < clip_threads; ++t)
        workers.emplace_back(clip_worker, std::ref(clip_q), std::ref(out_q),
                             std::cref(stubs5), std::cref(stubs3), min_length);

    // Producer (main thread): drain scan_buf first (buffered during detection),
    // then continue reading remaining records from the already-open reader.
    {
        FastqRecord rec;
        std::vector<FastqRecord> batch;
        batch.reserve(BATCH_SZ);
        uint64_t batch_id = 0;

        for (auto& r : scan_buf) {
            batch.push_back(std::move(r));
            if ((int)batch.size() == BATCH_SZ) {
                TrimBatch tb;
                tb.id      = batch_id++;
                tb.records = std::move(batch);
                clip_q.push(std::move(tb));
                batch.clear();
                batch.reserve(BATCH_SZ);
            }
        }
        scan_buf.clear();
        scan_buf.shrink_to_fit();

        while (rdr->read(rec)) {
            batch.push_back(std::move(rec));
            if ((int)batch.size() == BATCH_SZ) {
                TrimBatch tb;
                tb.id      = batch_id++;
                tb.records = std::move(batch);
                clip_q.push(std::move(tb));
                batch.clear();
                batch.reserve(BATCH_SZ);
            }
        }
        if (!batch.empty()) {
            TrimBatch tb;
            tb.id      = batch_id++;
            tb.records = std::move(batch);
            clip_q.push(std::move(tb));
        }
        total_batches = batch_id;
    }

    clip_q.set_done();
    for (auto& w : workers) w.join();
    out_q.set_done();
    writer_thread.join();

    std::cerr << "Reads in:         " << n_in      << "\n"
              << "5' clipped:       " << n_clipped5 << " ("
              << (n_in ? 100.0 * n_clipped5 / n_in : 0.0) << "%)\n"
              << "3' clipped:       " << n_clipped3 << " ("
              << (n_in ? 100.0 * n_clipped3 / n_in : 0.0) << "%)\n"
              << "Dropped (<" << min_length << " bp): " << n_dropped  << "\n"
              << "Reads out:        " << n_out      << "\n";

    return 0;
}
