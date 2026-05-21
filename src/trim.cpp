// fqdup trim — detect and remove 5'/3' adapter stubs from collapsed FASTQ
//
// Pass 1 (pre-scan, --scan-reads): hexamer frequency analysis → detect stubs
// Pass 2 (all reads):              clip matching stubs from seq+qual, write output
//
// Designed for post-fastp collapsed reads where fastp trims P7 from R1 3' but
// misses short P5-tail remnants (e.g. CTCTTC) at the read 5' end.

#include "fqdup/fastq_common.hpp"
#include "taph/frame_selector_decl.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/sample_damage_profile.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

static constexpr int LSD_L_MIN = 30;

static void usage() {
    std::cerr <<
        "Usage: fqdup trim -i FILE -o FILE [options]\n\n"
        "Detect and remove 5'/3' adapter stubs from collapsed FASTQ.\n"
        "Uses hexamer frequency analysis on the first --scan-reads reads\n"
        "to detect stubs, then trims all reads in a single streaming pass.\n\n"
        "Required:\n"
        "  -i FILE          Input FASTQ (.gz or plain)\n"
        "  -o FILE          Output FASTQ (.gz)\n\n"
        "Options:\n"
        "  -p N             Threads for output compression (default: all cores)\n"
        "  --scan-reads N   Reads sampled for stub detection (default: 1000000; 0=all)\n"
        "  --min-length N   Discard trimmed reads shorter than N bp (default: 15)\n"
        "  --stub5 SEQ      Force 5' stub (skip detection)\n"
        "  --stub3 SEQ      Force 3' stub (skip detection)\n"
        "  -h, --help       Show this help\n";
}

int trim_main(int argc, char** argv) {
    std::string in_path, out_path;
    int64_t scan_reads  = 1'000'000;
    int     min_length  = 15;
    int     n_threads   = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    std::string forced_stub5, forced_stub3;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      ((a == "-i" || a == "--input")   && i+1 < argc) { in_path = argv[++i]; }
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

    if (!forced_stub5.empty() || !forced_stub3.empty()) {
        if (!forced_stub5.empty()) stubs5.push_back(forced_stub5);
        if (!forced_stub3.empty()) stubs3.push_back(forced_stub3);
        std::cerr << "Stubs (manual):";
        if (!stubs5.empty()) std::cerr << " 5'=" << stubs5[0];
        if (!stubs3.empty()) std::cerr << " 3'=" << stubs3[0];
        std::cerr << "\n";
    } else {
        std::cerr << "Scanning " << (scan_reads ? std::to_string(scan_reads) : "all")
                  << " reads for adapter stubs...\n";
        taph::SampleDamageProfile scan_profile;
        std::array<uint32_t, 4096> hex3_terminal{};
        uint64_t n_hex3 = 0;
        int64_t n = 0;

        auto rdr = make_fastq_reader(in_path, 1);
        FastqRecord rec;
        while (rdr->read(rec) && (scan_reads == 0 || n < scan_reads)) {
            int L = static_cast<int>(rec.seq.size());
            if (L < LSD_L_MIN) continue;
            taph::FrameSelector::update_sample_profile(scan_profile, rec.seq);
            if (L >= 12) {
                int code = taph::encode_hex_at(rec.seq, L - 6);
                if (code >= 0) { ++hex3_terminal[code]; ++n_hex3; }
            }
            ++n;
        }

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

    // ---- trimming pass -----------------------------------------------------
    int64_t n_in = 0, n_out = 0, n_clipped5 = 0, n_clipped3 = 0, n_dropped = 0;
    bool compress = out_path.size() >= 3 &&
                    out_path.compare(out_path.size() - 3, 3, ".gz") == 0;
    int write_threads = std::min(n_threads, 16);

    auto rdr = make_fastq_reader(in_path, 1);
    FastqWriter writer(out_path, compress, write_threads);
    FastqRecord rec;

    while (rdr->read(rec)) {
        ++n_in;

        // 5' clip — one stub length (6 bp)
        if (!stubs5.empty() && static_cast<int>(rec.seq.size()) >= 6) {
            for (const auto& s : stubs5) {
                if (rec.seq.compare(0, 6, s) == 0) {
                    rec.seq.erase(0, 6);
                    rec.qual.erase(0, 6);
                    ++n_clipped5;
                    break;
                }
            }
        }

        // 3' clip — iterative until no match or read too short
        if (!stubs3.empty()) {
            bool trimmed;
            do {
                trimmed = false;
                int L = static_cast<int>(rec.seq.size());
                if (L < 12) break;
                for (const auto& s : stubs3) {
                    if (rec.seq.compare(L - 6, 6, s) == 0) {
                        rec.seq.erase(L - 6, 6);
                        rec.qual.erase(L - 6, 6);
                        trimmed = true;
                        ++n_clipped3;
                        break;
                    }
                }
            } while (trimmed);
        }

        if (static_cast<int>(rec.seq.size()) < min_length) {
            ++n_dropped;
            continue;
        }

        writer.write(rec);
        ++n_out;
    }

    std::cerr << "Reads in:       " << n_in     << "\n"
              << "5' clipped:     " << n_clipped5 << " ("
              << (n_in ? 100.0 * n_clipped5 / n_in : 0.0) << "%)\n"
              << "3' clipped:     " << n_clipped3 << " ("
              << (n_in ? 100.0 * n_clipped3 / n_in : 0.0) << "%)\n"
              << "Dropped (<" << min_length << " bp): " << n_dropped << "\n"
              << "Reads out:      " << n_out    << "\n";

    return 0;
}
