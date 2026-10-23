// fqdup profile: argument parsing and FASTQ input for the libtaph library profiler
// (taph::run_library_profile), which does all of the estimation and writes the outputs.
// `dart profile` runs the same profiler, so both give byte-identical results.

#include "fqdup/fastq_common.hpp"
#include "fqdup/damage_profile.hpp"
#include "fqdup/damage_model_io.hpp"
#include "fqdup/logger.hpp"
#include "fqdup/version.hpp"

#include "taph/library_profile.hpp"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: fqdup " << prog << " (-i FILE | -1 R1.fq -2 R2.fq) [options]\n"
        << "\nInput modes (mutually exclusive):\n"
        << "  -i FILE                    Single-end / merged input FASTQ (.gz or plain)\n"
        << "  -1 FILE  -2 FILE           Paired-end raw reads (un-merged); R1 contributes\n"
        << "                             5' counters, R2 (complement-mapped) contributes\n"
        << "                             3' counters. Read 3' ends ignored (adapter zone).\n"
        << "\nOptions:\n"
        << "  -p N                       Worker threads (default: all cores)\n"
        << "  --library-type auto|ds|ss  Library type for 3'-end interpretation (default: auto;\n"
        << "                             when omitted, inherits a merge-declared type from --subst-in)\n"
        << "  --mask-threshold FLOAT     Mask when excess P(deam) > T (default: 0.05)\n"
        << "  --tsv FILE                 Write per-position table as TSV\n"
        << "  --json FILE                Write full damage profile as JSON\n"
        << "  --json-level LEVEL         summary | standard | full (default: standard)\n"
        << "  --count-table-json FILE    Write Layer-0 stop-channel count tables as JSONL (golden-gate source of truth)\n"
        << "  --damage-summary FILE      Write the taph DamageSummary sidecar (.tdmg). This is the\n"
        << "                             profile DART consumes via --library-profile: it takes the fit\n"
        << "                             VERBATIM (pi.state, pi.point, both end shapes) instead of\n"
        << "                             re-deriving it, so the damage verdict is made once, here.\n"
        << "  --html FILE                Write interactive damage report as self-contained HTML\n"
        << "  --length-bins SPEC         Length-stratified damage: auto | off | N | e1,e2,... (default: auto)\n"
        << "  --adapter-scan-reads N     Reads sampled (single-thread) for adapter-stub detection\n"
        << "                             (default: 1000000; 0 = scan all reads)\n"
        << "  --min-length N             Floor for reads entering the estimator. N<30 lowers the\n"
        << "                             short-fragment floor to max(16,N), including 16-29 bp reads in a\n"
        << "                             JOINT length-isotonic re-fit (>=30 bulk delta/tau/preservation\n"
        << "                             then reflect ALL lengths, NOT purely additive). N>=30 keeps the\n"
        << "                             legacy >=30-only path (default: 0 = off)\n"
        << "  --max-length N             Drop reads longer than N from the estimator (default: 0 = off)\n"
        << "  --subst-in FILE[,FILE]     fqdup merge .bsubst files (comma-separated); also supplies the\n"
        << "                             library type when --library-type is not given\n"
        << "  --no-oxog                  Skip the 8-oxoG channels\n";
}

namespace {
class FastqSeqReader : public taph::SeqReader {
public:
    explicit FastqSeqReader(std::unique_ptr<FastqReaderBase> r) : reader_(std::move(r)) {}
    bool next(std::string& seq) override {
        if (!reader_->read(rec_)) return false;
        seq = std::move(rec_.seq);
        return true;
    }
private:
    std::unique_ptr<FastqReaderBase> reader_;
    FastqRecord rec_;
};
} // namespace

int profile_main(int argc, char** argv) {
    taph::LibraryProfileOptions o;
    o.n_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (o.n_threads < 1) o.n_threads = 1;
    std::string damage_json_out;   // scalar damage model for `split --damage-json`
    std::string model_bin_out;     // full split model for `split --model-bin`

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            o.in_path = argv[++i];
        } else if ((arg == "-1" || arg == "--read1") && i + 1 < argc) {
            o.r1_path = argv[++i];
        } else if ((arg == "-2" || arg == "--read2") && i + 1 < argc) {
            o.r2_path = argv[++i];
        } else if ((arg == "-p" || arg == "--threads") && i + 1 < argc) {
            o.n_threads = std::stoi(argv[++i]);
            if (o.n_threads < 1) o.n_threads = 1;
        } else if (arg == "--library-type" && i + 1 < argc) {
            std::string lt(argv[++i]);
            if (lt == "ss" || lt == "single-stranded")
                o.forced_lib = taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (lt == "ds" || lt == "double-stranded")
                o.forced_lib = taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (lt != "auto") {
                std::cerr << "Error: unknown --library-type: " << lt << "\n";
                return 1;
            }
        } else if (arg == "--mask-threshold" && i + 1 < argc) {
            o.mask_threshold = std::stod(argv[++i]);
        } else if (arg == "--tsv" && i + 1 < argc) {
            o.tsv_path = argv[++i];
        } else if (arg == "--json" && i + 1 < argc) {
            o.json_path = argv[++i];
        } else if (arg == "--damage-summary" && i + 1 < argc) {
            o.damage_summary_path = argv[++i];
        } else if (arg == "--damage-json-out" && i + 1 < argc) {
            damage_json_out = argv[++i];
        } else if (arg == "--model-bin-out" && i + 1 < argc) {
            model_bin_out = argv[++i];
        } else if (arg == "--json-level" && i + 1 < argc) {
            std::string lvl(argv[++i]);
            if (lvl == "summary")       o.json_level = 1;
            else if (lvl == "standard") o.json_level = 2;
            else if (lvl == "full")     o.json_level = 3;
            else {
                std::cerr << "Error: unknown --json-level: " << lvl
                          << " (choices: summary, standard, full)\n";
                return 1;
            }
        } else if (arg == "--count-table-json" && i + 1 < argc) {
            o.count_table_path = argv[++i];
        } else if (arg == "--html" && i + 1 < argc) {
            o.html_path = argv[++i];
        } else if (arg == "--subst-in" && i + 1 < argc) {
            std::stringstream ss_si(argv[++i]); std::string t_si;  // one flag, comma-separated list
            while (std::getline(ss_si, t_si, ',')) if (!t_si.empty()) o.subst_in_paths.push_back(t_si);
        } else if (arg == "--no-oxog") {
            o.run_oxog = false;
        } else if (arg == "--min-length" && i + 1 < argc) {
            o.min_length = std::stoi(argv[++i]);
        } else if (arg == "--max-length" && i + 1 < argc) {
            o.max_length = std::stoi(argv[++i]);
        } else if (arg == "--adapter-scan-reads" && i + 1 < argc) {
            long long v = std::stoll(argv[++i]);
            if (v < 0) {
                std::cerr << "Error: --adapter-scan-reads must be >= 0 (0 = scan all), got " << v << "\n";
                return 1;
            }
            o.adapter_scan_reads = static_cast<int64_t>(v);
        } else if (arg == "--length-bins" && i + 1 < argc) {
            bool ok = true;
            o.lb_opts = taph::parse_length_bins(argv[++i], ok);
            if (!ok) {
                std::cerr << "Error: invalid --length-bins: " << argv[i] << "\n";
                return 1;
            }
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Error: unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (o.r1_path.empty() != o.r2_path.empty()) {
        std::cerr << "Error: -1 and -2 must both be provided for paired input\n";
        print_usage(argv[0]);
        return 1;
    }
    if (o.in_path.empty() && o.r1_path.empty()) {
        std::cerr << "Error: -i FILE and/or -1/-2 required\n";
        print_usage(argv[0]);
        return 1;
    }
    o.producer = std::string("fqdup/") + FQDUP_VERSION;
    o.log_info = log_info;
    o.log_warn = log_warn;

    auto open_reader = [](const std::string& path, size_t threads) -> std::unique_ptr<taph::SeqReader> {
        return std::make_unique<FastqSeqReader>(make_fastq_reader(path, threads));
    };
    // Damage-model export for `fqdup split`: split loads these instead of re-deriving them.
    auto export_model = [&](const taph::ProfileModelExportInput& in) -> int {
        if (damage_json_out.empty() && model_bin_out.empty()) return 0;
        if (in.paired_mode) {
            std::cerr << "Warning: --damage-json-out/--model-bin-out require single-end (-i merged); "
                         "skipping in paired mode.\n";
            return 0;
        }
        const std::string in_digest = fqdup_model_io::compute_input_digest(in.in_path);
        if (!damage_json_out.empty()) {
            fqdup_model_io::write_damage_model_json(in.bulk_dp, damage_json_out, in.reads_scanned, in_digest);
            std::cout << "Damage model (JSON): " << damage_json_out << "\n";
        }
        if (!model_bin_out.empty()) {
            // Build the scorer the SAME way `split` does, so the binary IS what split would compute.
            LengthStratifiedDamageProfile sm_lsd = estimate_damage_split_model(
                in.in_path, in.bulk_dp, in.merged_lsd_hist.empty() ? nullptr : &in.merged_lsd_hist);
            DamageSplitModel sm = DamageSplitModel::build(sm_lsd, in.bulk_dp);
            fqdup_model_io::write_split_model_bin(sm, model_bin_out, in.reads_scanned, in_digest);
            std::cout << "Split model (binary): " << model_bin_out
                      << " (" << sm.bins.size() << " bins)\n";
        }
        return 0;
    };
    return taph::run_library_profile(o, open_reader, export_model);
}
