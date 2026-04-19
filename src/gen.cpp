// fqdup gen — synthetic FASTQ generator with configurable ancient DNA damage patterns
//
// Generates IID reads (controlled GC content) with:
//   Channel A: C→T deamination at 5' end (+ G→A at 3' for DS, + C→T at 3' for SS)
//   Channel D: 8-oxoG G→T transversions, uniform across read
//   Channel E: depurination — purine enrichment at pos 0
//
// Ancient and modern reads are mixed at a specified fraction.
// Background substitution rates apply to all reads (both classes).
//
// Usage: fqdup gen -o OUTPUT [options]

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <zlib.h>

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " gen -o OUTPUT [options]\n"
        << "\nOutput:\n"
        << "  -o FILE                    Output FASTQ (.gz or plain)\n"
        << "\nSequence:\n"
        << "  -n, --reads N              Number of reads (default: 1000000)\n"
        << "  --seed N                   PRNG seed (default: 1)\n"
        << "  --read-len N               Fixed read length (default: 60)\n"
        << "  --gc FLOAT                 GC content fraction (default: 0.45)\n"
        << "\nMixture:\n"
        << "  --f-damaged FLOAT          Fraction of ancient reads (default: 0.70)\n"
        << "  --lib-type ds|ss           Library type (default: ds)\n"
        << "\nDeamination (Channel A):\n"
        << "  --dmg5-max FLOAT           5' C→T amplitude at pos 0 (default: 0.18)\n"
        << "  --dmg5-lambda FLOAT        5' C→T decay rate (default: 0.35)\n"
        << "  --dmg3-ga-max FLOAT        3' G→A amplitude (default: same as dmg5-max)\n"
        << "  --dmg3-ga-lambda FLOAT     3' G→A decay rate (default: same as dmg5-lambda)\n"
        << "  --dmg3-ct-max FLOAT        3' C→T amplitude, SS only (default: 0 ds / dmg5-max ss)\n"
        << "  --dmg3-ct-lambda FLOAT     3' C→T decay rate (default: same as dmg5-lambda)\n"
        << "\n8-oxoG (Channel D):\n"
        << "  --oxog FLOAT               Per-G G→T probability (default: 0.0)\n"
        << "\nDepurination (Channel E):\n"
        << "  --depur FLOAT              Rate of pos-0 purine enrichment in ancient reads (default: 0.0)\n"
        << "\nBackground:\n"
        << "  --bg-ct FLOAT              Background C→T rate, all reads (default: 0.001)\n"
        << "  --bg-ga FLOAT              Background G→A rate, all reads (default: 0.001)\n"
        << "\nQuality:\n"
        << "  --q-score N                Constant Phred quality (default: 40)\n";
}

int gen_main(int argc, char** argv) {
    std::string out_path;
    uint64_t    n_reads      = 1000000;
    uint64_t    seed         = 1;
    int         read_len     = 60;
    double      gc           = 0.45;
    double      f_damaged    = 0.70;
    bool        is_ss        = false;
    double      dmg5_max     = 0.18;
    double      dmg5_lam     = 0.35;
    double      dmg3_ga_max  = -1.0;  // -1 = auto
    double      dmg3_ga_lam  = -1.0;
    double      dmg3_ct_max  = -1.0;
    double      dmg3_ct_lam  = -1.0;
    double      oxog_p_gt    = 0.0;
    double      depur_p      = 0.0;
    double      bg_ct        = 0.001;
    double      bg_ga        = 0.001;
    int         q_score      = 40;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        auto next_str = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Error: missing value for " << arg << "\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if (arg == "-o" || arg == "--output")      out_path      = next_str();
        else if (arg == "-n" || arg == "--reads") {
            auto s = next_str();
            long long v = std::stoll(s);
            if (v < 0) { std::cerr << "Error: --reads must be >= 0, got " << s << "\n"; return 1; }
            n_reads = static_cast<uint64_t>(v);
        }
        else if (arg == "--seed")                  seed          = std::stoull(next_str());
        else if (arg == "--read-len")              read_len      = std::stoi(next_str());
        else if (arg == "--gc")                    gc            = std::stod(next_str());
        else if (arg == "--f-damaged")             f_damaged     = std::stod(next_str());
        else if (arg == "--lib-type") {
            std::string lt = next_str();
            if (lt == "ss" || lt == "single-stranded") is_ss = true;
            else if (lt != "ds" && lt != "double-stranded") {
                std::cerr << "Error: unknown --lib-type: " << lt << "\n"; return 1;
            }
        }
        else if (arg == "--dmg5-max")              dmg5_max      = std::stod(next_str());
        else if (arg == "--dmg5-lambda")           dmg5_lam      = std::stod(next_str());
        else if (arg == "--dmg3-ga-max")           dmg3_ga_max   = std::stod(next_str());
        else if (arg == "--dmg3-ga-lambda")        dmg3_ga_lam   = std::stod(next_str());
        else if (arg == "--dmg3-ct-max")           dmg3_ct_max   = std::stod(next_str());
        else if (arg == "--dmg3-ct-lambda")        dmg3_ct_lam   = std::stod(next_str());
        else if (arg == "--oxog")                  oxog_p_gt     = std::stod(next_str());
        else if (arg == "--depur")                 depur_p       = std::stod(next_str());
        else if (arg == "--bg-ct")                 bg_ct         = std::stod(next_str());
        else if (arg == "--bg-ga")                 bg_ga         = std::stod(next_str());
        else if (arg == "--q-score")               q_score       = std::stoi(next_str());
        else if (arg == "-h" || arg == "--help")  { print_usage(argv[0]); return 0; }
        else { std::cerr << "Error: unknown argument: " << arg << "\n"; return 1; }
    }

    if (out_path.empty()) {
        std::cerr << "Error: -o OUTPUT required\n";
        print_usage(argv[0]);
        return 1;
    }
    if (read_len < 1) {
        std::cerr << "Error: --read-len must be >= 1, got " << read_len << "\n";
        return 1;
    }
    if (q_score < 0 || q_score > 93) {
        std::cerr << "Error: --q-score must be in [0, 93], got " << q_score << "\n";
        return 1;
    }

    // Validate probability parameters
    auto check01 = [&](double v, const char* name) -> bool {
        if (v < 0.0 || v > 1.0) {
            std::cerr << "Error: " << name << " must be in [0, 1], got " << v << "\n";
            return false;
        }
        return true;
    };
    if (!check01(gc,         "--gc"))         return 1;
    if (!check01(f_damaged,  "--f-damaged"))  return 1;
    if (!check01(bg_ct,      "--bg-ct"))      return 1;
    if (!check01(bg_ga,      "--bg-ga"))      return 1;
    if (!check01(oxog_p_gt,  "--oxog"))       return 1;
    if (!check01(depur_p,    "--depur"))      return 1;
    if (!check01(dmg5_max,   "--dmg5-max"))   return 1;
    if (dmg3_ga_max >= 0.0 && !check01(dmg3_ga_max, "--dmg3-ga-max")) return 1;
    if (dmg3_ct_max >= 0.0 && !check01(dmg3_ct_max, "--dmg3-ct-max")) return 1;
    if (dmg5_lam <= 0.0) {
        std::cerr << "Error: --dmg5-lambda must be > 0, got " << dmg5_lam << "\n";
        return 1;
    }

    // Resolve auto parameters
    if (dmg3_ga_max < 0.0) dmg3_ga_max = dmg5_max;
    if (dmg3_ga_lam < 0.0) dmg3_ga_lam = dmg5_lam;
    if (dmg3_ct_max < 0.0) dmg3_ct_max = is_ss ? dmg5_max : 0.0;
    if (dmg3_ct_lam < 0.0) dmg3_ct_lam = dmg5_lam;

    // Validate resolved lambdas (all must be > 0 after auto-resolution)
    if (dmg3_ga_lam <= 0.0) {
        std::cerr << "Error: --dmg3-ga-lambda must be > 0, got " << dmg3_ga_lam << "\n";
        return 1;
    }
    if (dmg3_ct_lam <= 0.0) {
        std::cerr << "Error: --dmg3-ct-lambda must be > 0, got " << dmg3_ct_lam << "\n";
        return 1;
    }

    // Base composition: pi_A = pi_T = (1-gc)/2,  pi_C = pi_G = gc/2
    double pi_at = (1.0 - gc) / 2.0;
    double pi_cg = gc / 2.0;
    // Cumulative boundaries for IID base draw: A, C, G, T
    double cum_a = pi_at;
    double cum_c = pi_at + pi_cg;
    double cum_g = pi_at + pi_cg + pi_cg;

    // Depurination: among purines {A, G}, choose A with probability pi_A / (pi_A + pi_G)
    double p_depur_to_a = pi_at / (pi_at + pi_cg);

    char q_char = static_cast<char>(q_score + 33);
    std::string qual(read_len, q_char);

    // Open output
    bool   gz_out = out_path.size() > 3 && out_path.substr(out_path.size() - 3) == ".gz";
    gzFile gzout  = nullptr;
    FILE*  fout   = nullptr;
    if (gz_out) {
        gzout = gzopen(out_path.c_str(), "wb6");
        if (!gzout) { std::cerr << "Error: cannot open " << out_path << "\n"; return 1; }
    } else {
        fout = std::fopen(out_path.c_str(), "w");
        if (!fout) { std::cerr << "Error: cannot open " << out_path << "\n"; return 1; }
    }

    auto write_buf = [&](const char* s, int n) {
        if (gz_out) {
            if (gzwrite(gzout, s, n) != n)
                throw std::runtime_error("gzwrite failed");
        } else {
            if (static_cast<int>(std::fwrite(s, 1, n, fout)) != n)
                throw std::runtime_error("fwrite failed");
        }
    };

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    std::string seq(read_len, 'A');
    char header_buf[64];

    for (uint64_t r = 0; r < n_reads; ++r) {
        bool damaged = (U(rng) < f_damaged);

        // Generate IID sequence
        for (int i = 0; i < read_len; ++i) {
            double u = U(rng);
            if      (u < cum_a) seq[i] = 'A';
            else if (u < cum_c) seq[i] = 'C';
            else if (u < cum_g) seq[i] = 'G';
            else                seq[i] = 'T';
        }

        // Depurination: with probability depur_p, replace pyrimidine at pos 0 with purine
        if (damaged && depur_p > 0.0 && (seq[0] == 'C' || seq[0] == 'T')) {
            if (U(rng) < depur_p)
                seq[0] = (U(rng) < p_depur_to_a) ? 'A' : 'G';
        }

        // Per-position damage
        for (int i = 0; i < read_len; ++i) {
            int d5 = i;
            int d3 = read_len - 1 - i;

            if (seq[i] == 'C') {
                if (U(rng) < bg_ct) { seq[i] = 'T'; continue; }
                if (damaged && dmg5_max > 0.0 && U(rng) < dmg5_max * std::exp(-dmg5_lam * d5))
                    { seq[i] = 'T'; continue; }
                if (damaged && is_ss && dmg3_ct_max > 0.0 && U(rng) < dmg3_ct_max * std::exp(-dmg3_ct_lam * d3))
                    { seq[i] = 'T'; continue; }
            } else if (seq[i] == 'G') {
                if (U(rng) < bg_ga) { seq[i] = 'A'; continue; }
                if (damaged && dmg3_ga_max > 0.0 && U(rng) < dmg3_ga_max * std::exp(-dmg3_ga_lam * d3))
                    { seq[i] = 'A'; continue; }
                if (damaged && oxog_p_gt > 0.0 && U(rng) < oxog_p_gt)
                    { seq[i] = 'T'; continue; }
            }
        }

        // Write FASTQ record
        int hlen = std::snprintf(header_buf, sizeof(header_buf),
                                 "@SYN:%llu:%s\n", (unsigned long long)r,
                                 damaged ? "damaged" : "undamaged");
        write_buf(header_buf, hlen);
        write_buf(seq.c_str(), read_len);
        write_buf("\n+\n", 3);
        write_buf(qual.c_str(), read_len);
        write_buf("\n", 1);
    }

    if (gz_out) {
        if (gzclose(gzout) != Z_OK)
            throw std::runtime_error("gzclose failed");
    } else {
        if (std::fclose(fout) != 0)
            throw std::runtime_error("fclose failed");
    }

    std::cout << "Generated " << n_reads << " reads → " << out_path << "\n"
              << "  f_damaged=" << f_damaged
              << "  lib=" << (is_ss ? "ss" : "ds")
              << "  d5_max=" << dmg5_max << "  lambda=" << dmg5_lam
              << "  d3_ga_max=" << dmg3_ga_max
              << "  oxog=" << oxog_p_gt
              << "  depur=" << depur_p
              << "\n";
    return 0;
}
