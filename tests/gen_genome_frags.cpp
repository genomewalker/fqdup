// gen_genome_frags — Simulate aDNA reads as overlapping fragments of a reference genome.
//
// Unlike gen_synthetic (isolated molecules), reads here are drawn from a shared
// reference so longer reads overlap shorter ones and provide k-mers that extend
// beyond the shorter reads' termini.  This is the correct model for testing
// fqdup extend.
//
// Usage:
//   gen_genome_frags [options] > out.fq
//
// Options:
//   --genome-len N     Reference genome length    (default: 2000)
//   --n-reads   N      Number of reads            (default: 50000)
//   --min-len   N      Min fragment length        (default: 30)
//   --max-len   N      Max fragment length        (default: 80)
//   --dmax5     F      C→T amplitude at 5' end    (default: 0.0)
//   --dmax3     F      G→A amplitude at 3' end    (default: 0.0)
//   --lambda5   F      5' decay rate              (default: 0.35)
//   --lambda3   F      3' decay rate              (default: 0.35)
//   --qual      N      Phred base quality         (default: 37)
//   --seed      N      RNG seed                   (default: 42)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// PCG-32 PRNG (same as gen_synthetic)
// ---------------------------------------------------------------------------
struct PCG32 {
    uint64_t state, inc;
    explicit PCG32(uint64_t seed = 42, uint64_t stream = 1)
        : state(0), inc((stream << 1u) | 1u) {
        next(); state += seed; next();
    }
    uint32_t next() {
        uint64_t old = state;
        state = old * 6364136223846793005ULL + inc;
        uint32_t xsh = uint32_t(((old >> 18u) ^ old) >> 27u);
        uint32_t rot = uint32_t(old >> 59u);
        return (xsh >> rot) | (xsh << ((-rot) & 31u));
    }
    double uniform() { return (next() >> 8) * (1.0 / (1u << 24)); }
    // Uniform integer in [0, n)
    uint32_t range(uint32_t n) { return next() % n; }
};

static constexpr char BASES[4] = {'A', 'C', 'G', 'T'};
static constexpr char COMP[4]  = {'T', 'G', 'C', 'A'};  // complement of A,C,G,T

static int base_index(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        default:            return 3;
    }
}

static std::string revcomp(const std::string& s) {
    std::string rc(s.size(), 'N');
    for (int i = 0; i < (int)s.size(); ++i)
        rc[s.size() - 1 - i] = COMP[base_index(s[i])];
    return rc;
}

static std::string apply_damage(const std::string& seq,
                                double dmax5, double lambda5,
                                double dmax3, double lambda3,
                                PCG32& rng) {
    std::string s = seq;
    int L = (int)s.size();
    for (int p = 0; p < std::min(15, L); ++p) {
        double prob = dmax5 * std::exp(-lambda5 * p);
        if (s[p] == 'C' && rng.uniform() < prob) s[p] = 'T';
    }
    for (int p = 0; p < std::min(15, L); ++p) {
        int dist = L - 1 - p;
        if (dist < 0) break;
        double prob = dmax3 * std::exp(-lambda3 * p);
        if (s[dist] == 'G' && rng.uniform() < prob) s[dist] = 'A';
    }
    return s;
}

int main(int argc, char** argv) {
    int      genome_len = 2000;
    int      n_reads    = 50000;
    int      min_len    = 30;
    int      max_len    = 80;
    double   dmax5      = 0.0;
    double   dmax3      = 0.0;
    double   lambda5    = 0.35;
    double   lambda3    = 0.35;
    int      qual       = 37;
    uint64_t seed       = 42;

    for (int i = 1; i < argc; ++i) {
        auto d = [&]{ return std::stod(argv[++i]); };
        auto n = [&]{ return std::stoi(argv[++i]); };
        const char* a = argv[i];
        if      (!strcmp(a, "--genome-len")) genome_len = n();
        else if (!strcmp(a, "--n-reads"))    n_reads    = n();
        else if (!strcmp(a, "--min-len"))    min_len    = n();
        else if (!strcmp(a, "--max-len"))    max_len    = n();
        else if (!strcmp(a, "--dmax5"))      dmax5      = d();
        else if (!strcmp(a, "--dmax3"))      dmax3      = d();
        else if (!strcmp(a, "--lambda5"))    lambda5    = d();
        else if (!strcmp(a, "--lambda3"))    lambda3    = d();
        else if (!strcmp(a, "--qual"))       qual       = n();
        else if (!strcmp(a, "--seed"))       seed       = (uint64_t)n();
        else { std::cerr << "Unknown option: " << a << '\n'; return 1; }
    }

    std::cerr << "gen_genome_frags:"
              << " genome_len=" << genome_len
              << " n_reads="    << n_reads
              << " min_len="    << min_len
              << " max_len="    << max_len
              << " dmax5="      << dmax5
              << " dmax3="      << dmax3
              << '\n';

    PCG32 rng(seed);

    // Generate reference genome
    std::string genome(genome_len, 'N');
    for (char& c : genome) c = BASES[rng.range(4)];

    const char qual_char = (char)(qual + 33);

    for (int i = 0; i < n_reads; ++i) {
        int flen = min_len + (int)rng.range(max_len - min_len + 1);
        int start = (int)rng.range(genome_len - flen + 1);
        std::string seq = genome.substr(start, flen);

        // Randomly orient (simulate both strands)
        if (rng.uniform() < 0.5) seq = revcomp(seq);

        // Apply ancient DNA damage
        if (dmax5 > 0.0 || dmax3 > 0.0)
            seq = apply_damage(seq, dmax5, lambda5, dmax3, lambda3, rng);

        std::cout << '@' << 'r' << (i + 1) << '\n'
                  << seq << '\n'
                  << '+' << '\n'
                  << std::string(seq.size(), qual_char) << '\n';
    }

    return 0;
}
