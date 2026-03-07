// gen_synthetic — Realistic FASTQ generator for fqdup testing.
//
// Models three stochastic processes:
//
//   1. ANCIENT-DNA DEAMINATION DAMAGE (Briggs et al. 2007 / DART model)
//      C→T at 5' end and G→A at 3' end with exponential decay:
//        P_CT(p) = d_max * exp(-lambda * p)     p = position from terminus
//      Applied independently to each read.
//
//   2. PCR GENEALOGY ERRORS (Pienaar et al. 2006 + Potapov & Ong 2017)
//
//      Full PCR tree model:
//        - At cycle k, the pool has 2^(k-1) template molecules.
//        - Each new strand synthesis has probability ε_eff per base of error.
//        - An error introduced at cycle k propagates to 2^(n-k) final copies.
//
//      Two error classes:
//        a) CHAINED errors (early cycles, k ≤ k*):
//             affect multiple reads from the same molecule.
//             k* = floor(log2(mean_coverage))  [≈ 6 for 50x coverage, 30 cycles]
//             Create "phantom sub-clones" — clusters with count > 1 that are
//             NOT independent singleton errors.  Error correction should NOT
//             absorb these (they look like real variants).
//
//        b) INDEPENDENT errors (late cycles, k > k*):
//             affect at most one read per original molecule.
//             Create singleton clusters.  Error correction SHOULD absorb these.
//
//      Error rate parameters (Potapov & Ong 2017):
//        ε_pol   = polymerase fidelity per base per doubling:
//                    Q5:      5.3e-7   (default — high-fidelity)
//                    Phusion: 3.9e-6
//                    KOD:     1.2e-5
//                    Taq:     1.5e-4
//        ε_thermo = thermocycling oxidative damage: 1.4e-6/base/cycle (97% C→T)
//                   Uniform distribution (not terminal-enriched).
//        ε_eff    = ε_pol  (thermocycling handled separately, see below)
//        D        = n_cycles * log2(1 + efficiency)   effective doublings
//        μ_total  = ε_pol * L * D                     expected total errors/read
//
//      Implementation of chained vs independent split:
//        For each molecule producing M reads (M = Poisson(coverage)):
//          Chained errors (cycle 1..k*):
//            n_chained ~ Poisson(ε_pol * L * k*)
//            Each chained error is at a random position; substituted base is
//            random non-original.  Applied to ALL M reads from this molecule.
//          Independent errors (cycle k*+1..n):
//            Per-read:  n_indep ~ Poisson(ε_pol * L * (D - k*))
//            Unique to each individual read.
//
//   3. DUPLICATE LIBRARY SAMPLING
//      Each read is drawn i.i.d. from the molecule pool (Poisson coverage).
//
// Usage:
//   gen_synthetic [options] > out.fq
//
// Options (all optional):
//   --n-unique  N      Distinct molecules           (default: 10000)
//   --n-reads   N      Total reads in estimation mode (default: 100000)
//   --read-len  N      Read length bp               (default: 75)
//   --len-sd    F      Std dev of fragment length distribution (default: 0 = fixed)
//   --min-len   N      Minimum fragment length when variable   (default: 30)
//   --max-len   N      Maximum fragment length when variable   (default: 150)
//   --seed      N      RNG seed                     (default: 42)
//
//   Damage:
//   --dmax5     F      5' C→T amplitude             (default: 0.25)
//   --dmax3     F      3' G→A amplitude             (default: 0.20)
//   --lambda5   F      5' decay rate                (default: 0.35)
//   --lambda3   F      3' decay rate                (default: 0.35)
//   --ss               Single-stranded: 5' damage only (dmax3 = 0)
//   --no-damage        Disable all deamination
//
//   PCR errors:
//   --pcr-cycles N     PCR cycles                   (default: 0 = no errors)
//   --pcr-eff    F     PCR efficiency 0–1           (default: 1.0)
//   --polymerase STR   q5|phusion|kod|taq           (default: q5)
//   --pcr-thermo       Include thermocyclic C→T damage (uniform, 1.4e-6/bp/cycle)
//   --pcr-rate   F     Direct per-base error probability (overrides polymerase model)
//
//   Oxidative damage (Channel C/D — 8-oxoG):
//   --ox-rate    F     Per-G-base G→T probability, uniform across read (default: 0)
//                      Models 8-oxoG oxidative damage.  EC should NOT absorb these
//                      (is_damage_sub protects G↔T differences).
//
//   Output mode:
//   --dup-pair         Each molecule → 2 independently-processed reads
//                      (dedup correctness test)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// PCG-32 PRNG
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
    int base4()      { return (int)(next() & 3u); }
    double normal01() {
        // Box-Muller transform
        double u1, u2;
        do { u1 = uniform(); } while (u1 <= 0.0);
        u2 = uniform();
        static const double M_PI_val = std::acos(-1.0);
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI_val * u2);
    }
    // Poisson variate (Knuth, exact for small μ; Gaussian approx for large)
    int poisson(double mu) {
        if (mu <= 0.0) return 0;
        if (mu > 30.0) {
            // Normal approx: mu ± sqrt(mu)
            double u = (uniform() + uniform() + uniform() +
                        uniform() + uniform() + uniform() - 3.0) * std::sqrt(mu);
            return std::max(0, (int)std::round(mu + u));
        }
        double L = std::exp(-mu), p = 1.0; int k = 0;
        do { ++k; p *= uniform(); } while (p > L);
        return k - 1;
    }
};

static constexpr char BASES[4] = {'A', 'C', 'G', 'T'};

static int base_index(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        default:            return 3; // T/t
    }
}

static std::string random_seq(int len, PCG32& rng) {
    std::string s(len, 'N');
    for (int i = 0; i < len; ++i) s[i] = BASES[rng.base4()];
    return s;
}

static int draw_length(PCG32& rng, int read_len, double len_sd, int min_len, int max_len) {
    if (len_sd <= 0.0) return read_len;
    int L = (int)std::round(read_len + len_sd * rng.normal01());
    return std::max(min_len, std::min(max_len, L));
}

// ---------------------------------------------------------------------------
// Model 1: Ancient-DNA deamination
// ---------------------------------------------------------------------------
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
        double prob = dmax3 * std::exp(-lambda3 * p);
        if (s[L - 1 - p] == 'G' && rng.uniform() < prob) s[L - 1 - p] = 'A';
    }
    return s;
}

// ---------------------------------------------------------------------------
// Model 2a: PCR chained errors
// error_list = vector of (position, substitute_base) to apply to all reads
// ---------------------------------------------------------------------------
struct Mutation { int pos; char base; };

static std::vector<Mutation> draw_chained_errors(int L,
                                                  double mu_chained,
                                                  PCG32& rng) {
    std::vector<Mutation> errors;
    int n = rng.poisson(mu_chained);
    for (int i = 0; i < n; ++i) {
        int pos = (int)(rng.next() % (unsigned)L);
        errors.push_back({pos, '\0'});  // base filled in at apply time
    }
    return errors;
}

static std::string apply_mutations(const std::string& seq,
                                   const std::vector<Mutation>& muts,
                                   PCG32& rng) {
    if (muts.empty()) return seq;
    std::string s = seq;
    for (auto& m : muts) {
        int orig = base_index(s[m.pos]);
        int sub;
        do { sub = rng.base4(); } while (sub == orig);
        s[m.pos] = BASES[sub];
    }
    return s;
}

// ---------------------------------------------------------------------------
// Model 2b: Independent late-cycle PCR errors
// Applied per-read, each base independently.
// ---------------------------------------------------------------------------
static std::string apply_indep_errors(const std::string& seq,
                                      double mu_indep,
                                      PCG32& rng) {
    if (mu_indep <= 0.0) return seq;
    std::string s = seq;
    int L = (int)s.size();
    int n = rng.poisson(mu_indep);  // total errors for this read
    for (int i = 0; i < n; ++i) {
        int pos  = (int)(rng.next() % (unsigned)L);
        int orig = base_index(s[pos]);
        int sub;
        do { sub = rng.base4(); } while (sub == orig);
        s[pos] = BASES[sub];
    }
    return s;
}

// ---------------------------------------------------------------------------
// Thermocyclic oxidative C→T (uniform across read, 1.4e-6/base/cycle, 97% C→T)
// Applied per-read.
// ---------------------------------------------------------------------------
static std::string apply_thermo(const std::string& seq,
                                double mu_thermo, PCG32& rng) {
    if (mu_thermo <= 0.0) return seq;
    std::string s = seq;
    for (char& c : s) {
        if (c == 'C' && rng.uniform() < mu_thermo) c = 'T';
    }
    return s;
}

// Channel C/D: 8-oxoG oxidative damage — G→T at any position, uniform.
// Each G is converted to T independently with probability ox_rate.
// EC must NOT absorb these reads (is_damage_sub protects G↔T pairs).
static std::string apply_oxidative(const std::string& seq,
                                   double ox_rate, PCG32& rng) {
    if (ox_rate <= 0.0) return seq;
    std::string s = seq;
    for (char& c : s)
        if (c == 'G' && rng.uniform() < ox_rate) c = 'T';
    return s;
}

static void write_read(const std::string& name, const std::string& seq) {
    std::cout << '@' << name  << '\n'
              << seq          << '\n'
              << '+'          << '\n'
              << std::string(seq.size(), 'I') << '\n';
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    int      n_unique   = 10000;
    int      n_reads    = 100000;
    int      read_len   = 75;
    double   len_sd     = 0.0;
    int      min_len    = 30;
    int      max_len    = 150;
    double   dmax5      = 0.25;
    double   dmax3      = 0.20;
    double   lambda5    = 0.35;
    double   lambda3    = 0.35;
    bool     ss_mode    = false;
    bool     no_damage  = false;
    int      pcr_cycles = 0;
    double   pcr_eff    = 1.0;
    double   epsilon    = 5.3e-7;   // Q5
    double   pcr_rate   = -1.0;     // <0 = use polymerase model; ≥0 = direct per-base prob
    double   ox_rate    = 0.0;      // per-G-base G→T probability (Channel C/D, 8-oxoG)
    bool     pcr_thermo = false;
    uint64_t seed       = 42;
    bool     dup_pair   = false;

    for (int i = 1; i < argc; ++i) {
        auto d = [&]{ return std::stod(argv[++i]); };
        auto n = [&]{ return std::stoi(argv[++i]); };
        const char* a = argv[i];
        if      (!strcmp(a,"--n-unique"))   n_unique   = n();
        else if (!strcmp(a,"--n-reads"))    n_reads    = n();
        else if (!strcmp(a,"--read-len"))   read_len   = n();
        else if (!strcmp(a,"--len-sd"))     len_sd     = d();
        else if (!strcmp(a,"--min-len"))    min_len    = n();
        else if (!strcmp(a,"--max-len"))    max_len    = n();
        else if (!strcmp(a,"--dmax5"))      dmax5      = d();
        else if (!strcmp(a,"--dmax3"))      dmax3      = d();
        else if (!strcmp(a,"--lambda5"))    lambda5    = d();
        else if (!strcmp(a,"--lambda3"))    lambda3    = d();
        else if (!strcmp(a,"--ss"))         ss_mode    = true;
        else if (!strcmp(a,"--no-damage"))  no_damage  = true;
        else if (!strcmp(a,"--pcr-cycles")) pcr_cycles = n();
        else if (!strcmp(a,"--pcr-eff"))    pcr_eff    = d();
        else if (!strcmp(a,"--pcr-thermo")) pcr_thermo = true;
        else if (!strcmp(a,"--pcr-rate"))   pcr_rate   = d();
        else if (!strcmp(a,"--ox-rate"))    ox_rate    = d();
        else if (!strcmp(a,"--seed"))       seed       = (uint64_t)n();
        else if (!strcmp(a,"--dup-pair"))   dup_pair   = true;
        else if (!strcmp(a,"--polymerase")) {
            const char* p = argv[++i];
            if      (!strcmp(p,"q5"))      epsilon = 5.3e-7;
            else if (!strcmp(p,"phusion")) epsilon = 3.9e-6;
            else if (!strcmp(p,"kod"))     epsilon = 1.2e-5;
            else if (!strcmp(p,"taq"))     epsilon = 1.5e-4;
            else { std::cerr << "Unknown polymerase: " << p << '\n'; return 1; }
        }
        else { std::cerr << "Unknown option: " << a << '\n'; return 1; }
    }

    if (ss_mode) dmax3 = 0.0;

    // -------------------------------------------------------------------------
    // PCR error rates (Pienaar 2006 + Potapov & Ong 2017)
    //
    // Per-molecule mu values scale with molecule length L.  We store per-base
    // rate scalars here and compute actual mu inside make_read / mol_chains.
    //
    // Effective doublings:  D = n_cycles * log2(1 + efficiency)
    // Total μ per read:     μ_total = ε_pol * L * D
    //
    // Chained vs independent split:
    //   Coverage per molecule ≈ n_reads / n_unique
    //   k* = max(0, floor(log2(coverage)))  — cycles where error affects >1 read
    //   rate_chained = ε_pol * k*           (per base, per molecule)
    //   rate_indep   = ε_pol * (D - k*)     (per base, per read)
    // -------------------------------------------------------------------------
    double rate_chained = 0.0;  // per-base rate for chained (early-cycle) errors
    double rate_indep   = 0.0;  // per-base rate for independent (late-cycle) errors
    double mu_thermo    = 0.0;  // per-C probability (not length-dependent)
    if (pcr_rate >= 0.0) {
        // Direct per-base error probability — bypasses polymerase/cycles model.
        // All errors are independent (no chaining).
        rate_indep = pcr_rate;
    } else if (pcr_cycles > 0) {
        double D    = pcr_cycles * std::log2(1.0 + pcr_eff);
        double cov  = (double)n_reads / n_unique;
        double k_star = std::max(0.0, std::floor(std::log2(std::max(1.0, cov))));
        rate_chained = epsilon * std::min(k_star, D);
        rate_indep   = epsilon * std::max(0.0, D - k_star);
        if (pcr_thermo)
            mu_thermo = 1.4e-6 * 0.97 * pcr_cycles;  // per C, per read
    }

    std::cerr << "gen_synthetic:"
              << " n_unique="  << n_unique
              << " n_reads="   << (dup_pair ? n_unique*2 : n_reads)
              << " L="         << read_len;
    if (len_sd > 0.0)
        std::cerr << " len_sd=" << len_sd
                  << " min_len=" << min_len << " max_len=" << max_len;
    std::cerr << " dmax5="     << dmax5 << " dmax3=" << dmax3
              << " lambda="    << lambda5;
    if (pcr_cycles > 0)
        std::cerr << " pcr=" << pcr_cycles << "cyc"
                  << " rate_chain=" << rate_chained
                  << " rate_indep=" << rate_indep;
    if (pcr_rate >= 0.0)
        std::cerr << " pcr_rate=" << pcr_rate;
    if (mu_thermo > 0.0)
        std::cerr << " mu_thermo=" << mu_thermo;
    if (ox_rate > 0.0)
        std::cerr << " ox_rate=" << ox_rate;
    std::cerr << '\n';

    PCG32 rng(seed);

    std::vector<std::string> molecules(n_unique);
    for (auto& m : molecules) m = random_seq(draw_length(rng, read_len, len_sd, min_len, max_len), rng);

    // -------------------------------------------------------------------------
    // Generate reads
    // -------------------------------------------------------------------------
    auto make_read = [&](const std::string& mol,
                         const std::vector<Mutation>& shared_muts) -> std::string {
        int L = (int)mol.size();
        // Apply shared chained PCR errors first (same for all reads from this molecule)
        std::string seq = apply_mutations(mol, shared_muts, rng);
        // Independent late-cycle PCR errors (unique to this read), scaled by length
        double mol_mu_indep = rate_indep * L;
        seq = apply_indep_errors(seq, mol_mu_indep, rng);
        // Thermocyclic oxidative C→T (uniform, independent per read)
        if (mu_thermo > 0.0) seq = apply_thermo(seq, mu_thermo, rng);
        // 8-oxoG oxidative G→T (Channel C/D, uniform, independent per read)
        if (ox_rate > 0.0) seq = apply_oxidative(seq, ox_rate, rng);
        // Ancient-DNA terminal deamination (independent per read)
        if (!no_damage)
            seq = apply_damage(seq, dmax5, lambda5, dmax3, lambda3, rng);
        return seq;
    };

    if (dup_pair) {
        // Each molecule → exactly 2 reads with independent processing.
        // Shared chained errors apply to both (they come from same template).
        for (int i = 0; i < n_unique; ++i) {
            int L = (int)molecules[i].size();
            double mol_mu_chained = rate_chained * L;
            auto shared = draw_chained_errors(L, mol_mu_chained, rng);
            for (int copy = 0; copy < 2; ++copy) {
                char name[64];
                std::snprintf(name, sizeof(name), "mol%05d_c%d", i, copy);
                write_read(name, make_read(molecules[i], shared));
            }
        }
        return 0;
    }

    // Estimation / dedup-rate mode: i.i.d. sampling from molecule pool.
    // Pre-draw shared errors per molecule (only needed once per molecule).
    // For simplicity in this mode: no per-molecule chaining (molecules are
    // rarely sampled >1 time at low coverage).  At high coverage
    // (n_reads >> n_unique), chained errors do matter and we pre-compute them.
    bool need_chain = (rate_chained > 0.0 && n_reads > n_unique);

    // Pre-compute chained mutations per molecule (only if coverage > 1)
    std::vector<std::vector<Mutation>> mol_chains;
    if (need_chain) {
        mol_chains.resize(n_unique);
        for (int i = 0; i < n_unique; ++i) {
            int L = (int)molecules[i].size();
            mol_chains[i] = draw_chained_errors(L, rate_chained * L, rng);
        }
    }

    static const std::vector<Mutation> empty_chain;
    for (int i = 0; i < n_reads; ++i) {
        int m = (int)(rng.next() % (unsigned)n_unique);
        const auto& chain = need_chain ? mol_chains[m] : empty_chain;
        char name[64];
        std::snprintf(name, sizeof(name), "r%07d_m%05d", i, m);
        write_read(name, make_read(molecules[m], chain));
    }
    return 0;
}
