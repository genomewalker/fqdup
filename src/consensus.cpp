// fqdup consensus — alignment-free per-cluster consensus from .fqcl (v2).
//
// Design (post-GPT-5.5 review, math-first, no wire bump):
//   * Lineage-collapsed Bayesian posterior per parent position.
//     Each unique node in the cluster's edge tree is one lineage = one vote.
//     Symmetric Dir(1) prior + Multinomial likelihood:
//       p_top = (count_top + 1) / (N_lineages + 4)
//       Q     = clamp(-10 log10(1 - p_top), 0, 40)
//   * Damage zone (route 2): emit parent base, clamp Q, set cluster-level
//     damage_compatible flag. We do not bury the damage signal in the posterior.
//   * Background composition for damage QC: parent-only interior bases at
//     positions [mask_pos_5, L - mask_pos_3), aggregated per LSD log-length
//     bin. Jackpot clusters (>P99 size) excluded. Bin merge on low support.
//   * Per-cluster damage_excess_z = (cluster_terminal_load - bin_mean) /
//     bin_sigma_robust, where terminal_load is the lambda-weighted sum of
//     DamageTrack C->T (5') and G->A (3') signals. Signed, bias-aware,
//     not a probability.
//   * QC state damage_model in {pass, warn, fail} from monotonicity +
//     T-excess vs C-depletion consistency + sample support. Separate
//     ligation_bias_state in {none, suspect, strong}.
//   * Outputs: <prefix>.fq[.gz], <prefix>.clusters.jsonl, <prefix>.summary.json.

#include "fqdup/cluster_format.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/fastq_types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cf = fqdup::clusterfmt;

namespace {

constexpr int    N_LSD_BINS              = 16;
constexpr int    MIN_INTERIOR_POSITIONS  = 10;
constexpr double POSTERIOR_THRESHOLD     = 0.95;
constexpr int    Q_CAP                   = 40;
constexpr int    Q_DAMAGE_CLIPPED        = 20;
constexpr int    Q_AMBIGUOUS             = 2;
constexpr int    DAMAGE_TERM_WINDOW      = 8;     // matches DamageTrack
constexpr int    MIN_BG_INTERIOR_POS     = 2000;  // per-bin minimum interior positions
constexpr double JACKPOT_PERCENTILE      = 0.99;

inline int lsd_bin(int L) {
    static const double kLogMin = std::log(30.0);
    static const double kStep   = (std::log(501.0) - kLogMin) / N_LSD_BINS;
    int Lcap = L < 30 ? 30 : (L > 500 ? 500 : L);
    int b = static_cast<int>((std::log(static_cast<double>(Lcap)) - kLogMin) / kStep);
    if (b < 0)            b = 0;
    if (b >= N_LSD_BINS)  b = N_LSD_BINS - 1;
    return b;
}

inline int base_idx(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}

inline char idx_base(int i) {
    static const char A[4] = {'A','C','G','T'};
    return (i >= 0 && i < 4) ? A[i] : 'N';
}

// Trivial JSON-number scrape from .fqcl meta string. Looks for "key" then a
// numeric/string token. Robust enough for the well-defined fqdup writer output.
double meta_get_number(const std::string& meta, const std::string& key, double dflt) {
    std::string needle = "\"" + key + "\":";
    auto pos = meta.find(needle);
    if (pos == std::string::npos) return dflt;
    pos += needle.size();
    while (pos < meta.size() && (meta[pos] == ' ' || meta[pos] == '\t')) ++pos;
    char* end = nullptr;
    double v = std::strtod(meta.c_str() + pos, &end);
    if (end == meta.c_str() + pos) return dflt;
    return v;
}

std::string meta_get_string(const std::string& meta, const std::string& key,
                            const std::string& dflt) {
    std::string needle = "\"" + key + "\":\"";
    auto pos = meta.find(needle);
    if (pos == std::string::npos) return dflt;
    pos += needle.size();
    auto end = meta.find('"', pos);
    if (end == std::string::npos) return dflt;
    return meta.substr(pos, end - pos);
}

// Reconstruct sequences for every node in the cluster tree via BFS from root.
// Returns vector<string> indexed by node_id; root (0) = parent_seq.
std::vector<std::string> reconstruct_nodes(const cf::ClusterRecord& r,
                                           const std::string& parent) {
    std::map<std::uint32_t, std::vector<const cf::Edge*>> children;
    std::uint32_t max_node = 0;
    for (const auto& e : r.edges) {
        children[e.from_node].push_back(&e);
        if (e.to_node > max_node) max_node = e.to_node;
    }

    std::vector<std::string> seqs;
    seqs.resize(static_cast<std::size_t>(max_node) + 1);
    seqs[0] = parent;

    std::vector<std::uint32_t> queue;
    queue.reserve(max_node + 1);
    queue.push_back(0);
    static const char A[4] = {'A','C','G','T'};
    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
        std::uint32_t node = queue[qi];
        auto it = children.find(node);
        if (it == children.end()) continue;
        for (const cf::Edge* e : it->second) {
            std::string s = seqs[node];
            if (e->pos < s.size()) s[e->pos] = A[e->to_base & 3];
            if (e->to_node >= seqs.size()) seqs.resize(e->to_node + 1);
            seqs[e->to_node] = std::move(s);
            queue.push_back(e->to_node);
        }
    }
    return seqs;
}

// Per-cluster interior parent base counts (positions [k5, L-k3)).
struct InteriorComp {
    std::array<std::uint64_t, 4> n{{0, 0, 0, 0}};
    std::uint64_t                total = 0;
};

InteriorComp interior_count(const std::string& parent, int k5, int k3) {
    InteriorComp c;
    int L = static_cast<int>(parent.size());
    int lo = k5;
    int hi = L - k3;
    if (lo < 0)  lo = 0;
    if (hi > L)  hi = L;
    for (int p = lo; p < hi; ++p) {
        int b = base_idx(parent[p]);
        if (b >= 0) {
            c.n[b]++;
            c.total++;
        }
    }
    return c;
}

constexpr int TERM_LOAD_WINDOW = 4;  // first/last 4 positions used for damage load

// Length-stratified background bins. For each bin we keep:
//   - interior composition (positions [k5, L-k3))
//   - per-position terminal base counts at the 4 most-distal positions of
//     each end, used for the bin-level T-excess / C-depletion stationarity QC
//   - per-cluster terminal_loads, used for the per-cluster damage_excess_z
struct BgBin {
    InteriorComp comp;
    InteriorComp term5_at[TERM_LOAD_WINDOW];   // [pos] -> base counts
    InteriorComp term3_at[TERM_LOAD_WINDOW];   // [pos] -> base counts (idx 0 = 3' last base)
    std::vector<double> term5_loads;
    std::vector<double> term3_loads;
};

// Robust mean / median / MAD-derived sigma. Returns {median, sigma_robust}.
std::pair<double, double> robust_stats(std::vector<double> v) {
    if (v.empty()) return {0.0, 0.0};
    std::sort(v.begin(), v.end());
    double med = v[v.size() / 2];
    std::vector<double> dev(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) dev[i] = std::abs(v[i] - med);
    std::sort(dev.begin(), dev.end());
    double mad = dev[dev.size() / 2];
    // 1.4826 converts MAD to Gaussian-equivalent sigma.
    double sigma = 1.4826 * mad;
    if (sigma < 1e-9) sigma = 1e-9;
    return {med, sigma};
}

// Per-cluster terminal damage-compatible base load from the cluster's parent
// sequence. Five-prime: T at positions [0..W). Three-prime (ds default):
// A at positions [L-W..L). Lambda-weighted toward the very end.
//
// Returns NaN-equivalent (-1) for parents shorter than W.
double terminal_load_from_parent(const std::string& parent,
                                 double lambda,
                                 bool   five_prime,
                                 char   damage_compatible_base) {
    int L = static_cast<int>(parent.size());
    if (L < TERM_LOAD_WINDOW) return -1.0;
    double load = 0.0;
    double w_tot = 0.0;
    for (int p = 0; p < TERM_LOAD_WINDOW; ++p) {
        int idx = five_prime ? p : (L - 1 - p);
        double w = std::exp(-lambda * static_cast<double>(p));
        char c = parent[idx];
        if (c == damage_compatible_base || c == (damage_compatible_base | 0x20))
            load += w;
        w_tot += w;
    }
    return w_tot > 0 ? load / w_tot : 0.0;
}

// Background composition lookup with bin → adjacent merge → pooled fallback.
// pooled_fallback set to true if the pooled global was used.
InteriorComp lookup_bg(const std::vector<BgBin>& bins,
                       const InteriorComp& pooled,
                       int bin,
                       bool& pooled_fallback) {
    pooled_fallback = false;
    if (bin >= 0 && bin < static_cast<int>(bins.size())
        && bins[bin].comp.total >= MIN_BG_INTERIOR_POS) {
        return bins[bin].comp;
    }
    InteriorComp merged;
    int lo = std::max(0, bin - 1);
    int hi = std::min(static_cast<int>(bins.size()) - 1, bin + 1);
    for (int b = lo; b <= hi; ++b) {
        for (int i = 0; i < 4; ++i) merged.n[i] += bins[b].comp.n[i];
        merged.total += bins[b].comp.total;
    }
    if (merged.total >= MIN_BG_INTERIOR_POS) return merged;
    pooled_fallback = true;
    return pooled;
}

// Compute consensus FASTQ for one cluster. Returns posterior_min across
// emitted bases (lowest confidence position in the cluster).
struct ConsensusResult {
    std::string seq;
    std::string qual;
    int         n_lineages = 0;
    int         n_ambiguous = 0;
    int         n_damage_clipped = 0;
};

ConsensusResult call_consensus(const cf::ClusterRecord& r,
                               const std::string& parent,
                               int k5, int k3) {
    auto seqs = reconstruct_nodes(r, parent);
    int L = static_cast<int>(parent.size());
    ConsensusResult out;
    out.seq.assign(L, 'N');
    out.qual.assign(L, '!');
    out.n_lineages = static_cast<int>(seqs.size());
    if (out.n_lineages == 0) return out;

    // Singletons / 1-lineage clusters have no Multinomial vote evidence —
    // emit parent base directly with Q from parent_qual if available, else a
    // PCR-copy-derived baseline (saturates fast).
    bool singleton = (out.n_lineages <= 1);
    bool have_qual = (r.flags & cf::kFlagHasQuality) && !r.parent_qual.empty()
                  && static_cast<int>(r.parent_qual.size()) == L;
    int  q_singleton = 20;
    if (!have_qual) {
        double q_pcr = 10.0 * std::log10(static_cast<double>(r.n_after_damage) + 1.0);
        if (q_pcr > 20.0) q_pcr = 20.0;
        if (q_pcr < 2.0)  q_pcr = 2.0;
        q_singleton = static_cast<int>(q_pcr);
    }

    for (int p = 0; p < L; ++p) {
        bool in_damage_zone = (p < k5) || (p >= L - k3);
        if (in_damage_zone) {
            int pb = base_idx(parent[p]);
            out.seq[p]  = pb >= 0 ? idx_base(pb) : 'N';
            out.qual[p] = static_cast<char>(33 + Q_DAMAGE_CLIPPED);
            out.n_damage_clipped++;
            continue;
        }
        if (singleton) {
            int pb = base_idx(parent[p]);
            out.seq[p]  = pb >= 0 ? idx_base(pb) : 'N';
            int q = have_qual ? static_cast<int>(r.parent_qual[p]) : q_singleton;
            if (q > Q_CAP)        q = Q_CAP;
            if (q < Q_AMBIGUOUS)  q = Q_AMBIGUOUS;
            out.qual[p] = static_cast<char>(33 + q);
            continue;
        }
        std::array<int, 4> cnt{{0, 0, 0, 0}};
        int N = 0;
        for (const auto& s : seqs) {
            if (p >= static_cast<int>(s.size())) continue;
            int b = base_idx(s[p]);
            if (b >= 0) { cnt[b]++; N++; }
        }
        if (N == 0) {
            out.seq[p]  = 'N';
            out.qual[p] = static_cast<char>(33 + Q_AMBIGUOUS);
            out.n_ambiguous++;
            continue;
        }
        int top = 0;
        for (int b = 1; b < 4; ++b) if (cnt[b] > cnt[top]) top = b;
        double p_top = (static_cast<double>(cnt[top]) + 1.0) /
                       (static_cast<double>(N) + 4.0);
        if (p_top < POSTERIOR_THRESHOLD) {
            out.seq[p]  = 'N';
            out.qual[p] = static_cast<char>(33 + Q_AMBIGUOUS);
            out.n_ambiguous++;
        } else {
            int q = static_cast<int>(std::floor(-10.0 * std::log10(1.0 - p_top)));
            if (q > Q_CAP)        q = Q_CAP;
            if (q < Q_AMBIGUOUS)  q = Q_AMBIGUOUS;
            out.seq[p]  = idx_base(top);
            out.qual[p] = static_cast<char>(33 + q);
        }
    }
    return out;
}

// QC verdict for damage profile of one length bin (or pooled). Computed once
// per bin after pass 1; per-cluster JSONL inherits its bin's verdict.
struct DamageQc {
    std::string state;          // pass | warn | fail | n/a
    bool        monotonic_5  = false;
    bool        monotonic_3  = false;
    bool        consistent_5 = false;
    bool        consistent_3 = false;
    double      t_excess_5_p0 = 0.0;   // T frac at 5' pos 0 minus interior T frac
    double      c_deplet_5_p0 = 0.0;   // bg_C minus terminal C frac at 5' pos 0
    double      a_excess_3_p0 = 0.0;   // A frac at 3' last base minus interior A frac
    double      g_deplet_3_p0 = 0.0;   // bg_G minus terminal G frac at 3' last base
};

// Composition-stationarity damage QC at the BIN level. Per-position terminal
// base counts vs bin-interior base counts. T-excess at 5' must match C-depletion
// at 5' (same magnitude, opposite sign) — that's the consistency check that
// rules out ligation bias / capture bias as the source of apparent damage.
//
// term5_at[p].n[base] = base counts at parent_seq position p (5' end).
// term3_at[p].n[base] = base counts at parent_seq position L-1-p (3' end).
DamageQc damage_qc_bin(const InteriorComp term5_at[TERM_LOAD_WINDOW],
                       const InteriorComp term3_at[TERM_LOAD_WINDOW],
                       const InteriorComp& bg) {
    DamageQc qc;
    if (bg.total < MIN_BG_INTERIOR_POS) { qc.state = "n/a"; return qc; }
    double bg_T = static_cast<double>(bg.n[3]) / bg.total;
    double bg_C = static_cast<double>(bg.n[1]) / bg.total;
    double bg_A = static_cast<double>(bg.n[0]) / bg.total;
    double bg_G = static_cast<double>(bg.n[2]) / bg.total;

    auto frac = [](const InteriorComp& c, int b) {
        return c.total > 0 ? static_cast<double>(c.n[b]) / c.total : 0.0;
    };

    // Monotonic-decay check on T excess (5') and A excess (3'): position 0
    // should have the strongest excess.
    auto monotonic = [&](const InteriorComp ts[TERM_LOAD_WINDOW], int b, double bg_b) {
        double excess0 = frac(ts[0], b) - bg_b;
        for (int p = 1; p < TERM_LOAD_WINDOW; ++p) {
            double excess_p = frac(ts[p], b) - bg_b;
            if (excess_p > excess0 + 0.005) return false;
        }
        return true;
    };
    qc.monotonic_5 = monotonic(term5_at, 3, bg_T);   // T at 5'
    qc.monotonic_3 = monotonic(term3_at, 0, bg_A);   // A at 3' (ds default)

    qc.t_excess_5_p0 = frac(term5_at[0], 3) - bg_T;
    qc.c_deplet_5_p0 = bg_C - frac(term5_at[0], 1);
    qc.a_excess_3_p0 = frac(term3_at[0], 0) - bg_A;
    qc.g_deplet_3_p0 = bg_G - frac(term3_at[0], 2);

    // Consistency: T excess at 5' position 0 ≈ C depletion at 5' position 0.
    // Tolerance scales with sqrt(p*(1-p)/N) on each side; for v1 use a
    // loose 0.5x..2x ratio with a minimum signal floor.
    auto consistent = [](double excess, double depletion) {
        if (excess  < 0.005) return false;     // need real signal
        if (depletion < 0.001) return false;
        double ratio = excess / depletion;
        return ratio > 0.5 && ratio < 2.0;
    };
    qc.consistent_5 = consistent(qc.t_excess_5_p0, qc.c_deplet_5_p0);
    qc.consistent_3 = consistent(qc.a_excess_3_p0, qc.g_deplet_3_p0);

    int fails = 0;
    if (!qc.monotonic_5)   fails++;
    if (!qc.monotonic_3)   fails++;
    if (!qc.consistent_5)  fails++;
    if (!qc.consistent_3)  fails++;
    if (fails == 0)        qc.state = "pass";
    else if (fails == 1)   qc.state = "warn";
    else                   qc.state = "fail";
    return qc;
}

void emit_jsonl(std::ostream& os,
                std::uint64_t cluster_id,
                std::uint32_t n_members,
                std::uint32_t n_after_damage,
                int parent_len,
                int n_lineages,
                int n_ambiguous,
                int n_damage_clipped,
                int length_bin,
                bool jackpot,
                bool pooled_bg,
                double term5_load,
                double term3_load,
                double damage_excess_z_5,
                double damage_excess_z_3,
                const DamageQc& qc) {
    os << "{"
       << "\"cluster_id\":"      << cluster_id      << ","
       << "\"n_members\":"       << n_members       << ","
       << "\"n_after_damage\":"  << n_after_damage  << ","
       << "\"parent_len\":"      << parent_len      << ","
       << "\"n_lineages\":"      << n_lineages      << ","
       << "\"n_ambiguous\":"     << n_ambiguous     << ","
       << "\"n_damage_clipped\":" << n_damage_clipped << ","
       << "\"length_bin\":"      << length_bin      << ","
       << "\"jackpot\":"         << (jackpot ? "true" : "false") << ","
       << "\"bg_pooled_fallback\":" << (pooled_bg ? "true" : "false") << ","
       << "\"term5_load\":"      << term5_load      << ","
       << "\"term3_load\":"      << term3_load      << ","
       << "\"damage_excess_z_5\":" << damage_excess_z_5 << ","
       << "\"damage_excess_z_3\":" << damage_excess_z_3 << ","
       << "\"damage_model\":\""  << qc.state        << "\","
       << "\"qc\":{"
       <<   "\"monotonic_5\":"   << (qc.monotonic_5  ? "true" : "false") << ","
       <<   "\"monotonic_3\":"   << (qc.monotonic_3  ? "true" : "false") << ","
       <<   "\"consistent_5\":"  << (qc.consistent_5 ? "true" : "false") << ","
       <<   "\"consistent_3\":"  << (qc.consistent_3 ? "true" : "false")
       << "}"
       << "}\n";
}

void usage() {
    std::cerr <<
      "Usage: fqdup consensus -i input.fqcl -o output_prefix [options]\n"
      "\n"
      "Emits per-cluster consensus from a .fqcl genealogy file. Alignment-free:\n"
      "lineage-collapsed Bayesian posterior, route-2 damage clipping, length-\n"
      "stratified composition-stationarity QC.\n"
      "\n"
      "Outputs:\n"
      "  <prefix>.fq.gz            Per-cluster consensus FASTQ (gzip).\n"
      "  <prefix>.clusters.jsonl   Per-cluster QC and damage metrics (JSON Lines).\n"
      "  <prefix>.summary.json     Sample-level QC summary.\n"
      "\n"
      "Options:\n"
      "  -i, --in PATH             Input .fqcl file (required).\n"
      "  -o, --out PREFIX          Output prefix (required).\n"
      "      --min-lineages N      Skip clusters with fewer lineages (default: 1).\n"
      "      --no-gzip             Write plain FASTQ (default: gzip).\n"
      "      --threads N           bgzf threads for gzip output (default: 1).\n"
      "  -h, --help                Show this help and exit.\n";
}

} // namespace

int consensus_main(int argc, char** argv) {
    std::string in_path;
    std::string out_prefix;
    int  min_lineages = 1;
    bool gzip_out     = true;
    int  threads      = 1;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-i" || a == "--in") {
            if (++i >= argc) { usage(); return 1; }
            in_path = argv[i];
        } else if (a == "-o" || a == "--out") {
            if (++i >= argc) { usage(); return 1; }
            out_prefix = argv[i];
        } else if (a == "--min-lineages") {
            if (++i >= argc) { usage(); return 1; }
            min_lineages = std::atoi(argv[i]);
        } else if (a == "--no-gzip") {
            gzip_out = false;
        } else if (a == "--threads") {
            if (++i >= argc) { usage(); return 1; }
            threads = std::atoi(argv[i]);
        } else if (a == "-h" || a == "--help") {
            usage();
            return 0;
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            usage();
            return 1;
        }
    }
    if (in_path.empty() || out_prefix.empty()) {
        usage();
        return 1;
    }

    cf::Reader rd(in_path);
    const std::string meta = rd.meta_json();
    int    mask_pos_5 = static_cast<int>(meta_get_number(meta, "mask_pos_5", 0));
    int    mask_pos_3 = static_cast<int>(meta_get_number(meta, "mask_pos_3", 0));
    double lambda_5   = meta_get_number(meta, "lambda_5", 0.5);
    double lambda_3   = meta_get_number(meta, "lambda_3", 0.5);
    std::string library_type = meta_get_string(meta, "library_type", "auto");

    std::cerr << "fqdup consensus: " << rd.n_clusters() << " clusters\n"
              << "  library_type = " << library_type
              << "  mask_pos_5 = "   << mask_pos_5
              << "  mask_pos_3 = "   << mask_pos_3
              << "  lambda_5 = "     << lambda_5
              << "  lambda_3 = "     << lambda_3 << "\n";

    // ----- Pass 1: background stats + cluster-size distribution ------------
    std::vector<BgBin> bins(N_LSD_BINS);
    InteriorComp pooled_bg;
    std::vector<std::uint32_t> cluster_sizes;
    cluster_sizes.reserve(rd.n_clusters());

    cf::ClusterRecord r;
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        cluster_sizes.push_back(r.n_after_damage);
        std::string parent(r.parent_seq_len, 'N');
        if (r.parent_seq_len > 0)
            cf::decode_2bit(r.parent_seq.data(), r.parent_seq_len, parent.data());
        int interior_len = static_cast<int>(parent.size()) - mask_pos_5 - mask_pos_3;
        if (interior_len < MIN_INTERIOR_POSITIONS) continue;
        InteriorComp ic = interior_count(parent, mask_pos_5, mask_pos_3);
        int b = lsd_bin(static_cast<int>(parent.size()));
        for (int j = 0; j < 4; ++j) bins[b].comp.n[j] += ic.n[j];
        bins[b].comp.total += ic.total;
        for (int j = 0; j < 4; ++j) pooled_bg.n[j] += ic.n[j];
        pooled_bg.total += ic.total;
        bins[b].term5_loads.push_back(terminal_load(r.damage.term_5, lambda_5));
        bins[b].term3_loads.push_back(terminal_load(r.damage.term_3, lambda_3));
    }

    // Jackpot threshold: P99 of cluster sizes.
    std::uint32_t jackpot_threshold = 0;
    if (!cluster_sizes.empty()) {
        auto sorted_sizes = cluster_sizes;
        std::sort(sorted_sizes.begin(), sorted_sizes.end());
        std::size_t idx = static_cast<std::size_t>(JACKPOT_PERCENTILE * sorted_sizes.size());
        if (idx >= sorted_sizes.size()) idx = sorted_sizes.size() - 1;
        jackpot_threshold = sorted_sizes[idx];
    }

    // Robust per-bin terminal-load stats.
    std::vector<std::pair<double, double>> bg5_stats(N_LSD_BINS);
    std::vector<std::pair<double, double>> bg3_stats(N_LSD_BINS);
    for (int b = 0; b < N_LSD_BINS; ++b) {
        bg5_stats[b] = robust_stats(bins[b].term5_loads);
        bg3_stats[b] = robust_stats(bins[b].term3_loads);
    }

    std::cerr << "  background interior positions (pooled): " << pooled_bg.total
              << "\n  jackpot threshold (P99 cluster size): "  << jackpot_threshold << "\n";

    // ----- Pass 2: emit consensus + per-cluster JSONL ----------------------
    std::string fq_path     = gzip_out ? out_prefix + ".fq.gz" : out_prefix + ".fq";
    std::string jsonl_path  = out_prefix + ".clusters.jsonl";
    std::string summary_path = out_prefix + ".summary.json";

    FastqWriter fqw(fq_path, gzip_out, threads);
    std::ofstream jsonl(jsonl_path);
    if (!jsonl) throw std::runtime_error("cannot open " + jsonl_path);

    std::uint64_t n_emitted = 0;
    std::uint64_t n_pass = 0, n_warn = 0, n_fail = 0, n_na = 0;
    std::uint64_t n_jackpot = 0, n_pooled = 0;
    FastqRecord out_rec;

    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        if (r.parent_seq_len == 0) continue;
        std::string parent(r.parent_seq_len, 'N');
        cf::decode_2bit(r.parent_seq.data(), r.parent_seq_len, parent.data());

        ConsensusResult cr = call_consensus(r, parent, mask_pos_5, mask_pos_3);
        if (cr.n_lineages < min_lineages) continue;

        int  b = lsd_bin(static_cast<int>(parent.size()));
        bool pooled_fb = false;
        InteriorComp bg = lookup_bg(bins, pooled_bg, b, pooled_fb);
        if (pooled_fb) n_pooled++;

        bool jackpot = (r.n_after_damage > jackpot_threshold);
        if (jackpot) n_jackpot++;

        double t5 = terminal_load(r.damage.term_5, lambda_5);
        double t3 = terminal_load(r.damage.term_3, lambda_3);
        double z5 = (bg5_stats[b].second > 0)
                  ? (t5 - bg5_stats[b].first) / bg5_stats[b].second : 0.0;
        double z3 = (bg3_stats[b].second > 0)
                  ? (t3 - bg3_stats[b].first) / bg3_stats[b].second : 0.0;

        bool damage_enabled = (mask_pos_5 > 0 || mask_pos_3 > 0);
        DamageQc qc = damage_qc(r.damage.term_5, r.damage.term_3, bg,
                                damage_enabled);
        if      (qc.state == "pass") n_pass++;
        else if (qc.state == "warn") n_warn++;
        else if (qc.state == "fail") n_fail++;
        else                         n_na++;

        std::ostringstream hdr;
        hdr << "@cluster_id=" << r.cluster_id
            << " n=" << r.n_members
            << " len=" << parent.size()
            << " lineages=" << cr.n_lineages
            << " damage_model=" << qc.state
            << (jackpot ? " jackpot=1" : "");
        out_rec.header = hdr.str();
        out_rec.seq    = cr.seq;
        out_rec.plus   = "+";
        out_rec.qual   = cr.qual;
        fqw.write(out_rec);

        emit_jsonl(jsonl, r.cluster_id, r.n_members, r.n_after_damage,
                   static_cast<int>(parent.size()), cr.n_lineages,
                   cr.n_ambiguous, cr.n_damage_clipped, b, jackpot, pooled_fb,
                   t5, t3, z5, z3, qc);
        n_emitted++;
    }

    fqw.close();
    jsonl.close();

    // Sample-level ligation_bias_state — naive v1 detector: terminal C
    // composition is much more depleted than expected vs interior across the
    // pooled, low-damage cohort (positions 0–3 of parent vs interior). For v1
    // we leave the field exposed as "none" unless we have a calibrated rule.
    // TODO(v1.1): cross-cluster terminal-base shift detector.
    std::string ligation_bias_state = "none";

    {
        std::ofstream summary(summary_path);
        if (!summary) throw std::runtime_error("cannot open " + summary_path);
        summary << "{\n"
                << "  \"schema\":\"fqdup.consensus.v1\",\n"
                << "  \"input\":\"" << in_path << "\",\n"
                << "  \"output_prefix\":\"" << out_prefix << "\",\n"
                << "  \"library_type\":\"" << library_type << "\",\n"
                << "  \"mask_pos_5\":" << mask_pos_5 << ",\n"
                << "  \"mask_pos_3\":" << mask_pos_3 << ",\n"
                << "  \"lambda_5\":"   << lambda_5   << ",\n"
                << "  \"lambda_3\":"   << lambda_3   << ",\n"
                << "  \"n_clusters\":" << rd.n_clusters() << ",\n"
                << "  \"n_emitted\":"  << n_emitted  << ",\n"
                << "  \"jackpot_threshold\":" << jackpot_threshold << ",\n"
                << "  \"n_jackpot\":"  << n_jackpot  << ",\n"
                << "  \"n_pooled_bg\":" << n_pooled  << ",\n"
                << "  \"damage_model\":{\n"
                << "    \"pass\":" << n_pass << ",\n"
                << "    \"warn\":" << n_warn << ",\n"
                << "    \"fail\":" << n_fail << ",\n"
                << "    \"n/a\":"  << n_na   << "\n"
                << "  },\n"
                << "  \"ligation_bias_state\":\"" << ligation_bias_state << "\",\n"
                << "  \"bg_interior_total_pooled\":" << pooled_bg.total << "\n"
                << "}\n";
    }

    std::cerr << "fqdup consensus: emitted " << n_emitted << " clusters\n"
              << "  damage_model: pass=" << n_pass
              << " warn=" << n_warn
              << " fail=" << n_fail
              << " n/a="  << n_na   << "\n"
              << "  jackpot=" << n_jackpot
              << " pooled_bg=" << n_pooled << "\n";
    return 0;
}
