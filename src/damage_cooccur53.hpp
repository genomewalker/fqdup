#pragma once
#include <cstdint>
#include <cmath>
#include <string>
#include <array>
#include <sstream>
#include <iomanip>

// Length-stratified per-read 5'x3' JOINT deamination estimator.
//
// A single ancient molecule tends to carry damage at BOTH ends (5' C->T and
// 3' G->A for double-stranded; 3' C->T for single-stranded). The bulk position
// estimator marginalises over reads, so as the ancient fraction pi -> 0 the
// terminal rate collapses to background and d_max reads 0. The WITHIN-molecule
// co-occurrence of the two ends does NOT dilute that way: its expectation is a
// property of the damaged molecules themselves; only the mixing weight pi (the
// COUNT of such molecules), not the conditional shape, shrinks with dilution.
//
// Estimand. Per length bin, from the joint of terminal indicators (U5,U3) and
// interior negative-control indicators (V5,V3):
//   delta_cov = Cov(U5,U3) - Cov(V5,V3) = pi(1-pi) * A5 * A3
// where A5=(1-bg5)*d5, A3=(1-bg3)*d3. The read's own compositional covariance
// (AT-rich reads elevate both ends) is present under BOTH pairs and CANCELS in
// the difference -- so the interior channel is a built-in composition null.
// The terminal MARGINALS are NOT clean (end-composition offset, the "negative
// terminal excess" prep fingerprint), so only the covariance difference is used
// for detection. pi and d are not separable from a single bin (one scalar,
// three unknowns): delta_cov is the identifiable, dilution-invariant detector;
// pi is reported only conditionally on an assumed d.
//
// Runs UNCONDITIONALLY in pass 1 (accumulate() is called from worker_ox_accumulate,
// NOT behind the bulk d_max>0.01 oxog gate), so it recovers damage in libraries
// the bulk fit calls absent.
namespace fqdup_j53 {

struct J53Acc {
    static constexpr int NBIN = 5;                 // [30,50) [50,70) [70,90) [90,120) [120,inf)
    static constexpr int LEN_LO[NBIN] = {30, 50, 70, 90, 120};
    static constexpr int LEN_HI[NBIN] = {50, 70, 90, 120, 1 << 30};
    // channel 0 = ds (3' proxy base 'A'), 1 = ss (3' proxy base 'T').
    // cell index = U5<<3 | U3<<2 | V5<<1 | V3  (terminal 5'/3', interior 5'/3').
    uint64_t j[NBIN][2][16] = {};

    static int bin(int L) {
        if (L < 50)  return 0;
        if (L < 70)  return 1;
        if (L < 90)  return 2;
        if (L < 120) return 3;
        return 4;
    }
    void merge(const J53Acc& o) {
        for (int b = 0; b < NBIN; ++b)
            for (int c = 0; c < 2; ++c)
                for (int k = 0; k < 16; ++k) j[b][c][k] += o.j[b][c][k];
    }
};

// Accumulate one read (L >= LSD_L_MIN == 30 guaranteed by caller).
inline void accumulate(J53Acc& a, const std::string& seq, int L) {
    if (L < 30) return;
    const int b  = J53Acc::bin(L);
    // Quarter-point interior anchors: deep enough that terminal deamination
    // decay (lambda~0.2 -> <2% by position 8) cannot leak into the null.
    const int i5 = L / 4;
    const int i3 = L - 1 - L / 4;
    if (i3 <= i5) return;                          // unreachable for L>=30
    const char c5  = seq[0]   & ~0x20;
    const char c3  = seq[L-1] & ~0x20;
    const char ci5 = seq[i5]  & ~0x20;
    const char ci3 = seq[i3]  & ~0x20;
    const int U5 = (c5  == 'T');
    const int V5 = (ci5 == 'T');
    ++a.j[b][0][(U5<<3) | ((c3 =='A')<<2) | (V5<<1) | (ci3 =='A')];  // ds: 3'->A
    ++a.j[b][1][(U5<<3) | ((c3 =='T')<<2) | (V5<<1) | (ci3 =='T')];  // ss: 3'->T
}

struct BinResult {
    int len_lo = 0, len_hi = 0;
    uint64_t n = 0;
    double p5 = 0, p3 = 0, pv5 = 0, pv3 = 0;
    double cov_term = 0, cov_int = 0, delta_cov = 0, se = 0, z = 0;
    double or_term = 0, or_int = 0;
};
struct J53Result {
    bool valid = false;
    bool is_ss = false;
    std::array<BinResult, J53Acc::NBIN> bins{};
    int  nbins = 0;                                // populated bins
    int  primary_bin = -1;
    double primary_delta_cov = 0, primary_z = 0;
    bool   authenticated = false;
    double pi_at_d[3] = {0, 0, 0};                 // conditional pi at d = 0.2, 0.4, 0.6
};

inline double odds_ratio(double n11, double n10, double n01, double n00) {
    return ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5));
}

// Detection gate parameters. z is trivially huge at 1e8 reads, so the real
// guard is CROSS-BIN STABILITY (delta_cov>0 in every populated bin) plus a
// conservative z floor (Bonferroni-safe over 5 bins).
inline J53Result compute(const J53Acc& a, bool is_ss) {
    J53Result r; r.is_ss = is_ss;
    const int ch = is_ss ? 1 : 0;
    const uint64_t N_MIN = 10000;
    int    best = -1; double best_dc = 0;
    int    n_pos_bins = 0, n_tested = 0;
    for (int b = 0; b < J53Acc::NBIN; ++b) {
        const uint64_t* c = a.j[b][ch];
        uint64_t n = 0; for (int k = 0; k < 16; ++k) n += c[k];
        BinResult br;
        br.len_lo = J53Acc::LEN_LO[b]; br.len_hi = J53Acc::LEN_HI[b]; br.n = n;
        if (n < N_MIN) { r.bins[b] = br; continue; }
        const double dn = static_cast<double>(n);
        // marginal & joint counts from the 16-cell joint
        double nU5 = 0, nU3 = 0, nV5 = 0, nV3 = 0, nU = 0, nV = 0;   // U = U5&U3, V = V5&V3
        // 2x2 tables for odds ratios
        double tu[2][2] = {{0,0},{0,0}}, tv[2][2] = {{0,0},{0,0}};
        for (int k = 0; k < 16; ++k) {
            const double m = static_cast<double>(c[k]);
            if (m == 0) continue;
            const int U5 = (k>>3)&1, U3 = (k>>2)&1, V5 = (k>>1)&1, V3 = k&1;
            nU5 += U5*m; nU3 += U3*m; nV5 += V5*m; nV3 += V3*m;
            nU  += (U5&U3)*m; nV += (V5&V3)*m;
            tu[U5][U3] += m; tv[V5][V3] += m;
        }
        const double p5 = nU5/dn, p3 = nU3/dn, pv5 = nV5/dn, pv3 = nV3/dn;
        const double cov_t = nU/dn - p5*p3;
        const double cov_v = nV/dn - pv5*pv3;
        const double dcov  = cov_t - cov_v;
        // Paired per-read SE: g = (U5-p5)(U3-p3) - (V5-pv5)(V3-pv3); computed on
        // the SAME reads, so the terminal/interior estimators are positively
        // dependent and this variance is smaller (not larger) than independent.
        double sg2 = 0;
        for (int k = 0; k < 16; ++k) {
            const double m = static_cast<double>(c[k]);
            if (m == 0) continue;
            const int U5 = (k>>3)&1, U3 = (k>>2)&1, V5 = (k>>1)&1, V3 = k&1;
            const double g = (U5-p5)*(U3-p3) - (V5-pv5)*(V3-pv3);
            sg2 += m * (g - dcov) * (g - dcov);
        }
        const double se = std::sqrt(sg2 / dn / dn);
        br.p5=p5; br.p3=p3; br.pv5=pv5; br.pv3=pv3;
        br.cov_term=cov_t; br.cov_int=cov_v; br.delta_cov=dcov;
        br.se=se; br.z = se>0 ? dcov/se : 0.0;
        br.or_term = odds_ratio(tu[1][1], tu[1][0], tu[0][1], tu[0][0]);
        br.or_int  = odds_ratio(tv[1][1], tv[1][0], tv[0][1], tv[0][0]);
        r.bins[b] = br;
        ++n_tested; if (dcov > 0) ++n_pos_bins;
        if (dcov > best_dc) { best_dc = dcov; best = b; }
    }
    r.nbins = n_tested;
    if (n_tested == 0) { r.valid = false; return r; }
    r.valid = true;
    r.primary_bin = best;
    if (best >= 0) { r.primary_delta_cov = r.bins[best].delta_cov; r.primary_z = r.bins[best].z; }
    // Authentication: positive in EVERY tested bin (stability) AND conservative z.
    r.authenticated = (best >= 0) && (n_pos_bins == n_tested)
                      && (r.primary_z > 6.0)
                      && (r.bins[best].or_term > r.bins[best].or_int);
    // Conditional pi at assumed d (d5=d3=d): A5A3 = (1-pv5)(1-pv3) d^2,
    // pi(1-pi) = delta_cov / A5A3, pi = smaller root in [0,0.5].
    if (best >= 0 && r.primary_delta_cov > 0) {
        const BinResult& pb = r.bins[best];
        const double ds[3] = {0.2, 0.4, 0.6};
        for (int i = 0; i < 3; ++i) {
            const double A5A3 = (1.0 - pb.pv5) * (1.0 - pb.pv3) * ds[i] * ds[i];
            const double q = A5A3 > 0 ? r.primary_delta_cov / A5A3 : 2.0;  // = pi(1-pi)
            r.pi_at_d[i] = (q <= 0.25 && q >= 0) ? 0.5 * (1.0 - std::sqrt(1.0 - 4.0*q))
                                                 : std::nan("");           // d too small for observed cov
        }
    }
    return r;
}

// Render a `"cooccur53": {...}` JSON fragment (no leading/trailing comma) for
// splicing into the profile --json via ProfileJsonInput::cooccur53_json.
inline std::string to_json_fragment(const J53Result& r) {
    auto jv = [](double v) -> std::string {
        if (std::isnan(v) || std::isinf(v)) return "null";
        std::ostringstream o; o << std::setprecision(8) << v; return o.str();
    };
    std::ostringstream o;
    // delta_cov = pi(1-pi)*A5*A3 identifies only the PRODUCT; pi and d are not
    // separable from one bin. The terminal marginal-excess route (pi = p5^2/p_both)
    // is defeated here by the end-prep composition fingerprint (terminal T is
    // DEPLETED below interior background), so pi is reported only conditionally on
    // an assumed d. Downstream must treat a bulk-cold cooccur53 detection as
    // "damage present, point magnitude not identifiable" -- never as a magnitude 0.
    o << "  \"cooccur53\": {\n"
      << "    \"channel\": \"" << (r.is_ss ? "ss" : "ds") << "\",\n"
      << "    \"valid\": " << (r.valid ? "true" : "false") << ",\n"
      << "    \"authenticated\": " << (r.authenticated ? "true" : "false") << ",\n"
      << "    \"identifiability\": \"NOT_IDENTIFIABLE_POINT\",\n"
      << "    \"point_damaged_fraction\": null,\n"
      << "    \"primary_delta_cov\": " << jv(r.primary_delta_cov) << ",\n"
      << "    \"primary_z\": " << jv(r.primary_z) << ",\n"
      << "    \"pi_at_d02\": " << jv(r.pi_at_d[0]) << ",\n"
      << "    \"pi_at_d04\": " << jv(r.pi_at_d[1]) << ",\n"
      << "    \"pi_at_d06\": " << jv(r.pi_at_d[2]) << ",\n"
      << "    \"bins\": [";
    bool first = true;
    for (int b = 0; b < J53Acc::NBIN; ++b) {
        const BinResult& br = r.bins[b];
        if (br.n < 10000) continue;
        if (!first) o << ","; first = false;
        o << "\n      {\"len_lo\":" << br.len_lo << ",\"len_hi\":" << br.len_hi
          << ",\"n\":" << br.n
          << ",\"delta_cov\":" << jv(br.delta_cov) << ",\"se\":" << jv(br.se)
          << ",\"z\":" << jv(br.z)
          << ",\"cov_term\":" << jv(br.cov_term) << ",\"cov_int\":" << jv(br.cov_int)
          << ",\"or_term\":" << jv(br.or_term) << ",\"or_int\":" << jv(br.or_int)
          << ",\"p5\":" << jv(br.p5) << ",\"p3\":" << jv(br.p3) << "}";
    }
    o << "]\n  }";
    return o.str();
}

// ---------------------------------------------------------------------------
// OVERHANG-IMMUNE variant: recover (pi, d) by modelling the overhang latent.
//
// The 3rd-order triple was FALSIFIED on real data (FLB16mNdAds1): a second
// latent -- single-strand overhang/tail length -- co-deaminates positions
// within a few bp of the SAME terminus, skewing every same-end moment. With
// only two ends, every triple carries a same-end pair, so cross-triple
// concordance (C=1.276<2) failed to catch it: all triples shared the same 5'
// overhang and agreed on the wrong pi (0.026 -> d3=1.03, impossible).
//
// Fix (approach b, ratified by fable5): make the overhang OBSERVABLE and
// calibrate it from the immune channel. Two latents: ancient w.p. pi; given
// ancient, geometric overhangs L5,L3 (survival s), INDEPENDENT across ends;
// covered_p <=> L>=p+1, so P(covered_p|anc)=s^{p+1}; a covered base deaminates
// w.p. beta. Interior-differencing removes composition, leaving:
//   cross-end (immune, opposite ends independent | ancient):
//     dcov(5'p,3'q) = pi(1-pi) beta5 beta3 s5^{p+1} s3^{q+1}
//   same-end 5' (overhang-contaminated, now explicitly modelled):
//     Cov5(0,b)     = pi beta5^2 s5^{b+1} (1 - pi s5)
// s5,s3 come from the IMMUNE cross-end decay; pi from
//   K = P'^2/(W5 W3) = (1-pi)^2 / [(1-pi s5)(1-pi s3)],   monotone in (0,1/2)
// with P'=P/(s5 s3), W5=Cov5(0,b)/s5^{b+1}, W3=Cov3(0,d)/s3^{d+1}. Then
// beta5=sqrt(W5/[pi(1-pi s5)]); the metaDMG-comparable terminal C->T rate is
// s5*beta5 (P(covered)*rate), NOT beta5.
//
// Two honesty guards (the part the triple lacked):
//  (1) s5_cross = [dcov(5'3,3'0)/dcov(5'0,3'0)]^{1/3}  (immune channel)
//      s5_shape = [Cov5(0,6)/Cov5(0,3)]^{1/3}          (contaminated, pi-free)
//      must agree: |ln s5_cross - ln s5_shape| > 3 se(Delta) -> REFUTE. The two
//      channels probe the SAME geometric s5; a mismatch is the overhang.
//  (2) K solved at two spacings: pi(K, b=3) vs pi(K, b=6); |Dpi| > 2 se -> REFUTE.
// Any guard trips, or beta not in (0,1], or pi not in (0,1/2) => detection-only.
//
// feature bits: 0=U5_0 1=U5_3 2=U5_6 3=U3_0 4=U3_3 5=U3_6
//               6=V5_0 7=V5_3 8=V5_6 9=V3_0 10=V3_3 11=V3_6
struct J53OhAcc {
    static constexpr int NBIN = J53Acc::NBIN;
    uint64_t t[NBIN][2][4096] = {};
    void merge(const J53OhAcc& o) {
        for (int b = 0; b < NBIN; ++b)
            for (int ch = 0; ch < 2; ++ch)
                for (int k = 0; k < 4096; ++k) t[b][ch][k] += o.t[b][ch][k];
    }
};

inline void accumulate_oh(J53OhAcc& a, const std::string& seq, int L) {
    if (L < 30) return;                            // L>=30 => terminal{0,3,6} & interior blocks disjoint
    const int b  = J53Acc::bin(L);
    const int i5 = L / 4;                           // 5' interior anchor
    const int i3 = L - 1 - L / 4;                   // 3' interior anchor
    auto is = [&](int idx, char base) { return (seq[idx] & ~0x20) == base ? 1 : 0; };
    const int u5_0 = is(0,'T'), u5_3 = is(3,'T'), u5_6 = is(6,'T');
    const int v5_0 = is(i5,'T'), v5_3 = is(i5+3,'T'), v5_6 = is(i5+6,'T');
    // ds channel: 3'->A ; ss channel: 3'->T (matches accumulate() convention)
    {   const int u3_0 = is(L-1,'A'), u3_3 = is(L-4,'A'), u3_6 = is(L-7,'A');
        const int v3_0 = is(i3,'A'), v3_3 = is(i3-3,'A'), v3_6 = is(i3-6,'A');
        ++a.t[b][0][u5_0 | (u5_3<<1) | (u5_6<<2) | (u3_0<<3) | (u3_3<<4) | (u3_6<<5)
                  | (v5_0<<6) | (v5_3<<7) | (v5_6<<8) | (v3_0<<9) | (v3_3<<10) | (v3_6<<11)]; }
    {   const int u3_0 = is(L-1,'T'), u3_3 = is(L-4,'T'), u3_6 = is(L-7,'T');
        const int v3_0 = is(i3,'T'), v3_3 = is(i3-3,'T'), v3_6 = is(i3-6,'T');
        ++a.t[b][1][u5_0 | (u5_3<<1) | (u5_6<<2) | (u3_0<<3) | (u3_3<<4) | (u3_6<<5)
                  | (v5_0<<6) | (v5_3<<7) | (v5_6<<8) | (v3_0<<9) | (v3_3<<10) | (v3_6<<11)]; }
}

struct J53OhResult {
    bool   valid = false, is_ss = false, identifiable = false;
    int    primary_bin = -1;
    const char* reason = "";                       // why not identifiable (empty if identified)
    double pi = std::nan(""), se_pi = std::nan("");
    double pi_k6 = std::nan(""), se_pi_diff = std::nan("");
    double beta5 = std::nan(""), beta3 = std::nan("");   // per-covered-base deamination rate
    double term5 = std::nan(""), term3 = std::nan("");   // s*beta = terminal C->T (compare metaDMG 0.41)
    double s5 = std::nan(""), s3 = std::nan("");          // overhang survival, immune (cross-end)
    double s5_shape = std::nan(""), s3_shape = std::nan("");
    double dln_s5 = std::nan(""), dln_s3 = std::nan("");
    double se_dln_s5 = std::nan(""), se_dln_s3 = std::nan("");
    double dcov[7] = {std::nan(""),std::nan(""),std::nan(""),std::nan(""),
                      std::nan(""),std::nan(""),std::nan("")};
};

inline J53OhResult compute_oh(const J53OhAcc& a, bool is_ss) {
    J53OhResult r; r.is_ss = is_ss;
    const int ch = is_ss ? 1 : 0;
    const uint64_t N_MIN = 50000;
    int pb = -1;
    for (int b = 0; b < J53OhAcc::NBIN; ++b) {     // shortest populated bin = strongest damage
        uint64_t n = 0; for (int k = 0; k < 4096; ++k) n += a.t[b][ch][k];
        if (n >= N_MIN) { pb = b; break; }
    }
    if (pb < 0) { r.reason = "no_bin"; return r; }
    r.valid = true; r.primary_bin = pb;
    const uint64_t* c = a.t[pb][ch];
    double n = 0; for (int k = 0; k < 4096; ++k) n += static_cast<double>(c[k]);

    double p[12] = {0};
    for (int k = 0; k < 4096; ++k) { double m = c[k]; if (!m) continue;
        for (int f = 0; f < 12; ++f) if ((k>>f)&1) p[f] += m; }
    for (int f = 0; f < 12; ++f) p[f] /= n;

    // 7 interior-differenced pair covariances (tx,ty terminal bits; vx,vy interior).
    // indices: 0=D_0030 1=D_5330 2=D_0033 3=C5_03 4=C5_06 5=C3_03 6=C3_06
    struct PR { int tx, ty, vx, vy; };
    const PR pr[7] = { {0,3, 6,9},  {1,3, 7,9},  {0,4, 6,10},
                       {0,1, 6,7},  {0,2, 6,8},  {3,4, 9,10}, {3,5, 9,11} };
    static thread_local double H[7][4096];         // per-cell h (pre-mean-subtraction)
    double& D0030 = r.dcov[0]; double& D5330 = r.dcov[1]; double& D0033 = r.dcov[2];
    double& C5_03 = r.dcov[3]; double& C5_06 = r.dcov[4];
    double& C3_03 = r.dcov[5]; double& C3_06 = r.dcov[6];
    for (int pi_ = 0; pi_ < 7; ++pi_) {
        const PR& q = pr[pi_];
        const double px=p[q.tx], py=p[q.ty], pvx=p[q.vx], pvy=p[q.vy];
        double s = 0;
        for (int k = 0; k < 4096; ++k) { double m = c[k]; if (!m) { H[pi_][k]=0; continue; }
            const double h = (((k>>q.tx)&1)-px)*(((k>>q.ty)&1)-py)
                           - (((k>>q.vx)&1)-pvx)*(((k>>q.vy)&1)-pvy);
            H[pi_][k] = h; s += m*h;
        }
        r.dcov[pi_] = s / n;
    }

    bool pos = true; for (int i = 0; i < 7; ++i) if (!(r.dcov[i] > 0)) pos = false;
    if (!pos) { r.reason = "nonpos_dcov"; return r; }

    // overhang survival: immune (cross-end decay) and contaminated (same-end shape)
    r.s5 = std::cbrt(D5330 / D0030);   r.s3 = std::cbrt(D0033 / D0030);
    r.s5_shape = std::cbrt(C5_06 / C5_03);   r.s3_shape = std::cbrt(C3_06 / C3_03);
    r.dln_s5 = std::log(r.s5) - std::log(r.s5_shape);
    r.dln_s3 = std::log(r.s3) - std::log(r.s3_shape);
    const double s5 = r.s5, s3 = r.s3;

    // per-cell influence phi_i = H_i - dcov_i; assemble linear log-combos, then se.
    auto se_combo = [&](const double w[7]) -> double {   // se of sum_i w_i * ln(dcov_i)
        double s2 = 0;
        for (int k = 0; k < 4096; ++k) { double m = c[k]; if (!m) continue;
            double phi = 0; for (int i = 0; i < 7; ++i) if (w[i] != 0.0)
                phi += w[i] * (H[i][k] - r.dcov[i]) / r.dcov[i];
            s2 += m * phi * phi;
        }
        return std::sqrt(s2) / n;
    };
    // Delta ln s5 = 1/3 (lnD5330 - lnD0030 - lnC5_06 + lnC5_03)
    { const double w[7] = {-1.0/3, 1.0/3, 0, 1.0/3, -1.0/3, 0, 0}; r.se_dln_s5 = se_combo(w); }
    { const double w[7] = {-1.0/3, 0, 1.0/3, 0, 0, 1.0/3, -1.0/3}; r.se_dln_s3 = se_combo(w); }

    // Newton solve K = (1-pi)^2 / [(1-pi s5)(1-pi s3)] for pi in (0,1/2)
    auto gprime = [&](double pi) { return -2.0/(1-pi) + s5/(1-pi*s5) + s3/(1-pi*s3); };
    auto solve_pi = [&](double lnK) -> double {
        double pi = 0.196;
        for (int it = 0; it < 60; ++it) {
            const double g = 2*std::log(1-pi) - std::log(1-pi*s5) - std::log(1-pi*s3) - lnK;
            const double gp = gprime(pi);
            double np = pi - g/gp;
            if (np < 1e-5) np = 1e-5; if (np > 0.4999) np = 0.4999;
            const double step = np - pi; pi = np;
            if (std::fabs(step) < 1e-10) break;
        }
        return pi;
    };
    // lnK3 = 2/3(lnD0030+lnD5330+lnD0033) - lnC5_03 - lnC3_03
    // lnK6 = -4/3 lnD0030 + 5/3 lnD5330 + 5/3 lnD0033 - lnC5_06 - lnC3_06
    const double lnK3 = (2.0/3)*(std::log(D0030)+std::log(D5330)+std::log(D0033))
                        - std::log(C5_03) - std::log(C3_03);
    const double lnK6 = (-4.0/3)*std::log(D0030) + (5.0/3)*std::log(D5330)
                        + (5.0/3)*std::log(D0033) - std::log(C5_06) - std::log(C3_06);
    r.pi    = solve_pi(lnK3);
    r.pi_k6 = solve_pi(lnK6);

    // se(pi3) and se(pi3-pi6) via stacked per-read influences. phi_pi = (1/g') phi_lnK.
    const double w3[7] = { 2.0/3, 2.0/3, 2.0/3, -1, 0, -1, 0 };   // d(lnK3)/d(ln dcov_i)
    const double w6[7] = { -4.0/3, 5.0/3, 5.0/3, 0, -1, 0, -1 };  // d(lnK6)/d(ln dcov_i)
    const double inv3 = 1.0 / gprime(r.pi), inv6 = 1.0 / gprime(r.pi_k6);
    double s2_3 = 0, s2_d = 0;
    for (int k = 0; k < 4096; ++k) { double m = c[k]; if (!m) continue;
        double lk3 = 0, lk6 = 0;
        for (int i = 0; i < 7; ++i) { const double e = (H[i][k]-r.dcov[i])/r.dcov[i];
            lk3 += w3[i]*e; lk6 += w6[i]*e; }
        const double phi3 = inv3*lk3, phi6 = inv6*lk6;
        s2_3 += m*phi3*phi3; s2_d += m*(phi3-phi6)*(phi3-phi6);
    }
    r.se_pi = std::sqrt(s2_3)/n; r.se_pi_diff = std::sqrt(s2_d)/n;

    // recover beta (per-covered rate) at b=3; terminal C->T = s*beta (metaDMG 0.41)
    const double W5 = C5_03 / std::pow(s5, 4), W3 = C3_03 / std::pow(s3, 4);
    r.beta5 = std::sqrt(W5 / (r.pi * (1 - r.pi*s5)));
    r.beta3 = std::sqrt(W3 / (r.pi * (1 - r.pi*s3)));
    r.term5 = s5 * r.beta5; r.term3 = s3 * r.beta3;

    // honesty guards
    const bool pi_ok    = std::isfinite(r.pi) && r.pi > 1e-4 && r.pi < 0.4999;
    const bool beta_ok  = std::isfinite(r.beta5) && r.beta5 > 0 && r.beta5 <= 1.0
                       && std::isfinite(r.beta3) && r.beta3 > 0 && r.beta3 <= 1.0;
    const bool s5_ok    = std::fabs(r.dln_s5) <= 3.0 * r.se_dln_s5;
    const bool s3_ok    = std::fabs(r.dln_s3) <= 3.0 * r.se_dln_s3;
    const bool conc_ok  = std::fabs(r.pi - r.pi_k6) <= 2.0 * r.se_pi_diff;
    r.identifiable = pi_ok && beta_ok && s5_ok && s3_ok && conc_ok;
    if      (!pi_ok)   r.reason = "pi_out_of_range";
    else if (!beta_ok) r.reason = "beta_out_of_range";
    else if (!s5_ok)   r.reason = "overhang_s5_mismatch";
    else if (!s3_ok)   r.reason = "overhang_s3_mismatch";
    else if (!conc_ok) r.reason = "pi_spacing_discordant";
    return r;
}

inline std::string to_json_oh_fragment(const J53OhResult& r) {
    auto jv = [](double v) -> std::string {
        if (std::isnan(v) || std::isinf(v)) return "null";
        std::ostringstream o; o << std::setprecision(8) << v; return o.str();
    };
    std::ostringstream o;
    o << "  \"cooccur53_ohim\": {\n"
      << "    \"channel\": \"" << (r.is_ss ? "ss" : "ds") << "\",\n"
      << "    \"valid\": " << (r.valid ? "true" : "false") << ",\n"
      << "    \"identifiable\": " << (r.identifiable ? "true" : "false") << ",\n"
      << "    \"reason\": \"" << r.reason << "\",\n"
      << "    \"primary_bin\": " << r.primary_bin << ",\n"
      << "    \"pi\": " << jv(r.pi) << ",\n"
      << "    \"se_pi\": " << jv(r.se_pi) << ",\n"
      << "    \"pi_k6\": " << jv(r.pi_k6) << ",\n"
      << "    \"se_pi_diff\": " << jv(r.se_pi_diff) << ",\n"
      << "    \"beta5\": " << jv(r.beta5) << ",\n"
      << "    \"beta3\": " << jv(r.beta3) << ",\n"
      << "    \"term5_metaDMG\": " << jv(r.term5) << ",\n"
      << "    \"term3_metaDMG\": " << jv(r.term3) << ",\n"
      << "    \"s5\": " << jv(r.s5) << ",\n"
      << "    \"s3\": " << jv(r.s3) << ",\n"
      << "    \"s5_shape\": " << jv(r.s5_shape) << ",\n"
      << "    \"s3_shape\": " << jv(r.s3_shape) << ",\n"
      << "    \"dln_s5\": " << jv(r.dln_s5) << ",\n"
      << "    \"se_dln_s5\": " << jv(r.se_dln_s5) << ",\n"
      << "    \"dln_s3\": " << jv(r.dln_s3) << ",\n"
      << "    \"se_dln_s3\": " << jv(r.se_dln_s3) << "\n"
      << "  }";
    return o.str();
}

} // namespace fqdup_j53
