#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace fqdup_tailhead {
// Fixed before evaluation in commit 292e039; never fit to a dataset.
inline constexpr double odds_threshold = 3276800000.0;
inline constexpr int max_alternatives = 4096;
inline int base_index(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3; default: return -1; }
}
struct EvidenceRead {
    std::string seq, qual;
    std::vector<bool> usable;
};
inline EvidenceRead reverse_read(const EvidenceRead& r) {
    EvidenceRead out = r;
    std::reverse(out.seq.begin(), out.seq.end());
    std::reverse(out.qual.begin(), out.qual.end());
    std::reverse(out.usable.begin(), out.usable.end());
    for (char& c : out.seq) {
        int b = base_index(c);
        c = b < 0 ? 'N' : "TGCA"[b];
    }
    return out;
}
template<class LowComplexity>
EvidenceRead evidence_read(const std::string& seq, const std::string& qual,
                           LowComplexity low_complexity) {
    EvidenceRead r{seq, qual, std::vector<bool>(seq.size(), true)};
    for (int i = 0; i + 11 <= static_cast<int>(seq.size()); ++i)
        if (low_complexity(seq.data() + i, 11))
            for (int j = i; j < i + 11; ++j) r.usable[j] = false;
    for (size_t i = 0; i < seq.size(); ++i)
        if (base_index(seq[i]) < 0 || static_cast<unsigned char>(qual[i]) < 53)
            r.usable[i] = false;
    return r;
}
inline double score(const EvidenceRead& a, int ai, const EvidenceRead& b, int bi,
                    int n, const std::array<double, 4>& composition) {
    static const std::array<double, 94> errors = [] {
        std::array<double, 94> out{};
        for (int q = 0; q < 94; ++q) out[q] = std::pow(10.0, -q / 10.0);
        return out;
    }();
    double llr = 0;
    for (int k = 0; k < n; ++k) {
        int i = ai + k, j = bi + k;
        if (!a.usable[i] || !b.usable[j]) continue;
        int x = base_index(a.seq[i]), y = base_index(b.seq[j]);
        int qa = std::min(93, static_cast<int>(static_cast<unsigned char>(a.qual[i])) - 33);
        int qb = std::min(93, static_cast<int>(static_cast<unsigned char>(b.qual[j])) - 33);
        double ea = errors[qa], eb = errors[qb];
        double shared = 0, pa = 0, pb = 0;
        for (int t = 0; t < 4; ++t) {
            double ca = x == t ? 1 - ea : ea / 3;
            double cb = y == t ? 1 - eb : eb / 3;
            shared += composition[t] * ca * cb;
            pa += composition[t] * ca;
            pb += composition[t] * cb;
        }
        // Disagreements contribute their full negative log evidence.
        llr += std::log(shared / (pa * pb));
    }
    return llr;
}
template<class LowComplexity>
bool unique(const std::string& r1, const std::string& q1,
            const std::string& rc2, const std::string& q2,
            int ov, int min_ov, LowComplexity low_complexity) {
    int l1 = static_cast<int>(r1.size()), l2 = static_cast<int>(rc2.size());
    if (q1.size() != r1.size() || q2.size() != rc2.size() ||
        ov < min_ov || ov > std::min(l1, l2)) return false;
    std::array<EvidenceRead, 4> reads;
    reads[0] = evidence_read(r1, q1, low_complexity);
    reads[1] = evidence_read(rc2, q2, low_complexity);
    reads[2] = reverse_read(reads[0]);
    reads[3] = reverse_read(reads[1]);
    std::array<double, 4> composition{{1, 1, 1, 1}};
    double total = 4;
    for (int s = 0; s < 2; ++s)
        for (char c : reads[s].seq) {
            int b = base_index(c);
            if (b >= 0) { ++composition[b]; ++total; }
        }
    for (double& p : composition) p /= total;
    const double bound = std::log(odds_threshold);
    const double best = score(reads[0], l1 - ov, reads[1], 0, ov, composition);
    if (!(best > bound)) return false;
    int alternatives = 0;
    auto inferior = [&](const EvidenceRead& a, int ai, const EvidenceRead& b,
                        int bi, int n) {
        if (++alternatives > max_alternatives) return false;
        return best - score(a, ai, b, bi, n, composition) > bound;
    };
    // All relative offsets, in both orientations, not only mismatch-gate survivors.
    for (int orientation : {1, 3})
        for (int shift = min_ov - l2; shift <= l1 - min_ov; ++shift) {
            if (orientation == 1 && shift == l1 - ov) continue;
            int ai = std::max(0, shift), bi = std::max(0, -shift);
            int n = std::min(l1 - ai, l2 - bi);
            if (!inferior(reads[0], ai, reads[orientation], bi, n)) return false;
        }
    // Either observed overlap segment must also be unique within the two reads.
    // Omit only the two windows supporting the selected physical placement.
    for (int query = 0; query < 2; ++query) {
        int qi = query == 0 ? l1 - ov : 0;
        for (int target = 0; target < 4; ++target)
            for (int j = 0; j + ov <= static_cast<int>(reads[target].seq.size()); ++j) {
                if ((target == 0 && j == l1 - ov) || (target == 1 && j == 0)) continue;
                if (!inferior(reads[query], qi, reads[target], j, ov)) return false;
            }
    }
    return true;
}
} // namespace fqdup_tailhead
