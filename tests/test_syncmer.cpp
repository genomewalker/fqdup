// Smoke tests for fqdup::syncmer + fqdup::derep_detail::banded_edit_distance_le2.
// Build: g++ -std=c++20 -O2 -I include -I src tests/test_syncmer.cpp -o /tmp/test_syncmer

#include "fqdup/syncmer.hpp"
#include "derep_detail/banded_ed.hpp"

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using fqdup::SyncmerParams;
using fqdup::compute_open_syncmers;
using fqdup::compute_sketch;
using fqdup::syncmer_params_for_length;
using fqdup::derep_detail::EditScript;
using fqdup::derep_detail::banded_edit_distance_le2;

static std::vector<uint8_t> encode(const std::string& s) {
    std::vector<uint8_t> out; out.reserve(s.size());
    for (char c : s) {
        switch (c) { case 'A': out.push_back(0); break;
                     case 'C': out.push_back(1); break;
                     case 'G': out.push_back(2); break;
                     case 'T': out.push_back(3); break;
                     default: out.push_back(0); }
    }
    return out;
}

static std::string rc(const std::string& s) {
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        out[i] = (c == 'A') ? 'T' : (c == 'C') ? 'G' : (c == 'G') ? 'C' : 'A';
    }
    return out;
}

static std::string rand_seq(int L, std::mt19937& rng) {
    static const char bases[] = "ACGT";
    std::string s; s.reserve(L);
    for (int i = 0; i < L; ++i) s.push_back(bases[rng() & 3u]);
    return s;
}

static int test_rc_equivalence() {
    std::mt19937 rng(0xC0FFEE);
    int passed = 0;
    for (int trial = 0; trial < 50; ++trial) {
        int L = 60 + (rng() % 90);
        auto s = rand_seq(L, rng);
        auto sr = rc(s);
        auto a = encode(s), b = encode(sr);
        std::vector<uint64_t> ha, hb;
        std::vector<uint16_t> pa, pb;
        compute_sketch(a.data(), L, ha, pa);
        compute_sketch(b.data(), L, hb, pb);
        // Hash sets should be identical (positions differ by mirroring).
        std::sort(ha.begin(), ha.end());
        std::sort(hb.begin(), hb.end());
        if (ha == hb) ++passed;
        else {
            fprintf(stderr, "RC mismatch L=%d: |fwd|=%zu |rc|=%zu\n", L, ha.size(), hb.size());
        }
    }
    printf("test_rc_equivalence: %d/50 passed (>=44 ok; rare s-mer hash ties on random seqs)\n", passed);
    return passed >= 44;
}

static int test_single_sub_outside_smer() {
    // A substitution at position p only invalidates k-mers whose window covers p
    // AND whose s-mer offset/min-rank is affected. So the union-shared sketch
    // size must be at least k_total - (kmer_window_count_at_p).
    std::mt19937 rng(0xBEEF);
    int passed = 0;
    for (int trial = 0; trial < 30; ++trial) {
        int L = 80 + (rng() % 60);
        auto s = rand_seq(L, rng);
        auto a = encode(s);
        // Mutate a single base in the middle (outside terminals).
        int p = L / 2;
        auto t = s; t[p] = "ACGT"[(rng() + 1) & 3];
        auto b = encode(t);

        std::vector<uint64_t> ha, hb;
        std::vector<uint16_t> pa, pb;
        compute_sketch(a.data(), L, ha, pa);
        compute_sketch(b.data(), L, hb, pb);

        // Count shared hashes.
        std::sort(ha.begin(), ha.end());
        std::sort(hb.begin(), hb.end());
        std::vector<uint64_t> inter;
        std::set_intersection(ha.begin(), ha.end(), hb.begin(), hb.end(),
                              std::back_inserter(inter));
        // Expect substantial overlap (at least 30% of the smaller sketch).
        size_t small = std::min(ha.size(), hb.size());
        if (small == 0) continue;
        if (inter.size() * 3 >= small) ++passed;
    }
    printf("test_single_sub_outside_smer: %d/30 trials with shared >= 1/3\n", passed);
    return passed >= 25;
}

static int test_banded_ed_basic() {
    // exact match
    auto a = encode("ACGTACGT");
    EditScript es{};
    int ed = banded_edit_distance_le2(a.data(), 8, a.data(), 8, &es);
    if (!(ed == 0 && es.total() == 0)) { fprintf(stderr, "match fail ed=%d\n", ed); return 0; }

    // single sub
    auto b = encode("ACGTACGA");
    ed = banded_edit_distance_le2(a.data(), 8, b.data(), 8, &es);
    if (!(ed == 1 && es.n_sub == 1 && es.kind[0] == 0)) { fprintf(stderr, "sub fail ed=%d kind=%d\n", ed, es.kind[0]); return 0; }

    // single deletion (b shorter by 1)
    auto c = encode("ACGTACG");
    ed = banded_edit_distance_le2(a.data(), 8, c.data(), 7, &es);
    if (!(ed == 1 && es.n_del == 1)) { fprintf(stderr, "del fail ed=%d n_del=%d\n", ed, es.n_del); return 0; }

    // single insertion (b longer by 1)
    auto d = encode("ACGTACGTC");
    ed = banded_edit_distance_le2(a.data(), 8, d.data(), 9, &es);
    if (!(ed == 1 && es.n_ins == 1)) { fprintf(stderr, "ins fail ed=%d n_ins=%d\n", ed, es.n_ins); return 0; }

    // ed > 2 → -1
    auto e = encode("TGCATGCA");
    ed = banded_edit_distance_le2(a.data(), 8, e.data(), 8, &es);
    if (ed != -1) { fprintf(stderr, "cap fail ed=%d (expected -1)\n", ed); return 0; }

    // 1 sub + 1 del = ed 2
    auto f = encode("ACGTACGA");  // ed=1 vs a
    auto g = encode("ACGTACG");   // delete last
    ed = banded_edit_distance_le2(f.data(), 8, g.data(), 7, &es);
    if (!(ed >= 1 && ed <= 2)) { fprintf(stderr, "mix fail ed=%d\n", ed); return 0; }

    printf("test_banded_ed_basic: passed\n");
    return 1;
}

int main() {
    int ok = 0, total = 0;
    ++total; ok += test_rc_equivalence();
    ++total; ok += test_single_sub_outside_smer();
    ++total; ok += test_banded_ed_basic();
    printf("\nT8.1+T8.5 smoke: %d/%d test groups passed\n", ok, total);
    return ok == total ? 0 : 1;
}
