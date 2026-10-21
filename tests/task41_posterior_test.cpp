#include "tailhead_posterior.hpp"
#include <cassert>
#include <iostream>

int main() {
    using namespace fqdup_tailhead;
    auto no_complexity = [](const char*, int) { return false; };
    auto all_complexity = [](const char*, int) { return true; };
    const std::array<double, 4> uniform{{.25, .25, .25, .25}};
    auto a = evidence_read("ACGTACGTACGT", "IIIIIIIIIIII", no_complexity);
    auto b = a;
    double d = 2e-4 - (4.0 / 3) * 1e-8;
    assert(std::abs(score(a, 0, b, 0, 12, uniform) - 12 * std::log(4 * (1-d))) < 1e-9);
    b.seq[0] = 'T';
    assert(std::abs(score(a, 0, b, 0, 12, uniform) -
        (11 * std::log(4 * (1-d)) + std::log(4*d/3))) < 1e-9);
    auto masked = evidence_read("NNNNNNNNNNNN", "IIIIIIIIIIII", no_complexity);
    assert(score(a, 0, masked, 0, 12, uniform) == 0);
    masked = evidence_read(a.seq, "444444444444", no_complexity); // Q19
    assert(score(a, 0, masked, 0, 12, uniform) == 0);
    masked = evidence_read(a.seq, a.qual, all_complexity);
    assert(score(a, 0, masked, 0, 12, uniform) == 0);
    assert(reverse_read(reverse_read(a)).seq == a.seq);
    assert(!unique(std::string(100, 'A'), std::string(100, 'I'),
                   std::string(100, 'A'), std::string(100, 'I'), 19, 11, no_complexity));
    // Deterministic independent flanks and one known shared 19-base segment.
    unsigned state = 139;
    auto sequence = [&](int n) {
        std::string s;
        while (n--) { state ^= state << 13; state ^= state >> 17; state ^= state << 5;
            s += "ACGT"[state & 3]; }
        return s;
    };
    int accepted = 0;
    for (int k = 0; k < 100; ++k) {
        auto overlap = sequence(19), r1 = sequence(60) + overlap;
        auto r2 = overlap + sequence(60);
        std::string q(79, 'I');
        bool pass = unique(r1, q, r2, q, 19, 11, no_complexity);
        accepted += pass;
        if (pass) {
            // A second exact occurrence of the segment inside either read vetoes it.
            r1.replace(0, 19, overlap);
            assert(!unique(r1, q, r2, q, 19, 11, no_complexity));
        }
    }
    assert(accepted > 0);
    std::cout << "posterior model, masking and repeat tests PASS; positive cases=" << accepted << '\n';
}
