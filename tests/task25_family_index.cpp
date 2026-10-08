#include "../src/merge_r8qkeep/family_index.hpp"
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
using fqdup::r8qkeep::FamilyPartIndex;
// Independent frozen enumeration: expand query Ns, then mutate <=2-ns
// non-N positions over the same five-code alphabet as ad1e1ca.
static std::set<std::string> old_candidates(const std::string& query) {
    std::vector<std::string> bases{query}; std::vector<int> normal;
    int ns = 0; const std::string alphabet = "ACGTN";
    for (int p = 0; p < int(query.size()); ++p) {
        if (query[p] != 'N') { normal.push_back(p); continue; }
        ++ns; std::vector<std::string> next;
        for (auto base : bases) for (char c : alphabet) { base[p] = c; next.push_back(base); }
        bases = std::move(next);
    }
    std::set<std::string> out;
    for (const auto& base : bases) {
        out.insert(base); if (ns == 2) continue;
        for (std::size_t a = 0; a < normal.size(); ++a) for (char c : alphabet) {
            const int p = normal[a]; if (c == query[p]) continue;
            auto first = base; first[p] = c; out.insert(first); if (ns == 1) continue;
            for (std::size_t b = a+1; b < normal.size(); ++b) for (char d : alphabet) {
                const int q = normal[b]; if (d == query[q]) continue;
                auto second = first; second[q] = d; out.insert(second);
            }
        }
    }
    return out;
}
int main() {
    try {
        std::mt19937 rng(2501); const std::string bases = "ACGT";
        for (int length : {13, 14}) for (int ns = 0; ns <= 2; ++ns) for (int trial = 0; trial < 12; ++trial) {
            std::string query(length, 'A'); for (char& c : query) c = bases[rng()%4];
            for (int n = 0; n < ns;) { auto& c = query[rng()%length]; if (c != 'N') { c = 'N'; ++n; } }
            const auto expected_values = old_candidates(query);
            std::vector<std::string> keys(expected_values.begin(), expected_values.end());
            // Every enumerated neighbor, duplicate part postings, and negatives.
            keys.push_back(keys.front());
            for (int i = 0; i < 100; ++i) {
                std::string k(length, 'A'); for (char& c : k) c = "ACGTN"[rng()%5]; keys.push_back(k);
            }
            std::vector<int> representatives(keys.size()); std::iota(representatives.begin(), representatives.end(), 0);
            FamilyPartIndex index; index.build(keys, representatives, 0, length);
            std::vector<int> actual, expected;
            std::uint64_t lookups = 0, candidates = 0;
            index.visit(FamilyPartIndex::pack(query, 0, length), [&](int i) { actual.push_back(i); return true; }, lookups, candidates);
            std::sort(actual.begin(), actual.end());
            for (int i = 0; i < int(keys.size()); ++i) if (expected_values.count(keys[i])) expected.push_back(i);
            if (actual != expected || lookups > 4) throw std::runtime_error("candidate set mismatch");
            // The exhaustive oracle includes mutations in the same sub-block
            // and separate sub-blocks, N/N differences, and duplicate members.
            int visits=0;lookups=candidates=0;
            index.visit(FamilyPartIndex::pack(query,0,length),[&](int){++visits;return false;},lookups,candidates);
            if(visits!=1)throw std::runtime_error("visitor early termination");
            std::string too_many_n(length,'A');too_many_n[0]=too_many_n[4]=too_many_n[length-1]='N';
            visits=0;
            index.visit(FamilyPartIndex::pack(too_many_n,0,length),[&](int){++visits;return true;},lookups,candidates);
            if(visits)throw std::runtime_error("three query Ns accepted");
        }
        std::cout << "Task 25 candidate sets: 72 exhaustive-neighborhood comparisons passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
