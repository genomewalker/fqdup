#include "../src/merge_r8qkeep/family_index.hpp"
#include "task25_family_index_ref.hpp"
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
static std::uint64_t expected_lookups(const std::string& query,std::size_t cells){
    const int length=query.size();const std::array<int,5> bounds{0,4,length==13?7:8,length==13?10:11,length};
    std::uint64_t clean=0;
    for(int b=0;b<4;++b)if(query.substr(bounds[b],bounds[b+1]-bounds[b]).find('N')==std::string::npos)++clean;
    return cells<=65536?clean*(clean-1)/2:clean;
}
static void boundary_and_rebuilds(){
    FamilyPartIndex index;fqdup::r8qkeep::FrozenFamilyPartIndex frozen;
    for(std::size_t count:{1023,1024,1025,65535,65536,65537,0,65536,1,65537}){
        std::vector<std::string> keys;keys.reserve(count);
        for(std::size_t i=0;i<count;++i){
            std::uint64_t word=(i*2654435761ULL)&((1ULL<<28)-1);std::string key(14,'A');
            for(char& c:key){c="ACGT"[word&3];word>>=2;}keys.push_back(key);
        }
        std::vector<int> representatives(count);std::iota(representatives.begin(),representatives.end(),0);
        index.build(keys,representatives,0,14);frozen.build(keys,representatives,0,14);
        if(index.cells()!=count||index.bytes()>frozen.bytes())throw std::runtime_error("boundary retained storage");
        std::vector<std::string> queries{std::string(14,'A'),std::string(14,'T')};
        if(count){queries.push_back(keys.back());queries.push_back(keys[count/2]);}
        for(auto query:queries)for(int ns=0;ns<=3;++ns){
            if(ns)query[4*(ns-1)]='N';
            for(std::size_t stop:{std::size_t(1),std::size_t(3),count+1}){
                std::vector<int> actual,expected;std::uint64_t lookups=0,candidates=0,old_lookups=0,old_candidates=0;
                const auto packed=FamilyPartIndex::pack(query,0,14);
                index.visit(packed,[&](int i){actual.push_back(i);return actual.size()<stop;},lookups,candidates);
                frozen.visit(packed,[&](int i){expected.push_back(i);return expected.size()<stop;},old_lookups,old_candidates);
                if(actual!=expected||candidates!=old_candidates||lookups!=expected_lookups(query,count))
                    throw std::runtime_error("boundary/rebuild/prefix oracle mismatch");
            }
        }
    }
}
int main() {
    try {
        boundary_and_rebuilds();
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
            fqdup::r8qkeep::FrozenFamilyPartIndex frozen;frozen.build(keys,representatives,0,length);
            std::vector<int> actual, expected;
            std::uint64_t lookups = 0, candidates = 0;
            index.visit(FamilyPartIndex::pack(query, 0, length), [&](int i) { actual.push_back(i); return true; }, lookups, candidates);
            std::vector<int> frozen_order;std::uint64_t frozen_lookups=0,frozen_candidates=0;
            frozen.visit(FamilyPartIndex::pack(query,0,length),[&](int i){frozen_order.push_back(i);return true;},frozen_lookups,frozen_candidates);
            if(actual!=frozen_order||lookups!=expected_lookups(query,index.cells())||candidates!=frozen_candidates||index.bytes()>frozen.bytes())
                throw std::runtime_error("frozen traversal/counter/storage mismatch");
            for(std::size_t stop:{std::size_t(1),std::size_t(3),frozen_order.size()}){
                if(!stop)continue;
                std::vector<int> prefix,old_prefix;lookups=candidates=frozen_lookups=frozen_candidates=0;
                index.visit(FamilyPartIndex::pack(query,0,length),[&](int i){prefix.push_back(i);return prefix.size()<stop;},lookups,candidates);
                frozen.visit(FamilyPartIndex::pack(query,0,length),[&](int i){old_prefix.push_back(i);return old_prefix.size()<stop;},frozen_lookups,frozen_candidates);
                if(prefix!=old_prefix||lookups!=expected_lookups(query,index.cells())||candidates!=frozen_candidates)throw std::runtime_error("frozen early-exit mismatch");
            }
            std::sort(actual.begin(), actual.end());
            for (int i = 0; i < int(keys.size()); ++i) if (expected_values.count(keys[i])) expected.push_back(i);
            if (actual != expected || lookups > 6) throw std::runtime_error("candidate set mismatch");
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
