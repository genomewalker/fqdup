#include "../src/merge_r8qkeep/canonical_pairs.hpp"
#include <array>
#include <iostream>
#include <random>

using namespace fqdup::r8qkeep;
using Row = std::array<std::string, 8>;
static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
static std::vector<Row> rows(const std::vector<FamilyRead>& r1,
                             const std::vector<FamilyRead>& r2) {
    std::vector<Row> out;
    for (std::size_t i = 0; i < r1.size(); ++i)
        out.push_back({r1[i].id, r1[i].h, r1[i].s, r1[i].q,
                       r2[i].h, r2[i].s, r2[i].q, r2[i].id});
    return out;
}
static void check(const std::vector<Row>& input,WorkerPool* pool=nullptr,int workers=1) {
    std::vector<FamilyRead> r1, r2;
    for (const auto& r : input) {
        r1.push_back({r[1], r[0], r[2], r[3]});
        r2.push_back({r[4], r[7], r[5], r[6]});
    }
    auto expected = input;
    std::sort(expected.begin(), expected.end());
    canonicalize_family_input(r1, r2,pool,workers);
    require(rows(r1, r2) == expected, "canonical output differs from row oracle");
    canonicalize_family_input(r1, r2,pool,workers);
    require(rows(r1, r2) == expected, "canonicalization is not idempotent");
}
int main() {
    check({});
    const Row base{"same", "@same/1", "AC", "II", "@same/2", "GT", "II", "same"};
    check({base});
    // Exhaust all permutations of records differing in each mate field.
    std::vector<Row> variants{base};
    for (std::size_t field = 1; field <= 6; ++field) {
        auto changed = base;
        changed[field] += "x";
        variants.push_back(std::move(changed));
    }
    std::sort(variants.begin(), variants.end());
    do { check(variants); } while (std::next_permutation(variants.begin(), variants.end()));
    variants.push_back(base); // Exact duplicates need no additional tie-break.
    auto other = base; other[0] = other[7] = "a";
    other[1] = "@a/1"; other[4] = "@a/2";
    variants.push_back(other);
    std::mt19937 rng(25000);
    for (int i = 0; i < 300; ++i) { std::shuffle(variants.begin(), variants.end(), rng); check(variants); }

    // Uneven merge ranges, tied IDs with different mate fields, exact duplicate
    // rows and already-sorted input exercise the parallel path and idempotence.
    std::vector<Row> large;
    for(int i=0;i<8193;++i){
        auto row=variants[rng()%variants.size()];
        row[0]=row[7]="id"+std::to_string(rng()%1024);large.push_back(std::move(row));
    }
    for(int workers:{2,3,8}){WorkerPool pool(workers);check(large,&pool,workers);}

    std::vector<FamilyRead> r1{{"@a/1", "a", "AC", "II"}, {"@b/1", "b", "GT", "II"}};
    std::vector<FamilyRead> r2{{"@b/2", "b", "TG", "II"}, {"@a/2", "a", "CA", "II"}};
    const auto original = rows(r1, r2);
    bool rejected = false;
    try { canonicalize_family_input(r1, r2); } catch (const std::runtime_error&) { rejected = true; }
    require(rejected && rows(r1, r2) == original, "mispaired inputs were sorted or changed");
    r2.pop_back(); rejected = false;
    try { canonicalize_family_input(r1, r2); } catch (const std::runtime_error&) { rejected = true; }
    require(rejected && r1.size() == 2 && r2.size() == 1 && r1[0].id == "a" && r2[0].id == "b",
            "unequal mate counts were accepted or changed");
    std::cout << "TASK25_CANONICAL_PAIRS_PASS\n";
}
