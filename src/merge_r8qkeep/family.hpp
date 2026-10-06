#pragma once
#include <array>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
struct FamilyRead { std::string h, id, s, q; };
struct FamilyPair { std::string id; std::array<FamilyRead, 2> r; };
struct FamilyResult {
    std::vector<FamilyPair> pairs;
    std::string edits_tsv, families_tsv, stats_json;
};
std::vector<FamilyRead> read_family_fastq(const std::string& path);
FamilyResult correct_family(std::vector<FamilyRead> r1, std::vector<FamilyRead> r2);
}
