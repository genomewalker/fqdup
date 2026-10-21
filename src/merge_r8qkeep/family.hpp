#pragma once
#include <array>
#include <memory>
#include <string>
#include <vector>
#include "family_key.hpp"
namespace fqdup::r8qkeep {
struct FamilyRead { std::string h, id, s, q; };
struct FamilyPair { std::string id; std::array<FamilyRead, 2> r; };
class WorkerPool;
struct FamilyGroup { std::vector<int> m; bool oversize=false; int diameter=0; };
struct FamilyPlan {
    std::size_t pairs=0;
    bool supported=true;
    std::vector<FamilyGroup> groups;
    std::shared_ptr<WorkerPool> pool;
};
struct FamilyResult {
    std::vector<FamilyPair> pairs;
    std::string edits_tsv, families_tsv, stats_json;
};
std::vector<FamilyRead> read_family_fastq(const std::string& path);
std::vector<FamilyKey> extract_family_keys(const std::vector<FamilyPair>& pairs);
FamilyPlan discover_family_keys(std::vector<FamilyKey> keys,int threads,WorkerPool* pool=nullptr);
void reorder_family_plan(FamilyPlan& plan,const std::vector<std::size_t>& order);
FamilyResult correct_family(std::vector<FamilyRead> r1, std::vector<FamilyRead> r2, int threads = 1,bool debug = true,WorkerPool* shared_pool=nullptr,FamilyPlan* plan=nullptr);
}
