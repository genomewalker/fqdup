#pragma once
#include <memory>
#include <atomic>
#include <cstdint>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
struct RepairShadow {
    std::atomic<uint64_t> skip_and_unchanged{0},skip_but_changed{0},noskip_unchanged{0},noskip_changed{0};
};
struct RoutingResult { std::string routed_tsv,audit_tsv;std::vector<bool> routed;std::vector<std::string> replacements;
    // One bound per native row, from this pipeline's router; -1 means full repair.
    std::vector<int64_t> repair_bounds; };
class RepairAlignmentCache {
    struct Impl;
    std::shared_ptr<Impl> impl_;
    friend class RepairPipeline;
public:
    RepairAlignmentCache(std::size_t start,std::size_t count);
};
class RepairPipeline {
    struct Impl;
    std::shared_ptr<Impl> impl_;
public:
    RepairPipeline(const std::string& nominal_prior,const std::string& full_prior);
    RoutingResult route(const std::string& native_tsv,bool debug=true,int threads=1,std::size_t start=0,RepairAlignmentCache* cache=nullptr) const;
    std::string repair(const std::string& input_tsv,int mode,double cutoff,std::size_t start=0,int threads=1,RepairShadow* shadow=nullptr,const std::vector<bool>* excluded=nullptr,const std::vector<int64_t>* repair_bounds=nullptr,RepairAlignmentCache* cache=nullptr) const;
};
}
