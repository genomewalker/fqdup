#pragma once
#include "stage_records.hpp"
#include "shared_alignment.hpp"
#include <memory>
#include <atomic>
#include <cstdint>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
struct RepairShadow {
    std::atomic<uint64_t> skip_and_unchanged{0},skip_but_changed{0},noskip_unchanged{0},noskip_changed{0};
    std::atomic<uint64_t> clean_and_unchanged{0},clean_but_changed{0},notclean{0};
};
struct RoutingResult { std::string routed_tsv,audit_tsv;std::vector<bool> routed;std::vector<std::string> replacements;
    // One bound per native row, from this pipeline's router; -1 means full repair.
    std::vector<int64_t> repair_bounds,clean_bounds;
    std::vector<int64_t> native_clean_bounds;
    std::size_t native_clean=0;
    std::vector<NativeRecord> record_replacements; };
class RepairAlignmentCache {
    struct Impl;
    std::shared_ptr<Impl> impl_;
    friend class RepairPipeline;
public:
    RepairAlignmentCache(std::size_t start,std::size_t count,const std::vector<std::shared_ptr<fqdup_shared_alignment::Pair>>* shared=nullptr,std::size_t offset=0,const std::vector<std::size_t>* selected=nullptr);
};
class RepairPipeline {
    struct Impl;
    std::shared_ptr<Impl> impl_;
public:
    RepairPipeline(const std::string& nominal_prior,const std::string& full_prior);
    std::vector<int64_t> native_clean_bounds(NativeRecords native,double cutoff) const;
    RoutingResult route_records(NativeRecords native,bool debug=false,std::size_t start=0,RepairAlignmentCache* cache=nullptr,bool shadow=false,double clean_cutoff=0,const std::vector<int64_t>* certified=nullptr) const;
    std::vector<RepairedRecord> repair_records(NativeRecords input,int mode,double cutoff,std::size_t start=0,RepairShadow* shadow=nullptr,const std::vector<bool>* excluded=nullptr,const std::vector<int64_t>* repair_bounds=nullptr,RepairAlignmentCache* cache=nullptr,const std::vector<int64_t>* clean_bounds=nullptr,const std::vector<bool>* selected=nullptr,const std::vector<NativeRecord>* replacements=nullptr,std::size_t mask_offset=0) const;
    // Explicit TSV compatibility adapters; production stages use records.
    RoutingResult route(const std::string& native_tsv,bool debug=true,int threads=1,std::size_t start=0,RepairAlignmentCache* cache=nullptr,bool shadow=false) const;
    std::string repair(const std::string& input_tsv,int mode,double cutoff,std::size_t start=0,int threads=1,RepairShadow* shadow=nullptr,const std::vector<bool>* excluded=nullptr,const std::vector<int64_t>* repair_bounds=nullptr,RepairAlignmentCache* cache=nullptr,const std::vector<int64_t>* clean_bounds=nullptr,const std::vector<bool>* selected=nullptr) const;
};
}
