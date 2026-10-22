#pragma once
#include "census.hpp"
#include "compact_family.hpp"
#include "mixture_fit.hpp"
#include "repair.hpp"
#include "worker_pool.hpp"
#include <array>
namespace fqdup::r8qkeep {
enum class OutcomesMode { full, summary, none };
struct LibraryPipeline {
    double merged_skip_cutoff=0.05;
    OutcomesMode outcomes_mode=OutcomesMode::full;
    FamilyResult family;
    std::shared_ptr<CompactPairStore> compact;
    std::size_t pair_count()const{return compact?compact->size():family.pairs.size();}
    NativeSettings native_settings;
    NativeResult native;
    CensusResult observations;
    std::string nominal_prior,full_prior;
    std::unique_ptr<RepairPipeline> repair;
    RoutingResult routing;
    std::shared_ptr<WorkerPool> pool;
};
struct FinalOutput { std::string merged,r1,r2,orphan,outcomes; std::array<std::size_t,5> outcome_counts{}; };
struct CutoffFloor { double floor; std::string diagnostics; };
CutoffFloor derive_cutoff_floor(const LibraryPipeline& library,int threads);
LibraryPipeline prepare_compact_library(std::shared_ptr<CompactPairStore> store,const NativeSettings& settings,int threads,bool debug,bool timings,bool counts,bool no_family=false);
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug=false,bool timings=false,bool defer_routing=false,bool counts=false,bool no_family=false);
void write_library_cutoffs(LibraryPipeline& library,const std::vector<double>& cutoffs,const std::vector<std::string>& directories,const std::vector<std::string>& labels,int threads,bool debug,bool timings,bool shadow,bool consume_records=false);
FinalOutput emit_library(const LibraryPipeline& library,double cutoff);
void write_library_outputs(const LibraryPipeline& library,double cutoff,const std::string& directory,int threads=1,RepairShadow* shadow=nullptr);
}
int merge_r8qkeep_main(int argc,char** argv);
