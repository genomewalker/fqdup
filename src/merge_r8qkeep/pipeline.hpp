#pragma once
#include "census.hpp"
#include "mixture_fit.hpp"
#include "repair.hpp"
namespace fqdup::r8qkeep {
struct LibraryPipeline {
    FamilyResult family;
    NativeResult native;
    CensusResult observations;
    std::string nominal_prior,full_prior;
    std::unique_ptr<RepairPipeline> repair;
    RoutingResult routing;
};
struct FinalOutput { std::string merged,r1,r2,orphan,outcomes; };
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug=false,bool timings=false,bool defer_routing=false);
void write_library_cutoffs(LibraryPipeline& library,const std::vector<double>& cutoffs,const std::vector<std::string>& directories,const std::vector<std::string>& labels,int threads,bool debug,bool timings,bool shadow);
FinalOutput emit_library(const LibraryPipeline& library,double cutoff);
void write_library_outputs(const LibraryPipeline& library,double cutoff,const std::string& directory,int threads=1,RepairShadow* shadow=nullptr);
}
int merge_r8qkeep_main(int argc,char** argv);
