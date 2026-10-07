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
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug=false,bool timings=false);
FinalOutput emit_library(const LibraryPipeline& library,double cutoff);
}
int merge_r8qkeep_main(int argc,char** argv);
