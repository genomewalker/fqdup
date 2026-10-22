#pragma once
#include "cutoff_grid.hpp"
#include "risk_diagnostics.hpp"
namespace fqdup::r8qkeep {
inline double calibrated_mate_risk(const OutputRisk& risk,const HoldoutBins& holdout){
    const auto bins=holdout.qualities();std::uint64_t total=0;double sum=0;
    for(int q=0;q<94;++q){const auto n=risk.qualities[1][q];total+=n;sum+=n*(bins[q].n>=calibration_min_bin?double(bins[q].disagree)/bins[q].n:std::pow(10.,-q/10.));}
    return total?sum/total:0;
}
inline std::size_t select_automatic_cutoff(const std::vector<OutputRisk>& risks,const HoldoutBins& holdout,double epsilon,double floor,bool calibration_ok,bool vocabulary_frozen){
    if(risks.size()!=cutoff_grid.size()||!std::isfinite(epsilon)||epsilon<0||epsilon>1)throw std::runtime_error("automatic cutoff inputs");
    const auto first=std::lower_bound(cutoff_grid.begin(),cutoff_grid.end(),floor);
    if(first==cutoff_grid.end())throw std::runtime_error("automatic cutoff floor");
    const auto fallback=std::size_t(first-cutoff_grid.begin());
    if(!calibration_ok||!vocabulary_frozen)return fallback;
    const auto reference=std::size_t(std::find(cutoff_grid.begin(),cutoff_grid.end(),.316)-cutoff_grid.begin());
    for(std::size_t i=risks.size();i-->fallback;)if(calibrated_mate_risk(risks[i],holdout)<=epsilon&&risks[i].adapter_bases<=risks[reference].adapter_bases)return i;
    return fallback;
}
}
