#pragma once
#include "native_quality.hpp"
namespace fqquality {
inline double retained_median_cutoff(const std::vector<double>& risk,const std::string& quality,double cut,double strength=3){
 if(quality.empty()||risk.size()!=quality.size()||!std::isfinite(cut)||cut<0||cut>1)throw std::invalid_argument("retained median input");
 const auto iv=fqalign::mott_segment(risk,risk.size(),cut);if(iv.end<=iv.begin)return cut;
 std::vector<double> raw;raw.reserve(iv.end-iv.begin);for(size_t i=iv.begin;i<iv.end;i++)raw.push_back(fqalign::phred_error(quality[i]));
 std::sort(raw.begin(),raw.end());const size_t n=raw.size();return (raw[(n-1)/2]+raw[n/2])/2<=cut/4?cut/strength:cut;
}
}

