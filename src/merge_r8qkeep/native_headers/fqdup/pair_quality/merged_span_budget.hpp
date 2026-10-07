#pragma once
#include "native_quality.hpp"
namespace fqquality {
// Longest contiguous span with mean estimated error <= cut. Fixed 1e-12 units
// make ties portable; ties choose lower total estimated error, then earliest.
inline fqalign::TrimInterval merged_span_budget(const std::vector<double>& risk,double cut) {
 if(!std::isfinite(cut)||cut<0||cut>1||risk.size()>512)throw std::invalid_argument("span budget input");
 constexpr int64_t scale=1000000000000LL;const int64_t budget=int64_t(std::floor(cut*scale+.5));
 std::vector<int64_t> sums(risk.size()+1,0);
 for(size_t i=0;i<risk.size();++i){if(!std::isfinite(risk[i])||risk[i]<0||risk[i]>1)throw std::invalid_argument("span budget risk");sums[i+1]=sums[i]+int64_t(std::floor(risk[i]*scale+.5));}
 fqalign::TrimInterval best{};int64_t best_error=0;
 for(size_t a=0;a<risk.size();++a)for(size_t b=a+1;b<=risk.size();++b){const int64_t e=sums[b]-sums[a];const size_t n=b-a,old=best.end-best.begin;
 if(e>budget*int64_t(n))continue;
 if(n>old||(n==old&&e<best_error)){best={a,b};best_error=e;}}
 return best;
}
}

