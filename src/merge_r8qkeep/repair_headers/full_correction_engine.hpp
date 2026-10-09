#pragma once
#include "full_correspondence.hpp"
namespace fqrepair {
using fqalign::Correspondence;
using FullResult=fqalign::FullResult;
// The alignment score, path order, delta and cap are identical to nominal full
// correspondence. Share its immutable cached summary, then apply this policy.
class FullCorrespondence {
 fqalign::FullCorrespondence engine;
 static FullResult correction_policy(FullResult result,const std::string& a,const std::string& b){
  if(!result.summary_complete)return result;
  auto& out=result.corr;out.eligible=false;
  int invariant_n=0;for(auto [i,j]:out.invariant)invariant_n+=a[i]!='N'&&b[j]!='N';
  if(result.min_informative<20||invariant_n<20){out.reason="insufficient_columns";return result;}
  if(result.max_fraction>.10){out.reason="excess_error_fraction";return result;}
  if(out.certified.empty()){out.reason="no_certified_anchor";return result;}
  out.eligible=true;out.reason="ok";return result;
 }
public:
 FullResult align(const std::string& a,const std::string& b){
  return correction_policy(engine.align(a,b),a,b);
 }
 FullResult align_uncached(const std::string& a,const std::string& b){
  return correction_policy(engine.align_uncached(a,b),a,b);
 }
};
}
