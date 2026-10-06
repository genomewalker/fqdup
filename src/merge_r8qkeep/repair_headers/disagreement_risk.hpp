#pragma once
#include "native_quality.hpp"
namespace fqquality {
struct RiskUpdate {std::vector<std::pair<int,int>> raised;};
// Research penalty under the nominal independent-substitution model.
// Not a calibrated error probability and not enabled in the production command.
inline RiskUpdate raise_anchored_disagreements(fqalign::PairRisk& risk,
                                             const std::string& raw1,const std::string& raw2) {
 RiskUpdate update;
 if(!risk.correspondence.eligible)return update;
 if(risk.r1.size()!=raw1.size()||risk.r2.size()!=raw2.size())throw std::invalid_argument("risk length mismatch");
 const auto rc2=fqalign::reverse_complement(raw2);
 const auto& columns=risk.correspondence.invariant;
 size_t begin=0;
 while(begin<columns.size()){
  size_t end=begin+1;
  while(end<columns.size()&&columns[end].first==columns[end-1].first+1
        &&columns[end].second==columns[end-1].second+1)++end;
  size_t total=0,before=0;
  for(size_t t=begin;t<end;++t){auto [i,j]=columns[t];total+=raw1.at(i)==rc2.at(j)&&raw1.at(i)!='N';}
  for(size_t t=begin;t<end;++t){auto [i,j]=columns[t];
   const bool agree=raw1.at(i)==rc2.at(j)&&raw1.at(i)!='N';
   if(agree){++before;continue;}
   if(raw1.at(i)=='N'||rc2.at(j)=='N'||before<5||total-before<5)continue;
   const int k=int(raw2.size())-1-j;
   const double e1=risk.r1.at(i),e2=risk.r2.at(k),d=e1+e2-4.0*e1*e2/3.0;
   if(!(d>0))continue;
   const double p1=std::clamp(e1*(1.0-e2/3.0)/d,0.0,1.0);
   const double p2=std::clamp(e2*(1.0-e1/3.0)/d,0.0,1.0);
   if(p1>e1||p2>e2)update.raised.emplace_back(i,k);
   risk.r1[i]=std::max(e1,p1);risk.r2[k]=std::max(e2,p2);
  }
  begin=end;
 }
 return update;
}
}

