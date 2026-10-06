#pragma once
#include "native_quality.hpp"
namespace fqquality {
// Experimental selection using the declared Mott cutoff, without new tuned
// thresholds. Nominal conditional risks remain assumption-dependent.
inline void raise_weaker_disagreements(fqalign::PairRisk& risk,
                                      const std::string& raw1,const std::string& raw2,double cutoff) {
 if(!risk.correspondence.eligible)return;
 if(risk.r1.size()!=raw1.size()||risk.r2.size()!=raw2.size())throw std::invalid_argument("risk length mismatch");
 const auto rc2=fqalign::reverse_complement(raw2);const auto& columns=risk.correspondence.invariant;
 for(size_t begin=0;begin<columns.size();){
  size_t end=begin+1;while(end<columns.size()&&columns[end].first==columns[end-1].first+1
        &&columns[end].second==columns[end-1].second+1)++end;
  size_t total=0,before=0;
  for(size_t t=begin;t<end;++t){auto [i,j]=columns[t];total+=raw1.at(i)==rc2.at(j)&&raw1.at(i)!='N';}
  for(size_t t=begin;t<end;++t){auto [i,j]=columns[t];const bool agree=raw1.at(i)==rc2.at(j)&&raw1.at(i)!='N';
   if(agree){++before;continue;}if(raw1.at(i)=='N'||rc2.at(j)=='N'||before<5||total-before<5)continue;
   const int k=int(raw2.size())-1-j;const double e1=risk.r1.at(i),e2=risk.r2.at(k),d=e1+e2-4.0*e1*e2/3.0;
   if(!(d>0))continue;const double p1=std::clamp(e1*(1.0-e2/3.0)/d,0.0,1.0),p2=std::clamp(e2*(1.0-e1/3.0)/d,0.0,1.0);
   if(p1<=cutoff&&p2>cutoff)risk.r2[k]=std::max(e2,p2);
   else if(p2<=cutoff&&p1>cutoff)risk.r1[i]=std::max(e1,p1);
  }
  begin=end;
 }
}
inline void protect_unique_coverage(fqalign::PairRisk& updated,
                                   const std::vector<double>& original1,const std::vector<double>& original2,
                                   size_t prefix1,size_t prefix2,double cutoff) {
 const std::array<fqalign::TrimInterval,2> old{
  fqalign::mott_segment(original1,prefix1,cutoff),fqalign::mott_segment(original2,prefix2,cutoff)};
 const std::array<fqalign::TrimInterval,2> proposed{
  fqalign::mott_segment(updated.r1,prefix1,cutoff),fqalign::mott_segment(updated.r2,prefix2,cutoff)};
 std::array<std::vector<int>,2> counterpart{std::vector<int>(original1.size(),-1),std::vector<int>(original2.size(),-1)};
 for(auto [i,j]:updated.correspondence.invariant){int k=int(original2.size())-1-j;counterpart[0].at(i)=k;counterpart[1].at(k)=i;}
 std::array<bool,2> keep{true,true};
 for(int mate=0;mate<2;++mate){const auto a=old[mate],b=proposed[mate];
  if(b.end>b.begin&&(b.begin<a.begin||b.end>a.end)){keep[mate]=false;continue;}
  for(size_t i=a.begin;i<a.end&&keep[mate];++i){if(i>=b.begin&&i<b.end)continue;
   int k=counterpart[mate].at(i);const auto o=old[1-mate],p=proposed[1-mate];
   if(k<0||size_t(k)<o.begin||size_t(k)>=o.end||size_t(k)<p.begin||size_t(k)>=p.end)keep[mate]=false;
  }
 }
 if(!keep[0])updated.r1=original1;if(!keep[1])updated.r2=original2;
 // Coverage protection does not identify the correct base at a disagreement.
}
}

