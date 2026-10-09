#pragma once
#include "native_correspondence.hpp"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>
namespace fqalign {

// Experimental unbanded boundary overlap. Score +1/-1, gap 3+g, delta2, cap256.
// Gates frozen before labels: all paths n>=30, invariant ACGT>=30, worst error fraction<=.30, certified anchor.
// Enumeration stops at cap+1: rejected counts are lower bounds, never exact claims.
struct FullResult {Correspondence corr;int min_informative=1000000000,max_gap=0;double max_fraction=0;};
class FullCorrespondence {
 static constexpr int neg=-1000000000;
 using Cell=fqalign::AlignmentCell;
 Cell* cells=nullptr;
 int n=0,m=0;
 bool inside(int i,int j) const {return i>=1&&i<=n&&j>=1&&j<=m;}
 Cell& at(int i,int j) {return cells[i*(m+1)+j];}
 static int sub(char a,char b) {return a==b&&a!='N'?1:-1;}
 template<bool next_i,bool next_j>
 void fill_suffix(int i,int j,const std::string& a,const std::string& b) {
    auto& c=at(i,j);c.through=0;int diag=neg;
    if constexpr(next_i&&next_j)diag=at(i+1,j+1).suffix[0]+sub(a[i],b[j]);
    c.suffix[0]=std::max((!next_i||!next_j)?0:neg,diag);
    if constexpr(next_i)c.suffix[0]=std::max(c.suffix[0],at(i+1,j).suffix[1]-4);
    if constexpr(next_j)c.suffix[0]=std::max(c.suffix[0],at(i,j+1).suffix[2]-4);
    c.suffix[1]=diag;c.suffix[2]=diag;
    if constexpr(next_i)c.suffix[1]=std::max(c.suffix[1],at(i+1,j).suffix[1]-1);
    if constexpr(next_j)c.suffix[2]=std::max(c.suffix[2],at(i,j+1).suffix[2]-1);
 }
 public:
 FullResult align(const std::string&a,const std::string&b) {
  FullResult result;auto& out=result.corr;
  if(a.size()>512||b.size()>512){out.reason="unsupported_read_length";return result;}
  if(a.empty()||b.empty()){out.reason="empty_read";return result;}
  if(!::fqdup_fast_bases::valid(a)||!::fqdup_fast_bases::valid(b)){out.reason="unsupported_base";return result;}
  n=a.size();m=b.size();
  int diagonal=0;
  if(::fqdup_fast_diagonal::full_offset(a,b,diagonal)){
   ::fqdup_fast_diagonal::fill(out,n,m,diagonal);
   result.min_informative=out.best;return result;
  }
  cells=fqalign::alignment_scratch((n+1)*(m+1));
  // Suffix after the current state. Terminal states must be aligned columns.
  // Preserve descending row-major updates; specialize only boundary predicates.
  fill_suffix<false,false>(n,m,a,b);
  for(int j=m-1;j>=1;--j)fill_suffix<false,true>(n,j,a,b);
  for(int i=n-1;i>=1;--i){
    fill_suffix<true,false>(i,m,a,b);
    for(int j=m-1;j>=1;--j)fill_suffix<true,true>(i,j,a,b);
  }
  // Original row-major boundary order: entire first row, then first column.
  for(int j=1;j<=m;++j)out.best=std::max(out.best,sub(a[0],b[j-1])+at(1,j).suffix[0]);
  for(int i=2;i<=n;++i)out.best=std::max(out.best,sub(a[i-1],b[0])+at(i,1).suffix[0]);
  if(out.best<=neg/2)return result;
  int threshold=out.best-2;
  int max_counted_score=neg;
  std::vector<std::pair<int,int>> path;int informative=0,mismatches=0,gaps=0;
  auto visit=[&](auto&& recurse,int i,int j,int state,int score)->void{
   if(out.paths>256||score+at(i,j).suffix[state]<threshold)return;
   bool info=state==0&&a[i-1]!='N'&&b[j-1]!='N';bool mm=info&&a[i-1]!=b[j-1];
   if(state==0)path.emplace_back(i,j);else ++gaps;informative+=info;mismatches+=mm;
   if(state==0&&(i==n||j==m)&&score>=threshold){
    result.min_informative=std::min(result.min_informative,informative);result.max_gap=std::max(result.max_gap,gaps);result.max_fraction=std::max(result.max_fraction,double(mismatches+gaps)/std::max(1,informative+gaps));
    ++out.paths;max_counted_score=std::max(max_counted_score,score);
    if(out.paths<=256)for(auto c:path)++at(c.first,c.second).through;
   }
   if(out.paths<=256){
    if(inside(i+1,j+1))recurse(recurse,i+1,j+1,0,score+sub(a[i],b[j]));
    if(state!=2&&inside(i+1,j))recurse(recurse,i+1,j,1,score-(state==1?1:4));
    if(state!=1&&inside(i,j+1))recurse(recurse,i,j+1,2,score-(state==2?1:4));
   }
   if(state==0)path.pop_back();else --gaps;informative-=info;mismatches-=mm;
  };
  for(int j=1;j<=m&&out.paths<=256;++j)visit(visit,1,j,0,sub(a[0],b[j-1]));
  for(int i=2;i<=n&&out.paths<=256;++i)visit(visit,i,1,0,sub(a[i-1],b[0]));
  if(out.paths>256){out.reason="cap_exceeded";return result;}
  if(out.paths<1||max_counted_score!=out.best){out.reason="internal_inconsistency";return result;}
  for(int i=1;i<=n;++i)for(int j=1;j<=m;++j)
   if(at(i,j).through==out.paths)out.invariant.emplace_back(i-1,j-1);
  // Maximal consecutive diagonal blocks; flanks count agreeing columns.
  for(size_t begin=0;begin<out.invariant.size();){
   size_t end=begin+1;
   while(end<out.invariant.size()&&out.invariant[end].first==out.invariant[end-1].first+1&&out.invariant[end].second==out.invariant[end-1].second+1)++end;
   int total=0,before=0;
   for(size_t k=begin;k<end;++k){auto c=out.invariant[k];total+=a[c.first]==b[c.second]&&a[c.first]!='N';}
   for(size_t k=begin;k<end;++k){auto c=out.invariant[k];bool agree=a[c.first]==b[c.second]&&a[c.first]!='N';
    if(agree&&before>=5&&total-before-1>=5)out.certified.push_back(c);
    before+=agree;
   }begin=end;
  }
  int invariant_n=0;for(auto [i,j]:out.invariant)invariant_n+=a[i]!='N'&&b[j]!='N';
  if(result.min_informative<30||invariant_n<30){out.reason="insufficient_columns";return result;}
  if(result.max_fraction>.30){out.reason="excess_error_fraction";return result;}
  if(out.certified.empty()){out.reason="no_certified_anchor";return result;}
  out.eligible=true;out.reason="ok";return result;
 }
};
}
