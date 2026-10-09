#pragma once
#include <algorithm>
#include <array>
#include <cstdlib>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include "alignment_cache.hpp"
namespace fqalign {
// Same alignment space/score/band/delta/cap as gapped_corr.py.
// Enumeration stops at cap+1: rejected counts are lower bounds, never exact claims.
struct Correspondence {
 int best=-1000000000, paths=0; bool eligible=false;
 std::string reason="no_alignment";
 std::vector<std::pair<int,int>> invariant, certified;
};
class BandedCorrespondence {
 static constexpr int neg=-1000000000;
 struct Cell {std::array<int,3> suffix{neg,neg,neg}; int through=0;};
 std::vector<Cell> cells;
 int n=0,m=0,d0=0;
 bool inside(int i,int j) const {return i>=1&&i<=n&&j>=1&&j<=m&&std::abs(j-i-d0)<=3;}
 Cell& at(int i,int j) {return cells[i*7+(j-i-d0+3)];}
 static int sub(char a,char b) {return a==b&&a!='N'?1:-1;}
 public:
 Correspondence align(const std::string&a,const std::string&b,int insert) {
  return alignment_cache::get<BandedCorrespondence>(a,b,insert,[&]{return align_uncached(a,b,insert);});
 }
 Correspondence align_uncached(const std::string&a,const std::string&b,int insert) {
  Correspondence out;
  if(a.size()>512||b.size()>512){out.reason="unsupported_read_length";return out;}
  if(a.empty()||b.empty()){out.reason="empty_read";return out;}
  for(char c:a+b)if(std::string("ACGTN").find(c)==std::string::npos){out.reason="unsupported_base";return out;}
  n=a.size();m=b.size();
  if(insert<0){out.reason="no_geometry";return out;}
  d0=m-insert;
  bool any=false;
  for(int i=1;i<=n;++i)if(std::max(1,i+d0-3)<=std::min(m,i+d0+3)){any=true;break;}
  if(!any){out.reason="empty_band";return out;}
  cells.assign((n+1)*7,Cell{});
  // Suffix after the current state. Terminal states must be aligned columns.
  for(int i=n;i>=1;--i) for(int j=std::min(m,i+d0+3);j>=std::max(1,i+d0-3);--j) {
    auto& c=at(i,j); int diag=neg;
    if(inside(i+1,j+1))diag=at(i+1,j+1).suffix[0]+sub(a[i],b[j]);
    c.suffix[0]=std::max((i==n||j==m)?0:neg,diag);
    if(inside(i+1,j))c.suffix[0]=std::max(c.suffix[0],at(i+1,j).suffix[1]-4);
    if(inside(i,j+1))c.suffix[0]=std::max(c.suffix[0],at(i,j+1).suffix[2]-4);
    c.suffix[1]=diag;c.suffix[2]=diag;
    if(inside(i+1,j))c.suffix[1]=std::max(c.suffix[1],at(i+1,j).suffix[1]-1);
    if(inside(i,j+1))c.suffix[2]=std::max(c.suffix[2],at(i,j+1).suffix[2]-1);
  }
  for(int i=1;i<=n;++i)for(int j=std::max(1,i+d0-3);j<=std::min(m,i+d0+3);++j)
   if(i==1||j==1)out.best=std::max(out.best,sub(a[i-1],b[j-1])+at(i,j).suffix[0]);
  if(out.best<=neg/2)return out;
  int threshold=out.best-2;
  int max_counted_score=neg;
  std::vector<std::pair<int,int>> path;
  auto visit=[&](auto&& recurse,int i,int j,int state,int score)->void{
   if(out.paths>256||score+at(i,j).suffix[state]<threshold)return;
   if(state==0)path.emplace_back(i,j);
   if(state==0&&(i==n||j==m)&&score>=threshold){
    ++out.paths;max_counted_score=std::max(max_counted_score,score);
    if(out.paths<=256)for(auto c:path)++at(c.first,c.second).through;
   }
   if(out.paths<=256){
    if(inside(i+1,j+1))recurse(recurse,i+1,j+1,0,score+sub(a[i],b[j]));
    if(state!=2&&inside(i+1,j))recurse(recurse,i+1,j,1,score-(state==1?1:4));
    if(state!=1&&inside(i,j+1))recurse(recurse,i,j+1,2,score-(state==2?1:4));
   }
   if(state==0)path.pop_back();
  };
  for(int i=1;i<=n&&out.paths<=256;++i)for(int j=std::max(1,i+d0-3);j<=std::min(m,i+d0+3)&&out.paths<=256;++j)
   if(i==1||j==1)visit(visit,i,j,0,sub(a[i-1],b[j-1]));
  if(out.paths>256){out.reason="cap_exceeded";return out;}
  if(out.paths<1||max_counted_score!=out.best){out.reason="internal_inconsistency";return out;}
  for(int i=1;i<=n;++i)for(int j=std::max(1,i+d0-3);j<=std::min(m,i+d0+3);++j)
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
  out.eligible=true;out.reason="ok";return out;
 }
};
}

