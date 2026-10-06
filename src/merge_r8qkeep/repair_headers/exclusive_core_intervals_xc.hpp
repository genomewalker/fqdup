#pragma once
#include "joint_pair_intervals.hpp"
#include <array>
#include <cstdint>
namespace fqquality {
inline bool exclusive_solve(const std::vector<double>& a,const std::vector<double>& b,
 size_t n,size_t m,const std::vector<JointColumn>& cols,double cut,
 std::array<fqalign::TrimInterval,2>& result,int excluded,const std::array<fqalign::TrimInterval,2>& core,int64_t& score,int& length) {
 if(cols.empty()||n>a.size()||m>b.size())return false;
 struct Event{int i,j;bool agree=false;};std::vector<Event> events;int i=0,j=0;
 for(const auto& c:cols){int k=int(m)-1-c.k;
  if(c.i<i||k<j||c.i>=int(n)||k>=int(m))return false;
  while(i<c.i)events.push_back({i++,-1});while(j<k)events.push_back({-1,j++});events.push_back({i++,j++,c.agree});
 }
 while(i<int(n))events.push_back({i++,-1});while(j<int(m))events.push_back({-1,j++});
 struct Cell{int64_t score=0;int len=0;bool valid=false;};std::array<Cell,9> dp{};dp[0].valid=true;
 std::vector<std::array<int,9>> history;history.reserve(events.size());
 for(const auto& ev:events){std::array<Cell,9> next{};std::array<int,9> back;back.fill(-1);
  for(int p=0;p<9;++p){if(!dp[p].valid)continue;int u=p/3,v=p%3;
   for(int x=u;x<3;++x)for(int y=v;y<3;++y){if((ev.i<0&&x!=u)||(ev.j<0&&y!=v))continue;
    bool s1=ev.i>=0&&x==1,s2=ev.j>=0&&y==1; if((s1&&excluded==1&&size_t(ev.i)>=core[0].begin&&size_t(ev.i)<core[0].end)||(s2&&excluded==2&&m-1-ev.j>=core[1].begin&&m-1-ev.j<core[1].end))continue;double reward=0;
    if(s1&&s2)reward=cut-(ev.agree?std::max(a.at(ev.i),b.at(m-1-ev.j)):1);
    else if(s1)reward=cut-a.at(ev.i);else if(s2)reward=cut-b.at(m-1-ev.j);
    const int z=x*3+y;Cell q{dp[p].score+int64_t(std::floor(reward*1e12+.5)),dp[p].len+int(s1)+int(s2),true};
    if(!next[z].valid||q.score>next[z].score||(q.score==next[z].score&&q.len>next[z].len)){next[z]=q;back[z]=p;}
   }
  }dp=next;history.push_back(back);
 }
 int state=0;for(int z=1;z<9;++z)if(dp[z].valid&&(!dp[state].valid||dp[z].score>dp[state].score||(dp[z].score==dp[state].score&&dp[z].len>dp[state].len)))state=z;
 score=dp[state].score;length=dp[state].len;
 std::array<std::vector<int>,2> selected;
 for(int t=int(events.size())-1;t>=0;--t){auto ev=events[t];if(ev.i>=0&&state/3==1)selected[0].push_back(ev.i);if(ev.j>=0&&state%3==1)selected[1].push_back(int(m)-1-ev.j);state=history[t][state];if(state<0)throw std::runtime_error("joint traceback");}
 for(int q=0;q<2;++q){auto& xs=selected[q];std::sort(xs.begin(),xs.end());result[q]={};if(xs.empty())continue;if(xs.back()-xs.front()+1!=int(xs.size()))throw std::runtime_error("joint noncontiguous interval");result[q]={size_t(xs.front()),size_t(xs.back()+1)};}
 return true;
}
inline bool exclusive_core_intervals(const std::vector<double>& a,const std::vector<double>& b,const std::vector<double>& raw1,const std::vector<double>& raw2,const std::vector<JointColumn>& cols,double cut,int fold,std::array<fqalign::TrimInterval,2>& result) {
 if(cols.empty()||a.empty()||b.empty()||raw1.size()!=a.size()||raw2.size()!=b.size())return false;
 std::array<fqalign::TrimInterval,2> core{{{a.size(),0},{b.size(),0}}};
 for(auto c:cols){if(c.i<0||c.k<0||size_t(c.i)>=a.size()||size_t(c.k)>=b.size())return false;core[0].begin=std::min(core[0].begin,size_t(c.i));core[0].end=std::max(core[0].end,size_t(c.i+1));core[1].begin=std::min(core[1].begin,size_t(c.k));core[1].end=std::max(core[1].end,size_t(c.k+1));}
 struct Candidate{std::array<fqalign::TrimInterval,2> iv;int supplier;int64_t mean,score;int length;};std::vector<Candidate> candidates;
 for(int supplier=1;supplier<=2;++supplier){Candidate x{};x.supplier=supplier;if(!exclusive_solve(a,b,a.size(),b.size(),cols,cut,x.iv,3-supplier,core,x.score,x.length))return false;const auto iv=x.iv[supplier-1],span=core[supplier-1];if(std::max(iv.begin,span.begin)>=std::min(iv.end,span.end))continue;const auto& raw=supplier==1?raw1:raw2;double sum=0;for(size_t i=span.begin;i<span.end;++i)sum+=std::max(1e-10,std::min(1.0,raw[i]));x.mean=int64_t(std::floor(sum/(span.end-span.begin)*1e12+.5));candidates.push_back(x);}
 if(candidates.empty())return false;
 std::sort(candidates.begin(),candidates.end(),[fold](const Candidate& x,const Candidate& y){if(x.mean!=y.mean)return x.mean<y.mean;if(x.score!=y.score)return x.score>y.score;if(x.length!=y.length)return x.length>y.length;return (x.supplier==1+fold)>(y.supplier==1+fold);});result=candidates[0].iv;return true;
}
}

