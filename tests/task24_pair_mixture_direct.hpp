// Frozen pre-Task-24 arithmetic oracle (6ace7cc).
#pragma once
#include "../src/merge_r8qkeep/repair_headers/native_quality.hpp"
#include <fstream>
#include <sstream>
namespace task24_direct {
struct PairMixture {
 static constexpr int grid=64,K=66;
 std::array<std::array<double,K>,2> weights{};
 std::array<bool,2> supported{};
 explicit PairMixture(std::istream& f) {
  std::string magic;int g;
  if(!(f>>magic>>g)||magic!="FQDUP_PAIR_MIXTURE_V1"||g!=grid)throw std::invalid_argument("invalid pair-mixture header");
  for(int row=0;row<2;++row){int fold;uint64_t pairs;
   if(!(f>>fold>>pairs)||fold!=row)throw std::invalid_argument("invalid pair-mixture fold");
   double sum=0;for(double& w:weights[fold]){if(!(f>>w)||!std::isfinite(w)||w<0||w>1)throw std::invalid_argument("invalid mixture weight");sum+=w;}
   supported[fold]=pairs>0;if((supported[fold]&&std::abs(sum-1)>1e-8)||(!supported[fold]&&sum!=0))throw std::invalid_argument("mixture weight mass");
  }
  std::string extra;if(f>>extra)throw std::invalid_argument("extra mixture data");
 }
 bool prior_floor(fqalign::PairRisk& risk,const std::string& s1,const std::string& q1,const std::string& s2,const std::string& q2)const {
  int training=1-fold(s1,s2);if(!supported[training])return false;
  double flat=0;for(int j=1;j<K;++j)flat+=weights[training][j]*component(j);
  for(size_t i=0;i<q1.size();++i)risk.r1.at(i)=std::max(risk.r1.at(i),weights[training][0]*fqalign::phred_error(q1[i])+flat);
  for(size_t i=0;i<q2.size();++i)risk.r2.at(i)=std::max(risk.r2.at(i),weights[training][0]*fqalign::phred_error(q2[i])+flat);
  return true;
 }
 static double component(int j){return .75*std::pow(double(j-1)/grid,2);}
 static double clipped_log(double p){return std::log(std::clamp(p,1e-15,1-1e-15));}
 static std::pair<double,double> conditional(double e1,double e2,bool agree){
  if(agree){double a=e1*e2/3,den=(1-e1)*(1-e2)+a;return den>0?std::pair<double,double>{a/den,a/den}:std::pair<double,double>{1,1};}
  double d=e1+e2-4*e1*e2/3;return d>0?std::pair<double,double>{e1*(1-e2/3)/d,e2*(1-e1/3)/d}:std::pair<double,double>{1,1};
 }
 struct Audit {int fold=-1;size_t n=0;double nominal=0,flat=0;bool supported=false;std::array<double,K> posterior{};};
 static int fold(const std::string& s1,const std::string& s2){uint64_t hash=14695981039346656037ULL;for(unsigned char c:s1){hash^=c;hash*=1099511628211ULL;}hash^='|';hash*=1099511628211ULL;for(unsigned char c:s2){hash^=c;hash*=1099511628211ULL;}return int(hash&1);}
 bool apply(fqalign::PairRisk& risk,const std::string& s1,const std::string& q1,
            const std::string& s2,const std::string& q2,int mode,Audit* audit=nullptr,bool full_trusted=false)const {
  if(audit)audit->fold=fold(s1,s2);
  if(!risk.correspondence.eligible)return false;

  int training=1-fold(s1,s2);if(!supported[training])return false;
  const auto rc=fqalign::reverse_complement(s2);std::vector<std::pair<int,int>> cols;for(auto c:risk.correspondence.invariant)if(s1.at(c.first)!='N'&&rc.at(c.second)!='N')cols.push_back(c);
  std::array<double,K> logpost{},post{};
  for(int j=0;j<K;++j)logpost[j]=weights[training][j]>0?std::log(weights[training][j]):-std::numeric_limits<double>::infinity();
  size_t n=0;for(auto [i,j]:cols){int k=int(s2.size())-1-j;if(s1.at(i)=='N'||rc.at(j)=='N')continue;++n;
   double e1=fqalign::phred_error(q1.at(i)),e2=fqalign::phred_error(q2.at(k));bool agree=s1.at(i)==rc.at(j);
   for(int t=0;t<K;++t){double e=component(t),d=t==0?e1+e2-4*e1*e2/3:2*e-4*e*e/3;logpost[t]+=clipped_log(agree?1-d:d);}
  }
  if(!n)return false;double maximum=*std::max_element(logpost.begin(),logpost.end()),mass=0;if(!std::isfinite(maximum))throw std::runtime_error("empty mixture posterior");
  for(int t=0;t<K;++t){post[t]=std::exp(logpost[t]-maximum);mass+=post[t];}for(double& p:post)p/=mass;
  double flat=0;for(int t=1;t<K;++t)flat+=post[t]*component(t);
  if(audit){audit->n=n;audit->nominal=post[0];audit->flat=flat;audit->supported=true;audit->posterior=post;}
  if((mode==2||mode==4)&&post[0]>=.5)return false;
  for(size_t i=0;i<q1.size();++i)risk.r1.at(i)=post[0]*fqalign::phred_error(q1[i])+flat;
  for(size_t i=0;i<q2.size();++i)risk.r2.at(i)=post[0]*fqalign::phred_error(q2[i])+flat;
  // All columns including N are kept when forming blocks, matching the native
  // verifier. N is neither an agreeing flank nor a conditional-risk update.
  for(size_t begin=0;begin<cols.size();){size_t end=begin+1;while(end<cols.size()&&cols[end].first==cols[end-1].first+1&&cols[end].second==cols[end-1].second+1)++end;
   size_t total=0,before=0;for(size_t t=begin;t<end;++t){auto [i,j]=cols[t];total+=s1.at(i)==rc.at(j)&&s1.at(i)!='N';}
   for(size_t t=begin;t<end;++t){auto [i,j]=cols[t];int k=int(s2.size())-1-j;bool agree=s1.at(i)==rc.at(j)&&s1.at(i)!='N';if(agree)++before;
    if(s1.at(i)=='N'||rc.at(j)=='N')continue;
    bool cert=std::find(risk.correspondence.certified.begin(),risk.correspondence.certified.end(),std::pair<int,int>{i,j})!=risk.correspondence.certified.end();
    if(!full_trusted&&!cert&&(agree||before<5||total-before<5))continue;
    double e1=fqalign::phred_error(q1.at(i)),e2=fqalign::phred_error(q2.at(k)),p1=0,p2=0;
    for(int z=0;z<K;++z){double e=component(z);auto p=conditional(z==0?e1:e,z==0?e2:e,agree);p1+=post[z]*p.first;p2+=post[z]*p.second;}
    risk.r1.at(i)=p1;risk.r2.at(k)=p2;
   }begin=end;
  }
  return true;
 }
};
inline void guard_pair_mixture(fqalign::PairRisk& risk,const std::vector<double>& base1,
 const std::vector<double>& base2,size_t prefix1,size_t prefix2,double cut) {
 std::array<std::vector<bool>,2> certified{std::vector<bool>(base1.size()),std::vector<bool>(base2.size())};
 std::array<std::vector<int>,2> counterpart{std::vector<int>(base1.size(),-1),std::vector<int>(base2.size(),-1)};
 for(auto [i,j]:risk.correspondence.invariant){int k=int(base2.size())-1-j;counterpart[0].at(i)=k;counterpart[1].at(k)=i;}
 for(auto [i,j]:risk.correspondence.certified){int k=int(base2.size())-1-j;certified[0].at(i)=true;certified[1].at(k)=true;}
 // Risk reductions may recover only certified agreement positions.
 for(size_t i=0;i<base1.size();++i)if(!certified[0][i])risk.r1[i]=std::max(risk.r1[i],base1[i]);
 for(size_t i=0;i<base2.size();++i)if(!certified[1][i])risk.r2[i]=std::max(risk.r2[i],base2[i]);
 const std::array<fqalign::TrimInterval,2> old{fqalign::mott_segment(base1,prefix1,cut),fqalign::mott_segment(base2,prefix2,cut)},proposed{fqalign::mott_segment(risk.r1,prefix1,cut),fqalign::mott_segment(risk.r2,prefix2,cut)};
 std::array<bool,2> keep{true,true};for(int m=0;m<2;++m){auto a=old[m],b=proposed[m];for(size_t i=b.begin;i<b.end&&keep[m];++i)if((i<a.begin||i>=a.end)&&!certified[m].at(i))keep[m]=false;
  for(size_t i=a.begin;i<a.end&&keep[m];++i){if(i>=b.begin&&i<b.end)continue;int k=counterpart[m].at(i);auto o=old[1-m],p=proposed[1-m];if(k<0||size_t(k)<o.begin||size_t(k)>=o.end||size_t(k)<p.begin||size_t(k)>=p.end)keep[m]=false;}
 }
 // Counterparts must survive both old and proposed intervals, so either revert
 // preserves the other mate's coverage check. Coverage is not correctness.
 if(!keep[0])risk.r1=base1;if(!keep[1])risk.r2=base2;
}
}

