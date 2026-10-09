#pragma once
#include "native_quality.hpp"
#include <cfenv>
#include <fstream>
#include <sstream>
namespace fqquality {
struct PairMixture {
 static constexpr int grid=64,K=66;
 std::array<std::array<double,K>,2> weights{};
 std::array<bool,2> supported{};
private:
 std::array<std::array<double,K>,2> original_weights{},original_logs{};
 int log_rounding=FE_TONEAREST;
public:
 // The production model is immutable; legacy callers may still edit weights.
 // Reuse only logs of the same weights computed in the same rounding mode.
 std::array<double,K> initial_log_weights(int training)const {
  if(std::fegetround()==log_rounding&&weights[training]==original_weights[training])return original_logs[training];
  std::array<double,K> out{};
  for(int j=0;j<K;++j)out[j]=weights[training][j]>0?std::log(weights[training][j]):-std::numeric_limits<double>::infinity();
  return out;
 }
 explicit PairMixture(std::istream& f) {
  std::string magic;int g;
  if(!(f>>magic>>g)||magic!="FQDUP_PAIR_MIXTURE_V1"||g!=grid)throw std::invalid_argument("invalid pair-mixture header");
  for(int row=0;row<2;++row){int fold;uint64_t pairs;
   if(!(f>>fold>>pairs)||fold!=row)throw std::invalid_argument("invalid pair-mixture fold");
   double sum=0;for(double& w:weights[fold]){if(!(f>>w)||!std::isfinite(w)||w<0||w>1)throw std::invalid_argument("invalid mixture weight");sum+=w;}
   supported[fold]=pairs>0;if((supported[fold]&&std::abs(sum-1)>1e-8)||(!supported[fold]&&sum!=0))throw std::invalid_argument("mixture weight mass");
  }
  std::string extra;if(f>>extra)throw std::invalid_argument("extra mixture data");
  original_weights=weights;log_rounding=std::fegetround();
  for(int fold=0;fold<2;++fold)for(int j=0;j<K;++j)
   original_logs[fold][j]=weights[fold][j]>0?std::log(weights[fold][j]):-std::numeric_limits<double>::infinity();
 }
 bool prior_floor(fqalign::PairRisk& risk,const std::string& s1,const std::string& q1,const std::string& s2,const std::string& q2)const {
  int training=1-fold(s1,s2);if(!supported[training])return false;
  double flat=0;for(int j=1;j<K;++j)flat+=weights[training][j]*component(j);
  for(size_t i=0;i<q1.size();++i)risk.r1.at(i)=std::max(risk.r1.at(i),weights[training][0]*fqalign::phred_error(q1[i])+flat);
  for(size_t i=0;i<q2.size();++i)risk.r2.at(i)=std::max(risk.r2.at(i),weights[training][0]*fqalign::phred_error(q2[i])+flat);
  return true;
 }
 static double component_direct(int j){return .75*std::pow(double(j-1)/grid,2);}
 static double clipped_log(double p){return std::log(std::clamp(p,1e-15,1-1e-15));}
 static std::pair<double,double> conditional(double e1,double e2,bool agree){
  if(agree){double a=e1*e2/3,den=(1-e1)*(1-e2)+a;return den>0?std::pair<double,double>{a/den,a/den}:std::pair<double,double>{1,1};}
  double d=e1+e2-4*e1*e2/3;return d>0?std::pair<double,double>{e1*(1-e2/3)/d,e2*(1-e1/3)/d}:std::pair<double,double>{1,1};
 }
 // Tables retain the exact scalar expressions and ordered accumulation below.
 // Function-local initialization is thread safe; no per-row allocations.
 struct Entry {double log;std::pair<double,double> risk;};
 struct Tables {
  std::array<double,K> components{};
  std::array<std::array<Entry,2>,K> fixed{};
  std::array<std::array<double,K>,2> fixed_log{};
  std::array<std::array<Entry,2>,94*94> nominal{};
 };
 static const Tables& tables(){
  static const Tables value=[] {
   Tables out;
   for(int t=0;t<K;++t){double e=component_direct(t);out.components[t]=e;
    if(t)for(int a=0;a<2;++a){double d=2*e-4*e*e/3;out.fixed[t][a]={clipped_log(a?1-d:d),conditional(e,e,a)};out.fixed_log[a][t]=out.fixed[t][a].log;}
   }
   for(int q1=33;q1<=126;++q1)for(int q2=33;q2<=126;++q2){
    double e1=fqalign::phred_error(char(q1)),e2=fqalign::phred_error(char(q2)),d=e1+e2-4*e1*e2/3;
    for(int a=0;a<2;++a)out.nominal[(q1-33)*94+q2-33][a]={clipped_log(a?1-d:d),conditional(e1,e2,a)};
   }
   return out;
  }();
  return value;
 }
 static double component(int t){return t>=0&&t<K?tables().components[t]:component_direct(t);}
 static Entry lookup(int t,char q1,char q2,bool agree){
  const unsigned a=static_cast<unsigned char>(q1),b=static_cast<unsigned char>(q2);
  if(t>=0&&t<K&&a>=33&&a<=126&&b>=33&&b<=126)
   return t?tables().fixed[t][agree]:tables().nominal[(a-33)*94+b-33][agree];
  // Preserve the original Phred validation and scalar path outside the domain.
  double e1=fqalign::phred_error(q1),e2=fqalign::phred_error(q2),e=component_direct(t);
  double d=t==0?e1+e2-4*e1*e2/3:2*e-4*e*e/3;
  return {clipped_log(agree?1-d:d),conditional(t==0?e1:e,t==0?e2:e,agree)};
 }
 // Bound to one immutable posterior. Collisions only evict an entry; the
 // full quality/agreement key is checked before reusing the ordered sums.
 struct ConditionalProducts {
  struct Entry {unsigned key=0;std::pair<double,double> risk{};};
  const std::array<double,K>& posterior;
  std::array<Entry,64> entries{};
  std::array<int,K> active{};
  int active_size=1;
  explicit ConditionalProducts(const std::array<double,K>& p):posterior(p){
   // Keep zero for quality validation; positive terms retain component order.
   for(int z=1;z<K;++z)if(p[z]!=0)active[active_size++]=z;
  }
  std::pair<double,double> get(char q1,char q2,bool agree){
   const unsigned key=1+static_cast<unsigned char>(q1)+
       (unsigned(static_cast<unsigned char>(q2))<<8)+(unsigned(agree)<<16);
   auto& entry=entries[(key*2654435761U)>>26];
   if(entry.key!=key){
    double p1=0,p2=0;
    for(int a=0;a<active_size;++a){const int z=active[a];auto p=lookup(z,q1,q2,agree).risk;p1+=posterior[z]*p.first;p2+=posterior[z]*p.second;}
    entry={key,{p1,p2}};
   }
   return entry.risk;
  }
 };
 struct Audit {int fold=-1;size_t n=0;double nominal=0,flat=0;bool supported=false;std::array<double,K> posterior{};};
 struct CachedApply {
  const PairMixture* model;
  int mode;bool full_trusted;
  std::string s1,q1,s2,q2;
  std::vector<std::pair<int,int>> invariant,certified;
  std::vector<double> r1,r2;
  Audit audit;
 };
 inline static thread_local std::vector<CachedApply>* active_cache=nullptr;
 struct CacheScope {
  std::vector<CachedApply>* previous;
  explicit CacheScope(std::vector<CachedApply>* cache):previous(active_cache){active_cache=cache;}
  ~CacheScope(){active_cache=previous;}
 };
 static int fold(const std::string& s1,const std::string& s2){uint64_t hash=14695981039346656037ULL;for(unsigned char c:s1){hash^=c;hash*=1099511628211ULL;}hash^='|';hash*=1099511628211ULL;for(unsigned char c:s2){hash^=c;hash*=1099511628211ULL;}return int(hash&1);}
 bool apply(fqalign::PairRisk& risk,const std::string& s1,const std::string& q1,
            const std::string& s2,const std::string& q2,int mode,Audit* audit=nullptr,bool full_trusted=false)const {
  Audit cached_audit;
  if(active_cache&&!audit)audit=&cached_audit;
  if(audit)audit->fold=fold(s1,s2);
  if(!risk.correspondence.eligible)return false;

  int training=1-fold(s1,s2);if(!supported[training])return false;
  // Only 2/4 have the nominal-posterior early return. All other modes
  // execute identical arithmetic, including router 1 and repair 9.
  const int cache_mode=(mode==2||mode==4)?2:1;
  if(active_cache)for(const auto& entry:*active_cache){
   if(entry.model==this&&entry.mode==cache_mode&&entry.full_trusted==full_trusted&&
      entry.s1==s1&&entry.q1==q1&&entry.s2==s2&&entry.q2==q2&&
      entry.invariant==risk.correspondence.invariant&&entry.certified==risk.correspondence.certified&&
      entry.r1.size()==risk.r1.size()&&entry.r2.size()==risk.r2.size()){
    risk.r1=entry.r1;risk.r2=entry.r2;*audit=entry.audit;return true;
   }
  }
  // apply does not recurse; no scratch references escape the call or cache.
  struct Scratch {std::string rc;std::vector<std::pair<int,int>> cols;};
  static thread_local Scratch scratch;auto& rc=scratch.rc;auto& cols=scratch.cols;
  fqalign::reverse_complement_into(s2,rc);cols.clear();
  for(auto c:risk.correspondence.invariant)if(s1.at(c.first)!='N'&&rc.at(c.second)!='N')cols.push_back(c);
  auto logpost=initial_log_weights(training);std::array<double,K> post{};
  std::array<int,K> active{};int active_size=0;
  for(int t=1;t<K;++t)if(weights[training][t]>0)active[active_size++]=t;
  // Zero-weight log posteriors stay -infinity under every finite update.
  // Use the original dense loop for broad support to retain vectorization.
  const bool sparse=active_size<16;
  size_t n=0;for(auto [i,j]:cols){int k=int(s2.size())-1-j;if(s1.at(i)=='N'||rc.at(j)=='N')continue;++n;
   bool agree=s1.at(i)==rc.at(j);
   // Component zero validates both qualities before the fixed table reads.
   logpost[0]+=lookup(0,q1.at(i),q2.at(k),agree).log;
   const auto& fixed_log=tables().fixed_log[agree];
   if(sparse){for(int a=0;a<active_size;++a){const int t=active[a];logpost[t]+=fixed_log[t];}}
   else for(int t=1;t<K;++t)logpost[t]+=fixed_log[t];
  }
  if(!n)return false;double maximum=*std::max_element(logpost.begin(),logpost.end()),mass=0;if(!std::isfinite(maximum))throw std::runtime_error("empty mixture posterior");
  for(int t=0;t<K;++t){post[t]=std::exp(logpost[t]-maximum);mass+=post[t];}for(double& p:post)p/=mass;
  double flat=0;for(int t=1;t<K;++t)flat+=post[t]*component(t);
  if(audit){audit->n=n;audit->nominal=post[0];audit->flat=flat;audit->supported=true;audit->posterior=post;}
  if((mode==2||mode==4)&&post[0]>=.5)return false;
  for(size_t i=0;i<q1.size();++i)risk.r1.at(i)=post[0]*fqalign::phred_error(q1[i])+flat;
  for(size_t i=0;i<q2.size();++i)risk.r2.at(i)=post[0]*fqalign::phred_error(q2[i])+flat;
  ConditionalProducts conditional_products(post);
  // All columns including N are kept when forming blocks, matching the native
  // verifier. N is neither an agreeing flank nor a conditional-risk update.
  for(size_t begin=0;begin<cols.size();){size_t end=begin+1;while(end<cols.size()&&cols[end].first==cols[end-1].first+1&&cols[end].second==cols[end-1].second+1)++end;
   size_t total=0,before=0;for(size_t t=begin;t<end;++t){auto [i,j]=cols[t];total+=s1.at(i)==rc.at(j)&&s1.at(i)!='N';}
   for(size_t t=begin;t<end;++t){auto [i,j]=cols[t];int k=int(s2.size())-1-j;bool agree=s1.at(i)==rc.at(j)&&s1.at(i)!='N';if(agree)++before;
    if(s1.at(i)=='N'||rc.at(j)=='N')continue;
    bool cert=std::find(risk.correspondence.certified.begin(),risk.correspondence.certified.end(),std::pair<int,int>{i,j})!=risk.correspondence.certified.end();
    if(!full_trusted&&!cert&&(agree||before<5||total-before<5))continue;
    const auto conditional=conditional_products.get(q1.at(i),q2.at(k),agree);
    risk.r1.at(i)=conditional.first;risk.r2.at(k)=conditional.second;
   }begin=end;
  }
  if(active_cache)active_cache->push_back({this,cache_mode,full_trusted,s1,q1,s2,q2,
      risk.correspondence.invariant,risk.correspondence.certified,risk.r1,risk.r2,*audit});
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
