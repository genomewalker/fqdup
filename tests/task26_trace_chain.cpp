#include "../src/merge_r8qkeep/repair.cpp"
namespace fqdup_r8qkeep_repair_internal {
#include "../tests/task25_full_scan_ref.hpp"
}
#include <random>
#include <chrono>
namespace internal=fqdup_r8qkeep_repair_internal;
template<class A,class B>void same(const A&a,const B&b,int trial){
 auto&x=a.corr;auto&y=b.corr;
 if(x.best!=y.best||x.paths!=y.paths||x.eligible!=y.eligible||x.reason!=y.reason||x.invariant!=y.invariant||x.certified!=y.certified||a.min_informative!=b.min_informative||a.max_gap!=b.max_gap||a.max_fraction!=b.max_fraction||a.summary_complete!=b.summary_complete){std::cerr<<"mismatch trial="<<trial<<" best="<<x.best<<","<<y.best<<" paths="<<x.paths<<","<<y.paths<<" reasons="<<x.reason<<","<<y.reason<<'\n';throw std::runtime_error("score band oracle");}
}
int main(int argc,char**){
 std::mt19937 rng(261040);internal::fqalign::FullCorrespondence old;internal::frozen_full_scan::FullCorrespondence oracle;
 std::vector<std::pair<std::string,std::string>> bench;
 for(int trial=0;trial<3000;++trial){
  int n=1+rng()%512,m=1+rng()%512;if(trial%3)n=m=150;
  std::string a(n,'A'),b(m,'A');for(auto&c:a)c="ACGTN"[rng()%(trial%4?4:5)];for(auto&c:b)c="ACGTN"[rng()%5];
  if(trial%5){b=a;if(trial%5==1)b=std::string(trial%31,'G')+b;else if(trial%5==2)b.erase(0,std::min<int>(trial%31,n-1));else if(trial%5==3)b.insert(b.size()/2,1,'T');
   for(int i=0;i<int(b.size());++i)if(rng()%100<1+trial%15)b[i]="ACGTN"[rng()%5];
   b.resize(std::min<std::size_t>(512,b.size()));
  }
  if(trial%31==0){a.assign(n,'A');b.assign(m,'A');}
  const auto expected=oracle.align_uncached(a,b);same(old.align_uncached(a,b),expected,trial);
  if(bench.size()<200&&a.size()==150&&b.size()<=180)bench.emplace_back(a,b);
 }
 std::cout<<"TRACE_CHAIN_PASS cases=3000\n";
}
