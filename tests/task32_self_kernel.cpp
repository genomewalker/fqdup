#include "../src/merge_r8qkeep/repair.cpp"
#include <random>
#include <fstream>
namespace internal=fqdup_r8qkeep_repair_internal;
using Engine=internal::fqalign::SelfCorrespondence;
static std::size_t checks=0;
static void compare(const std::string& a,const std::string& b){
 Engine scalar,optimized;
 const auto x=scalar.align_scalar(a,b),y=optimized.align_uncached(a,b);
 const auto& p=x.corr;const auto& q=y.corr;
 if(p.best!=q.best||p.paths!=q.paths||p.eligible!=q.eligible||p.reason!=q.reason||p.invariant!=q.invariant||p.certified!=q.certified||x.min_informative!=y.min_informative||x.max_gap!=y.max_gap||x.max_fraction!=y.max_fraction||x.summary_complete!=y.summary_complete)
  throw std::runtime_error("Self kernel mismatch lengths "+std::to_string(a.size())+","+std::to_string(b.size()));
 ++checks;
}
int main(int argc,char** argv){
 unsetenv("FQDUP_SCALAR_KERNELS");
 std::mt19937 rng(32);const char* bases="ACGTN";
 if(argc==1){
 for(int n=0;n<=512;++n){
  std::string a(n,'A'),b(rng()%513,'A');
  for(char& c:a)c=bases[rng()%5];for(char& c:b)c=bases[rng()%5];
  compare(a,b);compare(b,a);compare(a,a);
  for(int k=0;k<n;++k)a[k]="ACGT"[k%4];
  compare(a,a);compare(a,a.substr(n?1:0));
  if(n<=150){
   for(int pos=0;pos<n;++pos){auto unknown=a;unknown[pos]='N';compare(a,unknown);compare(unknown,a);}
   compare(std::string(n,'N'),a);compare(std::string(n,'A'),std::string(n,'A'));
  }
 }
 for(int n=0;n<=12;++n)for(int m=0;m<=12;++m){
  compare(std::string(n,'A'),std::string(m,'A'));compare(std::string(n,'N'),std::string(m,'N'));
 }
 compare("acgt","ACGT");compare("X","A");compare(std::string(513,'A'),"A");
 }
 if(argc==3){
  std::ifstream f1(argv[1]),f2(argv[2]);if(!f1||!f2)throw std::runtime_error("prefix input open");
  std::string h,a,p,q,h2,b,p2,q2;
  for(int row=0;row<256;++row){
   if(!std::getline(f1,h)||!std::getline(f1,a)||!std::getline(f1,p)||!std::getline(f1,q)||!std::getline(f2,h2)||!std::getline(f2,b)||!std::getline(f2,p2)||!std::getline(f2,q2))throw std::runtime_error("short real prefix");
   compare(a,a);compare(b,b);compare(a,b);
   setenv("FQDUP_SCALAR_KERNELS","1",1);
   Engine runtime;const auto old=runtime.align_uncached(a,b);
   unsetenv("FQDUP_SCALAR_KERNELS");const auto next=runtime.align_uncached(a,b);
   if(old.corr.invariant!=next.corr.invariant||old.corr.certified!=next.corr.certified||old.corr.best!=next.corr.best||old.corr.reason!=next.corr.reason)throw std::runtime_error("runtime prefix oracle");
  }
 }
 std::cout<<"TASK32_SELF_ORACLE_PASS checks="<<checks<<"\n";
}
