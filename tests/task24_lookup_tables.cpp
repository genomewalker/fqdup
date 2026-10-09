#include "../src/merge_r8qkeep/repair.cpp"
#include <cstring>
#include <random>
#include <limits>
namespace fqalign = fqdup_r8qkeep_repair_internal::fqalign;
#include "task24_pair_mixture_direct.hpp"
using Model = fqdup_r8qkeep_repair_internal::fqquality::PairMixture;
using Direct = task24_direct::PairMixture;
static void require(bool ok){if(!ok)throw std::runtime_error("Task 24 lookup parity failure");}
static void bits(double a,double b){require(std::memcmp(&a,&b,sizeof(a))==0);}
static void vector_bits(const std::vector<double>& a,const std::vector<double>& b){require(a.size()==b.size());for(size_t i=0;i<a.size();++i)bits(a[i],b[i]);}
int main(){
 for(int t=-1;t<=Model::K;++t)bits(Model::component(t),Direct::component(t));
 for(int t=1;t<Model::K;++t)for(bool agree:{false,true}){
  double e=Direct::component(t),d=2*e-4*e*e/3;
  const auto x=Model::tables().fixed[t][agree];const auto p=Direct::conditional(e,e,agree);
  bits(x.log,Direct::clipped_log(agree?1-d:d));bits(x.risk.first,p.first);bits(x.risk.second,p.second);
  bits(x.log,Model::tables().fixed_log[agree][t]);
 }
 for(int a=33;a<=126;++a)for(int b=33;b<=126;++b)for(bool agree:{false,true}){
  double e1=fqalign::phred_error(char(a)),e2=fqalign::phred_error(char(b)),d=e1+e2-4*e1*e2/3;
  const auto x=Model::lookup(0,char(a),char(b),agree);const auto p=Direct::conditional(e1,e2,agree);
  bits(x.log,Direct::clipped_log(agree?1-d:d));bits(x.risk.first,p.first);bits(x.risk.second,p.second);
 }
 for(int t:{-1,66})for(bool agree:{false,true}){
  double e=Direct::component(t),d=2*e-4*e*e/3;auto x=Model::lookup(t,'I','J',agree);auto p=Direct::conditional(e,e,agree);
  bits(x.log,Direct::clipped_log(agree?1-d:d));bits(x.risk.first,p.first);bits(x.risk.second,p.second);
 }
 for(int q:{0,32,127,255})for(int t:{0,1,65})for(bool first:{false,true}){
  bool threw=false;try{Model::lookup(t,first?char(q):'I',first?'I':char(q),true);}catch(const std::invalid_argument&){threw=true;}require(threw);
 }
 // Sparse conditional products retain bitwise reductions under every rounding
 // mode, including zero nominal mass, underflowed terms and signed zeros.
 const int original_rounding=std::fegetround();
 for(int rounding:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
  require(std::fesetround(rounding)==0);
  for(int variant=0;variant<8;++variant){
   std::array<double,Model::K> posterior{};
   posterior[variant%2?18:0]=.75;
   posterior[variant%3?65:1]=.25;
   posterior[33]=variant%2?std::numeric_limits<double>::denorm_min():-0.;
   Model::ConditionalProducts products(posterior);
   for(char q1:{'!','I','~'})for(char q2:{'!','J','~'})for(bool agree:{false,true}){
    double p1=0,p2=0;
    for(int z=0;z<Model::K;++z){const auto p=Model::lookup(z,q1,q2,agree).risk;p1+=posterior[z]*p.first;p2+=posterior[z]*p.second;}
    for(int repeat=0;repeat<2;++repeat){const auto p=products.get(q1,q2,agree);bits(p.first,p1);bits(p.second,p2);}
   }
   for(int q:{0,32,127,255}){
    bool threw=false;try{products.get(char(q),'I',true);}catch(const std::invalid_argument&){threw=true;}require(threw);
   }
  }
 }
 require(std::fesetround(original_rounding)==0);
 std::mt19937 rng(24);
 for(int row=0;row<350;++row){
  std::ostringstream prior;prior.precision(17);prior<<"FQDUP_PAIR_MIXTURE_V1 64\n";
  for(int f=0;f<2;++f){std::array<double,66>w{};double total=0;for(int t=0;t<66;++t){w[t]=row%5==0?(t==row%66):row%5==1?(t==0||t==18?double(1+rng()%100):0):row%5==2?(t==1||t==65?double(1+rng()%100):0):double(rng()%100);total+=w[t];}
   prior<<f<<' '<<(row%17!=0);for(double v:w)prior<<' '<<(row%17==0?0:v/total);prior<<'\n';}
  std::istringstream a(prior.str()),b(prior.str());Model model(a);Direct direct(b);
  // Public weight mutation must bypass constructor logs. Signed zero has the
  // same nonpositive-weight result; swapped positive components change logs.
  if(row%3==0)for(int f=0;f<2;++f){std::swap(model.weights[f][1],model.weights[f][33]);std::swap(direct.weights[f][1],direct.weights[f][33]);}
  const int saved_rounding=std::fegetround();
  for(int rounding:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
   require(std::fesetround(rounding)==0);
   for(int f=0;f<2;++f){auto logs=model.initial_log_weights(f);for(int t=0;t<Model::K;++t)
    bits(logs[t],model.weights[f][t]>0?std::log(model.weights[f][t]):-std::numeric_limits<double>::infinity());}
  }
  require(std::fesetround(saved_rounding)==0);
  const size_t n=11+rng()%100;std::string s1(n,'A'),s2(n,'T'),q1(n,'I'),q2(n,'I');
  for(size_t i=0;i<n;++i){s1[i]="ACGTN"[rng()%5];q1[i]=char(33+rng()%94);q2[i]=char(33+rng()%94);}
  s2=fqalign::reverse_complement(s1);for(char& c:s2)if(rng()%8==0)c="ACGTN"[rng()%5];
  fqalign::PairRisk seed;seed.r1.resize(n,.2);seed.r2.resize(n,.3);seed.correspondence.eligible=row%13!=0;
  for(size_t i=0;i<n;++i){if(rng()%12==0)continue;seed.correspondence.invariant.emplace_back(i,i);if(rng()%3==0)seed.correspondence.certified.emplace_back(i,i);}
  for(int mode:{1,2,4,9})for(bool trusted:{false,true}){
   auto x=seed,y=seed;Model::Audit ax;Direct::Audit ay;
   require(model.apply(x,s1,q1,s2,q2,mode,&ax,trusted)==direct.apply(y,s1,q1,s2,q2,mode,&ay,trusted));
   vector_bits(x.r1,y.r1);vector_bits(x.r2,y.r2);require(ax.fold==ay.fold&&ax.n==ay.n&&ax.supported==ay.supported);
   bits(ax.nominal,ay.nominal);bits(ax.flat,ay.flat);for(int t=0;t<66;++t)bits(ax.posterior[t],ay.posterior[t]);
  }
  // Router mode 1 and repair mode 9 share arithmetic; 2/4 form a distinct
  // early-return class. Compare real cache hits against the uncached oracle.
  for(bool trusted:{false,true}){
   std::vector<Model::CachedApply> cache;Model::CacheScope scope(&cache);
   bool first_applied=false;
   for(int mode:{1,9,2,4,7,1}){
    const auto before=cache.size();auto x=seed,y=seed;Model::Audit ax;Direct::Audit ay;
    const bool applied=model.apply(x,s1,q1,s2,q2,mode,&ax,trusted);
    require(applied==direct.apply(y,s1,q1,s2,q2,mode,&ay,trusted));
    vector_bits(x.r1,y.r1);vector_bits(x.r2,y.r2);
    require(ax.fold==ay.fold&&ax.n==ay.n&&ax.supported==ay.supported);
    bits(ax.nominal,ay.nominal);bits(ax.flat,ay.flat);for(int t=0;t<66;++t)bits(ax.posterior[t],ay.posterior[t]);
    if(mode==1)first_applied=applied;
    if((mode==9||mode==7)&&first_applied)require(cache.size()==before);
    if(!applied)require(cache.size()==before);
   }
   require(cache.size()<=2);
  }
  auto x=seed,y=seed;require(model.prior_floor(x,s1,q1,s2,q2)==direct.prior_floor(y,s1,q1,s2,q2));vector_bits(x.r1,y.r1);vector_bits(x.r2,y.r2);
 }
 std::cout<<"TASK24_LOOKUP_TABLES_PASS\n";
}
