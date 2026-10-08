#include "../src/merge_r8qkeep/repair.cpp"
namespace fqdup_r8qkeep_repair_internal {
namespace fqalign {
#include "task24_original_quality.hpp"
}
#include "task24_original_guard.hpp"
}
#include <random>
#include <cstring>
using namespace fqdup_r8qkeep_repair_internal;
static void require(bool b){if(!b)throw std::runtime_error("guard scratch differential failure");}
static std::string prior(int component){std::ostringstream o;o<<"FQDUP_PAIR_MIXTURE_V1 64\n";for(int f=0;f<2;++f){o<<f<<' '<<(component>=0);for(int z=0;z<66;++z)o<<' '<<(z==component?1:0);o<<'\n';}return o.str();}
static void same(const fqalign::PairRisk& a,const fqalign::PairRisk& b){
 require(a.r1==b.r1&&a.r2==b.r2&&a.discounted==b.discounted);
 const auto& x=a.correspondence;const auto& y=b.correspondence;
 require(x.best==y.best&&x.paths==y.paths&&x.eligible==y.eligible&&x.reason==y.reason&&x.invariant==y.invariant&&x.certified==y.certified);
}
int main(){
 std::mt19937 rng(2403);NativeQualityScratch<Record> scratch;fqalign::QualityProcessor processor;fqalign::Task24OriginalQualityProcessor original_processor;
 for(int component:{-1,0,40,65}){
  auto mixture=model(prior(component));
  for(int row=0;row<10;++row){
   const size_t n=row==9?7:35+row*4;std::string seq(n,'A'),q(n,'I');
   for(size_t i=0;i<n;++i){seq[i]="ACGT"[rng()%4];q[i]=char(33+rng()%61);}
   Record r1{"@same",seq,q},r2{"@same",fqalign::reverse_complement(seq),q};
   if(row%3==1&&n){r2.seq[8]=r2.seq[8]=='A'?'C':'A';}
   for(bool enabled:{true,false})for(int insert:{int(n),-1}){
    auto expected=original_processor.risks(r1.seq,r1.qual,r2.seq,r2.qual,insert,enabled);
    processor.risks_into(scratch.pair,r1.seq,r1.qual,r2.seq,r2.qual,insert,enabled);same(expected,scratch.pair);
   }
   for(int mode=1;mode<=9;++mode)for(double cutoff:{-1.0,0.0,.01,.2,1.0})for(int kind=0;kind<5;++kind)for(int policy=0;policy<4;++policy){
    Batch before,after;MergeRecord emitted;emitted.merged=r1;emitted.unmerged2=r2;
    if(n&&row%2){emitted.merged.seq.resize(n-3);emitted.merged.qual.resize(n-3);emitted.unmerged2.seq.resize(n-1);emitted.unmerged2.qual.resize(n-1);}
    emitted.is_merged=kind==0;emitted.is_orphan=kind==2||kind==3;emitted.orphan_r1=kind==2;
    if(kind==3)emitted.merged=emitted.unmerged2;
    if(emitted.is_merged||emitted.is_orphan)emitted.unmerged2={};
    auto fill=[&](Batch& b){if(kind==4)return;b.records.push_back(emitted);b.n_merged=kind==0;b.n_unmerged=kind==1;b.n_orphan=kind==2||kind==3;b.n_orphan_r1=kind==2;b.n_orphan_r2=kind==3;};
    const auto* m=policy==3?nullptr:&mixture;
    {Task24OriginalGuard<Batch,Record> guard(before,r1,r2,original_processor,cutoff,policy==1,policy==2?2:0,m,mode,m,policy==1);fill(before);}
    {NativeQualityGuard<Batch,Record> guard(after,r1,r2,processor,cutoff,policy==1,policy==2?2:0,m,mode,m,policy==1,&scratch);fill(after);}
    require(fqdup::r8qkeep::serialize_repair(before,23)==fqdup::r8qkeep::serialize_repair(after,23));
    require(before.n_merged==after.n_merged&&before.n_unmerged==after.n_unmerged&&before.n_orphan==after.n_orphan&&before.n_orphan_r1==after.n_orphan_r1&&before.n_orphan_r2==after.n_orphan_r2&&before.n_dropped==after.n_dropped);
   }
  }
 }
 // Repeated identical shapes keep the primary risk and raw snapshot storage.
 const auto capacity=scratch.pair.r1.capacity();processor.risks_into(scratch.pair,"A","I","T","I",-1);require(scratch.pair.r1.capacity()==capacity);
 std::cout<<"TASK24_GUARD_SCRATCH_PASS\n";
}
