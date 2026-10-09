#include "../src/merge_r8qkeep/repair.cpp"
#include <random>
#include <cstring>
namespace fqdup_r8qkeep_repair_internal {
#include "task25_full_correction_ref.hpp"
#include "task25_full_scan_ref.hpp"
#include "task25_geometry_ref.hpp"
#include "task25_native_correspondence_ref.hpp"
#include "task25_self_correspondence_ref.hpp"
}
using namespace fqdup::r8qkeep;
namespace internal=fqdup_r8qkeep_repair_internal;
static void require(bool value){if(!value)throw std::runtime_error("alignment cache differential failure");}
static std::string prior(int component){std::ostringstream out;out<<"FQDUP_PAIR_MIXTURE_V1 64\n";for(int f=0;f<2;++f){out<<f<<' '<<(component>=0);for(int z=0;z<66;++z)out<<' '<<(z==component?1:0);out<<'\n';}return out.str();}
static void same(const internal::fqalign::Correspondence& a,const internal::fqalign::Correspondence& b){require(a.best==b.best&&a.paths==b.paths&&a.eligible==b.eligible&&a.reason==b.reason&&a.invariant==b.invariant&&a.certified==b.certified);}
template<class A,class B> static void same_full(const A& a,const B& b){
    same(a.corr,b.corr);require(a.min_informative==b.min_informative&&a.max_gap==b.max_gap&&a.max_fraction==b.max_fraction);
}
static void geometry_parity(const std::string& a,const std::string& b){
    const auto expected=internal::frozen_geometry::geometry(a,b);
    const auto actual=internal::fqalign::quality_independent_geometry_uncached(a,b);
    require(expected.insert==actual.insert&&expected.informative==actual.informative&&
            expected.mismatches==actual.mismatches&&expected.fraction==actual.fraction&&
            expected.second==actual.second&&expected.eligible==actual.eligible);
}
static void suffix_cells(){
#if defined(__AVX2__)
 using Cell=internal::fqalign::AlignmentCell;constexpr int neg=-1000000000;
 const auto scalar=[](std::vector<Cell>& cells,int n,int m,int i,int j,const std::string& a,const std::string& b){
  const auto at=[&](int x,int y)->Cell&{return cells[x*(m+1)+y];};
  auto& c=at(i,j);c.through=0;int diag=neg;
  if(i<n&&j<m)diag=at(i+1,j+1).suffix[0]+(a[i]==b[j]&&a[i]!='N'?1:-1);
  c.suffix[0]=std::max(i==n||j==m?0:neg,diag);
  if(i<n)c.suffix[0]=std::max(c.suffix[0],at(i+1,j).suffix[1]-4);
  if(j<m)c.suffix[0]=std::max(c.suffix[0],at(i,j+1).suffix[2]-4);
  c.suffix[1]=diag;c.suffix[2]=diag;
  if(i<n)c.suffix[1]=std::max(c.suffix[1],at(i+1,j).suffix[1]-1);
  if(j<m)c.suffix[2]=std::max(c.suffix[2],at(i,j+1).suffix[2]-1);
 };
 std::mt19937 rng(2535);std::vector<Cell> expected,actual;
 for(int n:{1,2,7,8,9,15,16,17,30,31,32,33,100,150,511,512})
  for(int m:{1,2,7,8,9,15,16,17,30,31,32,33,100,150,511,512})
   for(int pattern=0;pattern<4;++pattern){
    std::string a(n,'A'),b(m,'A');
    if(pattern==1){a.assign(n,'N');b.assign(m,'N');}
    if(pattern>=2){for(char& c:a)c="ACGTN"[rng()%5];for(char& c:b)c="ACGTN"[rng()%5];}
    if(pattern==3){b=a;b.resize(m,'A');}
    expected.assign((n+1)*(m+1),Cell{});for(auto& c:expected)c.through=157;actual=expected;
    for(int i=n;i>=1;--i)for(int j=m;j>=1;--j)scalar(expected,n,m,i,j,a,b);
    for(int j=m;j>=1;--j)scalar(actual,n,m,n,j,a,b);
    for(int i=n-1;i>=1;--i){
     scalar(actual,n,m,i,m,a,b);int j=m-1;
     for(;j>=8;j-=8)fqdup_fast_dp::suffix_eight(&actual[i*(m+1)+j-7],&actual[(i+1)*(m+1)+j-7],a[i],b.data()+j-7);
     for(;j>=1;--j)scalar(actual,n,m,i,j,a,b);
    }
    require(std::memcmp(expected.data(),actual.data(),expected.size()*sizeof(Cell))==0);
   }
 std::cout<<"TASK25_SUFFIX_CELLS_PASS cases=1024\n";
#endif
}
int main(){
    suffix_cells();
    // All tail lengths and unaligned starts, including arbitrary byte equality.
    std::mt19937 count_rng(2528);
    std::array<char,576> left{},right{};
    for(int variant=0;variant<4;++variant) {
        for(std::size_t i=0;i<left.size();++i) {
            left[i]=variant<2?"ACGTN"[count_rng()%5]:char(count_rng()%256);
            right[i]=variant==0||variant==2?left[i]:
                variant==1?"ACGTN"[count_rng()%5]:char(count_rng()%256);
        }
        for(int offset=0;offset<32;++offset)for(int n=0;n<=512;++n) {
            const char* a=left.data()+offset;
            const char* b=right.data()+31-offset;
            int known=0,mismatch=0;
            for(int i=0;i<n;++i)if(a[i]!='N'&&b[i]!='N') {
                ++known;if(a[i]!=b[i])++mismatch;
            }
            require(geometry_counts(a,b,n)==std::make_pair(known,mismatch));
        }
    }
    std::mt19937 geometry_rng(2512);
    for(int n:{0,1,10,11,12,29,30,31,32,33,63,64,65,100,151,511,512,513})
        for(int m:{0,1,10,11,12,29,30,31,32,33,63,64,65,100,151,511,512,513})
            for(int variant=0;variant<4;++variant){
                std::string a(n,'A'),b(m,'T');
                if(variant==1){a.assign(n,'N');b.assign(m,'N');}
                if(variant>=2){
                    for(auto& c:a)c="ACGTN"[geometry_rng()%(variant==2?4:5)];
                    for(auto& c:b)c="ACGTN"[geometry_rng()%(variant==2?4:5)];
                }
                geometry_parity(a,b);
            }
    for(auto [a,b]:std::vector<std::pair<std::string,std::string>>{
            {"X","A"},{"A","X"},{"","X"},{std::string(513,'A'),"X"}}){
        bool expected=false,actual=false;
        try{internal::frozen_geometry::geometry(a,b);}catch(const std::invalid_argument&){expected=true;}
        try{internal::fqalign::quality_independent_geometry_uncached(a,b);}catch(const std::invalid_argument&){actual=true;}
        require(expected==actual);
    }
    // Degenerate first-row/column scans, rectangular inputs, caps and maximum size.
    for(int n:{0,1,2,7,31,100,512,513})for(int m:{0,1,2,7,31,100,512,513})
        for(int variant=0;variant<3;++variant){
            std::string a(n,variant==2?'N':'A'),b(m,variant==0?'A':variant==1?'C':'N');
            internal::fqalign::FullCorrespondence nominal;
            internal::frozen_full_scan::FullCorrespondence old_nominal;
            const auto actual=nominal.align_uncached(a,b);const auto expected=old_nominal.align_uncached(a,b);
            same_full(actual,expected);require(actual.summary_complete==expected.summary_complete);
            internal::fqrepair::FullCorrespondence correction;
            internal::frozen_full_correction::FullCorrespondence old_correction;
            same_full(correction.align_uncached(a,b),old_correction.align_uncached(a,b));
        }
    std::mt19937 policy_rng(2511);bool relaxed_length=false,stricter_error=false;
    std::vector<int> lengths;for(int length=20;length<=100;length+=5)lengths.push_back(length);
    lengths.insert(lengths.end(),{512,31,511,20}); // grow, shrink, reuse across engines
    for(int length:lengths)for(int variant=0;variant<8;++variant){
        std::string a(length,'A');for(auto& base:a)base="ACGT"[policy_rng()%4];auto b=a;
        if(variant)for(int i=7;i<length;i+=std::max(5,20-variant*2))b[i]=b[i]=='A'?'C':'A';
        if(variant==6)b.insert(b.begin()+length/2,'G');
        if(variant==7){a[3]='N';b[length/2]='N';}
        internal::fqalign::FullCorrespondence nominal;
        internal::fqrepair::FullCorrespondence correction;
        internal::frozen_full_correction::FullCorrespondence oracle;
        internal::fqalign::SelfCorrespondence self;
        internal::frozen_traceback::SelfCorrespondence old_self;
        same_full(self.align_uncached(a,a),old_self.align_uncached(a,a));
        same_full(self.align_uncached(a,b),old_self.align_uncached(a,b));
        internal::fqalign::BandedCorrespondence banded;
        internal::frozen_traceback::BandedCorrespondence old_banded;
        for(int insert:{length/2,length,length+5})
            same(banded.align_uncached(a,b,insert),old_banded.align_uncached(a,b,insert));
        const auto expected=oracle.align_uncached(a,b);
        const auto normal=nominal.align_uncached(a,b);
        internal::frozen_full_scan::FullCorrespondence old_nominal;
        same_full(normal,old_nominal.align_uncached(a,b));
        relaxed_length|=expected.corr.eligible&&!normal.corr.eligible;
        stricter_error|=!expected.corr.eligible&&normal.corr.eligible;
        internal::alignment_cache::Pair cache;internal::alignment_cache::Scope scope(&cache);
        if(variant%2)same_full(correction.align(a,b),expected);
        same_full(nominal.align(a,b),normal);
        same_full(correction.align(a,b),expected);
        same_full(correction.align_uncached(a,b),expected);
        require(cache.misses==1);
    }
    require(relaxed_length&&stricter_error);

    std::mt19937 rng(24);std::string input;
    for(int row=0;row<12;++row){std::string a(60,'A'),q(60,'I');for(char& c:a)c="ACGT"[rng()%4];auto b=internal::fqalign::reverse_complement(a);if(row%3==0)b[20]=b[20]=='A'?'C':'A';
        const auto raw="@same\t"+a+"\t"+q+"\t@same\t"+b+"\t"+q+"\t";
        input+=raw+"merged\t@same\t"+a+"\t"+q+"\t\t\t\n";
        input+=raw+"unmerged\t@same\t"+a.substr(0,50)+"\t"+q.substr(0,50)+"\t@same\t"+b.substr(0,53)+"\t"+q.substr(0,53)+"\n";
        input+=raw+"orphan_r1\t@same\t"+a+"\t"+q+"\t\t\t\n";
        input+=raw+"dropped\t\t\t\t\t\t\n";
        internal::alignment_cache::Pair cache;internal::alignment_cache::Scope scope(&cache);
        internal::fqalign::FullCorrespondence full;internal::fqrepair::FullCorrespondence repair;internal::fqalign::SelfCorrespondence self;internal::fqalign::BandedCorrespondence banded;
        same(full.align(a,a).corr,full.align_uncached(a,a).corr);same(repair.align(a,a).corr,repair.align_uncached(a,a).corr);same(self.align(a,a).corr,self.align_uncached(a,a).corr);same(banded.align(a,a,60),banded.align_uncached(a,a,60));
        require(cache.misses==3);auto modified=full.align(a,a);modified.corr.invariant.clear();same(full.align(a,a).corr,full.align_uncached(a,a).corr);self.align(a,a);repair.align(a,a);banded.align(a,a,60);require(cache.misses==3);
        full.align(a.substr(0,49),a);banded.align(a,a,59);require(cache.misses==5);
        const auto g=internal::fqalign::quality_independent_geometry(a,b);const auto direct=internal::fqalign::quality_independent_geometry_uncached(a,b);require(g.insert==direct.insert&&g.informative==direct.informative&&g.mismatches==direct.mismatches&&g.fraction==direct.fraction&&g.second==direct.second&&g.eligible==direct.eligible);
        internal::fqalign::quality_independent_geometry(a,b);require(cache.misses==6);
    }
    const std::size_t count=48,start=17;
    for(int component:{-1,0,65}){const auto p=prior(component);RepairPipeline pipeline(p,p);RepairAlignmentCache cache(start,count);
        const auto direct=pipeline.route(input,true,1,start),cached=pipeline.route(input,true,1,start,&cache);
        require(direct.routed_tsv==cached.routed_tsv&&direct.audit_tsv==cached.audit_tsv&&direct.routed==cached.routed&&direct.repair_bounds==cached.repair_bounds&&direct.replacements==cached.replacements);
        std::vector<bool> excluded(start,false);excluded.insert(excluded.end(),direct.routed.begin(),direct.routed.end());std::vector<int64_t> bounds(start,-1);bounds.insert(bounds.end(),direct.repair_bounds.begin(),direct.repair_bounds.end());
        for(double cut:{0.0,.0001,.01,.1,.5,1.0})for(int mode:{7,9}){
            RepairShadow a,b;const auto expected=pipeline.repair(input,mode,cut,start,1,&a,&excluded,&bounds);const auto actual=pipeline.repair(input,mode,cut,start,1,&b,&excluded,&bounds,&cache);require(expected==actual);require(a.skip_and_unchanged==b.skip_and_unchanged.load()&&a.skip_but_changed==b.skip_but_changed.load()&&a.noskip_unchanged==b.noskip_unchanged.load()&&a.noskip_changed==b.noskip_changed.load());
            require(pipeline.repair(direct.routed_tsv,mode,cut,start)==pipeline.repair(direct.routed_tsv,mode,cut,start,1,nullptr,nullptr,nullptr,&cache));
        }
    }
    std::cout<<"TASK24_ALIGNMENT_CACHE_PASS\n";
}
