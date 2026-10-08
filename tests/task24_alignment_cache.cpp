#include "../src/merge_r8qkeep/repair.cpp"
#include <random>
using namespace fqdup::r8qkeep;
namespace internal=fqdup_r8qkeep_repair_internal;
static void require(bool value){if(!value)throw std::runtime_error("alignment cache differential failure");}
static std::string prior(int component){std::ostringstream out;out<<"FQDUP_PAIR_MIXTURE_V1 64\n";for(int f=0;f<2;++f){out<<f<<' '<<(component>=0);for(int z=0;z<66;++z)out<<' '<<(z==component?1:0);out<<'\n';}return out.str();}
static void same(const internal::fqalign::Correspondence& a,const internal::fqalign::Correspondence& b){require(a.best==b.best&&a.paths==b.paths&&a.eligible==b.eligible&&a.reason==b.reason&&a.invariant==b.invariant&&a.certified==b.certified);}
int main(){
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
        require(cache.misses==4);auto modified=full.align(a,a);modified.corr.invariant.clear();same(full.align(a,a).corr,full.align_uncached(a,a).corr);self.align(a,a);repair.align(a,a);banded.align(a,a,60);require(cache.misses==4);
        full.align(a.substr(0,49),a);banded.align(a,a,59);require(cache.misses==6);
        const auto g=internal::fqalign::quality_independent_geometry(a,b);const auto direct=internal::fqalign::quality_independent_geometry_uncached(a,b);require(g.insert==direct.insert&&g.informative==direct.informative&&g.mismatches==direct.mismatches&&g.fraction==direct.fraction&&g.second==direct.second&&g.eligible==direct.eligible);
        internal::fqalign::quality_independent_geometry(a,b);require(cache.misses==7);
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
