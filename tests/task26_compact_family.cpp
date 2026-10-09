#include "../src/merge_r8qkeep/compact_family.hpp"
#include "../src/merge_r8qkeep/canonical_pairs.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
using namespace fqdup::r8qkeep;
static void require(bool yes){if(!yes)throw std::runtime_error("compact family mismatch");}
int main(){
    std::mt19937 rng(260101);std::vector<FamilyRead> first,second;
    const unsigned qs[]={0,2,10,20,30,40};
    for(unsigned group=0;group<90;++group){
        std::array<std::string,2> source;
        for(auto& sequence:source)for(unsigned j=0;j<18+group%134;++j)sequence+="ACGT"[rng()%4];
        for(unsigned member=0;member<2+group%15;++member){
            const std::string id="group"+std::to_string(group)+"-member"+std::to_string(member);
            for(unsigned mate=0;mate<2;++mate){
                FamilyRead r;r.id=id;r.h="@"+id+"/"+std::to_string(mate+1)+" comment";r.s=source[mate];
                r.q.assign(r.s.size(),char(33+qs[(group+member+mate)%6]));
                for(unsigned j=20;j<r.s.size();++j)if(rng()%29==0)r.s[j]="ACGTN"[rng()%5];
                if(member%5==0&&r.s.size()>24){r.s.resize(r.s.size()-4);r.q.resize(r.s.size());}
                (mate?second:first).push_back(std::move(r));
            }
        }
    }
    for(int threads:{1,8})for(bool canonical:{false,true})for(bool debug:{false,true}){
        auto a=first,b=second;CompactPairStore store;
        for(std::size_t i=0;i<a.size();++i)store.append(a[i],b[i]);
        if(canonical){canonicalize_family_input(a,b);require(store.canonicalize());}
        const auto ref=correct_family(std::move(a),std::move(b),threads,debug);
        const auto compact=correct_compact_family(store,threads,debug,nullptr);
        if(debug){require(ref.edits_tsv==compact.edits_tsv);require(ref.families_tsv==compact.families_tsv);}require(ref.stats_json==compact.stats_json);
        FamilyPair pair;CompactPairStore::Decoder decoder;
        for(std::size_t i=0;i<store.size();++i){store.read(i,pair,decoder);require(pair.id==ref.pairs[i].id);
            for(int mate=0;mate<2;++mate){const auto& x=pair.r[mate];const auto& y=ref.pairs[i].r[mate];require(x.h==y.h&&x.id==y.id&&x.s==y.s&&x.q==y.q);}}
    }
    for(int threads:{1,8}){CompactPairStore empty;const auto compact=correct_compact_family(empty,threads,true,nullptr);const auto ref=correct_family({},{},threads,true);require(compact.edits_tsv==ref.edits_tsv&&compact.families_tsv==ref.families_tsv&&compact.stats_json==ref.stats_json);}
    std::cout<<"TASK26_COMPACT_FAMILY_PASS\n";
}
