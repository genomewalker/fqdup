#include "../src/merge_r8qkeep/compact_store.hpp"
#include "../src/merge_r8qkeep/canonical_pairs.hpp"
#include <random>
#include <iostream>
using namespace fqdup::r8qkeep;
int main(){
    std::mt19937 rng(2701);WorkerPool pool(8);
    for(unsigned shared=0;shared<=65;++shared){
        std::vector<FamilyRead> first,second;
        const unsigned count=shared==65?12017:517;
        for(unsigned i=0;i<count;++i){
            std::string id(shared,'P');
            if(shared==65)id="A00706:736:HMGMNDSX5:1:"+std::to_string(1101+rng()%16)+":"+std::to_string(rng()%30000)+":"+std::to_string(rng()%30000);
            else if(i%5)for(unsigned j=0,n=rng()%48;j<n;++j)id+=char('A'+rng()%26);
            if(i%7==0&&!first.empty())id=first.back().id;
            if(i%13==0)id+=std::string("\0Z",2);
            for(unsigned m=0;m<2;++m){FamilyRead r;r.id=id;
                r.h="@"+id+(id.empty()?std::string():"/"+std::to_string(m+1))+(i%2?" 1:N:0:TAG FQT:i:27":" comment");
                r.s=std::string(rng()%25,"ACGTN"[rng()%5]);r.q=std::string(r.s.size(),char(33+rng()%80));
                (m?second:first).push_back(std::move(r));}
        }
        auto a=first,b=second;canonicalize_family_input(a,b);
        bool unique=true;for(unsigned i=1;i<count;++i)if(a[i-1].id==a[i].id)unique=false;
        for(int workers:{1,8}){
            CompactPairStore store;for(unsigned i=0;i<count;++i)store.append(first[i],second[i]);
            if(store.canonicalize(&pool,workers)!=unique)throw std::runtime_error("uniqueness");
            FamilyPair pair;CompactPairStore::Decoder decoder;
            for(unsigned i=0;i<count;++i){store.read(i,pair,decoder);
                for(unsigned m=0;m<2;++m){const auto& expected=m?b[i]:a[i];const auto& got=pair.r[m];
                    if(std::tie(got.id,got.h,got.s,got.q)!=std::tie(expected.id,expected.h,expected.s,expected.q))throw std::runtime_error("ordering");}}
        }
    }
    std::cout<<"TASK27_CANONICAL_PASS shared=0..65 Illumina variable-length duplicate comments FQT binary exact-tuple p1/p8\n";
}
