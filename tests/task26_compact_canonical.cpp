#include "../src/merge_r8qkeep/compact_store.hpp"
#include "../src/merge_r8qkeep/canonical_pairs.hpp"
#include <iostream>
#include <random>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("compact canonical mismatch");}
static void equal(const FamilyPair& observed,const FamilyRead& a,const FamilyRead& b){
    require(observed.id==a.id);
    for(unsigned m=0;m<2;++m){const auto& expected=m?b:a;const auto& got=observed.r[m];
        require(std::tie(got.h,got.id,got.s,got.q)==std::tie(expected.h,expected.id,expected.s,expected.q));}
}
int main(){
    WorkerPool pool(8);std::mt19937 rng(2622);
    for(unsigned mode=0;mode<7;++mode){
        std::vector<FamilyRead> first,second;
        for(unsigned i=0;i<8219;++i){
            auto id=mode==0?std::string():mode==1?std::string("same"):mode==2?std::to_string(i):
                mode==3?"run-A-long-prefix-"+std::to_string(i%417):
                "run-"+std::string(i%2?"A-common-":"B-common-")+std::to_string(rng()%99999999);
            if(mode==4&&i%17==0)id+=std::string("\0suffix",7);
            if(mode==5)id=std::string(127,'p')+std::to_string(8219-i);
            if(mode==6){id.clear();for(unsigned j=0;j<i%17+1;++j)id+=char(128+rng()%128);}
            for(unsigned m=0;m<2;++m){FamilyRead r;
                r.h="@"+id+(id.empty()?std::string():"/"+std::to_string(m+1))+" comment-"+std::to_string(i%3);
                r.id=id;r.s=std::string(i%19,"ACGTN"[i%5]);r.q=std::string(i%19,char(33+i%80));
                (m?second:first).push_back(std::move(r));}
        }
        auto a=first,b=second;canonicalize_family_input(a,b);
        bool unique=true;for(std::size_t i=1;i<a.size();++i)if(a[i-1].id==a[i].id)unique=false;
        for(int workers:{1,2,3,8}){
            CompactPairStore store;
            for(unsigned begin=0;begin<first.size();begin+=4096){CompactPairStore part;
                for(auto i=begin;i<std::min(first.size(),std::size_t(begin+4096));++i)part.append(first[i],second[i]);
                store.append_store(part);}
            require(store.canonicalize(&pool,workers)==unique);
            FamilyPair pair;CompactPairStore::Decoder decoder;
            for(std::size_t i=0;i<a.size();++i){store.read(i,pair,decoder);equal(pair,a[i],b[i]);}
            store.input_order();decoder={};
            for(std::size_t i=0;i<first.size();++i){store.read(i,pair,decoder);equal(pair,first[i],second[i]);}
        }
    }
    for(unsigned count:{0,1,15,16,31,4095,4096}){
        CompactPairStore store;for(unsigned i=0;i<count;++i){const auto h="@"+std::to_string(i);store.append(h,"","",h,"","");}
        require(store.canonicalize(&pool,8));
    }
    WorkerPool stopped(2);stopped.stop();CompactPairStore failed;
    for(unsigned i=0;i<4096;++i){const auto h="@"+std::to_string(i);failed.append(h,"","",h,"","");}
    bool rejected=false;try{failed.canonicalize(&stopped,2);}catch(const std::runtime_error&){rejected=true;}require(rejected);
    std::cout<<"TASK26_COMPACT_CANONICAL_PASS scalar=p1=p2=p3=p8 joined chunks binary IDs duplicate tuples failure drain\n";
}
