#include "../src/merge_r8qkeep/compact_store.hpp"
#include <iostream>
#include <random>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("compact join mismatch");}
int main(){
    std::mt19937 rng(2604);
    for(unsigned symbols:{1,2,3,5,9,16,17,94,256}){
        CompactPairStore joined,serial;std::vector<FamilyPair> truth;
        for(unsigned block=0;block<9;++block){
            CompactPairStore part;
            for(unsigned i=0;i<(block==8?7U:32U);++i){
                FamilyPair pair;pair.id="join-"+std::to_string(truth.size());
                for(int m=0;m<2;++m){auto&r=pair.r[m];r.id=pair.id;r.h="@"+pair.id+"/"+std::to_string(m+1)+" comment";
                    for(unsigned j=0,n=rng()%513;j<n;++j){r.s+="ACGTNry"[rng()%7];r.q+=char((block+(rng()%symbols))%256);}}
                part.append(pair.r[0],pair.r[1]);serial.append(pair.r[0],pair.r[1]);truth.push_back(pair);
            }
            joined.append_store(part);
        }
        require(joined.size()==truth.size()&&joined.quality_width()==serial.quality_width());
        FamilyPair observed;CompactPairStore::Decoder cursor;
        for(std::size_t i=0;i<truth.size();++i){joined.read(i,observed,cursor);require(observed.id==truth[i].id);
            for(int m=0;m<2;++m)require(observed.r[m].h==truth[i].r[m].h&&observed.r[m].s==truth[i].r[m].s&&observed.r[m].q==truth[i].r[m].q);}
        for(std::size_t i=truth.size();i>0;--i){joined.read(i-1,observed,cursor);require(observed.r[0].h==truth[i-1].r[0].h);}
    }
    std::cout<<"encoded arena joins: offsets, headers, N bytes and quality palettes PASS\n";
}
