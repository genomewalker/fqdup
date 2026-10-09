#include "../src/merge_r8qkeep/compact_store.hpp"
#include "../src/merge_r8qkeep/canonical_pairs.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
using namespace fqdup::r8qkeep;
static void require(bool yes){if(!yes)throw std::runtime_error("compact store mismatch");}
static void equal(const FamilyPair& a,const FamilyPair& b){
    require(a.id==b.id);
    for(int m=0;m<2;++m)require(a.r[m].h==b.r[m].h&&a.r[m].id==b.r[m].id&&a.r[m].s==b.r[m].s&&a.r[m].q==b.r[m].q);
}
int main(){
    std::mt19937 rng(2601);
    for(unsigned alphabet:{1,2,3,5,9,16,17,94,256}){
        CompactPairStore store;std::vector<FamilyPair> reference;
        for(unsigned length=0;length<=512;++length){
            FamilyPair pair;pair.id="library-"+std::string(193, 'p')+":"+std::to_string(length);
            for(int m=0;m<2;++m){auto& r=pair.r[m];r.id=pair.id;r.h="@"+pair.id+"/"+std::to_string(m+1)+(m?"\tcomment FQT=verbatim":" comment");
                for(unsigned j=0;j<length;++j){r.s+="ACGTNry"[rng()%7];r.q+=char(1+rng()%alphabet);}}
            store.append(pair.r[0],pair.r[1]);reference.push_back(std::move(pair));
        }
        const unsigned expected=alphabet>16?8:alphabet>8?4:alphabet>4?3:alphabet>2?2:alphabet>1?1:0;
        require(store.quality_width()==expected);
        CompactPairStore::Decoder decoder;FamilyPair observed;
        for(unsigned i=0;i<reference.size();++i){store.read(i,observed,decoder);equal(observed,reference[i]);std::string id;store.id(i,id,decoder);require(id==reference[i].id);}
        for(unsigned n=0;n<2000;++n){const auto i=rng()%reference.size();store.read(i,observed,decoder);equal(observed,reference[i]);std::string id;store.id(i,id,decoder);require(id==reference[i].id);}
        // Sparse correction after parsing must expand the palette without
        // changing any previous quality, including non-printable byte values.
        for(unsigned n=0;n<600;++n){const auto i=1+rng()%(reference.size()-1),m=rng()%2,pos=rng()%reference[i].r[m].s.size();
            const char base="ACGTN"[rng()%5],quality=char(rng()%256);store.edit(i,m,pos,base,quality);
            reference[i].r[m].s[pos]=base;reference[i].r[m].q[pos]=quality;}
        for(unsigned i=0;i<reference.size();++i){store.read(i,observed,decoder);equal(observed,reference[i]);std::string id;store.id(i,id,decoder);require(id==reference[i].id);}
        // Long shared prefixes deliberately force the canonical comparator to
        // decode ties after its eight-byte prefix acceleration.
        std::vector<FamilyRead> first,second;
        for(const auto& p:reference){first.push_back(p.r[0]);second.push_back(p.r[1]);}
        canonicalize_family_input(first,second);store.canonicalize();
        for(unsigned i=0;i<reference.size();++i){store.read(i,observed,decoder);equal(observed,FamilyPair{first[i].id,{first[i],second[i]}});}
    }
    // Distinct IDs can collide in the eight-byte sorting prefix even after
    // removing the whole-library common prefix. Only headers break these ties.
    CompactPairStore collisions;std::vector<FamilyRead> ca,cb;
    for(unsigned i=0;i<160;++i){const auto id=std::string(i%2?"lane-A-common-prefix-":"lane-B-common-prefix-")+std::to_string(160-i);
        FamilyRead a{"@"+id+"/1 comment",id,"ACGTN","I!#?I"},b{"@"+id+"/2 comment",id,"TN","!I"};
        ca.push_back(a);cb.push_back(b);collisions.append(a,b);}
    canonicalize_family_input(ca,cb);require(collisions.canonicalize());
    for(std::size_t i=0;i<ca.size();++i){FamilyPair p;collisions.read(i,p);equal(p,FamilyPair{ca[i].id,{ca[i],cb[i]}});}
    CompactPairStore empty;require(empty.empty());empty.canonicalize();
    // Duplicate normalized IDs must retain all seven scalar tie-break fields.
    CompactPairStore ties;std::vector<FamilyRead> a,b;
    for(unsigned i=0;i<100;++i){FamilyRead x{"@id/1 x","id",std::string(i%7,'A'),std::string(i%7,'I')};
        FamilyRead y{"@id/2 y","id",std::string(i%9,'C'),std::string(i%9,'!')};a.push_back(x);b.push_back(y);ties.append(x,y);}
    canonicalize_family_input(a,b);ties.canonicalize();
    for(unsigned i=0;i<a.size();++i){FamilyPair p;ties.read(i,p);equal(p,FamilyPair{"id",{a[i],b[i]}});}
    bool rejected=false;try{CompactPairStore bad;bad.append({"@x","","A","I"},{"@y","","A","I"});}catch(const std::runtime_error&){rejected=true;}require(rejected);
    std::cout<<"TASK26_COMPACT_STORE_PASS\n";
}
