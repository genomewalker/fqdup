#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_set>
namespace fqdup::r8qkeep {
// B4-ADAPTER-DEF: frozen before scoring; no reverse complements or output mining.
class AdapterVocabulary {
    struct Seeds {std::unordered_set<std::uint32_t> exact,near;};
    std::array<Seeds,2> mates;
    std::size_t minimum;
    static int base(char c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return -1;}}
    static void add(Seeds& seeds,const std::string& word){
        for(std::size_t p=0;p+12<=word.size();++p){
            std::uint32_t key=0;bool valid=true;
            for(std::size_t j=0;j<12;++j){const int b=base(word[p+j]);if(b<0){valid=false;break;}key=(key<<2)|b;}
            if(!valid||!seeds.exact.insert(key).second)continue;
            seeds.near.insert(key);
            for(int j=0;j<12;++j)for(unsigned b=0;b<4;++b)seeds.near.insert((key&~(3u<<(2*j)))|(b<<(2*j)));
        }
    }
public:
    AdapterVocabulary(const std::string& first,const std::string& second,std::size_t min_length):minimum(min_length){
        add(mates[0],first);add(mates[1],second);
        for(auto& mate:mates)add(mate,"AGATCGGAAGAG");
        // The r8qkeep wrapper does not call detect_params: its only frozen
        // constructs are NativeSettings.adapter1/2 (including the declared ss adapter).
    }
    bool frozen()const{return !mates[0].exact.empty()&&!mates[1].exact.empty();}
    std::uint64_t residual(std::string_view read,int mate)const{
        if(mate<1||mate>2)return 0;
        const auto& seeds=mates[mate-1];std::uint32_t key=0,invalid=0;
        for(std::size_t i=0;i<read.size();++i){
            const int b=base(read[i]);key=((key<<2)|unsigned(b<0?0:b))&0xffffffu;invalid=((invalid<<1)|unsigned(b<0))&0xfffu;
            if(i+1<12||i+1-12<minimum)continue;
            bool match=false;
            if(!invalid)match=seeds.near.count(key);
            else if((invalid&(invalid-1))==0){const unsigned shift=2*__builtin_ctz(invalid);for(unsigned x=0;x<4;++x)if(seeds.exact.count((key&~(3u<<shift))|(x<<shift))){match=true;break;}}
            if(match)return read.size()-(i+1-12);
        }
        return 0;
    }
};
}
