#pragma once
#include <cstdint>
#include <string>
#include <string_view>

namespace fqdup::r8qkeep {
// Two halves, each containing 20 two-bit bases and a 20-bit N mask.
// The remaining bits mark eligibility and an unsupported input alphabet.
// No read data, IDs, heap allocation, or pointers are retained.
struct FamilyKey {
    std::uint64_t lo=0,hi=0;
    static constexpr std::uint64_t present=1ULL<<63;
    FamilyKey()=default;
    FamilyKey(const std::string& text):FamilyKey(std::string_view(text)){}
    explicit FamilyKey(std::string_view text){
        if(text.empty())return;
        lo=present;
        if(text.size()!=40){hi=present;return;}
        for(int i=0;i<40;++i){
            auto& word=i<20?lo:hi;const int p=i%20;
            unsigned base=0;
            switch(text[i]){
                case 'A':break;case 'C':base=1;break;case 'G':base=2;break;case 'T':base=3;break;
                case 'N':word|=1ULL<<(40+p);break;
                default:hi|=present;break;
            }
            word|=std::uint64_t(base)<<(2*p);
        }
    }
    bool empty()const{return !(lo&present);}
    bool supported()const{return !(hi&present);}
    std::size_t size()const{return empty()?0:40;}
    unsigned ns()const{return __builtin_popcountll(lo>>40&0xfffffULL)
        +__builtin_popcountll(hi>>40&0xfffffULL);}
    char operator[](std::size_t i)const{
        const auto word=i<20?lo:hi;const auto p=i%20;
        return word&(1ULL<<(40+p))?'N':"ACGT"[(word>>(2*p))&3];
    }
    std::string substr(std::size_t begin,std::size_t count)const{
        std::string text;for(std::size_t i=begin;i<begin+count&&i<size();++i)text.push_back((*this)[i]);return text;
    }
    bool operator==(const FamilyKey& other)const{return lo==other.lo&&hi==other.hi;}
    std::size_t hash()const{
        auto mix=[](std::uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);};
        return mix(lo)^mix(hi+0x9e3779b97f4a7c15ULL);
    }
};
static_assert(sizeof(FamilyKey)==16);
}
