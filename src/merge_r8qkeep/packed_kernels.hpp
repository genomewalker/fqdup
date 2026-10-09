#pragma once
#include "family_key.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string_view>
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#endif
namespace fqdup::r8qkeep {
namespace packed {
inline constexpr auto code_table=[](){std::array<unsigned char,256> codes{};for(auto& c:codes)c=4;codes['A']=0;codes['C']=1;codes['G']=2;codes['T']=3;return codes;}();
inline std::uint64_t spread(std::uint32_t mask){
    std::uint64_t x=mask;x=(x|(x<<16))&0x0000ffff0000ffffULL;
    x=(x|(x<<8))&0x00ff00ff00ff00ffULL;x=(x|(x<<4))&0x0f0f0f0f0f0f0f0fULL;
    x=(x|(x<<2))&0x3333333333333333ULL;return (x|(x<<1))&0x5555555555555555ULL;
}
struct Block {std::uint64_t codes=0;std::uint32_t exceptions=0;bool supported=true;};
inline Block scalar_block(const char* input){
    Block out;
    for(unsigned i=0;i<32;++i){const auto c=static_cast<unsigned char>(input[i]),code=code_table[c];
        out.codes|=std::uint64_t(code&3)<<(2*i);
        if(code==4){out.exceptions|=std::uint32_t{1}<<i;if(c!='N')out.supported=false;}}
    return out;
}
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("avx2"))) inline Block avx2_block(const char* input){
    const auto bytes=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input));
    const auto valid=_mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi8(bytes,_mm256_set1_epi8('A')),_mm256_cmpeq_epi8(bytes,_mm256_set1_epi8('C'))),
                                    _mm256_or_si256(_mm256_cmpeq_epi8(bytes,_mm256_set1_epi8('G')),_mm256_cmpeq_epi8(bytes,_mm256_set1_epi8('T'))));
    const auto values=_mm256_and_si256(valid,_mm256_and_si256(_mm256_xor_si256(_mm256_srli_epi16(bytes,1),_mm256_srli_epi16(bytes,2)),_mm256_set1_epi8(3)));
    const auto nibbles=_mm256_maddubs_epi16(values,_mm256_set1_epi16(0x0401));
    const auto words=_mm256_madd_epi16(nibbles,_mm256_set1_epi32(0x00100001));
    const auto shuffle=_mm_setr_epi8(0,4,8,12,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
    const auto lo=std::uint32_t(_mm_cvtsi128_si32(_mm_shuffle_epi8(_mm256_castsi256_si128(words),shuffle)));
    const auto hi=std::uint32_t(_mm_cvtsi128_si32(_mm_shuffle_epi8(_mm256_extracti128_si256(words,1),shuffle)));
    const auto five=_mm256_or_si256(valid,_mm256_cmpeq_epi8(bytes,_mm256_set1_epi8('N')));
    return {std::uint64_t(lo)|(std::uint64_t(hi)<<32),~std::uint32_t(_mm256_movemask_epi8(valid)),std::uint32_t(_mm256_movemask_epi8(five))==0xffffffffU};
}
#endif
inline Block block(const char* input){
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    static const auto encode=__builtin_cpu_supports("avx2")?avx2_block:scalar_block;
    return encode(input);
#else
    return scalar_block(input);
#endif
}
inline constexpr auto decoded_bytes=[](){
    std::array<std::array<char,4>,256> table{};
    for(unsigned b=0;b<256;++b)for(unsigned i=0;i<4;++i)table[b][i]="ACGT"[(b>>(2*i))&3];
    return table;
}();
inline void store32(unsigned char* bytes,unsigned shift,std::uint64_t codes){
    const auto low=(codes<<shift)|(bytes[0]&((1U<<shift)-1));
    for(unsigned i=0;i<8;++i)bytes[i]=static_cast<unsigned char>(low>>(8*i));
    if(shift)bytes[8]=static_cast<unsigned char>(codes>>(64-shift));
}
inline void decode32(const unsigned char* bytes,unsigned shift,char* output){
    for(unsigned i=0;i<8;++i){
        const unsigned code=(bytes[i]>>shift)|(shift?unsigned(bytes[i+1])<<(8-shift):0);
        std::memcpy(output+4*i,decoded_bytes[code&255].data(),4);
    }
}
enum class NRule {ignore,always_mismatch,equal};
class Bases {
    std::array<std::uint64_t,17> codes_{},ns_{};
    std::size_t length_=0;bool supported_=true;
    static std::uint64_t window(const std::array<std::uint64_t,17>& words,std::size_t offset){
        const auto index=offset/32,shift=2*(offset%32);
        return (words[index]>>shift)|(shift?words[index+1]<<(64-shift):0);
    }
public:
    explicit Bases(std::string_view text):length_(text.size()){
        if(length_>512){supported_=false;return;}
        std::size_t i=0;
        for(;i+32<=length_;i+=32){const auto value=block(text.data()+i);codes_[i/32]=value.codes;ns_[i/32]=spread(value.exceptions);supported_&=value.supported;}
        for(;i<length_;++i){const auto c=static_cast<unsigned char>(text[i]),code=code_table[c];codes_[i/32]|=std::uint64_t(code&3)<<(2*(i%32));
            if(code==4){ns_[i/32]|=std::uint64_t{1}<<(2*(i%32));if(c!='N')supported_=false;}}
    }
    bool supported()const{return supported_;}
    std::size_t size()const{return length_;}
    static int mismatches(const Bases& a,std::size_t first,const Bases& b,std::size_t second,std::size_t length,NRule rule,int limit=std::numeric_limits<int>::max()){
        if(!a.supported_||!b.supported_||first>a.length_||length>a.length_-first||second>b.length_||length>b.length_-second)throw std::invalid_argument("packed mismatch range");
        int total=0;
        for(std::size_t i=0;i<length;i+=32){
            const auto x=window(a.codes_,first+i)^window(b.codes_,second+i);
            auto mismatch=(x|(x>>1))&0x5555555555555555ULL;
            const auto na=window(a.ns_,first+i),nb=window(b.ns_,second+i);
            if(rule==NRule::ignore)mismatch&=~(na|nb);
            else if(rule==NRule::always_mismatch)mismatch|=na|nb;
            else mismatch=(mismatch&~(na|nb))|(na^nb);
            const auto n=std::min<std::size_t>(32,length-i);if(n<32)mismatch&=(std::uint64_t{1}<<(2*n))-1;
            total+=__builtin_popcountll(mismatch);if(total>limit)return total;
        }
        return total;
    }
};
}
// Direct distance on the stored two-bit family halves. Part values retain the
// original three-bit index alphabet, including the explicit N code of four.
struct TwoBitFamilyKey {
    std::uint64_t lo,hi;
    explicit TwoBitFamilyKey(const FamilyKey& key):lo(key.lo),hi(key.hi){if(key.empty()||!key.supported())throw std::invalid_argument("packed family alphabet");}
    explicit TwoBitFamilyKey(const std::string& key):TwoBitFamilyKey(FamilyKey(key)){}
    std::uint64_t part(int p)const{
        if(p<0||p>2)throw std::out_of_range("family part");
        const unsigned begin=p==0?0:p==1?13:26,end=p==2?40:begin+13;
        std::uint64_t value=0;
        for(unsigned i=begin;i<end;++i){const auto word=i<20?lo:hi;const auto pos=i%20;
            const auto code=word&(1ULL<<(40+pos))?4:(word>>(2*pos))&3;value|=code<<(3*(i-begin));}
        return value;
    }
    static std::uint64_t differences(std::uint64_t a,std::uint64_t b){
        const auto x=a^b;auto bits=(x|(x>>1))&0x5555555555ULL;
        const auto ns=std::uint32_t((a|b)>>40)&0xfffffU;if(ns)bits|=packed::spread(ns);
        return bits;
    }
};
static_assert(sizeof(TwoBitFamilyKey)==16);
}
