#include "../src/merge_r8qkeep/packed_kernels.hpp"
#include "../src/merge_r8qkeep/family_index.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("packed/scalar mismatch");}
static int scalar(std::string_view a,std::string_view b,packed::NRule rule){
    int n=0;for(std::size_t i=0;i<a.size();++i){const bool any=a[i]=='N'||b[i]=='N';
        n+=rule==packed::NRule::ignore?(a[i]!=b[i]&&!any):rule==packed::NRule::always_mismatch?(a[i]!=b[i]||any):a[i]!=b[i];}return n;
}
int main(){
    std::mt19937 rng(2602);
    for(unsigned rep=0;rep<4096;++rep){
        std::string text(32,'A');for(char& c:text)c=rep<256?char(rep):"ACGTNry"[rng()%7];
        const auto a=packed::scalar_block(text.data()),b=packed::block(text.data());
        require(a.codes==b.codes&&a.exceptions==b.exceptions&&a.supported==b.supported);
        for(unsigned shift:{0U,2U,4U,6U}){
            std::array<unsigned char,9> bytes;bytes.fill(255);packed::store32(bytes.data(),shift,b.codes);
            require((bytes[0]&((1U<<shift)-1))==((1U<<shift)-1));
            std::string decoded(32,'?');packed::decode32(bytes.data(),shift,decoded.data());
            for(unsigned i=0;i<32;++i)require(decoded[i]=="ACGT"[packed::code_table[static_cast<unsigned char>(text[i])]&3]);
        }
    }
    for(unsigned length=0;length<=512;++length){
        for(unsigned trial=0;trial<12;++trial){
            std::string a(length,'A'),b(length,'A');for(char& c:a)c="ACGTN"[rng()%5];b=a;
            for(char& c:b)if(rng()%3==0)c="ACGTN"[rng()%5];
            const packed::Bases pa(a),pb(b);require(pa.supported()&&pb.supported());
            for(auto rule:{packed::NRule::ignore,packed::NRule::always_mismatch,packed::NRule::equal}){
                for(unsigned k=0;k<35;++k){
                    const auto x=k==0?0U:rng()%(length+1),y=k==0?0U:rng()%(length+1);
                    const auto n=k==0?length:rng()%(std::min(length-x,length-y)+1);
                    const int expected=scalar(std::string_view(a).substr(x,n),std::string_view(b).substr(y,n),rule);
                    require(packed::Bases::mismatches(pa,x,pb,y,n,rule)==expected);
                    for(int limit:{0,1,2,5,31,32,33,511}){
                        const int got=packed::Bases::mismatches(pa,x,pb,y,n,rule,limit);
                        require((got>limit)==(expected>limit));if(expected<=limit)require(got==expected);
                    }
                }
            }
        }
    }
    for(unsigned rep=0;rep<20000;++rep){
        std::string a(40,'A'),b(40,'A');for(char& c:a)c="ACGTN"[rng()%5];b=a;
        for(char& c:b)if(rng()%3==0)c="ACGTN"[rng()%5];
        const TwoBitFamilyKey pa(a),pb(b);
        const auto d0=TwoBitFamilyKey::differences(pa.lo,pb.lo),d1=TwoBitFamilyKey::differences(pa.hi,pb.hi);
        require(__builtin_popcountll(d0)+__builtin_popcountll(d1)==scalar(a,b,packed::NRule::always_mismatch));
        require(__builtin_popcountll(d0&((1ULL<<26)-1))==scalar(std::string_view(a).substr(0,13),std::string_view(b).substr(0,13),packed::NRule::always_mismatch));
        require(__builtin_popcountll(d0>>26)+__builtin_popcountll(d1&((1ULL<<12)-1))==scalar(std::string_view(a).substr(13,13),std::string_view(b).substr(13,13),packed::NRule::always_mismatch));
        for(int p=0;p<3;++p)require(pa.part(p)==FamilyPartIndex::pack(a,p==0?0:p==1?13:26,p==0?13:p==1?26:40));
    }
    require(!packed::Bases(std::string(513,'A')).supported());require(!packed::Bases("ACGTry").supported());
    std::cout<<"packed/scalar blocks, all lengths 0..512, offsets, N rules, thresholds and family parts PASS\n";
}
