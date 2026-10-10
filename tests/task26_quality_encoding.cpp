#include "../src/merge_r8qkeep/quality_encoding.hpp"
#include "task26_quality_encoding_oracle.hpp"
#include <random>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
int main(){std::mt19937 rng(260113);std::array<unsigned char,256> alphabet;for(unsigned i=0;i<256;++i)alphabet[i]=i;std::shuffle(alphabet.begin(),alphabet.end(),rng);std::uint64_t checked=0;
 for(unsigned count:{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,94,256})for(unsigned trial=1;trial<=640;++trial){const unsigned length=trial<=512?trial:trial==640?4096:trial==639?4095:513+rng()%(4096-512);std::string text(length,'!');for(unsigned i=0;i<length;){const auto c=char(alphabet[count ? rng()%count : (rng()%100<98?0:1+rng()%2)]);const unsigned end=std::min<unsigned>(length,i+(trial%2?1:1+rng()%128));std::fill(text.begin()+i,text.begin()+end,c);i=end;}
   std::vector<unsigned char> a(17,0x5a),b=a;frozen_quality_encoding::encode(text,a);fqdup::r8qkeep::quality_encoding::encode(text,b);if(a!=b)throw std::runtime_error("packed quality bytes");++checked;
 }std::cout<<"COMBINED_QUALITY_PASS blocks="<<checked<<" byte-identical prefixes alphabet1..256 widths0..4/8 raw/RLE tails\n";return 0;
}
