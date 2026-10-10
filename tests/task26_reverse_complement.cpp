#include "../src/merge_r8qkeep/reverse_complement.hpp"
#include <random>
#include <iostream>
#include <chrono>
inline char scalar_complement(char c){switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';case 'N':return 'N';default:throw std::invalid_argument("unsupported base");}}
inline void scalar(std::string_view s,std::string&out){out.resize(s.size());for(std::size_t i=0;i<s.size();++i)out[s.size()-1-i]=scalar_complement(s[i]);}
int main(){
 std::mt19937 rng(261041);std::size_t cases=0;std::string a,b;
 const auto check=[&](const std::string&s){a.assign(s.size(),'x');b=a;bool ta=false,tb=false;try{scalar(s,a);}catch(const std::invalid_argument&){ta=true;}try{fqdup_fast_bases::reverse_complement_into(s,b);}catch(const std::invalid_argument&){tb=true;}if(ta!=tb||a!=b)throw std::runtime_error("reverse complement mismatch");++cases;};
 for(int n=0;n<=1024;++n){std::string s(n,'A');for(auto&c:s)c="ACGTN"[rng()%5];check(s);if(n){for(int byte=0;byte<256;++byte){auto bad=s;bad[rng()%n]=char(byte);check(bad);}}}
 for(int pos=0;pos<97;++pos)for(int byte=0;byte<256;++byte){std::string s(97,'C');s[pos]=char(byte);check(s);}
 std::cout<<"TASK26_REVERSE_COMPLEMENT_PASS cases="<<cases<<'\n';
}
