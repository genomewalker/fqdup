#include "adapter_residual.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
using fqdup::r8qkeep::AdapterVocabulary;
std::uint64_t reference(const std::string& s,const std::string& adapter,int minimum){
    std::set<std::string> seeds;
    for(const auto& word:{adapter,std::string("AGATCGGAAGAG")})for(std::size_t p=0;p+12<=word.size();++p){auto seed=word.substr(p,12);if(seed.find_first_not_of("ACGT")==std::string::npos)seeds.insert(seed);}
    for(std::size_t p=minimum;p+12<=s.size();++p)for(const auto& seed:seeds){int differences=0;for(int j=0;j<12;++j)differences+=s[p+j]!=seed[j];if(differences<=1)return s.size()-p;}
    return 0;
}
int main(){
    const std::string a="AGATCGGAAGAGCACACGTCTGAACTCCAG",b="GGAAGAGCGTCGTGTAGGGAAAGAGTGTA";
    std::mt19937 random(47);const std::string alphabet="ACGTNX";
    for(int minimum:{1,10,30}){
        AdapterVocabulary vocabulary(a,b,minimum);
        for(int trial=0;trial<10000;++trial){
            std::string read(random()%120,'A');for(char& c:read)c=alphabet[random()%alphabet.size()];
            const int mate=1+random()%2;const auto& adapter=mate==1?a:b;
            if(read.size()>=12&&trial%2==0){const auto p=random()%(read.size()-11);const auto start=random()%(adapter.size()-11);read.replace(p,12,adapter.substr(start,12));for(int n=0;n<trial%4;++n)read[p+random()%12]=alphabet[random()%alphabet.size()];}
            if(vocabulary.residual(read,mate)!=reference(read,adapter,minimum))throw std::runtime_error("adapter residual mismatch");
        }
    }
    std::cout<<"adapter residual: 30000 brute-force cases pass\n";
}
