#include "../src/merge_r8qkeep/family_index.hpp"
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
using fqdup::r8qkeep::FamilyPartIndex;
static int distance(const std::string& a,const std::string& b){
    int d=0;for(std::size_t i=0;i<a.size();++i)d+=a[i]=='N'||a[i]!=b[i];return d;
}
int main(){
    std::mt19937 rng(2026101430);
    for(int length:{13,14}){
        std::vector<std::string> keys;
        // Local neighborhoods exercise every query N position, stored Ns,
        // two mismatches, duplicate postings, and bitmap suffix/prefix boundaries.
        std::string base(length,'A');for(char& c:base)c="ACGT"[rng()%4];
        keys.push_back(base);keys.push_back(base);
        for(int a=0;a<length;++a)for(char c:std::string("ACGTN")){
            auto first=base;first[a]=c;keys.push_back(first);
            for(int b=a+1;b<length;++b)for(char d:std::string("CGTN")){
                auto second=first;second[b]=d;keys.push_back(second);
            }
        }
        for(int i=0;i<4200;++i){auto key=base;for(char& c:key)c="ACGT"[rng()%4];keys.push_back(key);}
        std::vector<int> ids(keys.size());std::iota(ids.begin(),ids.end(),0);
        std::vector<std::string> queries{base};
        for(int a=0;a<length;++a){auto q=base;q[a]='N';queries.push_back(q);}
        auto q=base;q[0]=q[length-1]='N';queries.push_back(q);
        for(int radius:{1,2})for(bool dense:{false,true}){
            FamilyPartIndex index(radius);index.build(keys,ids,0,length,dense);
            for(const auto& query:queries){
                std::set<int> expected;for(int i:ids)if(distance(query,keys[i])<=radius)expected.insert(i);
                for(bool ordered:{false,true}){
                    std::set<int> actual;std::uint64_t lookups=0,candidates=0;
                    auto visitor=[&](int i){if(!actual.insert(i).second)throw std::runtime_error("duplicate candidate");return true;};
                    if(ordered)index.visit(FamilyPartIndex::pack(query,0,length),visitor,lookups,candidates);
                    else index.visit_unordered(FamilyPartIndex::pack(query,0,length),visitor,lookups,candidates);
                    if(actual!=expected)throw std::runtime_error("radius oracle mismatch");
                }
            }
        }
    }
    std::cout<<"TASK30_RADIUS_ORACLE_PASS\n";
}
