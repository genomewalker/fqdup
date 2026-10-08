#include "merge_r8qkeep/prior_sample.hpp"
#include <cassert>
#include <iostream>
using namespace fqdup::r8qkeep;
int main(){
    std::vector<FamilyPair> pairs(10000);for(std::size_t i=0;i<pairs.size();++i)pairs[i].id="fixture_"+std::to_string(i);
    auto ids=[&](std::size_t n){std::vector<std::string> out;for(auto i:prior_sample_indices(pairs,n))out.push_back(pairs[i].id);return out;};
    const auto sample=ids(2500);assert(!sample.empty()&&sample.size()<pairs.size());
    std::reverse(pairs.begin(),pairs.end());assert(ids(2500)==sample);
    for(auto& p:pairs){p.r[0].s="NNNN";p.r[0].q="!!!!";}assert(ids(2500)==sample);
    const auto all=ids(20000);assert(all.size()==pairs.size());
    bool rejected=false;try{prior_sample_indices(pairs,0);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    pairs[0].id=pairs[1].id;rejected=false;try{prior_sample_indices(pairs,20000);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    std::cout<<"TASK22_ID_SELECTION_PASS selected="<<sample.size()<<" seed="<<prior_sample_seed<<'\n';
}
