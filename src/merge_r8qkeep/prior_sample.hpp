#pragma once
#include <xxhash.h>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include "family.hpp"
namespace fqdup::r8qkeep {
// ID-only Bernoulli selection, followed by deterministic ID ordering for fitting.
inline constexpr std::uint64_t prior_sample_seed=0x66716475703232ULL;
inline std::vector<std::size_t> prior_sample_indices(const std::vector<FamilyPair>& pairs,std::size_t requested) {
    if(requested==0)throw std::invalid_argument("prior subsample must be positive");
    std::vector<std::size_t> selected;
    if(pairs.empty())return selected;
    const auto threshold=(static_cast<unsigned __int128>(std::min(requested,pairs.size()))<<64)/pairs.size();
    for(std::size_t i=0;i<pairs.size();++i){const auto& id=pairs[i].id;const auto hash=XXH64(id.data(),id.size(),prior_sample_seed);if(static_cast<unsigned __int128>(hash)<threshold)selected.push_back(i);}
    std::sort(selected.begin(),selected.end(),[&](auto a,auto b){return pairs[a].id<pairs[b].id;});
    // Ambiguous duplicate IDs cannot define an ID-only deterministic fit order.
    for(std::size_t i=1;i<selected.size();++i)if(pairs[selected[i-1]].id==pairs[selected[i]].id)throw std::invalid_argument("subsample prior fit requires unique selected read IDs");
    return selected;
}
}
