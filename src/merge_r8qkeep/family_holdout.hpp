#pragma once
#include "family.hpp"
#include <algorithm>
#include <stdexcept>
namespace fqdup::r8qkeep {
// Raw g0, before any family edits or member exclusions (B3-POPULATION).
inline void record_family_holdout(const std::vector<FamilyPair>& pairs,const std::vector<int>& group,HoldoutBins& bins){
    if(group.size()<3||group.size()>50)return;
    const std::string alphabet="ACGT";
    for(int side=0;side<2;++side){
        std::size_t length=0;for(int i:group)length=std::max(length,pairs[i].r[side].s.size());
        for(std::size_t pos=0;pos<length;++pos){
            std::array<int,256> counts{};
            for(int i:group){const auto& r=pairs[i].r[side];if(pos<r.s.size())++counts[static_cast<unsigned char>(r.s[pos])];}
            for(int i:group){
                const auto& r=pairs[i].r[side];if(pos>=r.s.size())continue;
                const unsigned char b=r.s[pos];const auto base=alphabet.find(b);if(base==std::string::npos)continue;
                int best=-1,top=0,second=0;
                for(int x=0;x<256;++x)if(counts[x]){const int n=counts[x]-(x==b);if(n>top){second=top;top=n;best=x;}else second=std::max(second,n);}
                if(top<2||top==second||best=='N')continue;
                const int q=int(static_cast<unsigned char>(r.q.at(pos)))-33;if(q<0||q>=94)throw std::runtime_error("holdout quality out of range");
                const int window=pos<10?0:pos+10>=r.s.size()?2:1;
                auto& cell=bins.cells[HoldoutBins::index(side,window,int(base),q)];++cell.n;cell.disagree+=best!=b;
            }
        }
    }
}
}
