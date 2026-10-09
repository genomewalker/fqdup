#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace fqdup::r8qkeep {
// Compact postings over distinct 13/14-base parts, not over all reads.
class FamilyPartIndex {
    struct Part { std::uint64_t value; int begin, end; };
    std::vector<int> members_;
    std::vector<Part> parts_;
    std::array<std::vector<int>, 4> blocks_;
    std::array<std::vector<std::uint32_t>, 4> offsets_;
    std::array<int, 5> bounds_{};
    int length_ = 0;
    static int code(char c) {
        switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 4; }
    }
    std::uint32_t block(std::uint64_t value, int b) const {
        return (value >> (3 * bounds_[b])) & ((std::uint64_t(1) << (3 * (bounds_[b+1] - bounds_[b]))) - 1);
    }
    bool has_n(std::uint64_t value, int b) const {
        for (int p = bounds_[b]; p < bounds_[b+1]; ++p) if (((value >> (3*p)) & 7) == 4) return true;
        return false;
    }
    static constexpr std::array<int,6> first_block{0,0,0,1,1,2};
    static constexpr std::array<int,6> second_block{1,2,3,2,3,3};
    std::uint32_t pair_key(std::uint64_t value,int stream) const {
        return (block(value,first_block[stream])<<12)|block(value,second_block[stream]);
    }
    int posting(int stream,std::size_t position) const {
        return (std::uint32_t(blocks_[stream/2][position])>>(16*(stream%2)))&65535;
    }
    template<class Visitor>
    void visit_pairs(std::uint64_t query,Visitor visitor,std::uint64_t& lookups,std::uint64_t& candidates) const {
        std::array<std::size_t,6> at{},end{};
        for(int b=0;b<6;++b){
            if(has_n(query,first_block[b])||has_n(query,second_block[b]))continue;
            ++lookups;const auto target=pair_key(query,b);
            const auto key_at=[&](std::size_t pos){return pair_key(parts_[posting(b,pos)].value,b);};
            std::size_t lo=0,hi=parts_.size();
            while(lo<hi){const auto mid=lo+(hi-lo)/2;if(key_at(mid)<target)lo=mid+1;else hi=mid;}
            at[b]=end[b]=lo;
            if(lo==parts_.size()||key_at(lo)!=target)continue;
            // Upper-bound ties are usually nearby; bracket before binary search.
            const auto begin=lo,remaining=parts_.size()-begin;std::size_t span=1;
            while(span<remaining&&key_at(begin+span)==target)span*=2;
            lo=begin+(span==1?1:span/2+1);hi=begin+std::min(span+1,remaining);
            while(lo<hi){const auto mid=lo+(hi-lo)/2;if(key_at(mid)<=target)lo=mid+1;else hi=mid;}
            end[b]=lo;
        }
        for(;;){
            int minimum=-1;
            for(int b=0;b<6;++b)if(at[b]<end[b]){
                const int id=posting(b,at[b]);if(minimum<0||id<minimum)minimum=id;
            }
            if(minimum<0)return;
            ++candidates;const auto& part=parts_[minimum];const auto x=query^part.value;
            const auto different=((x|(x>>1)|(x>>2))|(query>>2))&0x9249249249ULL;
            auto remaining=different&(different-1);remaining&=remaining-1;
            if(remaining==0)for(int i=part.begin;i<part.end;++i)if(!visitor(members_[i]))return;
            for(int b=0;b<6;++b)if(at[b]<end[b]&&posting(b,at[b])==minimum)++at[b];
        }
    }
public:
    static std::uint64_t pack(const std::string& key, int begin, int end) {
        std::uint64_t value = 0;
        for (int p = begin; p < end; ++p) value |= std::uint64_t(code(key[p])) << (3 * (p-begin));
        return value;
    }
    void build(const std::vector<std::string>& keys, const std::vector<int>& representatives, int begin, int end) {
        length_ = end - begin;
        bounds_ = {0, 4, length_ == 13 ? 7 : 8, length_ == 13 ? 10 : 11, length_};
        members_.resize(representatives.size());
        std::iota(members_.begin(), members_.end(), 0);
        // Pack once instead of re-reading 13/14 bases on every sort comparison.
        std::vector<std::uint64_t> packed;packed.reserve(representatives.size());
        for(int i:representatives)packed.push_back(pack(keys[i],begin,end));
        const auto value = [&](int i) { return packed[i]; };
        parts_.clear();
        std::sort(members_.begin(), members_.end(), [&](int a, int b) {
            const auto x = value(a), y = value(b); return x < y || (x == y && a < b);
        });
        std::size_t count = 0;
        for (std::size_t i = 0; i < members_.size(); ++i) if (!i || value(members_[i]) != value(members_[i-1])) ++count;
        parts_.reserve(count);
        for (std::size_t i = 0; i < members_.size();) {
            const auto v = value(members_[i]); const int first = i++;
            while (i < members_.size() && value(members_[i]) == v) ++i;
            parts_.push_back({v, first, int(i)});
        }
        std::vector<std::uint64_t>().swap(packed);
        if(parts_.size()<=65536){
            // Six 16-bit ID streams occupy three existing 32-bit vectors.
            // Pair keys are reconstructed from parts_, so no key arrays persist.
            static_assert(sizeof(int)==sizeof(std::uint32_t));
            for(auto& offsets:offsets_)std::vector<std::uint32_t>().swap(offsets);
            std::vector<int>().swap(blocks_[3]);
            for(int b=0;b<3;++b)blocks_[b].resize(parts_.size());
            std::vector<std::uint16_t> ids(parts_.size());
            for(int b=0;b<6;++b){
                std::iota(ids.begin(),ids.end(),std::uint16_t(0));
                std::sort(ids.begin(),ids.end(),[&](int a,int c){
                    const auto x=pair_key(parts_[a].value,b),y=pair_key(parts_[c].value,b);
                    return x<y||(x==y&&a<c);
                });
                const int shift=16*(b%2);
                for(std::size_t i=0;i<ids.size();++i){
                    auto& entry=blocks_[b/2][i];
                    const std::uint32_t word=(std::uint32_t(entry)&~(std::uint32_t(65535)<<shift))|(std::uint32_t(ids[i])<<shift);
                    std::memcpy(&entry,&word,sizeof(word));
                }
            }
            return;
        }
        for (int b = 0; b < 4; ++b) {
            const std::size_t domain=std::size_t(1) << (3*(bounds_[b+1]-bounds_[b]));
            auto& offsets=offsets_[b];offsets.assign(domain+1,0);
            for(const auto& part:parts_)++offsets[block(part.value,b)+1];
            std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
            blocks_[b].resize(parts_.size());
            auto cursor=offsets;
            // Ascending part IDs give sorted postings within each block value.
            for(std::size_t i=0;i<parts_.size();++i)blocks_[b][cursor[block(parts_[i].value,b)]++]=i;
        }
    }
    // Two differences (including query N/N) leave two exact N-free blocks.
    // Merge four sorted posting ranges with fixed query scratch. A singleton
    // ID can jump to the next range's minimum because no earlier ID can recur.
    template<class Visitor>
    void visit(std::uint64_t query, Visitor visitor, std::uint64_t& lookups, std::uint64_t& candidates) const {
        if(offsets_[0].empty()){visit_pairs(query,visitor,lookups,candidates);return;}
        std::array<std::size_t,4> at{},end{};
        for(int b=0;b<4;++b){
            if(has_n(query,b))continue;
            const auto target=block(query,b);++lookups;
            at[b]=offsets_[b][target];end[b]=offsets_[b][target+1];
        }
        for(;;){
            int minimum=-1,second=-1,first=-1,hits=0,active=0;
            for(int b=0;b<4;++b)if(at[b]<end[b]){
                ++active;const int id=blocks_[b][at[b]];
                if(minimum<0||id<minimum){second=minimum;minimum=id;first=b;hits=1;}
                else if(id==minimum)++hits;
                else if(second<0||id<second)second=id;
            }
            if(active<2)return;
            if(hits==1){
                const auto& entries=blocks_[first];
                const auto begin=at[first]+1,limit=end[first];
                if(begin==limit||entries[begin]>=second)at[first]=begin;
                else {
                    // The destination is often nearby. Exponentially bracket
                    // the same lower_bound without scanning the whole tail.
                    std::size_t span=1;const auto remaining=limit-begin;
                    while(span<remaining&&entries[begin+span]<second)span*=2;
                    const auto low=begin+span/2+1;
                    const auto high=begin+std::min(span+1,remaining);
                    at[first]=std::lower_bound(entries.begin()+low,entries.begin()+high,second)-entries.begin();
                }
                continue;
            }
            ++candidates;
            const auto& part=parts_[minimum];
            const auto x=query^part.value;
            // Collapse each 3-bit code to one difference bit; query code 4
            // adds its N bit even when the stored code is also N.
            const auto different=((x|(x>>1)|(x>>2))|(query>>2))&0x9249249249ULL;
            // Unsigned zero stays zero; clearing twice leaves zero iff <=2
            // difference bits were set. Candidate order and N rules are intact.
            auto remaining=different&(different-1);
            remaining&=remaining-1;
            if(remaining==0)
                for(int i=part.begin;i<part.end;++i)if(!visitor(members_[i]))return;
            for(int b=0;b<4;++b)if(at[b]<end[b]&&blocks_[b][at[b]]==minimum)++at[b];
        }
    }
    std::size_t cells() const { return parts_.size(); }
    std::size_t max_bucket() const {
        std::size_t n = 0; for (const auto& p : parts_) n = std::max(n, std::size_t(p.end-p.begin)); return n;
    }
    std::size_t bytes() const {
        std::size_t n = members_.capacity()*sizeof(int) + parts_.capacity()*sizeof(Part);
        for (const auto& b : blocks_) n += b.capacity()*sizeof(int);
        for (const auto& offsets : offsets_) n += offsets.capacity()*sizeof(std::uint32_t);
        return n;
    }
};
}
