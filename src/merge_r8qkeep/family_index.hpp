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
    std::array<std::vector<std::uint32_t>,6> pair_offsets_;
    std::array<int,6> pair_shift_{};
    std::array<int, 5> bounds_{};
    int length_ = 0;
    bool paired_=true;
    bool compact_pairs_=false;
    template<class Id,class Key>
    static void radix_ids(std::vector<Id>& ids,Key key,int bits) {
        // Stable digits preserve the ascending-ID tie break from initial iota.
        // Keys are at most 42 bits (part) or 24 bits (two sub-blocks).
        std::vector<Id> scratch(ids.size());
        std::array<std::size_t,4096> positions{};
        for(int shift=0;shift<bits;shift+=12){
            positions.fill(0);
            for(auto id:ids)++positions[(key(id)>>shift)&4095];
            std::size_t offset=0;
            for(auto& count:positions){const auto next=offset+count;count=offset;offset=next;}
            for(auto id:ids)scratch[positions[(key(id)>>shift)&4095]++]=id;
            ids.swap(scratch);
        }
    }
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
    int pair_bits(int stream) const {
        const int a=first_block[stream],b=second_block[stream];
        return 2*(bounds_[a+1]-bounds_[a]+bounds_[b+1]-bounds_[b]);
    }
    static std::uint32_t compact_block(std::uint32_t value){
        return (value&3)|((value>>1)&12)|((value>>2)&48)|((value>>3)&192);
    }
    std::uint32_t compact_pair_key(std::uint64_t value,int stream) const {
        const int a=first_block[stream],b=second_block[stream];
        const auto first=block(value,a),second=block(value,b);
        if((first|second)&0x924U)return std::uint32_t(1)<<pair_bits(stream);
        return (compact_block(first)<<(2*(bounds_[b+1]-bounds_[b])))|compact_block(second);
    }
    int posting(int stream,std::size_t position) const {
        return (std::uint32_t(blocks_[stream/2][position])>>(16*(stream%2)))&65535;
    }
    template<class Visitor>
    void visit_pairs(std::uint64_t query,Visitor visitor,std::uint64_t& lookups,std::uint64_t& candidates) const {
        std::array<std::size_t,6> at{},end{};
        constexpr int exhausted=65536; // one past every possible 16-bit ID
        std::array<int,6> heads;heads.fill(exhausted);
        std::array<std::uint32_t,4> query_blocks;
        std::array<bool,4> n_free;
        for(int b=0;b<4;++b){query_blocks[b]=block(query,b);n_free[b]=!has_n(query,b);}
        for(int b=0;b<6;++b){
            if(!n_free[first_block[b]]||!n_free[second_block[b]])continue;
            ++lookups;const auto target=query_blocks[second_block[b]];
            const auto key_at=[&](std::size_t pos){return block(parts_[posting(b,pos)].value,second_block[b]);};
            const auto& directory=compact_pairs_?pair_offsets_[b]:offsets_[first_block[b]];
            const auto prefix=compact_pairs_?(compact_pair_key(query,b)>>pair_shift_[b]):query_blocks[first_block[b]];
            std::size_t lo=directory.empty()?0:directory[prefix];
            const std::size_t limit=directory.empty()?parts_.size():directory[prefix+1];
            std::size_t hi=limit;
            while(lo<hi){const auto mid=lo+(hi-lo)/2;if(key_at(mid)<target)lo=mid+1;else hi=mid;}
            at[b]=end[b]=lo;
            if(lo==limit||key_at(lo)!=target)continue;
            // Upper-bound ties are usually nearby; bracket before binary search.
            const auto begin=lo,remaining=limit-begin;std::size_t span=1;
            while(span<remaining&&key_at(begin+span)==target)span*=2;
            lo=begin+(span==1?1:span/2+1);hi=begin+std::min(span+1,remaining);
            while(lo<hi){const auto mid=lo+(hi-lo)/2;if(key_at(mid)<=target)lo=mid+1;else hi=mid;}
            end[b]=lo;
            heads[b]=posting(b,at[b]);
        }
        for(;;){
            const int minimum=std::min(std::min(heads[0],heads[1]),
                std::min(std::min(heads[2],heads[3]),std::min(heads[4],heads[5])));
            if(minimum==exhausted)return;
            ++candidates;const auto& part=parts_[minimum];const auto x=query^part.value;
            const auto different=((x|(x>>1)|(x>>2))|(query>>2))&0x9249249249ULL;
            auto remaining=different&(different-1);remaining&=remaining-1;
            if(remaining==0)for(int i=part.begin;i<part.end;++i)if(!visitor(members_[i]))return;
            for(int b=0;b<6;++b)if(heads[b]==minimum){
                ++at[b];heads[b]=at[b]<end[b]?posting(b,at[b]):exhausted;
            }
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
        if(members_.size()>=1024)radix_ids(members_,value,3*length_);
        else std::sort(members_.begin(), members_.end(), [&](int a, int b) {
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
        paired_=parts_.size()<=65536;
        compact_pairs_=paired_&&parts_.size()>=4096;
        for(auto& directory:pair_offsets_)std::vector<std::uint32_t>().swap(directory);
        if(paired_){
            // Six 16-bit ID streams occupy three existing 32-bit vectors.
            // Pair keys are reconstructed from parts_, so no key arrays persist.
            static_assert(sizeof(int)==sizeof(std::uint32_t));
            for(auto& offsets:offsets_)std::vector<std::uint32_t>().swap(offsets);
            if(compact_pairs_){
                // Six domains each <=P/8 cost <=3P+24 bytes. Every prefix
                // includes the entire first block (at most8 compact bits).
                int prefix_bits=0;while((std::size_t(1)<<(prefix_bits+1))<=parts_.size()/8)++prefix_bits;
                for(int b=0;b<6;++b){
                    const int bits=pair_bits(b),prefix=std::min(bits,prefix_bits);
                    pair_shift_[b]=bits-prefix;
                    const std::size_t domain=std::size_t(1)<<prefix;
                    auto& directory=pair_offsets_[b];directory.assign(domain+1,0);
                    for(const auto& part:parts_){const auto key=compact_pair_key(part.value,b);
                        if(key<(std::uint32_t(1)<<bits))++directory[(key>>pair_shift_[b])+1];}
                    std::partial_sum(directory.begin(),directory.end(),directory.begin());
                }
            }else{
                // Small indexes keep the shared first-block directories.
                for(int b=0;b<3;++b){
                    const std::size_t domain=std::size_t(1)<<(3*(bounds_[b+1]-bounds_[b]));
                    auto& directory=offsets_[b];directory.assign(domain+1,0);
                    for(const auto& part:parts_)++directory[block(part.value,b)+1];
                    std::partial_sum(directory.begin(),directory.end(),directory.begin());
                }
            }
            std::vector<int>().swap(blocks_[3]);
            for(int b=0;b<3;++b)blocks_[b].resize(parts_.size());
            std::vector<std::uint16_t> ids(parts_.size());
            for(int b=0;b<6;++b){
                std::iota(ids.begin(),ids.end(),std::uint16_t(0));
                const auto key=[&](int id){return compact_pairs_?compact_pair_key(parts_[id].value,b):pair_key(parts_[id].value,b);};
                if(ids.size()>=1024)radix_ids(ids,key,24);
                else std::sort(ids.begin(),ids.end(),[&](int a,int c){
                    const auto x=key(a),y=key(c);
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
        if(paired_){visit_pairs(query,visitor,lookups,candidates);return;}
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
        for (const auto& offsets : pair_offsets_) n += offsets.capacity()*sizeof(std::uint32_t);
        return n;
    }
};
}
