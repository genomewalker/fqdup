#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
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
                at[first]=std::lower_bound(entries.begin()+at[first]+1,entries.begin()+end[first],second)-entries.begin();
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
