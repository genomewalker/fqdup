#pragma once
#include "packed_kernels.hpp"
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
    struct Part { int begin, end; };
    std::vector<int> members_;
    // Contiguous postings share adjacent endpoints; retain each boundary once.
    class PartRanges {
        std::vector<int> offsets_;
    public:
        void clear(){offsets_.clear();}
        void reserve(std::size_t count){offsets_.reserve(count+1);}
        std::size_t size()const{return offsets_.empty()?0:offsets_.size()-1;}
        void push_back(Part p){if(offsets_.empty())offsets_.push_back(p.begin);else if(offsets_.back()!=p.begin)throw std::logic_error("noncontiguous family postings");offsets_.push_back(p.end);}
        Part operator[](std::size_t i)const{return {offsets_[i],offsets_[i+1]};}
        std::size_t owned_bytes()const{return offsets_.capacity()*sizeof(int);}
    };
    PartRanges parts_;

    std::vector<std::uint64_t> values_;
    std::array<std::vector<int>, 4> blocks_;
    std::array<std::vector<std::uint32_t>, 4> offsets_;
    std::array<std::vector<std::uint32_t>,6> pair_offsets_;
    std::array<int,6> pair_shift_{};
    std::array<int, 5> bounds_{};
    int length_ = 0;
    bool paired_=true;
    bool compact_pairs_=false;
    // Full ACGT domain: one word covers the low three bases. Ranks map
    // matching bits back to the original sorted part IDs.
    std::vector<std::uint64_t> dense_words_;
    std::vector<std::uint32_t> dense_ranks_;
    std::vector<int> dense_ids_;
    bool sparse_=false;
    std::vector<std::uint64_t> sparse_words_,sparse_presence_;
    std::vector<std::uint32_t> sparse_ranks_,sparse_offsets_;
    void clear_sparse(){std::vector<std::uint64_t>().swap(sparse_words_);std::vector<std::uint64_t>().swap(sparse_presence_);std::vector<std::uint32_t>().swap(sparse_ranks_);std::vector<std::uint32_t>().swap(sparse_offsets_);}
    std::uint32_t sparse_slot(std::uint32_t row)const{
        const std::uint64_t present=sparse_presence_[row>>6],bit=1ULL<<(row&63);
        // Slot zero is an immutable empty word. Masking avoids an unpredictable
        // occupied-row branch before the usual rare-candidate branch.
        const auto slot=sparse_offsets_[row>>6]+__builtin_popcountll(present&(bit-1))+1;
        return slot&(0U-std::uint32_t((present&bit)!=0));
    }


    struct NNode { std::array<int,5> child{{-1,-1,-1,-1,-1}}; int id=-1; };
    std::vector<NNode> n_nodes_;
    static const auto& suffix_masks() {
        static const auto table=[] {
            std::array<std::array<std::uint64_t,3>,512> result{};
            for(int q=0;q<512;++q)for(int s=0;s<64;++s){
                int distance=0;
                for(int p=0;p<3;++p){const int c=(q>>(3*p))&7;
                    distance+=c>=4||c!=((s>>(2*p))&3);}
                for(int budget=distance;budget<=2;++budget)result[q][budget]|=1ULL<<s;
            }
            return result;
        }();
        return table;
    }
    struct PrefixDeltas {
        std::array<std::uint32_t,33> one{};
        std::array<std::uint32_t,495> two{};
        int singles=0,doubles=0;
    };
    static const auto& prefix_deltas(){
        static const auto table=[] {
            std::array<PrefixDeltas,2> result{};
            for(int t=0;t<2;++t){auto& out=result[t];
                for(int a=0;a<10+t;++a)for(std::uint32_t x=1;x<=3;++x){
                    out.one[out.singles++]=x<<(2*a);
                    for(int b=0;b<a;++b)for(std::uint32_t y=1;y<=3;++y)
                        out.two[out.doubles++]=(x<<(2*a))|(y<<(2*b));
                }
            }
            return result;
        }();
        return table;
    }
    void build_dense() {
        for(auto& v:blocks_)std::vector<int>().swap(v);
        for(auto& v:offsets_)std::vector<std::uint32_t>().swap(v);
        const std::size_t rows=std::size_t(1)<<(2*(length_-3));
        // Upper-bound row occupancy by the number of distinct parts. Quarter
        // occupancy limits extra rank-directory work to genuinely sparse data.
        sparse_=parts_.size()<=rows/4;clear_sparse();
        if(sparse_){std::vector<std::uint64_t>().swap(dense_words_);std::vector<std::uint32_t>().swap(dense_ranks_);
            sparse_presence_.assign(rows/64,0);sparse_offsets_.assign(rows/64,0);sparse_words_.reserve(parts_.size()+1);sparse_ranks_.reserve(parts_.size()+1);sparse_words_.push_back(0);sparse_ranks_.push_back(0);
        }else{dense_words_.assign(rows,0);dense_ranks_.assign(rows,0);}
        std::uint32_t previous_row=UINT32_MAX,clean_rank=0;

        dense_ids_.clear();n_nodes_.clear();
        for(std::size_t id=0;id<parts_.size();++id){
            const auto value=values_[id];
            if(value&0x24924924924ULL){
                if(n_nodes_.empty()){
                    // Materialize earlier identity IDs only when the first N
                    // actually requires a rank-to-part translation.
                    for(std::size_t previous=0;previous<id;++previous)dense_ids_.push_back(previous);
                    n_nodes_.emplace_back();
                }

                int node=0;
                for(int p=length_-1;p>=0;--p){
                    const int c=(value>>(3*p))&7;
                    int next=n_nodes_[node].child[c];
                    if(next<0){next=n_nodes_.size();n_nodes_[node].child[c]=next;n_nodes_.emplace_back();}
                    node=next;
                }
                n_nodes_[node].id=id;
            }else{
                std::uint32_t word=0;
                for(int p=0;p<length_;++p)word|=((value>>(3*p))&3)<<(2*p);
                if(sparse_){const auto row=word>>6;
                    if(row!=previous_row){sparse_presence_[row>>6]|=1ULL<<(row&63);sparse_words_.push_back(0);sparse_ranks_.push_back(clean_rank);previous_row=row;}
                    sparse_words_.back()|=1ULL<<(word&63);++clean_rank;
                }else dense_words_[word>>6]|=1ULL<<(word&63);

                if(!n_nodes_.empty())dense_ids_.push_back(id);
            }
        }
        // With no N-bearing parts, dense rank already equals the part ID.
        if(dense_ids_.empty())std::vector<int>().swap(dense_ids_);
        std::uint32_t rank=0;
        if(sparse_)for(std::size_t i=0;i<sparse_presence_.size();++i){sparse_offsets_[i]=rank;rank+=__builtin_popcountll(sparse_presence_[i]);}
        else for(std::size_t i=0;i<rows;++i){dense_ranks_[i]=rank;rank+=__builtin_popcountll(dense_words_[i]);}

        (void)suffix_masks();(void)prefix_deltas();
    }
    template<bool Ordered=true,bool Sparse=false,class Visitor>
    void visit_dense(std::uint64_t query,Visitor visitor,std::uint64_t& lookups,std::uint64_t& candidates) const {
        // Radius two over fourteen five-letter positions has at most
        // 1+14*4+C(14,2)*16=1513 values, independent of library size.
        std::array<int,1513> hits;std::size_t count=0;
        int ns=0;for(int p=0;p<length_;++p)ns+=((query>>(3*p))&7)==4;
        if(ns>2)return;
        std::uint32_t prefix=0;
        std::array<int,11> normal{};int normals=0;
        std::array<int,2> unknown{};int unknowns=0;
        for(int p=3;p<length_;++p){
            const int c=(query>>(3*p))&7,shift=2*(p-3);
            if(c==4)unknown[unknowns++]=shift;
            else{prefix|=std::uint32_t(c)<<shift;normal[normals++]=shift;}
        }
        const auto& masks=suffix_masks()[query&511];
        const auto collect=[&](std::uint32_t row,int budget){
            std::uint64_t present;
            if constexpr(Sparse)present=sparse_words_[sparse_slot(row)];
            else present=dense_words_[row];
            auto matches=present&masks[budget];

            while(matches){
                hits[count++]=(row<<6)|__builtin_ctzll(matches);
                matches&=matches-1;}
        };
        const int budget=2-unknowns,expansions=1<<(2*unknowns);
        lookups+=expansions*(1+(budget>=1?3*normals:0)+(budget>=2?9*normals*(normals-1)/2:0));
        if(unknowns==0){
            const auto& deltas=prefix_deltas()[length_-13];
            collect(prefix,2);
            for(int i=0;i<deltas.singles;++i)collect(prefix^deltas.one[i],1);
            for(int i=0;i<deltas.doubles;++i)collect(prefix^deltas.two[i],0);
        }else for(int expansion=0;expansion<expansions;++expansion){
            auto base=prefix;
            for(int i=0;i<unknowns;++i)base|=((expansion>>(2*i))&3)<<unknown[i];
            collect(base,budget);
            if(budget<1)continue;
            for(int a=0;a<normals;++a)for(std::uint32_t x=1;x<=3;++x){
                const auto first=base^(x<<normal[a]);collect(first,budget-1);
                if(budget<2)continue;
                for(int b=0;b<a;++b)for(std::uint32_t y=1;y<=3;++y)collect(first^(y<<normal[b]),0);
            }
        }
        // Compact ACGT codes have the same order as packed part values.
        // Separate dependent table loads from the enumeration hot loop.
        if constexpr(Ordered)std::sort(hits.begin(),hits.begin()+count);
        for(std::size_t i=0;i<count;++i){const auto code=std::uint32_t(hits[i]);
            const auto row=code>>6;const auto bit=std::uint64_t(1)<<(code&63);
            if constexpr(Sparse){const auto slot=sparse_slot(row);hits[i]=sparse_ranks_[slot]+__builtin_popcountll(sparse_words_[slot]&(bit-1));}
            else hits[i]=dense_ranks_[row]+__builtin_popcountll(dense_words_[row]&(bit-1));}

        if(!dense_ids_.empty())for(std::size_t i=0;i<count;++i)hits[i]=dense_ids_[hits[i]];
        if(!n_nodes_.empty()){
            const auto walk=[&](auto&& self,int node,int p,int remaining)->void {
                if(p<0){hits[count++]=n_nodes_[node].id;return;}
                const int q=(query>>(3*p))&7;
                for(int c=0;c<5;++c){const int next=n_nodes_[node].child[c];
                    const int cost=q==4||q!=c;
                    if(next>=0&&cost<=remaining){++lookups;self(self,next,p-1,remaining-cost);}}
            };
            walk(walk,0,length_-1,2);
            if constexpr(Ordered)std::sort(hits.begin(),hits.begin()+count);
        }
        for(std::size_t i=0;i<count;++i){++candidates;const auto& part=parts_[hits[i]];
            for(int j=part.begin;j<part.end;++j)if(!visitor(members_[j]))return;}
    }
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
        for(int b=0;b<4;++b){query_blocks[b]=block(query,b);n_free[b]=(query_blocks[b]&0x924U)==0;}
        for(int b=0;b<6;++b){
            if(!n_free[first_block[b]]||!n_free[second_block[b]])continue;
            ++lookups;const auto target=query_blocks[second_block[b]];
            const auto key_at=[&](std::size_t pos){return block(values_[posting(b,pos)],second_block[b]);};
            const auto& directory=compact_pairs_?pair_offsets_[b]:offsets_[first_block[b]];
            const auto prefix=compact_pairs_?(compact_pair_key(query,b)>>pair_shift_[b]):query_blocks[first_block[b]];
            std::size_t lo=directory.empty()?0:directory[prefix];
            const std::size_t limit=directory.empty()?parts_.size():directory[prefix+1];
            // A complete compact key already gives the exact equal range.
            if(compact_pairs_&&pair_shift_[b]==0){
                at[b]=lo;end[b]=limit;if(lo<limit)heads[b]=posting(b,lo);
                continue;
            }
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
            ++candidates;const auto x=query^values_[minimum];
            const auto different=((x|(x>>1)|(x>>2))|(query>>2))&0x9249249249ULL;
            auto remaining=different&(different-1);remaining&=remaining-1;
            if(remaining==0){const auto& part=parts_[minimum];
                for(int i=part.begin;i<part.end;++i)if(!visitor(members_[i]))return;}
            for(int b=0;b<6;++b)if(heads[b]==minimum){
                ++at[b];heads[b]=at[b]<end[b]?posting(b,at[b]):exhausted;
            }
        }
    }
public:
    // One entry per representative, sorted by this part and then original ID.
    const std::vector<int>& ordered_queries()const{return members_;}
    template<class Key>
    static std::uint64_t pack(const Key& key, int begin, int end) {
        std::uint64_t value = 0;
        for (int p = begin; p < end; ++p) value |= std::uint64_t(code(key[p])) << (3 * (p-begin));
        return value;
    }
    static std::uint64_t pack(const FamilyKey& key,int begin,int end){
        if(!key.empty()&&key.supported()){
            if(begin==0&&end==13)return TwoBitFamilyKey(key).part(0);
            if(begin==13&&end==26)return TwoBitFamilyKey(key).part(1);
            if(begin==26&&end==40)return TwoBitFamilyKey(key).part(2);
        }
        return pack<FamilyKey>(key,begin,end);
    }
    template<class Key>
    void build(const std::vector<Key>& keys, const std::vector<int>& representatives, int begin, int end, bool force_dense=false) {
        length_ = end - begin;
        bounds_ = {0, 4, length_ == 13 ? 7 : 8, length_ == 13 ? 10 : 11, length_};
        members_.resize(representatives.size());
        std::iota(members_.begin(), members_.end(), 0);
        // Pack once instead of re-reading 13/14 bases on every sort comparison.
        std::vector<std::uint64_t> packed;packed.reserve(representatives.size());
        for(int i:representatives)packed.push_back(pack(keys[i],begin,end));
        const auto value = [&](int i) { return packed[i]; };
        parts_.clear();values_.clear();
        if(members_.size()>=1024)radix_ids(members_,value,3*length_);
        else std::sort(members_.begin(), members_.end(), [&](int a, int b) {
            const auto x = value(a), y = value(b); return x < y || (x == y && a < b);
        });
        std::size_t count = 0;
        for (std::size_t i = 0; i < members_.size(); ++i) if (!i || value(members_[i]) != value(members_[i-1])) ++count;
        parts_.reserve(count);values_.reserve(count);
        for (std::size_t i = 0; i < members_.size();) {
            const auto v = value(members_[i]); const int first = i++;
            while (i < members_.size() && value(members_[i]) == v) ++i;
            parts_.push_back({first, int(i)});values_.push_back(v);
        }
        std::vector<std::uint64_t>().swap(packed);
        paired_=!force_dense&&parts_.size()<=65536;
        compact_pairs_=paired_&&parts_.size()>=4096;
        for(auto& directory:pair_offsets_)std::vector<std::uint32_t>().swap(directory);
        if(!paired_){
            build_dense();
            // Dense queries use bitmap ranks and the N trie, never values_.
            // Keep values only for the small paired-posting representation.
            std::vector<std::uint64_t>().swap(values_);
            return;
        }
        sparse_=false;clear_sparse();
        std::vector<std::uint64_t>().swap(dense_words_);
        std::vector<std::uint32_t>().swap(dense_ranks_);
        std::vector<int>().swap(dense_ids_);std::vector<NNode>().swap(n_nodes_);
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
                    for(const auto value:values_){const auto key=compact_pair_key(value,b);
                        if(key<(std::uint32_t(1)<<bits))++directory[(key>>pair_shift_[b])+1];}
                    std::partial_sum(directory.begin(),directory.end(),directory.begin());
                }
            }else{
                // Small indexes keep the shared first-block directories.
                for(int b=0;b<3;++b){
                    const std::size_t domain=std::size_t(1)<<(3*(bounds_[b+1]-bounds_[b]));
                    auto& directory=offsets_[b];directory.assign(domain+1,0);
                    for(const auto value:values_)++directory[block(value,b)+1];
                    std::partial_sum(directory.begin(),directory.end(),directory.begin());
                }
            }
            std::vector<int>().swap(blocks_[3]);
            for(int b=0;b<3;++b)blocks_[b].resize(parts_.size());
            std::vector<std::uint16_t> ids(parts_.size());
            for(int b=0;b<6;++b){
                std::iota(ids.begin(),ids.end(),std::uint16_t(0));
                const auto key=[&](int id){return compact_pairs_?compact_pair_key(values_[id],b):pair_key(values_[id],b);};
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
    }
    // Both representations return sorted part IDs and sorted member IDs.
    template<class Visitor>
    void visit(std::uint64_t query, Visitor visitor, std::uint64_t& lookups, std::uint64_t& candidates) const {
        if(paired_){visit_pairs(query,visitor,lookups,candidates);return;}
        if(sparse_)visit_dense<true,true>(query,visitor,lookups,candidates);else visit_dense<true,false>(query,visitor,lookups,candidates);
    }
    // Identical unique candidates in unspecified order, for commutative
    // component unions. The ordered API above retains its existing contract.
    template<class Visitor>
    void visit_unordered(std::uint64_t query,Visitor visitor,std::uint64_t& lookups,std::uint64_t& candidates)const {
        if(paired_){visit_pairs(query,visitor,lookups,candidates);return;}
        if(sparse_)visit_dense<false,true>(query,visitor,lookups,candidates);else visit_dense<false,false>(query,visitor,lookups,candidates);
    }
    // Return original read IDs while preserving posting and query order.
    template<class Key>
    void build_source(const std::vector<Key>& keys,const std::vector<int>& representatives,int begin,int end,bool force_dense=false){
        build(keys,representatives,begin,end,force_dense);
        for(auto& id:members_)id=representatives[id];
    }
    std::size_t cells() const { return parts_.size(); }
    std::size_t max_bucket() const {
        std::size_t n = 0; for (std::size_t i=0;i<parts_.size();++i){const auto p=parts_[i];n=std::max(n,std::size_t(p.end-p.begin));} return n;
    }
    std::size_t bytes() const {
        std::size_t n = members_.capacity()*sizeof(int) + parts_.owned_bytes() + values_.capacity()*sizeof(std::uint64_t);
        for (const auto& b : blocks_) n += b.capacity()*sizeof(int);
        for (const auto& offsets : offsets_) n += offsets.capacity()*sizeof(std::uint32_t);
        for (const auto& offsets : pair_offsets_) n += offsets.capacity()*sizeof(std::uint32_t);
        n+=dense_words_.capacity()*sizeof(std::uint64_t)+dense_ranks_.capacity()*sizeof(std::uint32_t);
        n+=(sparse_words_.capacity()+sparse_presence_.capacity())*sizeof(std::uint64_t)+(sparse_ranks_.capacity()+sparse_offsets_.capacity())*sizeof(std::uint32_t);
        n+=dense_ids_.capacity()*sizeof(int)+n_nodes_.capacity()*sizeof(NNode);
        return n;
    }
};
}
