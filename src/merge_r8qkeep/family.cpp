// Native family (duplicate) correction pre-pass: C++ port of famcorrect_r4.mjs (FAMILY_PREREG.md Revision 4).
// Must be byte-identical to the JS reference on R1.fastq, R2.fastq, family-edits.tsv and families.tsv.
// Usage: famcorrect R1.fastq R2.fastq OUTDIR     (plain FASTQ; 4-line records)
#include "family.hpp"
#include "atomic_family_union.hpp"
#include <atomic>
#include "compact_family.hpp"
#include "family_index.hpp"
#include "family_key.hpp"
#include "packed_kernels.hpp"
#include "ordered_chunks.hpp"
#include "memory_status.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <array>
#include <cstdint>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <numeric>
#include <type_traits>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <sys/stat.h>

namespace {
constexpr int KEYLEN = 20, MAXMIS = 8, MAXN = 2, BLOCKS = 9, MAXCOMP = 12, MAXDIAM = 16, ADSEED = 10, ADZONE = 2 * KEYLEN;
constexpr double DUST_MAX = 0.10, MAXDIS = 0.5;
constexpr bool g_rev8 = true;
// --phred-unresolved (arm r8qkeep): an unresolved family column sets each copy to -10log10 P(copy wrong), with
// P(allele a true) proportional to prod_{i in a}(1-e_i) * prod_{j not in a} e_j/3 from the copies own qualities, clamped to Q2..Q40 (the data floor; Q0 makes risk 1, rejected downstream), instead of Q2
constexpr bool g_phred = true;
const std::array<std::string, 2> ADAPTERS{"AGATCGGAAGAGCACACGTCTGAACTCCAG", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAG"};

using Read = fqdup::r8qkeep::FamilyRead;
using Pair = fqdup::r8qkeep::FamilyPair;

std::string base_id(const std::string& h) {
    std::string x = h.substr(1, h.find_first_of(" \t") == std::string::npos ? std::string::npos : h.find_first_of(" \t") - 1);
    if (x.size() > 2 && x[x.size() - 2] == '/' && (x.back() == '1' || x.back() == '2')) x.resize(x.size() - 2);
    return x;
}
std::vector<Read> read_fastq(const std::string& path) {
    std::ifstream f(path); if (!f) throw std::runtime_error("cannot open " + path);
    std::vector<Read> out; std::string h, s, p, q;
    while (std::getline(f, h)) { if (!std::getline(f, s) || !std::getline(f, p) || !std::getline(f, q)) throw std::runtime_error("truncated " + path); out.push_back({h, base_id(h), s, q}); }
    return out;
}
// native full_complexity::overlap_is_low_complexity (window 64)
bool low_complexity(std::string_view seq) {
    const int len = int(seq.size()); if (len < 5) return false;
    auto b2 = [](char c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; } return -1; };
    const int nt = len - 2, wt = std::min(62, nt); if (wt < 2) return false;
    std::array<int,62> tri{};
    const auto triplet=[&](int i){int a=b2(seq[i]),b=b2(seq[i+1]),c=b2(seq[i+2]);return a<0||b<0||c<0?-1:a*16+b*4+c;};
    int cnt[64] = {}; long valid = 0, score = 0;
    auto add = [&](int t) { if (t >= 0) { score += cnt[t]; ++cnt[t]; ++valid; } };
    auto del = [&](int t) { if (t >= 0) { --cnt[t]; score -= cnt[t]; --valid; } };
    auto exceeds = [&]() { return valid >= 2 && double(score) / (double(valid) * (valid - 1) / 2) > DUST_MAX; };
    for (int i = 0; i < wt; ++i) add(tri[i]=triplet(i)); if (exceeds()) return true;
    int cursor=0;
    for (int i = wt; i < nt; ++i) {auto& slot=tri[cursor];del(slot);add(slot=triplet(i));if(exceeds())return true;if(++cursor==wt)cursor=0;}
    return false;
}
template<class Key>
int mismatches(const Key& a, const Key& b) { int m = 0; for (size_t i = 0; i < a.size(); ++i) if (a[i] == 'N' || b[i] == 'N' || a[i] != b[i]) ++m; return m; }

// FamilyKey stores the exact N masks used by the scalar character comparison.
// Group diameter requires the full distance, including an N/N mismatch.
int mismatches(const fqdup::r8qkeep::FamilyKey& a,const fqdup::r8qkeep::FamilyKey& b){
    if(a.empty())return 0;
    using fqdup::r8qkeep::TwoBitFamilyKey;
    return __builtin_popcountll(TwoBitFamilyKey::differences(a.lo,b.lo))
        +__builtin_popcountll(TwoBitFamilyKey::differences(a.hi,b.hi));
}

using Group=fqdup::r8qkeep::FamilyGroup;
// Masks index the original group order; views never own or reorder read IDs.
struct MemberMask {
    const std::vector<int>* group;
    std::uint16_t bits;
    struct Iterator {
        const std::vector<int>* group;
        std::uint16_t bits;
        int operator*() const {return (*group)[__builtin_ctz(unsigned(bits))];}
        Iterator& operator++(){bits&=bits-1;return *this;}
        bool operator!=(const Iterator& other) const {return bits!=other.bits;}
    };
    Iterator begin() const {return {group,bits};}
    Iterator end() const {return {group,0};}
    std::size_t size() const {return __builtin_popcount(unsigned(bits));}
    int operator[](std::size_t index) const {
        auto remaining=bits;while(index--)remaining&=remaining-1;
        return (*group)[__builtin_ctz(unsigned(remaining))];
    }
};
static_assert(MAXCOMP<=255,"family vote count fits one byte");
struct ColumnVotes {
    const std::vector<int>* group=nullptr;
    std::array<char,MAXCOMP> bases{};
    std::array<std::uint16_t,MAXCOMP> masks{};
    std::uint8_t count=0;
    struct Vote {char first;MemberMask second;};
    Vote operator[](std::size_t index) const {return {bases[index],{group,masks[index]}};}
    std::size_t size() const {return count;}
    struct Iterator {
        const ColumnVotes* votes;
        std::size_t index;
        Vote operator*() const {return (*votes)[index];}
        Iterator& operator++(){++index;return *this;}
        bool operator!=(const Iterator& other) const {return index!=other.index;}
    };
    Iterator begin() const {return {this,0};}
    Iterator end() const {return {this,count};}
};
struct Column {int n=0,top=-1;ColumnVotes votes;}; // first-seen base order

// Forty three-bit bases fit in two sixty-bit words. This avoids cache-line
// straddles in candidate records while retaining N as an explicit code.
struct PackedFamilyKey {
    std::uint64_t lo,hi;
    template<class Key>
    explicit PackedFamilyKey(const Key& key):
        lo(fqdup::r8qkeep::FamilyPartIndex::pack(key,0,20)),
        hi(fqdup::r8qkeep::FamilyPartIndex::pack(key,20,40)){}
    std::uint64_t part(int p) const {
        if(p==0)return lo&((1ULL<<39)-1);
        if(p==1)return (lo>>39)|((hi&((1ULL<<18)-1))<<21);
        return hi>>18;
    }
    static std::uint64_t differences(std::uint64_t a,std::uint64_t b){
        const auto x=a^b;
        return ((x|(x>>1)|(x>>2))|(a>>2))&0x9249249249249249ULL;
    }
};
static_assert(sizeof(PackedFamilyKey)==16);

// Integer slots reference immutable full keys; hash collisions never imply equality.
template<class Key>
class ExactFamilyKeys {
    const std::vector<Key>& keys_;
    std::vector<int> slots_;
public:
    ExactFamilyKeys(const std::vector<Key>& keys,std::size_t eligible):keys_(keys){
        std::size_t capacity=1;
        if(eligible>slots_.max_size()/2)throw std::length_error("family key table size");
        while(capacity<2*eligible)capacity*=2;
        slots_.assign(capacity,-1);
    }
    int insert(int id,std::size_t hash){
        const auto mask=slots_.size()-1;auto pos=hash&mask;
        while(slots_[pos]>=0){const int old=slots_[pos];
            if(keys_[old]==keys_[id])return old;
            pos=(pos+1)&mask;
        }
        slots_[pos]=id;return id;
    }
};

// Radius8 over40 bases implies radius2 in at least one13/13/14 block.
// Pigeonhole candidates with the frozen N/N mismatch rule.
// Union roots are not output: groups/members are emitted by ascending read id.
bool indexable_key(const fqdup::r8qkeep::FamilyKey& key){return key.supported()&&key.ns()<=MAXN;}
bool indexable_key(const std::string& key){
    const auto code=[](char c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;case 'N':return 4;default:return -1;}};
    if(key.size()!=40||std::count(key.begin(),key.end(),'N')>MAXN)return false;
    for(char c:key)if(code(c)<0)return false;return true;
}
std::size_t key_hash(const std::string& key){return std::hash<std::string_view>{}(key);}
std::size_t key_hash(const fqdup::r8qkeep::FamilyKey& key){return key.hash();}
std::size_t key_bytes(const std::vector<std::string>& keys){return fqdup::r8qkeep::strings_bytes(keys);}
std::size_t key_bytes(const std::vector<fqdup::r8qkeep::FamilyKey>& keys){return fqdup::r8qkeep::vector_bytes(keys);}
template<class Key>
bool neighbor_families_impl(const std::vector<Key>& keys,std::vector<int>& parent,int threads=1,std::size_t edge_limit=64,fqdup::r8qkeep::WorkerPool* shared_pool=nullptr){
    std::size_t eligible_keys=0;
    for(const auto& key:keys)if(!key.empty()){++eligible_keys;if(!indexable_key(key))return false;}
    const auto started=std::chrono::steady_clock::now();
    const auto root=[&](int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
    const auto unite=[&](int i,int j){i=root(i);j=root(j);if(i!=j)parent[i]=j;};
    std::vector<int> representative;
    {
        ExactFamilyKeys exact(keys,eligible_keys);
        for(std::size_t i=0;i<keys.size();++i)if(!keys[i].empty()){
            const int first=exact.insert(i,key_hash(keys[i]));
            if(first==int(i))representative.push_back(i);else unite(i,first);
        }
    } // Release duplicate lookup storage before constructing part indexes.
    std::array<fqdup::r8qkeep::FamilyPartIndex,3> index; const int bounds[]={0,13,26,40};
    const int workers=std::max(1,std::min(threads,8));
    std::unique_ptr<fqdup::r8qkeep::WorkerPool> owned_pool;
    auto* pool=workers>1?shared_pool:nullptr;
    if(workers>1&&!pool){owned_pool=std::make_unique<fqdup::r8qkeep::WorkerPool>(workers);pool=owned_pool.get();}
    const auto build_part=[&](int part){if constexpr(std::is_same_v<Key,fqdup::r8qkeep::FamilyKey>)index[part].build_source(keys,representative,bounds[part],bounds[part+1]);else index[part].build(keys,representative,bounds[part],bounds[part+1]);};
    if(pool){
        std::array<std::future<void>,3> builds;
        try{for(int part=0;part<3;++part)builds[part]=pool->submit([&,part]{build_part(part);});for(auto& build:builds)build.get();}
        catch(...){pool->stop();throw;}
    }else for(int part=0;part<3;++part)build_part(part);
    std::size_t cells=0,max_bucket=0,index_bytes=0;
    for(int part=0;part<3;++part){cells+=index[part].cells();max_bucket=std::max(max_bucket,index[part].max_bucket());index_bytes+=index[part].bytes();}
    const auto unique_keys=representative.size();
    // Compact postings already hold original read IDs. The immutable key array
    // is the verification cache, so its former ID translation can be released.
    if constexpr(std::is_same_v<Key,fqdup::r8qkeep::FamilyKey>)std::vector<int>().swap(representative);
    std::vector<fqdup::r8qkeep::TwoBitFamilyKey> packed_keys;

    if constexpr(!std::is_same_v<Key,fqdup::r8qkeep::FamilyKey>){packed_keys.reserve(representative.size());for(int i:representative)packed_keys.emplace_back(keys[i]);}
    const auto key_at=[&](int id)->const auto&{if constexpr(std::is_same_v<Key,fqdup::r8qkeep::FamilyKey>)return keys[id];else return packed_keys[id];};
    const auto source_id=[&](int id){if constexpr(std::is_same_v<Key,fqdup::r8qkeep::FamilyKey>)return id;else return representative[id];};
    if(fqdup::r8qkeep::memory_enabled()){
        fqdup::r8qkeep::MemoryUsage m;
        m.family=index_bytes+fqdup::r8qkeep::vector_bytes(representative)+fqdup::r8qkeep::vector_bytes(packed_keys)+fqdup::r8qkeep::vector_bytes(parent)+key_bytes(keys);
        m.report("family_index",keys.size(),"index_only_excludes_raw");
    }
    const auto generated=std::chrono::steady_clock::now();
    #include "parallel_family.inc"
    if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS"))
{const auto finished=std::chrono::steady_clock::now();std::cerr<<"FAMILY_INDEX pairs="<<keys.size()<<" unique_keys="<<unique_keys<<" index_bytes="<<index_bytes<<" cells="<<cells<<" max_bucket="<<max_bucket<<" lookups="<<lookups<<" candidates="<<candidates<<" verified="<<verified<<" workers="<<workers<<" scratch_bytes="<<scratch_bytes<<" overflow_queries="<<overflows<<" generation_seconds="<<std::chrono::duration<double>(generated-started).count()<<" verification_seconds="<<std::chrono::duration<double>(finished-generated).count()<<'\n';}
    return true;
}

bool neighbor_families(const std::vector<std::string>& keys,std::vector<int>& parent,int threads=1,std::size_t edge_limit=64,fqdup::r8qkeep::WorkerPool* pool=nullptr){return neighbor_families_impl(keys,parent,threads,edge_limit,pool);}
template<class Key=std::string>
std::vector<Key> keys_of(const std::vector<Pair>& P,fqdup::r8qkeep::WorkerPool* pool=nullptr,int threads=1) {
    std::vector<Key> keys(P.size());
    const auto collect=[&](std::size_t begin,std::size_t count){
    for (size_t i = begin; i < begin+count; ++i) {
        const auto &a = P[i].r[0].s, &b = P[i].r[1].s; if (int(a.size()) < KEYLEN || int(b.size()) < KEYLEN) continue;
        std::array<char,2*KEYLEN> bytes;
        std::copy_n(a.data(),KEYLEN,bytes.data());std::copy_n(b.data(),KEYLEN,bytes.data()+KEYLEN);
        const std::string_view k(bytes.data(),bytes.size());
        if (std::count(k.begin(), k.end(), 'N') > MAXN) continue;
        if (!g_rev8) {
            if (low_complexity(a) || low_complexity(b)) continue;
            bool adapter = false; for (int s = 0; s < 2; ++s) { auto at = P[i].r[s].s.find(ADAPTERS[s].substr(0, ADSEED)); if (at != std::string::npos && int(at) < ADZONE) adapter = true; }
            if (!adapter) keys[i] = Key(k);
        } else { // Revision 8: complexity on the biological prefix (before adapter seed or terminal G>=10 run); adapter must start at >= KEYLEN
            bool ok = true;
            for (int s = 0; s < 2 && ok; ++s) { const auto& r = P[i].r[s].s; auto at = r.find(ADAPTERS[s].substr(0, ADSEED)); if (at != std::string::npos && int(at) < KEYLEN) { ok = false; break; }
                size_t e = r.size(); if (at != std::string::npos) e = std::min(e, at); size_t g = r.size(); while (g > 0 && r[g - 1] == 'G') --g; if (r.size() - g >= 10) e = std::min(e, g);
                if (low_complexity(std::string_view(r).substr(0, e))) ok = false; }
            if (ok) keys[i] = Key(k);
        }
    }
    };
    if(pool)fqdup::r8qkeep::ordered_index_chunks(*pool,P.size(),threads,0,
        [&](std::size_t begin,std::size_t count){collect(begin,count);return 0;},
        [](int){},"family_keys");
    else collect(0,P.size());
    return keys;
}
template<class Key>
std::vector<Group> families(const std::vector<Pair>& P, std::vector<Key>& keys, int threads,fqdup::r8qkeep::WorkerPool* shared_pool,bool supplied=false) {
    auto last=std::chrono::steady_clock::now();
    const auto phase=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS"))std::cerr<<"FAMILY_PHASE "<<name<<" seconds="<<std::chrono::duration<double>(now-last).count()<<'\n';last=now;};
    const int workers=std::max(1,std::min(threads,8));
    std::unique_ptr<fqdup::r8qkeep::WorkerPool> owned_pool;
    auto* pool=workers>1?shared_pool:nullptr;
    if(workers>1&&!pool){owned_pool=std::make_unique<fqdup::r8qkeep::WorkerPool>(workers);pool=owned_pool.get();}
    if(!supplied)keys=keys_of<Key>(P,pool,threads);
    const int n=int(keys.size());phase("keys");
    if constexpr(std::is_same_v<Key,fqdup::r8qkeep::FamilyKey>){
        // Unusual alphabets retain the complete frozen string fallback.
        for(const auto& key:keys)if(!key.supported()){
            if(supplied)throw std::logic_error("unsupported precomputed family key");
            std::vector<std::string> original;return families(P,original,threads,pool);
        }
    }
    if(fqdup::r8qkeep::memory_enabled()){
        fqdup::r8qkeep::MemoryUsage m;m.pairs(P);m.family=key_bytes(keys);m.report("family_keys",n);
    }
    std::vector<int> parent(n); std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    if(!neighbor_families_impl(keys,parent,threads,64,pool)){
    std::array<int, BLOCKS + 1> bounds{}; for (int b = 0; b <= BLOCKS; ++b) bounds[b] = int(std::lround(double(b) * 2 * KEYLEN / BLOCKS));
    std::unordered_map<std::string, std::vector<int>> buckets;
    for (int i = 0; i < n; ++i) { if (keys[i].empty()) continue; for (int b = 0; b < BLOCKS; ++b) { auto s = keys[i].substr(bounds[b], bounds[b + 1] - bounds[b]); if (s.find('N') != std::string::npos) continue; buckets[std::to_string(b) + ":" + s].push_back(i); } }
    // Reproduce the original unordered-bucket traversal and first pair encounter
    // without retaining an O(candidate-pairs) tree. A pair is new iff it shares
    // no earlier bucket; ranks require only BLOCKS integers per input pair.
    std::vector<std::array<int,BLOCKS>> rank(n);for(auto& r:rank)r.fill(-1);
    int next_rank=0;
    for(const auto& [k,ids]:buckets){const int block=k[0]-'0';for(int i:ids)rank[i][block]=next_rank;++next_rank;}
    int current_rank=0;std::uint64_t candidates=0,unique=0;
    for (auto& [k, ids] : buckets) {
        for (size_t x = 0; x < ids.size(); ++x) for (size_t y = x + 1; y < ids.size(); ++y) {
            int i=ids[x],j=ids[y];++candidates;bool earlier=false;
            for(int block=0;block<BLOCKS;++block)if(rank[i][block]>=0&&rank[i][block]<current_rank&&rank[i][block]==rank[j][block]){earlier=true;break;}
            if(earlier)continue;++unique;
            if(mismatches(keys[i],keys[j])<=MAXMIS){int a=find(i),b=find(j);if(a!=b)parent[a]=b;}
        }
        ++current_rank;
    }
    if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS"))std::cerr<<"FAMILY_CANDIDATES pairs="<<n<<" buckets="<<buckets.size()<<" candidates="<<candidates<<" unique="<<unique<<" rank_bytes="<<rank.capacity()*sizeof(rank[0])<<'\n';
    }
    phase("graph");
    // Count roots first: singleton groups never enter the correction pre-pass.
    // Negative slots then encode first-seen group indices, preserving order.
    std::vector<int> slot(n,0);std::vector<Group> groups;
    for(int i=0;i<n;++i)if(!keys[i].empty())++slot[find(i)];
    for(int i=0;i<n;++i)if(!keys[i].empty()){
        auto& at=slot[find(i)];if(at==1)continue;
        if(at>1){const int members=at;at=-int(groups.size())-1;groups.push_back({});groups.back().m.reserve(members);}
        groups[-at-1].m.push_back(i);
    }
    std::vector<Group> out;
    for (auto& g : groups) { if (g.m.size() < 2) continue; int diam = 0;
        if (int(g.m.size()) <= MAXCOMP) for (size_t x = 0; x < g.m.size(); ++x) for (size_t y = x + 1; y < g.m.size(); ++y) diam = std::max(diam, mismatches(keys[g.m[x]], keys[g.m[y]]));
        g.oversize = int(g.m.size()) > MAXCOMP || diam > MAXDIAM; g.diameter = diam; out.push_back(std::move(g)); }
    phase("groups");
    return out;
}
const std::vector<Column>& columns(const std::vector<Pair>& P, const std::vector<int>& g, int side) {
    if(g.size()>MAXCOMP)throw std::logic_error("unbounded family column group");
    // Cache immutable views once; reuse scratch and initialize only live vote slots.
    std::array<std::string_view,MAXCOMP> sequences;size_t L=0;for(std::size_t member=0;member<g.size();++member){sequences[member]=P[g[member]].r[side].s;L=std::max(L,sequences[member].size());}
    static thread_local std::array<std::vector<Column>,2> scratch;
    auto& cols=scratch.at(side);cols.resize(L);
    for (size_t c = 0; c < L; ++c) { auto& col=cols[c];col.n=0;col.votes.count=0;col.votes.group=&g;
        for(std::size_t member=0;member<g.size();++member){
            const auto& s=sequences[member];
            if(c>=s.size()||s[c]=='N')continue;
            ++col.n;const char base=s[c];std::size_t slot=0;
            while(slot<col.votes.count&&col.votes.bases[slot]!=base)++slot;
            if(slot==col.votes.count){col.votes.bases[slot]=base;col.votes.masks[slot]=0;++col.votes.count;}
            col.votes.masks[slot]|=std::uint16_t(1U<<member);
        }
        int top = -1; for (size_t v = 0; v < col.votes.size(); ++v) if (top < 0 || col.votes[v].second.size() > col.votes[top].second.size()) top = int(v);
        bool tie = false; if (top >= 0) for (size_t v = 0; v < col.votes.size(); ++v) if (int(v) != top && col.votes[v].second.size() == col.votes[top].second.size()) tie = true;
        col.top = tie ? -1 : top; }
    return cols;
}
bool second_molecule(const std::vector<Pair>& P, const std::vector<int>& g,
                     std::array<const std::vector<Column>*,2>* prepared=nullptr) {
    if(g.size()<2)return false;
    const std::array<const std::vector<Column>*,2> all_columns{&columns(P,g,0),&columns(P,g,1)};
    if(prepared)*prepared=all_columns;
    std::array<std::array<int,MAXCOMP>,MAXCOMP> shared{}; long multi = 0, cover = 0;
    for (int side = 0; side < 2; ++side) for (const auto& col : *all_columns[side]) if (col.n >= 2) { ++cover; if (col.votes.size() > 1) ++multi; }
    if (cover && double(multi) / cover > MAXDIS) return true;
    for (int side = 0; side < 2; ++side) for (const auto& col : *all_columns[side]) { if (col.top < 0) continue;
        for (size_t v = 0; v < col.votes.size(); ++v) { if (int(v) == col.top || col.votes[v].second.size() < 2) continue;
            const auto mask=col.votes[v].second.bits;
            for(auto x=mask;x;x&=x-1)for(auto y=std::uint16_t(x&(x-1));y;y&=y-1)
                if(++shared[__builtin_ctz(unsigned(x))][__builtin_ctz(unsigned(y))]>=2)return true;
        }
    }
    return false;
}
int qv(const Read& r, size_t c) { return r.q[c] - 33; }
std::array<int,MAXCOMP> high_quality_counts(const std::vector<Pair>& P,const std::vector<int>& g,
                                         std::array<const std::vector<Column>*,2>& prepared) {
    std::array<int,MAXCOMP> hiq{};
    if(g.size()<2)return hiq;
    for(int side=0;side<2;++side){
        if(!prepared[side])prepared[side]=&columns(P,g,side);
        const auto& cols=*prepared[side];
        for(size_t c=0;c<cols.size();++c){const auto& col=cols[c];if(col.top<0)continue;
            int topq=-1;for(int i:col.votes[col.top].second)topq=std::max(topq,qv(P[i].r[side],c));
            if(topq<30)continue;
            for(size_t v=0;v<col.votes.size();++v)if(int(v)!=col.top)
                for(auto bits=col.votes[v].second.bits;bits;bits&=bits-1){
                    const auto member=__builtin_ctz(unsigned(bits));
                    if(qv(P[g[member]].r[side],c)>=30)++hiq[member];
                }
        }
    }
    return hiq;
}
std::array<double,256> family_quality_errors() {
    std::array<double,256> out{};
    // Preserve char signedness and the original qv/pow expression for every byte.
    for(int byte=0;byte<256;++byte){const int q=char(byte)-33;out[byte]=std::pow(10.0,-q/10.0);}
    return out;
}
// --damage-safe (DAMAGE_PREREG.md, validated): never change a C<->T / G<->A base, nor lower its Q, at read position < 10
// or within 10 positions before the first adapter-seed occurrence in that read.
constexpr int DSZONE = 10;
bool transition(char a, char b) { return (a == 'C' && b == 'T') || (a == 'T' && b == 'C') || (a == 'G' && b == 'A') || (a == 'A' && b == 'G'); }
bool is_protected(const Read& r, size_t c, int side) { if (int(c) < DSZONE) return true; auto at = r.s.find(ADAPTERS[side].substr(0, ADSEED)); return at != std::string::npos && c < at && int(c) >= int(at) - DSZONE; }
} // namespace

namespace fqdup::r8qkeep {
std::vector<FamilyKey> extract_family_keys(const std::vector<FamilyPair>& pairs){return keys_of<FamilyKey>(pairs);}
FamilyPlan discover_family_keys(std::vector<FamilyKey> keys,int threads,WorkerPool* pool){
    FamilyPlan plan;plan.pairs=keys.size();
    if(keys.size()>std::size_t(std::numeric_limits<int>::max()))throw std::length_error("family pair count");
    for(const auto& key:keys)if(!key.supported()){plan.supported=false;return plan;}
    const std::vector<FamilyPair> no_records;
    plan.groups=families(no_records,keys,threads,pool,true);
    return plan;
}
void reorder_family_plan(FamilyPlan& plan,const std::vector<std::size_t>& order){
    if(!plan.supported)return;
    if(order.size()!=plan.pairs)throw std::logic_error("family plan order size");
    std::vector<int> rank(order.size());
    for(std::size_t i=0;i<order.size();++i)rank.at(order[i])=int(i);
    for(auto& group:plan.groups){
        for(auto& member:group.m)member=rank.at(member);
        std::sort(group.m.begin(),group.m.end());
    }
    std::sort(plan.groups.begin(),plan.groups.end(),[](const auto& a,const auto& b){return a.m.front()<b.m.front();});
}
FamilyResult correct_family(std::vector<FamilyRead> R1, std::vector<FamilyRead> R2, int threads,bool debug,WorkerPool* shared_pool,FamilyPlan* plan) {
        std::unique_ptr<WorkerPool> owned_pool;
        if(!shared_pool){owned_pool=std::make_unique<WorkerPool>(std::max(1,std::min(threads,8)));shared_pool=owned_pool.get();}
        auto& pool=*shared_pool;
        auto last=std::chrono::steady_clock::now();
        const auto phase=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS"))std::cerr<<"FAMILY_PHASE "<<name<<" seconds="<<std::chrono::duration<double>(now-last).count()<<'\n';last=now;};
        constexpr bool damage_safe = false;
        if (R1.size() != R2.size()) throw std::runtime_error("pair count");
        std::vector<Pair> P; P.reserve(R1.size());
        for (size_t i = 0; i < R1.size(); ++i) { if (R1[i].id != R2[i].id) throw std::runtime_error("pair order " + R1[i].id); P.push_back({R1[i].id, {std::move(R1[i]), std::move(R2[i])}}); }
        std::vector<FamilyRead>().swap(R1);std::vector<FamilyRead>().swap(R2);
        phase("pair_storage");
        if(memory_enabled()){MemoryUsage m;m.pairs(P);m.report("pair_storage",P.size());}
        if(plan&&plan->pairs!=P.size())throw std::runtime_error("input pair count changed between passes");
        std::vector<FamilyKey> keys;
        auto groups=plan&&plan->supported?std::move(plan->groups):families(P,keys,threads,shared_pool);
        phase("families_total");
        if(memory_enabled()){
            MemoryUsage m;m.pairs(P);m.family=key_bytes(keys)+vector_bytes(groups);
            for(const auto& group:groups)m.family+=vector_bytes(group.m);
            m.report("family_groups",P.size());
        }
        std::vector<FamilyKey>().swap(keys);
        std::vector<char> disband(groups.size());
        const auto quality_errors=family_quality_errors();
        std::ostringstream ledger; long prot = 0, rows = 0, fam = 0, dis = 0, members = 0, corrected = 0, confirmed = 0, unresolved = 0, excluded = 0;
        struct CorrectedChunk {std::string ledger;std::array<long,9> counts;};
        // Groups partition P; each worker changes only its own group's reads.
        // Ordered consumption preserves ledger rows and all integer statistics.
        ordered_index_chunks(pool,groups.size(),threads,0,[&](std::size_t begin,std::size_t n){
        std::ostringstream ledger; long prot = 0, rows = 0, fam = 0, dis = 0, members = 0, corrected = 0, confirmed = 0, unresolved = 0, excluded = 0;
        for (size_t gi = begin; gi < begin+n; ++gi) {
            const auto& g0 = groups[gi].m;
            std::array<const std::vector<Column>*,2> prepared{};
            disband[gi]=groups[gi].oversize||second_molecule(P,g0,&prepared);
            if(disband[gi]){++dis;continue;}
            const auto hiq=high_quality_counts(P,g0,prepared);
            std::vector<int> g;for(size_t member=0;member<g0.size();++member)if(hiq[member]<2)g.push_back(g0[member]);
            excluded += long(g0.size() - g.size()); if (g.size() < 2) { ++dis; continue; }
            ++fam; members += long(g.size());
            // With no exclusions, masks still refer to the identical ordered
            // g0 members. Other-side edits cannot invalidate these columns.
            for (int side = 0; side < 2; ++side) { const auto& cols = g.size()==g0.size()?*prepared[side]:columns(P, g, side);
                for (size_t c = 0; c < cols.size(); ++c) { auto& col = cols[c];
                    auto set = [&](int i, char base, int q) { auto& r = P[i].r[side]; if (r.s[c] == base && qv(r, c) == q) return; char ob = r.s[c]; int oq = qv(r, c);
                        r.s[c] = base; r.q[c] = char(q + 33); if(debug)ledger << P[i].id << '\t' << side + 1 << '\t' << c << '\t' << ob << '\t' << oq << '\t' << base << '\t' << q << '\n'; ++rows; };
                    if (col.n >= 3 && col.votes.size() == 1) { ++confirmed; for (int i : col.votes[0].second) set(i, col.votes[0].first, std::max(qv(P[i].r[side], c), 30)); }
                    else if (col.n >= 3 && col.top >= 0 && col.votes[col.top].second.size() >= 3 && col.votes[col.top].second.size() * 3 >= size_t(2 * col.n)) {
                        ++corrected; int q = 0; for (int i : col.votes[col.top].second) q = std::max(q, qv(P[i].r[side], c)); q = std::min(40, q); char top = col.votes[col.top].first;
                        for (auto v : col.votes) for (int i : v.second) { if (damage_safe && v.first != top && transition(v.first, top) && is_protected(P[i].r[side], c, side)) { ++prot; continue; } set(i, top, q); } }
                    else if (col.n >= 2 && col.votes.size() > 1) { ++unresolved; std::array<double,MAXCOMP> w;w.fill(1.0); double Z = 0; for (size_t a = 0; a < col.votes.size(); ++a) { for (size_t b = 0; b < col.votes.size(); ++b) for (int i : col.votes[b].second) { double e = quality_errors[static_cast<unsigned char>(P[i].r[side].q[c])]; w[a] *= a == b ? 1 - e : e / 3; } Z += w[a]; }
                        for (size_t a = 0; a < col.votes.size(); ++a) { auto v = col.votes[a]; for (int i : v.second) { bool tr = false; for (auto w : col.votes) if (w.first != v.first && transition(v.first, w.first)) tr = true; if (damage_safe && tr && is_protected(P[i].r[side], c, side)) { ++prot; continue; } set(i, P[i].r[side].s[c], g_phred ? std::clamp(int(std::lround(-10 * std::log10(std::max(1 - w[a] / Z, 1e-4)))), 2, 40) : 2); } } } } }
        }
        return CorrectedChunk{ledger.str(),{prot,rows,fam,dis,members,corrected,confirmed,unresolved,excluded}};
        },[&](CorrectedChunk part){
            if(debug)ledger<<part.ledger;const auto& c=part.counts;
            prot+=c[0];rows+=c[1];fam+=c[2];dis+=c[3];members+=c[4];
            corrected+=c[5];confirmed+=c[6];unresolved+=c[7];excluded+=c[8];
        },"family_correct");
        phase("correction");
        std::ostringstream f, stats;
        if(debug){f << "family\tsize\tdiameter\tdisbanded\tmembers\n";
          for (size_t k = 0; k < groups.size(); ++k) { f << k << '\t' << groups[k].m.size() << '\t' << groups[k].diameter << '\t' << (disband[k] ? 1 : 0) << '\t';
              for (size_t x = 0; x < groups[k].m.size(); ++x) f << (x ? "," : "") << P[groups[k].m[x]].id; f << '\n'; }}        stats << "{\"pairs\":" << P.size() << ",\"groups\":" << groups.size() << ",\"families\":" << fam << ",\"disbanded\":" << dis << ",\"members\":" << members
                  << ",\"corrected\":" << corrected << ",\"confirmed\":" << confirmed << ",\"unresolved\":" << unresolved << ",\"excluded\":" << excluded << (damage_safe ? ",\"protected\":" + std::to_string(prot) : std::string()) << ",\"ledger_rows\":" << rows << "}\n";
        FamilyResult result{std::move(P), debug?"pair_id\tmate\tpos\told_base\told_Q\tbase\tQ\n" + ledger.str():std::string(), f.str(), stats.str()};
        phase("serialization");return result;
}
#include "compact_family.inc"
std::vector<FamilyRead> read_family_fastq(const std::string& path) {
    const auto start=std::chrono::steady_clock::now();auto reads=read_fastq(path);
    if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS"))std::cerr<<"FAMILY_READ path="<<path<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    return reads;
}
}
