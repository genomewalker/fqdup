#pragma once
#include <string>
#include <vector>
#include <memory>
#include <utility>
// Pair-owned immutable alignment products shared across native/router/repair;
// keys include the engine, exact sequence bytes and insert, never DP scratch.
namespace fqdup_shared_alignment {
struct BandedEngine {};
struct Correspondence {
 int best=-1000000000, paths=0; bool eligible=false;
 std::string reason="no_alignment";
 std::vector<std::pair<int,int>> invariant, certified;
};
struct Geometry {
 int insert=-1, informative=0, mismatches=0;
 double fraction=1.0, second=9.0; bool eligible=false;
};
struct Pair {
    struct Bucket { const void* tag; explicit Bucket(const void* t):tag(t){} virtual ~Bucket()=default; };
    template<class T> struct Values : Bucket {
        struct Entry { std::string a,b; int extra; T result; };
        std::vector<Entry> entries;
        explicit Values(const void* t):Bucket(t){}
    };
    std::vector<Bucket*> buckets;
    std::size_t misses=0;
    Pair()=default;
    Pair(const Pair&)=delete;
    Pair& operator=(const Pair&)=delete;
    ~Pair(){for(auto* bucket:buckets)delete bucket;}
};
inline thread_local Pair* active=nullptr;
struct Scope {
    Pair* previous;
    explicit Scope(Pair* pair):previous(active){active=pair;}
    ~Scope(){active=previous;}
};
template<class Tag> inline const char engine_token=0;
template<class Tag,class Compute>
auto get(const std::string& a,const std::string& b,int extra,Compute compute) {
    if(!active)return compute();
    using Result=decltype(compute());
    using Values=Pair::Values<Result>;
    Values* values=nullptr;
    for(auto* bucket:active->buckets)if(bucket->tag==&engine_token<Tag>){values=static_cast<Values*>(bucket);break;}
    if(!values){values=new Values(&engine_token<Tag>);try{active->buckets.push_back(values);}catch(...){delete values;throw;}}
    for(const auto& entry:values->entries)if(entry.extra==extra&&entry.a==a&&entry.b==b)return entry.result;
    auto result=compute();
    values->entries.push_back({a,b,extra,result});++active->misses;
    return result;
}
}
