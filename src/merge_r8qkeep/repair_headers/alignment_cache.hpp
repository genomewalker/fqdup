#pragma once
// Each namespace copy owns its cache. Entries contain only immutable products,
// keyed by engine type and exact input bytes; DP scratch is never retained here.
namespace alignment_cache {
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
template<class Tag,class Compute>
auto get(const std::string& a,const std::string& b,int extra,Compute compute) {
    if(!active)return compute();
    using Result=decltype(compute());
    static const char token=0;
    using Values=Pair::Values<Result>;
    Values* values=nullptr;
    for(auto* bucket:active->buckets)if(bucket->tag==&token){values=static_cast<Values*>(bucket);break;}
    if(!values){values=new Values(&token);try{active->buckets.push_back(values);}catch(...){delete values;throw;}}
    for(const auto& entry:values->entries)if(entry.extra==extra&&entry.a==a&&entry.b==b)return entry.result;
    auto result=compute();
    values->entries.push_back({a,b,extra,result});++active->misses;
    return result;
}
}
