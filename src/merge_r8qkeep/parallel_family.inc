// Included inside neighbor_families_impl after immutable indexes are built.
// Component roots are internal; families() emits members and groups by ascending
// read ID. Atomically union every verified edge, then join before that traversal.
// Each part's existing sorted posting order also orders its queries. Visiting
// one part at a time keeps related bitmap neighborhoods together in cache.
struct QueryStats {std::uint64_t lookups=0,candidates=0,verified=0;};
fqdup::r8qkeep::AtomicFamilyUnion components(parent);
const auto enumerate=[&](int part,std::size_t query,QueryStats& stats){
    const fqdup::r8qkeep::TwoBitFamilyKey key(key_at(query));
    index[part].visit_unordered(key.part(part),[&](int other){
        if(std::size_t(other)>=query)return true;
        const auto& candidate=key_at(other);
        const auto d0=fqdup::r8qkeep::TwoBitFamilyKey::differences(key.lo,candidate.lo);
        const auto d1=fqdup::r8qkeep::TwoBitFamilyKey::differences(key.hi,candidate.hi);
        if(part>0&&__builtin_popcountll(d0&((1ULL<<26)-1))<=2)return true;
        if(part>1&&__builtin_popcountll(d0>>26)+__builtin_popcountll(d1&((1ULL<<12)-1))<=2)return true;
        ++stats.verified;
        if(__builtin_popcountll(d0)+__builtin_popcountll(d1)<=MAXMIS)components.unite(source_id(query),source_id(other));
        return true;
    },stats.lookups,stats.candidates);
};
std::uint64_t lookups=0,candidates=0,verified=0,overflows=0;
constexpr std::size_t scratch_bytes=sizeof(QueryStats)*8;
const auto consume=[&](QueryStats stats){lookups+=stats.lookups;candidates+=stats.candidates;verified+=stats.verified;};
for(int part=0;part<3;++part){
    const auto& queries=index[part].ordered_queries();
    if(!pool){QueryStats stats;for(int query:queries)enumerate(part,query,stats);consume(stats);}
    else{
        // Long-lived workers claim bounded ranges without per-query futures,
        // result journals, ordered consumption, or overflow recomputation.
        std::atomic<std::size_t> next{0};std::vector<std::future<QueryStats>> tasks;tasks.reserve(workers);
        try{
            for(int worker=0;worker<workers;++worker)tasks.push_back(pool->submit([&]{
                QueryStats stats;
                for(;;){const auto begin=next.fetch_add(128,std::memory_order_relaxed);if(begin>=queries.size())break;
                    const auto end=std::min(queries.size(),begin+128);for(auto q=begin;q<end;++q)enumerate(part,queries[q],stats);}
                return stats;
            }));
            for(auto& task:tasks)consume(task.get());
        }catch(...){for(auto& task:tasks)if(task.valid())task.wait();throw;}
    }
}
(void)edge_limit; // the concurrent path needs no bounded edge journal
