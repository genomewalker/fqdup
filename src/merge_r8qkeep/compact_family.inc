// Included by family.cpp inside fqdup::r8qkeep, reusing its scalar consensus
// primitives. Families partition the input. Workers decode at most MAXCOMP
// members and retain their string capacities across groups.
namespace {
FamilyKey compact_family_key(const FamilyPair& pair) {
    const auto& a=pair.r[0].s;const auto& b=pair.r[1].s;
    if(a.size()<KEYLEN||b.size()<KEYLEN)return {};
    std::array<char,2*KEYLEN> key;
    std::copy_n(a.data(),KEYLEN,key.data());std::copy_n(b.data(),KEYLEN,key.data()+KEYLEN);
    if(std::count(key.begin(),key.end(),'N')>MAXN)return {};
    for(int side=0;side<2;++side){
        const auto& read=pair.r[side].s;
        const auto at=read.find(ADAPTERS[side].data(),0,ADSEED);
        if(at!=std::string::npos&&at<KEYLEN)return {};
        auto end=std::min(read.size(),at),g=read.size();
        while(g&&read[g-1]=='G')--g;
        if(read.size()-g>=10)end=std::min(end,g);
        if(low_complexity(std::string_view(read).substr(0,end)))return {};
    }
    return FamilyKey(std::string_view(key.data(),key.size()));
}
struct CompactEdit {std::uint32_t pair,position;unsigned char mate;char base,quality;bool base_changed;};
static_assert(sizeof(CompactEdit)==12); // Flag uses the existing padding byte.
}
FamilyResult correct_compact_family(CompactPairStore& store,int threads,bool debug,WorkerPool* shared_pool) {
    std::unique_ptr<WorkerPool> owned;
    if(!shared_pool){owned=std::make_unique<WorkerPool>(std::max(1,std::min(threads,8)));shared_pool=owned.get();}
    auto& pool=*shared_pool;
    std::vector<FamilyKey> keys(store.size());
    ordered_index_chunks(pool,store.size(),threads,0,[&](std::size_t begin,std::size_t n){
        FamilyPair pair;CompactPairStore::Decoder decoder;
        for(std::size_t i=begin;i<begin+n;++i){store.sequences(i,pair);keys[i]=compact_family_key(pair);}
        return 0;
    },[](int){},"compact_family_keys");
    auto plan=discover_family_keys(std::move(keys),threads,shared_pool);
    if(!plan.supported)throw UnsupportedCompactAlphabet{};
    auto groups=std::move(plan.groups);
    auto phase_start=std::chrono::steady_clock::now();
    const auto phase=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS"))std::cerr<<"COMPACT_FAMILY_PHASE "<<name<<" seconds="<<std::chrono::duration<double>(now-phase_start).count()<<'\n';phase_start=now;};
    constexpr bool damage_safe=false;
    std::vector<CompactEdit> all_edits;
        std::vector<char> disband(groups.size());
        const auto quality_errors=family_quality_errors();
        std::ostringstream ledger; long prot = 0, rows = 0, fam = 0, dis = 0, members = 0, corrected = 0, confirmed = 0, unresolved = 0, excluded = 0;
        struct CorrectedChunk {std::string ledger;std::array<long,9> counts;std::vector<CompactEdit> edits;};
        // Groups partition P; each worker changes only its own group's reads.
        // Ordered consumption preserves ledger rows and all integer statistics.
        ordered_index_chunks(pool,groups.size(),threads,0,[&](std::size_t begin,std::size_t n){
        static thread_local std::vector<Pair> P(MAXCOMP);
        CompactPairStore::Decoder decoder;std::vector<CompactEdit> edits;
        std::ostringstream ledger; long prot = 0, rows = 0, fam = 0, dis = 0, members = 0, corrected = 0, confirmed = 0, unresolved = 0, excluded = 0;
        for (size_t gi = begin; gi < begin+n; ++gi) {
            if(groups[gi].oversize){disband[gi]=1;++dis;continue;}
            const auto& original=groups[gi].m;std::vector<int> g0(original.size());
            std::iota(g0.begin(),g0.end(),0);
            for(std::size_t member=0;member<original.size();++member){
                if(debug)store.read(original[member],P[member],decoder);
                else store.content(original[member],P[member],decoder);
            }
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
                        r.s[c] = base; r.q[c] = char(q + 33); edits.push_back({std::uint32_t(original[i]),std::uint32_t(c),static_cast<unsigned char>(side),base,char(q+33),base!=ob}); if(debug)ledger << P[i].id << '\t' << side + 1 << '\t' << c << '\t' << ob << '\t' << oq << '\t' << base << '\t' << q << '\n'; ++rows; };
                    if (col.n >= 3 && col.votes.size() == 1) { ++confirmed; for (int i : col.votes[0].second) set(i, col.votes[0].first, std::max(qv(P[i].r[side], c), 30)); }
                    else if (col.n >= 3 && col.top >= 0 && col.votes[col.top].second.size() >= 3 && col.votes[col.top].second.size() * 3 >= size_t(2 * col.n)) {
                        ++corrected; int q = 0; for (int i : col.votes[col.top].second) q = std::max(q, qv(P[i].r[side], c)); q = std::min(40, q); char top = col.votes[col.top].first;
                        for (auto v : col.votes) for (int i : v.second) { if (damage_safe && v.first != top && transition(v.first, top) && is_protected(P[i].r[side], c, side)) { ++prot; continue; } set(i, top, q); } }
                    else if (col.n >= 2 && col.votes.size() > 1) { ++unresolved; std::array<double,MAXCOMP> w;w.fill(1.0); double Z = 0; for (size_t a = 0; a < col.votes.size(); ++a) { for (size_t b = 0; b < col.votes.size(); ++b) for (int i : col.votes[b].second) { double e = quality_errors[static_cast<unsigned char>(P[i].r[side].q[c])]; w[a] *= a == b ? 1 - e : e / 3; } Z += w[a]; }
                        for (size_t a = 0; a < col.votes.size(); ++a) { auto v = col.votes[a]; for (int i : v.second) { bool tr = false; for (auto w : col.votes) if (w.first != v.first && transition(v.first, w.first)) tr = true; if (damage_safe && tr && is_protected(P[i].r[side], c, side)) { ++prot; continue; } set(i, P[i].r[side].s[c], g_phred ? std::clamp(int(std::lround(-10 * std::log10(std::max(1 - w[a] / Z, 1e-4)))), 2, 40) : 2); } } } } }
        }
        return CorrectedChunk{ledger.str(),{prot,rows,fam,dis,members,corrected,confirmed,unresolved,excluded},std::move(edits)};
        },[&](CorrectedChunk part){
            all_edits.insert(all_edits.end(),part.edits.begin(),part.edits.end());
            if(debug)ledger<<part.ledger;const auto& c=part.counts;
            prot+=c[0];rows+=c[1];fam+=c[2];dis+=c[3];members+=c[4];
            corrected+=c[5];confirmed+=c[6];unresolved+=c[7];excluded+=c[8];
        },"family_correct");
        phase("correction");
        const auto group_count=groups.size();
        if(memory_enabled()){MemoryUsage memory;memory.family=vector_bytes(all_edits);memory.report("family_edit_journal",store.size(),"journal_only_excludes_store_and_groups");}
        if(!debug){std::vector<FamilyGroup>().swap(groups);std::vector<char>().swap(disband);}
        store.apply_edits(std::move(all_edits),&pool,std::max(1,std::min(threads,8)),[](const CompactEdit& edit){return edit.base_changed;});
        phase("apply_edits");
        std::ostringstream f, stats;CompactPairStore::Decoder names;std::string id;
        if(debug){f << "family\tsize\tdiameter\tdisbanded\tmembers\n";
          for (size_t k = 0; k < groups.size(); ++k) { f << k << '\t' << groups[k].m.size() << '\t' << groups[k].diameter << '\t' << (disband[k] ? 1 : 0) << '\t';
              for(size_t x=0;x<groups[k].m.size();++x){store.id(groups[k].m[x],id,names);f<<(x?",":"")<<id;} f << '\n'; }}        stats << "{\"pairs\":" << store.size() << ",\"groups\":" << group_count << ",\"families\":" << fam << ",\"disbanded\":" << dis << ",\"members\":" << members
                  << ",\"corrected\":" << corrected << ",\"confirmed\":" << confirmed << ",\"unresolved\":" << unresolved << ",\"excluded\":" << excluded << (damage_safe ? ",\"protected\":" + std::to_string(prot) : std::string()) << ",\"ledger_rows\":" << rows << "}\n";
        return FamilyResult{{},debug?"pair_id\tmate\tpos\told_base\told_Q\tbase\tQ\n"+ledger.str():std::string(),f.str(),stats.str()};
}
