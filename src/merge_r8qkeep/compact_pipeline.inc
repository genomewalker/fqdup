// Included by pipeline.cpp inside fqdup::r8qkeep. Fitting consumes exactly the
// same ordered observations; only the source of decoded pair scratch changes.
LibraryPipeline prepare_compact_library(std::shared_ptr<CompactPairStore> store,const NativeSettings& settings,int threads,bool debug,bool timings,bool counts,bool no_family) {
    auto start=std::chrono::steady_clock::now();
    const auto stage=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(timings||debug)std::cerr<<"MERGE_STAGE "<<name<<" seconds="<<std::chrono::duration<double>(now-start).count()<<memory_status()<<'\n';start=now;};
    LibraryPipeline out;out.native_settings=settings;out.compact=store;
    const int workers=fqdup::r8qkeep::bounded_worker_count(threads);out.pool=std::make_shared<WorkerPool>(workers);
    const bool memory=memory_enabled();
    const auto memory_stage=[&](const char* name){if(memory){library_memory(out).report(name,store->size());const auto bytes=store->capacity_bytes();std::cerr<<"MERGE_COMPACT stage="<<name<<" pairs="<<store->size()<<" bases_bytes="<<bytes[0]<<" quality_bytes="<<bytes[1]<<" header_bytes="<<bytes[2]<<" offset_bytes="<<bytes[3]<<" exception_bytes="<<bytes[4]<<" quality_width="<<store->quality_width()<<" alphabet="<<store->quality_alphabet_size()<<'\n';}};
    memory_stage("compact_input");
    std::vector<std::pair<std::string,std::size_t>> selected;
    const bool unique=store->canonicalize(out.pool.get(),workers,counts?&selected:nullptr,
        [](std::string_view id)->std::uint64_t{return XXH64(id.data(),id.size(),prior_sample_seed);});
    if(!counts)store->input_order();
    stage("canonical_input");
    // Full census merges repeated IDs in one scalar accumulator. Preserve that
    // unusual-input path until a compact global accumulator is available.
    const auto legacy=[&](){
        std::vector<FamilyRead> first(store->size()),second(store->size());FamilyPair pair;CompactPairStore::Decoder decoder;
        for(std::size_t i=0;i<store->size();++i){store->read(i,pair,decoder);first[i]=std::move(pair.r[0]);second[i]=std::move(pair.r[1]);}
        out.compact.reset();store.reset();
        return prepare_library(std::move(first),std::move(second),settings,threads,debug,timings,true,counts,no_family);
    };
    if(!counts&&!unique)return legacy();
    if(!no_family)try{out.family=correct_compact_family(*store,threads,debug,out.pool.get());}
    catch(const UnsupportedCompactAlphabet&){return legacy();}
    stage("family");memory_stage("family");
    if(counts){
        constexpr std::size_t requested=25000;
        CompactPairStore::Decoder decoder;std::string id;
        if(std::getenv("FQDUP_PRIOR_SAMPLE_VERIFY")){
            std::vector<std::pair<std::string,std::size_t>> reference;
            if(store->size()){
            const auto threshold=(static_cast<unsigned __int128>(std::min(requested,store->size()))<<64)/store->size();
            for(std::size_t i=0;i<store->size();++i){store->id(i,id,decoder);if(static_cast<unsigned __int128>(XXH64(id.data(),id.size(),prior_sample_seed))<threshold)reference.emplace_back(id,i);}
            }
            std::sort(reference.begin(),reference.end(),[](const auto& a,const auto& b){return a.first<b.first;});
            std::sort(selected.begin(),selected.end(),[](const auto& a,const auto& b){return a.first<b.first;});
            if(const char* prefix=std::getenv("FQDUP_PRIOR_SAMPLE_DUMP")){
                const auto dump=[&](const char* suffix,const auto& rows){
                    std::ofstream file(std::string(prefix)+suffix);
                    if(!file)throw std::runtime_error("cannot open prior sample dump");
                    for(const auto& row:rows)file<<row.first<<'\t'<<row.second<<'\n';
                    file.close();if(!file)throw std::runtime_error("cannot write prior sample dump");
                };
                dump(".before.tsv",reference);dump(".after.tsv",selected);
            }
            if(selected!=reference)throw std::logic_error("canonical prior sample differs from reference scan");

            std::cerr<<"PRIOR_SAMPLE_VERIFY pairs="<<reference.size()<<" identical=1\n";
        }
        std::sort(selected.begin(),selected.end(),[](const auto& a,const auto& b){return a.first<b.first;});
        for(std::size_t i=1;i<selected.size();++i)if(selected[i-1].first==selected[i].first)throw std::invalid_argument("subsample prior fit requires unique selected read IDs");
        std::vector<FamilyPair> sample(selected.size());decoder={};
        for(std::size_t i=0;i<selected.size();++i)store->read(selected[i].second,sample[i],decoder);
        if(timings||debug)std::cerr<<"PRIOR_SAMPLE requested="<<requested<<" selected="<<sample.size()<<" total="<<store->size()<<" seed="<<prior_sample_seed<<'\n';
        const auto native=merge_native(sample,settings,timings||debug,debug,threads,out.pool.get());stage("sample_native");
        out.observations=census(sample,native,debug,threads,0,true,out.pool.get());stage("sample_census");
    }else{
        struct Part {CensusResult observations;std::string native;};bool first=true;
        ordered_index_chunks(*out.pool,store->size(),threads,0,[&](std::size_t begin,std::size_t n){
            static thread_local std::vector<FamilyPair> raw;raw.resize(n);CompactPairStore::Decoder decoder;
            for(std::size_t i=0;i<n;++i)store->read(begin+i,raw[i],decoder);
            auto native=merge_native_range(raw,0,n,settings,false,debug);
            auto observations=census(raw,native,debug,1,begin,false);
            return Part{std::move(observations),std::move(native.native_tsv)};
        },[&](Part part){
            const auto append=[&](std::string& into,const std::string& text){into+=first?text:text.substr(text.find('\n')+1);};
            if(debug){append(out.observations.nominal_tsv,part.observations.nominal_tsv);append(out.observations.full_tsv,part.observations.full_tsv);out.native.native_tsv+=part.native;}
            out.observations.proposed_edits_tsv+=part.observations.proposed_edits_tsv;
            for(auto& p:part.observations.nominal_pairs)out.observations.nominal_pairs.push_back(std::move(p));
            for(auto& p:part.observations.full_pairs)out.observations.full_pairs.push_back(std::move(p));
            first=false;
        },"compact_full_census",native_batch_size);
        if(store->empty()){std::vector<FamilyPair> empty;NativeResult native;out.observations=census(empty,native,debug,1,0,false);}
        stage("census");memory_stage("census");
    }
    std::array<std::vector<PairLikelihood>,2> pairs{std::move(out.observations.nominal_pairs),std::move(out.observations.full_pairs)};
    std::array<std::array<std::optional<MixtureFit>,2>,2> fits;
    // Run the two folds of one mode at once, sharing the whole thread budget.
    // This keeps nominal workers busy even when the full fits are tiny.
    const int concurrent=std::max(1,std::min(2,workers));
    // At most two fit coordinators occupy pool workers. Nested counts work
    // uses only the remaining capacity; at p1/p2 no nested jobs are submitted.
    const int auxiliary=counts&&concurrent>1?workers-concurrent:0;
    const int per_fit=counts?auxiliary+1:std::max(1,threads/concurrent);
    for(int begin=0;begin<4;begin+=concurrent){std::array<std::future<std::optional<MixtureFit>>,4> futures;
      try{
        for(int k=begin;k<std::min(4,begin+concurrent);++k){const int mode=k/2,f=k%2;std::vector<const PairLikelihood*> training;for(const auto& p:pairs[mode])if(p.fold==f)training.push_back(&p);
            if(concurrent>1){
                auto fit=[training=std::move(training),per_fit,counts,pool=out.pool.get()](){return counts?fit_mixture_counts(training,per_fit,pool):fit_mixture_views(training,per_fit);};
                futures[k]=counts?out.pool->submit(std::move(fit)):std::async(std::launch::async,std::move(fit));
            }else fits[mode][f]=counts?fit_mixture_counts(training,per_fit):fit_mixture_views(training,per_fit);}
        if(concurrent>1)for(int k=begin;k<std::min(4,begin+concurrent);++k)fits[k/2][k%2]=futures[k].get();
      }catch(...){
        // Packaged-task futures do not join on destruction. Keep training
        // rows alive until every sibling has finished, even when one fails.
        for(auto& future:futures)if(future.valid())future.wait();
        throw;
      }
    }
    out.nominal_prior=prior_text(fits[0]);out.full_prior=prior_text(fits[1]);stage("prior_fit");
    if(memory){auto m=library_memory(out);for(const auto& rows:pairs)m.census+=likelihood_bytes(rows);m.report("prior_fit",store->size());}
    if(timings)for(int mode=0;mode<2;++mode)for(int f=0;f<2;++f)if(fits[mode][f]){const auto& fit=*fits[mode][f];std::cerr<<"FIT_CACHE mode="<<mode<<" fold="<<f<<" pairs="<<fit.training_pairs<<" keys="<<fit.unique_likelihoods<<" calls="<<fit.exp_calls<<" hits="<<fit.exp_cache_hits<<'\n';std::cerr<<"FIT_PHASE mode="<<mode<<" fold="<<f<<" em_iterations="<<fit.iterations<<" em_seconds="<<fit.em_seconds<<" polish_iterations="<<fit.polish_iterations<<" polish_seconds="<<fit.polish_seconds<<" last_em_zero_cells="<<fit.last_em_zero_cells<<" last_em_subnormal_cells="<<fit.last_em_subnormal_cells<<" kkt_gap="<<fit.kkt_gap<<" newton_reentries="<<fit.newton_reentries<<" em_cells="<<fit.unique_likelihoods*mixture_components<<'\n';}
    for(auto& observations:pairs)std::vector<PairLikelihood>().swap(observations);
    if(!debug){std::string().swap(out.family.edits_tsv);std::string().swap(out.family.families_tsv);std::string().swap(out.observations.proposed_edits_tsv);}
    out.repair=std::make_unique<RepairPipeline>(out.nominal_prior,out.full_prior);memory_stage("prepared");return out;
}
