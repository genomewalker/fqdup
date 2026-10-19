#include "content_memo.hpp"
#include "deferred_adapter_clip.hpp"
#include "pipeline.hpp"
#include "canonical_pairs.hpp"
#include "prior_sample.hpp"
#include "ordered_chunks.hpp"
#include "ordered_writer.hpp"
#include "streaming_output.hpp"
#include "parallel_output.hpp"
#include "decoded_pair_pool.hpp"
#include "memory_status.hpp"
#include "string_input.hpp"
#include <algorithm>
#include <string_view>
#include <array>
#include <future>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
namespace fqdup::r8qkeep {
namespace {
std::size_t stage_read_bytes(const StageRead& r){return string_heap(r.header)+string_heap(r.seq)+string_heap(r.qual);}
template<class Record> std::size_t record_bytes(const std::vector<Record>& rows){
    std::size_t bytes=vector_bytes(rows);for(const auto& r:rows)bytes+=string_heap(r.outcome)+stage_read_bytes(r.first)+stage_read_bytes(r.second);return bytes;
}
std::size_t native_bytes(const NativeResult& n){
    MemoryUsage raw;if(n.raw_owner)raw.pairs(*n.raw_owner);
    return raw.raw+raw.ids+record_bytes(n.records)+string_heap(n.merged_fq)+string_heap(n.r1_fq)+string_heap(n.r2_fq)+string_heap(n.orphan_fq)+string_heap(n.outcomes_jsonl)+string_heap(n.native_tsv)+string_heap(n.native_outcomes_json);
}
std::size_t likelihood_bytes(const std::vector<PairLikelihood>& rows){
    std::size_t bytes=vector_bytes(rows);for(const auto& p:rows)bytes+=string_heap(p.id);return bytes;
}
std::size_t routing_bytes(const RoutingResult& r){
    return string_heap(r.routed_tsv)+string_heap(r.audit_tsv)+vector_bytes(r.routed)+strings_bytes(r.replacements)+vector_bytes(r.repair_bounds)+vector_bytes(r.clean_bounds)+vector_bytes(r.native_clean_bounds)+record_bytes(r.record_replacements);
}
MemoryUsage library_memory(const LibraryPipeline& lib){
    MemoryUsage m;m.pairs(lib.family.pairs);if(lib.compact){const auto bytes=lib.compact->capacity_bytes();m.ids+=bytes[2];m.raw+=lib.compact->owned_bytes()-bytes[2];}m.native=native_bytes(lib.native);
    m.census=likelihood_bytes(lib.observations.nominal_pairs)+likelihood_bytes(lib.observations.full_pairs);
    m.routing=routing_bytes(lib.routing);
    m.diagnostics=string_heap(lib.family.edits_tsv)+string_heap(lib.family.families_tsv)+string_heap(lib.family.stats_json)+string_heap(lib.observations.nominal_tsv)+string_heap(lib.observations.full_tsv)+string_heap(lib.observations.proposed_edits_tsv)+string_heap(lib.nominal_prior)+string_heap(lib.full_prior);
    return m;
}
void json_string(std::string& out,const std::string& s){static const char hex[]="0123456789abcdef";out+='"';for(unsigned char c:s){switch(c){case '"':out+="\\\"";break;case '\\':out+="\\\\";break;case '\b':out+="\\b";break;case '\f':out+="\\f";break;case '\n':out+="\\n";break;case '\r':out+="\\r";break;case '\t':out+="\\t";break;default:if(c<32){out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}else out+=char(c);}}out+='"';}
}
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug,bool timings,bool defer_routing,bool counts,bool no_family){
    auto start=std::chrono::steady_clock::now();
    const auto stage=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(timings||debug)std::cerr<<"MERGE_STAGE "<<name<<" seconds="<<std::chrono::duration<double>(now-start).count()<<memory_status()<<'\n';start=now;};
    LibraryPipeline out;out.native_settings=settings;
    const bool memory=memory_enabled();
    const auto memory_stage=[&](const char* name){if(memory)library_memory(out).report(name,out.family.pairs.size());};
    if(memory){MemoryUsage m;m.reads(r1);m.reads(r2);m.report("input",r1.size());}
    const int workers=fqdup::r8qkeep::bounded_worker_count(threads);out.pool=std::make_shared<WorkerPool>(workers);
    if(counts){canonicalize_family_input(r1,r2,out.pool.get(),workers);stage("canonical_input");}
    if(no_family){
        if(r1.size()!=r2.size())throw std::runtime_error("pair count");
        out.family.pairs.reserve(r1.size());
        for(std::size_t i=0;i<r1.size();++i){
            if(r1[i].id!=r2[i].id)throw std::runtime_error("pair order "+r1[i].id);
            out.family.pairs.push_back({r1[i].id,{std::move(r1[i]),std::move(r2[i])}});
        }
        std::vector<FamilyRead>().swap(r1);std::vector<FamilyRead>().swap(r2);
    }else out.family=correct_family(std::move(r1),std::move(r2),threads,debug,out.pool.get());
    stage("family");

    memory_stage("family");
    if(counts){
        constexpr std::size_t requested=25000;
        const auto selected=prior_sample_indices(out.family.pairs,requested);
        std::vector<FamilyPair> sample;sample.reserve(selected.size());
        for(auto i:selected)sample.push_back(std::move(out.family.pairs[i]));
        if(timings||debug)std::cerr<<"PRIOR_SAMPLE requested="<<requested<<" selected="<<sample.size()<<" total="<<out.family.pairs.size()<<" seed="<<prior_sample_seed<<'\n';
        const auto native=merge_native(sample,settings,timings||debug,debug,threads,out.pool.get());stage("sample_native");
        if(memory){auto m=library_memory(out);m.pairs(sample);m.native+=native_bytes(native);m.census+=vector_bytes(selected);m.report("sample_native",out.family.pairs.size());}
        out.observations=census(sample,native,debug,threads,0,true,out.pool.get());
        if(memory){auto m=library_memory(out);m.pairs(sample);m.native+=native_bytes(native);m.census+=vector_bytes(selected);m.report("sample_census",out.family.pairs.size());}
        // Both kernels have joined; observations own their IDs and statistics.
        // Restore the selected slots before fitting or streaming the library.
        for(std::size_t i=0;i<selected.size();++i)out.family.pairs[selected[i]]=std::move(sample[i]);
        stage("sample_census");
    }else{
        out.native=merge_native(out.family.pairs,settings,timings||debug,debug,threads,out.pool.get());stage("native");
        memory_stage("native");
        out.observations=census(out.family.pairs,out.native,debug,threads,0,false,out.pool.get());stage("census");
        memory_stage("census");
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
    if(memory){auto m=library_memory(out);for(const auto& rows:pairs)m.census+=likelihood_bytes(rows);m.report("prior_fit",out.family.pairs.size());}
    if(timings)for(int mode=0;mode<2;++mode)for(int f=0;f<2;++f)if(fits[mode][f]){const auto& fit=*fits[mode][f];std::cerr<<"FIT_CACHE mode="<<mode<<" fold="<<f<<" pairs="<<fit.training_pairs<<" keys="<<fit.unique_likelihoods<<" calls="<<fit.exp_calls<<" hits="<<fit.exp_cache_hits<<'\n';std::cerr<<"FIT_PHASE mode="<<mode<<" fold="<<f<<" em_iterations="<<fit.iterations<<" em_seconds="<<fit.em_seconds<<" polish_iterations="<<fit.polish_iterations<<" polish_seconds="<<fit.polish_seconds<<" last_em_zero_cells="<<fit.last_em_zero_cells<<" last_em_subnormal_cells="<<fit.last_em_subnormal_cells<<" kkt_gap="<<fit.kkt_gap<<" newton_reentries="<<fit.newton_reentries<<" em_cells="<<fit.unique_likelihoods*mixture_components<<'\n';}
    for(auto& observations:pairs)std::vector<PairLikelihood>().swap(observations);
    if(counts&&!defer_routing){out.native=merge_native(out.family.pairs,settings,timings||debug,debug,threads,out.pool.get());stage("native");}
    if(defer_routing){NativeResult empty;empty.native_tsv=std::move(out.native.native_tsv);out.native=std::move(empty);}
    if(!debug){std::string().swap(out.family.edits_tsv);std::string().swap(out.family.families_tsv);std::string().swap(out.observations.proposed_edits_tsv);}
    out.repair=std::make_unique<RepairPipeline>(out.nominal_prior,out.full_prior);if(!defer_routing){out.routing=out.repair->route_records(out.native.records,debug);stage("router");}memory_stage("prepared");return out;
}
#include "compact_pipeline.inc"
namespace {
std::string_view read_id(const std::string& h){
    if(h.empty()||h[0]!='@')throw std::runtime_error("record identity");
    const auto end=h.find_first_of(" \t");
    std::string_view view(h.data()+1,(end==std::string::npos?h.size():end)-1);
    if(view.size()>2&&view[view.size()-2]=='/'&&(view.back()=='1'||view.back()=='2'))view.remove_suffix(2);
    return view;
}
void append_record(FinalOutput& out,const RepairedRecord& record,const FamilyPair& pair,OutcomesMode mode){
    const auto& kind=record.native?record.native->outcome:record.outcome;
    const auto& first=record.native?record.native->first:record.first;
    const auto& second=record.native?record.native->second:record.second;
    const auto add=[&](const StageRead& r,std::string& text,const char* channel,int mate){
        if(read_id(r.header)!=std::string_view(pair.id))throw std::runtime_error("record identity");
        text+=r.header;
        if(record.native&&!r.seq.empty()){
            const auto length=std::to_string(r.seq.size());
            text+=" FQT=v1;ch=";text+=channel;text+=";m=";text+=char('0'+mate);
            text+=";olen=";text+=length;text+=";s=0;e=";text+=length;
        }
        text+='\n';text+=r.seq;text+="\n+\n";text+=r.qual;text+='\n';
    };
    if(kind=="merged")add(first,out.merged,"merged",0);
    else if(kind=="unmerged"){add(first,out.r1,"r1",1);add(second,out.r2,"r2",2);}
    else if(kind=="orphan_r1"||kind=="orphan_r2")add(first,out.orphan,"orphan",kind=="orphan_r1"?1:2);
    else if(kind!="dropped")throw std::runtime_error("chunk outcome");
    if(mode==OutcomesMode::none)return;
    if(mode==OutcomesMode::summary){
        const std::size_t index=kind=="merged"?0:kind=="unmerged"?1:kind=="orphan_r1"?2:kind=="orphan_r2"?3:4;
        ++out.outcome_counts[index];return;
    }
    out.outcomes+="{\"r1_header\":";json_string(out.outcomes,pair.r[0].h);
    out.outcomes+=",\"r2_header\":";json_string(out.outcomes,pair.r[1].h);
    out.outcomes+=",\"outcome\":";json_string(out.outcomes,kind);out.outcomes+="}\n";
}
struct OutputChunk {
    std::vector<FinalOutput> outputs;
    std::vector<double> seconds;
    double router_seconds=0,native_seconds=0;
    std::string routed,audit,native;
    std::size_t count=0;
    std::size_t native_clean=0,native_direct=0;
    std::vector<std::size_t> native_clean_pass;
};
template<class Consume>
void process_records(const LibraryPipeline& lib,const std::vector<double>& cuts,int threads,bool debug,RepairShadow* counters,Consume consume){
    using Clock=std::chrono::steady_clock;
    const auto elapsed=[](auto start){return std::chrono::duration<double>(Clock::now()-start).count();};
    const int workers=fqdup::r8qkeep::bounded_worker_count(threads);
    const std::size_t chunk_rows=128;
    // Deeper single-cut queues hide uneven pair work and output stalls. Multi-
    // cutoff jobs retain at most 512 output rows per worker.
    const std::size_t window=workers*(cuts.size()==1?16:4);
    const bool memory=memory_enabled();
    std::atomic<std::size_t> max_native{0},max_routing{0},max_output{0};
    const auto maximum=[](std::atomic<std::size_t>& dst,std::size_t value){auto old=dst.load();while(old<value&&!dst.compare_exchange_weak(old,value)){};};
    std::unique_ptr<WorkerPool> owned_pool;
    if(!lib.pool)owned_pool=std::make_unique<WorkerPool>(workers);
    auto& pool=lib.pool?*lib.pool:*owned_pool;const auto pipeline_begin=Clock::now();
    const bool native_debug=debug&&lib.native.native_tsv.empty();
    const auto structural_switch=std::getenv("FQDUP_STRUCTURAL_SKIP");
    const bool structural=!debug&&!counters&&(!structural_switch||std::string_view(structural_switch)!="0")
        &&!cuts.empty()&&std::all_of(cuts.begin(),cuts.end(),[](double cut){return std::isfinite(cut)&&cut>=0&&cut<=1;});
    const auto process=[&](NativeRecords native,std::size_t start,double native_seconds,std::vector<MemoValue>* collected=nullptr,const std::vector<std::shared_ptr<fqdup_shared_alignment::Pair>>* shared=nullptr,std::size_t shared_offset=0){
        const std::size_t n=native.size();if(collected)collected->resize(n);
        OutputChunk chunk;chunk.count=n;chunk.native_seconds=native_seconds;chunk.outputs.resize(cuts.size());chunk.seconds.resize(cuts.size());
        if(native_debug)for(const auto& record:native)chunk.native+=native_record_tsv(record);
        const auto route_begin=Clock::now();const double minimum=*std::min_element(cuts.begin(),cuts.end());
        std::vector<int64_t> certified,unresolved_bounds;std::vector<std::size_t> unresolved;
        NativeRecords work=native;std::size_t work_start=start;
        if(structural){
            certified=lib.repair->native_clean_bounds(native,minimum);
            unresolved.reserve(n);unresolved_bounds.reserve(n);
            const auto budget=int64_t(std::floor(minimum*1000000000000LL+.5));
            for(std::size_t i=0;i<n;++i){
                if(certified[i]>=0&&certified[i]<=budget*int64_t(native[i][8].size()))continue;
                unresolved.push_back(i);unresolved_bounds.push_back(certified[i]);
            }
            work=native.select(unresolved);work_start=0;chunk.native_direct=n-work.size();
        }
        // Compact slots only key pair-local scratch here. Diagnostic paths keep
        // original indices and run the complete reference router/repair path.
        RepairAlignmentCache cache(work_start,work.size(),shared,shared_offset,structural?&unresolved:nullptr);
        const auto routing=lib.repair->route_records(work,debug,work_start,&cache,counters!=nullptr,minimum,structural?&unresolved_bounds:nullptr);
        const auto& bounds=structural?certified:routing.native_clean_bounds;
        chunk.native_clean=std::count_if(bounds.begin(),bounds.end(),[](int64_t bound){return bound>=0;});
        chunk.native_clean_pass.resize(cuts.size());
        for(std::size_t c=0;c<cuts.size();++c)for(std::size_t i=0;i<n;++i){
            const auto bound=bounds[i];
            if(bound>=0&&bound<=int64_t(std::floor(cuts[c]*1000000000000LL+.5))*int64_t(native[i][8].size()))++chunk.native_clean_pass[c];
        }
        chunk.router_seconds=elapsed(route_begin);
        const bool alternative=std::find(routing.routed.begin(),routing.routed.end(),true)!=routing.routed.end();
        if(debug){chunk.routed=routing.routed_tsv;chunk.audit=routing.audit_tsv;}
        for(std::size_t c=0;c<cuts.size();++c){
            const auto begin=Clock::now();
            const auto base=lib.repair->repair_records(work,9,cuts[c],work_start,counters?&counters[c]:nullptr,&routing.routed,&routing.repair_bounds,&cache,&routing.clean_bounds);
            const auto changed=alternative?lib.repair->repair_records(work,7,cuts[c],work_start,nullptr,nullptr,nullptr,&cache,nullptr,&routing.routed,&routing.record_replacements):std::vector<RepairedRecord>{};
            if(memory)maximum(max_routing,routing_bytes(routing)+record_bytes(base)+record_bytes(changed)+vector_bytes(certified)+vector_bytes(unresolved)+vector_bytes(unresolved_bounds));
            if(base.size()!=work.size()||(alternative&&changed.size()!=work.size()))throw std::runtime_error("chunk cardinality");
            const auto emit=[&](std::size_t i,const RepairedRecord& selected){
                RepairedRecord clipped;
                const auto& record = native[i].pending_adapter_clip >= 0
                    ? (clipped = final_adapter_clip(selected, native[i].pending_adapter_clip))
                    : selected;
                if(collected){const auto bound=bounds[i];
                    (*collected)[i]=memo_owned(record,bound>=0,bound>=0&&bound<=int64_t(std::floor(cuts[c]*1000000000000LL+.5))*int64_t(native[i][8].size()));}
                else append_record(chunk.outputs[c],record,*native[i].raw,lib.outcomes_mode);
            };
            std::size_t slot=0;RepairedRecord direct;
            for(std::size_t i=0;i<n;++i){
                if(structural&&(slot==unresolved.size()||unresolved[slot]!=i)){
                    direct.native=&native[i];emit(i,direct);
                }else{const auto at=structural?slot++:i;emit(i,routing.routed[at]?changed[at]:base[at]);}
            }
            chunk.seconds[c]=elapsed(begin);
        }
        if(memory){
            std::size_t bytes=vector_bytes(chunk.outputs)+vector_bytes(chunk.seconds)+vector_bytes(chunk.native_clean_pass)+string_heap(chunk.routed)+string_heap(chunk.audit)+string_heap(chunk.native);
            for(const auto& out:chunk.outputs)bytes+=string_heap(out.merged)+string_heap(out.r1)+string_heap(out.r2)+string_heap(out.orphan)+string_heap(out.outcomes);
            maximum(max_output,bytes);
        }
        return chunk;
    };
    const auto content_memo_switch=std::getenv("FQDUP_CONTENT_MEMO");
    const bool force_memo=content_memo_switch&&std::string_view(content_memo_switch)=="1";
    #include "memo_pipeline.inc"
    const auto memo_switch=content_memo_switch?content_memo_switch:std::getenv("FQDUP_MEMO");
    const bool use_memo=lib.compact&&cuts.size()==1&&!debug&&!counters&&(!memo_switch||std::string_view(memo_switch)!="0");
    OrderedWriter<OutputChunk> writer(window,std::move(consume));
    std::size_t output_rows=0,output_chunks=0;
    const auto dispatch=[&](NativeRecords records,std::shared_ptr<const NativeResult> owner,std::size_t offset,double native_seconds){
        for(std::size_t pos=0;pos<records.size();pos+=chunk_rows){
            const auto n=std::min(chunk_rows,records.size()-pos),start=offset+pos;
            const auto native=records.slice(pos,n);const double seconds=pos==0?native_seconds:0;
            // Native views stay valid even after the coordinator advances to
            // another batch. Results own their emitted strings before returning.
            writer.submit(pool.submit([&,native,start,seconds,pos,keep_alive=owner]{
                return process(native,start,seconds,nullptr,keep_alive?&keep_alive->alignments:nullptr,pos);
            }));
            output_rows+=n;++output_chunks;
        }
    };
    try{
        if(use_memo){
            for(std::size_t start=0;start<lib.pair_count();start+=512){const auto n=std::min<std::size_t>(512,lib.pair_count()-start);
                writer.submit(pool.submit([&,start,n]{return memo_process(start,n);}));output_rows+=n;++output_chunks;}
        }else if(!lib.native.records.empty())dispatch(lib.native.records,{},0,0);
        else{
            struct NativeChunk {NativeResult result;std::size_t start;double seconds;};
            auto raw_pool=std::make_shared<DecodedPairPool>(native_batch_size);
            std::deque<std::future<NativeChunk>> pending;std::size_t next=0;
            const auto launch=[&]{
                const auto start=next,n=std::min(native_batch_size,lib.pair_count()-next);next+=n;
                pending.push_back(pool.submit([&,start,n]{const auto begin=Clock::now();NativeResult result;
                    if(lib.compact){
                        auto raw=raw_pool->acquire();CompactPairStore::Decoder decoder;
                        for(std::size_t i=0;i<n;++i)lib.compact->read(start+i,(*raw)[i],decoder);
                        result=merge_native_range(*raw,0,n,lib.native_settings,false,false);result.raw_owner=std::move(raw);
                    }else result=merge_native_range(lib.family.pairs,start,n,lib.native_settings,false,false);return NativeChunk{std::move(result),start,elapsed(begin)};}));
            };
            while(next<lib.pair_count()&&pending.size()<std::size_t(workers))launch();
            while(!pending.empty()){
                auto chunk=pending.front().get();pending.pop_front();
                if(next<lib.pair_count())launch();
                auto owner=std::make_shared<NativeResult>(std::move(chunk.result));
                if(memory)maximum(max_native,native_bytes(*owner));
                dispatch(owner->records,owner,chunk.start,chunk.seconds);
            }
        }
        writer.finish();
        std::cerr<<"MERGE_CHUNK rows="<<output_rows<<" chunks="<<output_chunks<<" chunk_rows="<<(use_memo?512:chunk_rows)<<" workers="<<workers<<" max_inflight="<<window<<'\n';
        std::cerr<<"MERGE_MEMO enabled="<<use_memo<<" forced="<<force_memo<<" pilot_rows="<<memo_pilot_rows<<" rows="<<memo_rows.load()<<" hits="<<memo_hits.load()<<" bypassed="<<memo_bypassed.load()<<" hit_rate="<<(memo_rows.load()?double(memo_hits.load())/memo_rows.load():0)<<" cache_capacity_bound_per_worker="<<ContentMemo::capacity_bound()<<'\n';
    }catch(...){pool.stop();throw;}
    if(memory){
        MemoryUsage m;m.native=max_native.load();m.routing=max_routing.load();m.output=max_output.load();
        m.report("streaming_chunk_max",lib.pair_count(),"independent_chunk_maxima_excludes_cache_and_dp_scratch");
        // Native batches may be referenced by queued outputs after the producer
        // advances. This is a conservative capacity bound, not simultaneous use.
        m.native*=use_memo?workers:workers+window+1;m.routing*=workers;m.output*=window+workers+1;
        m.report("streaming_buffer_bound",lib.pair_count(),"capacity_bound_excludes_cache_and_dp_scratch");
        library_memory(lib).report("streaming_finished",lib.pair_count());
    }
    if(std::getenv("FQDUP_CHUNK_PROFILE"))std::fprintf(stderr,
        "CHUNK_PROFILE stage=streaming_pipeline workers=%d chunk_rows=%zu max_inflight=%zu rows=%zu chunks=%zu wall=%.6f timing_basis=wall_including_native_and_consume\n",
        workers,chunk_rows,window,output_rows,output_chunks,elapsed(pipeline_begin));
}
using OutputFiles=std::array<StreamingOutput,5>;
void open_file(StreamingOutput& file,const std::string& directory,std::size_t index){
    static const char* names[]={"merged.fq","r1.fq","r2.fq","orphan.fq","outcomes.jsonl"};
    file.open(std::filesystem::path(directory)/names[index]);
    if(!file)throw std::runtime_error("open output");
}
void open_files(OutputFiles& files,const std::string& directory){
    for(std::size_t i=0;i<files.size();++i)open_file(files[i],directory,i);
}
void write_files(OutputFiles& files,const FinalOutput& out){
    files[0].append(out.merged);files[1].append(out.r1);files[2].append(out.r2);files[3].append(out.orphan);if(!out.outcomes.empty())files[4].append(out.outcomes);
}
void close_files(OutputFiles& files){for(auto& file:files){file.close();if(!file)throw std::runtime_error("write output");}}
}
FinalOutput emit_library(const LibraryPipeline& lib,double cutoff){
    FinalOutput result;std::size_t count=0;
    process_records(lib,{cutoff},1,false,nullptr,[&](OutputChunk chunk){
        count+=chunk.count;const auto& out=chunk.outputs[0];
        result.merged+=out.merged;result.r1+=out.r1;result.r2+=out.r2;result.orphan+=out.orphan;result.outcomes+=out.outcomes;
    });
    if(count!=lib.pair_count())throw std::runtime_error("final cardinality");
    return result;
}
void write_library_outputs(const LibraryPipeline& lib,double cutoff,const std::string& directory,int threads,RepairShadow* shadow){
    OutputFiles files;open_files(files,directory);std::size_t count=0;
    process_records(lib,{cutoff},threads,false,shadow,[&](OutputChunk chunk){count+=chunk.count;write_files(files,chunk.outputs[0]);});
    if(count!=lib.pair_count())throw std::runtime_error("final cardinality");
    close_files(files);
}
void write_library_cutoffs(LibraryPipeline& lib,const std::vector<double>& cuts,const std::vector<std::string>& directories,const std::vector<std::string>& labels,int threads,bool debug,bool timings,bool shadow,bool consume_records){
    if(cuts.size()!=directories.size()||cuts.size()!=labels.size())throw std::runtime_error("cutoff cardinality");
    using Clock=std::chrono::steady_clock;
    const bool measure=timings||debug;const auto stream_begin=Clock::now();
    const auto elapsed=[](auto begin){return std::chrono::duration<double>(Clock::now()-begin).count();};
    std::vector<OutputFiles> files(cuts.size());
    std::vector<std::array<std::size_t,5>> outcome_counts(cuts.size());
    const auto each_file=[&](auto work,const char* label){
        if(threads>1&&lib.pool){
            ordered_index_chunks(*lib.pool,files.size()*5,threads,0,
                [&](std::size_t begin,std::size_t n){for(std::size_t i=begin;i<begin+n;++i)if(i%5!=4||lib.outcomes_mode!=OutcomesMode::none)work(i/5,i%5);return 0;},
                [](int){},label,1);
        }else for(std::size_t c=0;c<files.size();++c)for(std::size_t f=0;f<5;++f)if(f!=4||lib.outcomes_mode!=OutcomesMode::none)work(c,f);
    };
    each_file([&](std::size_t c,std::size_t f){if(f==4&&lib.outcomes_mode==OutcomesMode::summary){files[c][f].open(std::filesystem::path(directories[c])/"outcomes.summary.tsv");if(!files[c][f])throw std::runtime_error("open output");}else open_file(files[c][f],directories[c],f);},"output_open");
    const double open_seconds=measure?elapsed(stream_begin):0;double write_seconds=0;
    std::unique_ptr<ParallelOutputFiles> file_writers;
    if(threads>=8&&cuts.size()==1)file_writers=std::make_unique<ParallelOutputFiles>(files[0],measure);
    std::vector<RepairShadow> counters(cuts.size());std::vector<double> seconds(cuts.size());
    double router_seconds=0,native_seconds=0;std::size_t count=0,native_clean=0,native_direct=0;
    std::vector<std::size_t> native_clean_pass(cuts.size());
    process_records(lib,cuts,threads,debug,shadow?counters.data():nullptr,[&](OutputChunk chunk){
        count+=chunk.count;router_seconds+=chunk.router_seconds;native_seconds+=chunk.native_seconds;
        native_clean+=chunk.native_clean;native_direct+=chunk.native_direct;
        for(std::size_t c=0;c<cuts.size();++c)native_clean_pass[c]+=chunk.native_clean_pass[c];
        if(debug){lib.routing.routed_tsv+=chunk.routed;lib.routing.audit_tsv+=chunk.audit;lib.native.native_tsv+=chunk.native;}
        const auto write_begin=measure?Clock::now():Clock::time_point{};
        for(std::size_t c=0;c<cuts.size();++c){seconds[c]+=chunk.seconds[c];auto& out=chunk.outputs[c];
            if(lib.outcomes_mode==OutcomesMode::summary)for(std::size_t k=0;k<5;++k)outcome_counts[c][k]+=out.outcome_counts[k];
            if(file_writers)file_writers->append({&out.merged,&out.r1,&out.r2,&out.orphan,&out.outcomes});
            else write_files(files[c],out);}
        if(measure)write_seconds+=elapsed(write_begin);
    });
    double enqueue_seconds=0,drain_seconds=0;
    if(file_writers){
        const auto drain_begin=measure?Clock::now():Clock::time_point{};
        file_writers->finish();drain_seconds=measure?elapsed(drain_begin):0;
        enqueue_seconds=write_seconds;write_seconds=file_writers->seconds();
        if(memory_enabled()){MemoryUsage memory;memory.output=file_writers->text_capacity_bound();memory.report("parallel_output_text_bound",lib.pair_count(),"queued_text_capacity_only");}
    }
    if(count!=lib.pair_count())throw std::runtime_error("final cardinality");
    std::future<double> reclaimed;double release_seconds=0;
    if(consume_records){
        // process_records joined every reader and writer. The task owns all
        // released storage, including on output-close failure. Debug text stays.
        auto release=[compact=std::move(lib.compact),pairs=std::move(lib.family.pairs),native=std::move(lib.native.records),
                      replacements=std::move(lib.routing.record_replacements),measure]() mutable {
            const auto begin=measure?Clock::now():Clock::time_point{};
            std::vector<NativeRecord>().swap(replacements);
            std::vector<NativeRecord>().swap(native);
            std::vector<FamilyPair>().swap(pairs);compact.reset();
            return measure?std::chrono::duration<double>(Clock::now()-begin).count():0;
        };
        if(threads>1&&lib.pool)reclaimed=lib.pool->submit(std::move(release));
        else release_seconds=release();
    }
    const auto close_begin=measure?Clock::now():Clock::time_point{};
    if(lib.outcomes_mode==OutcomesMode::summary){
        static const char* kinds[]={"merged","unmerged","orphan_r1","orphan_r2","dropped"};
        for(std::size_t c=0;c<cuts.size();++c){
            std::string text="outcome\tcount\n";std::size_t total=0;
            for(std::size_t k=0;k<5;++k){text+=kinds[k];text+='\t';text+=std::to_string(outcome_counts[c][k]);text+='\n';total+=outcome_counts[c][k];}
            text+="total\t"+std::to_string(total)+"\n";files[c][4].append(text);
        }
    }
    each_file([&](std::size_t c,std::size_t f){files[c][f].close();if(!files[c][f])throw std::runtime_error("write output");},"output_close");
    const double close_seconds=measure?elapsed(close_begin):0;
    if(reclaimed.valid())release_seconds=reclaimed.get();
    if(measure){
        if(consume_records)std::cerr<<"MERGE_STAGE release_input seconds="<<release_seconds<<" timing_basis=worker_elapsed_overlaps_output_close\n";
        std::uint64_t writeback_bytes=0;std::size_t disabled=0;double hint_seconds=0;
        for(const auto& cut:files)for(const auto& file:cut){writeback_bytes+=file.writeback_bytes();hint_seconds+=file.writeback_seconds();disabled+=!file.writeback_enabled();}
        std::cerr<<"MERGE_WRITEBACK submitted_bytes="<<writeback_bytes<<" disabled_streams="<<disabled<<" hint_seconds="<<hint_seconds<<" timing_basis=caller_elapsed"<<'\n';
        const double output_seconds=open_seconds+write_seconds+close_seconds;
        const double stream_seconds=elapsed(stream_begin);
        std::cerr<<"MERGE_STAGE output_open seconds="<<open_seconds<<" timing_basis=wall\n";
        std::cerr<<"MERGE_STAGE output_write seconds="<<write_seconds<<" timing_basis="<<(file_writers?"file_writer_elapsed_sum":"writer_elapsed")<<'\n';
        if(file_writers){
            std::cerr<<"MERGE_STAGE output_enqueue seconds="<<enqueue_seconds<<" timing_basis=ordered_writer_elapsed\n";
            std::cerr<<"MERGE_STAGE output_drain seconds="<<drain_seconds<<" timing_basis=wall_after_producer_join\n";
            std::cerr<<"MERGE_OUTPUT_WRITERS workers=5 queue_limit=16 text_capacity_bound="<<file_writers->text_capacity_bound()<<'\n';
        }
        std::cerr<<"MERGE_STAGE output_close seconds="<<close_seconds<<" timing_basis=wall\n";
        std::cerr<<"MERGE_STAGE output_io seconds="<<output_seconds<<" timing_basis="<<(file_writers?"file_writer_elapsed_sum_plus_open_close_wall":"writer_elapsed_plus_open_close_wall")<<memory_status()<<'\n';
        std::cerr<<"MERGE_STAGE streaming_wall seconds="<<stream_seconds<<" timing_basis=wall_including_output"<<memory_status()<<'\n';
    }
    if(timings||debug)std::cerr<<"MERGE_STAGE streaming_native seconds="<<native_seconds<<" timing_basis=worker_elapsed"<<memory_status()<<'\n';
    if(timings||debug)std::cerr<<"MERGE_STAGE router seconds="<<router_seconds<<" timing_basis=worker_elapsed"<<memory_status()<<'\n';
    for(std::size_t c=0;c<cuts.size();++c){
        if(timings||debug)std::cerr<<"CLEAN_FIRST cutoff="<<labels[c]<<" native_certified="<<native_clean<<" bound_pass="<<native_clean_pass[c]<<" pairs="<<count<<" reference_checks="<<(debug||shadow)<<'\n';
        if(timings||debug)std::cerr<<"STRUCTURAL_SKIP cutoff="<<labels[c]<<" direct_clean_rows="<<native_direct<<" pairs="<<count<<" fraction="<<(count?double(native_direct)/count:0)<<" reference_checks="<<(debug||shadow)<<'\n';
        if(shadow)std::cerr<<"REPAIR_SHADOW cutoff="<<labels[c]<<" skip_and_unchanged="<<counters[c].skip_and_unchanged.load()<<" skip_but_changed="<<counters[c].skip_but_changed.load()<<" noskip_unchanged="<<counters[c].noskip_unchanged.load()<<" noskip_changed="<<counters[c].noskip_changed.load()<<" clean_and_unchanged="<<counters[c].clean_and_unchanged.load()<<" clean_but_changed="<<counters[c].clean_but_changed.load()<<" notclean="<<counters[c].notclean.load()<<'\n';
        if(timings||debug)std::cerr<<"MERGE_STAGE repair_assemble_write cutoff="<<labels[c]<<" seconds="<<seconds[c]<<" timing_basis=worker_elapsed_excluding_write"<<memory_status()<<'\n';
    }
}
}
