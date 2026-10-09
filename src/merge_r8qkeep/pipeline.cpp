#include "pipeline.hpp"
#include "canonical_pairs.hpp"
#include "prior_sample.hpp"
#include "ordered_chunks.hpp"
#include "memory_status.hpp"
#include "string_input.hpp"
#include <algorithm>
#include <string_view>
#include <array>
#include <future>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
namespace fqdup::r8qkeep {
namespace {
void json_string(std::string& out,const std::string& s){static const char hex[]="0123456789abcdef";out+='"';for(unsigned char c:s){switch(c){case '"':out+="\\\"";break;case '\\':out+="\\\\";break;case '\b':out+="\\b";break;case '\f':out+="\\f";break;case '\n':out+="\\n";break;case '\r':out+="\\r";break;case '\t':out+="\\t";break;default:if(c<32){out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}else out+=char(c);}}out+='"';}
}
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug,bool timings,bool defer_routing,bool counts){
    auto start=std::chrono::steady_clock::now();
    const auto stage=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(timings||debug)std::cerr<<"MERGE_STAGE "<<name<<" seconds="<<std::chrono::duration<double>(now-start).count()<<memory_status()<<'\n';start=now;};
    if(counts){canonicalize_family_input(r1,r2);stage("canonical_input");}
    LibraryPipeline out;out.native_settings=settings;out.family=correct_family(std::move(r1),std::move(r2),threads,debug);stage("family");
    if(counts){
        constexpr std::size_t requested=25000;
        const auto selected=prior_sample_indices(out.family.pairs,requested);
        std::vector<FamilyPair> sample;sample.reserve(selected.size());
        for(auto i:selected)sample.push_back(out.family.pairs[i]);
        if(timings||debug)std::cerr<<"PRIOR_SAMPLE requested="<<requested<<" selected="<<sample.size()<<" total="<<out.family.pairs.size()<<" seed="<<prior_sample_seed<<'\n';
        const auto native=merge_native(sample,settings,timings||debug,debug);stage("sample_native");
        out.observations=census(sample,native,debug,threads);stage("sample_census");
    }else{
        out.native=merge_native(out.family.pairs,settings,timings||debug,debug);stage("native");
        out.observations=census(out.family.pairs,out.native,debug,threads);stage("census");
    }

    std::array<std::vector<PairLikelihood>,2> pairs{std::move(out.observations.nominal_pairs),std::move(out.observations.full_pairs)};
    std::array<std::array<std::optional<MixtureFit>,2>,2> fits;
    // Counts fits use two concurrent folds and split the thread budget between
    // their persistent workers. Default full-fit scheduling is unchanged.
    const int concurrent=std::max(1,std::min(counts?2:4,threads));
    const int per_fit=std::max(1,threads/concurrent);
    for(int begin=0;begin<4;begin+=concurrent){std::array<std::future<std::optional<MixtureFit>>,4> futures;
        for(int k=begin;k<std::min(4,begin+concurrent);++k){const int mode=k/2,f=k%2;std::vector<const PairLikelihood*> training;for(const auto& p:pairs[mode])if(p.fold==f)training.push_back(&p);
            if(concurrent>1)futures[k]=std::async(std::launch::async,[training=std::move(training),per_fit,counts](){return counts?fit_mixture_counts(training,per_fit):fit_mixture_views(training,per_fit);});else fits[mode][f]=counts?fit_mixture_counts(training,per_fit):fit_mixture_views(training,per_fit);}
        if(concurrent>1)for(int k=begin;k<std::min(4,begin+concurrent);++k)fits[k/2][k%2]=futures[k].get();
    }
    out.nominal_prior=prior_text(fits[0]);out.full_prior=prior_text(fits[1]);stage("prior_fit");
    if(timings)for(int mode=0;mode<2;++mode)for(int f=0;f<2;++f)if(fits[mode][f]){const auto& fit=*fits[mode][f];std::cerr<<"FIT_CACHE mode="<<mode<<" fold="<<f<<" pairs="<<fit.training_pairs<<" keys="<<fit.unique_likelihoods<<" calls="<<fit.exp_calls<<" hits="<<fit.exp_cache_hits<<'\n';std::cerr<<"FIT_PHASE mode="<<mode<<" fold="<<f<<" em_iterations="<<fit.iterations<<" em_seconds="<<fit.em_seconds<<" polish_iterations="<<fit.polish_iterations<<" polish_seconds="<<fit.polish_seconds<<" last_em_zero_cells="<<fit.last_em_zero_cells<<" last_em_subnormal_cells="<<fit.last_em_subnormal_cells<<" em_cells="<<fit.unique_likelihoods*mixture_components<<'\n';}
    for(auto& observations:pairs)std::vector<PairLikelihood>().swap(observations);
    if(counts&&!defer_routing){out.native=merge_native(out.family.pairs,settings,timings||debug,debug);stage("native");}
    if(defer_routing){NativeResult empty;empty.native_tsv=std::move(out.native.native_tsv);out.native=std::move(empty);}
    if(!debug){std::string().swap(out.family.edits_tsv);std::string().swap(out.family.families_tsv);std::string().swap(out.observations.proposed_edits_tsv);}
    out.repair=std::make_unique<RepairPipeline>(out.nominal_prior,out.full_prior);if(!defer_routing){out.routing=out.repair->route_records(out.native.records,debug);stage("router");}return out;
}
namespace {
std::string_view read_id(const std::string& h){
    if(h.empty()||h[0]!='@')throw std::runtime_error("record identity");
    const auto end=h.find_first_of(" \t");
    std::string_view view(h.data()+1,(end==std::string::npos?h.size():end)-1);
    if(view.size()>2&&view[view.size()-2]=='/'&&(view.back()=='1'||view.back()=='2'))view.remove_suffix(2);
    return view;
}
void append_record(FinalOutput& out,const RepairedRecord& record,const FamilyPair& pair){
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
};
template<class Consume>
void process_records(const LibraryPipeline& lib,const std::vector<double>& cuts,int threads,bool debug,RepairShadow* counters,Consume consume){
    using Clock=std::chrono::steady_clock;
    const auto elapsed=[](auto start){return std::chrono::duration<double>(Clock::now()-start).count();};
    WorkerPool pool(std::max(1,std::min(threads,8)));
    const bool native_debug=debug&&lib.native.native_tsv.empty();
    const auto process=[&](NativeRecords records,std::size_t offset,double native_seconds){
    bool first=true;
    ordered_index_chunks(pool,records.size(),threads,offset,[&](std::size_t start,std::size_t n){
        OutputChunk chunk;chunk.count=n;chunk.outputs.resize(cuts.size());chunk.seconds.resize(cuts.size());
        const auto native=records.slice(start-offset,n);
        if(native_debug)for(const auto& record:native)chunk.native+=native_record_tsv(record);
        RepairAlignmentCache cache(start,n);const auto route_begin=Clock::now();
        const auto routing=lib.repair->route_records(native,debug,start,&cache,counters!=nullptr);
        chunk.router_seconds=elapsed(route_begin);
        const bool alternative=std::find(routing.routed.begin(),routing.routed.end(),true)!=routing.routed.end();
        if(debug){chunk.routed=routing.routed_tsv;chunk.audit=routing.audit_tsv;}
        for(std::size_t c=0;c<cuts.size();++c){
            const auto begin=Clock::now();
            const auto base=lib.repair->repair_records(native,9,cuts[c],start,counters?&counters[c]:nullptr,&routing.routed,&routing.repair_bounds,&cache,&routing.clean_bounds);
            const auto changed=alternative?lib.repair->repair_records(native,7,cuts[c],start,nullptr,nullptr,nullptr,&cache,nullptr,&routing.routed,&routing.record_replacements):std::vector<RepairedRecord>{};
            if(base.size()!=n||(alternative&&changed.size()!=n))throw std::runtime_error("chunk cardinality");
            for(std::size_t i=0;i<n;++i)append_record(chunk.outputs[c],routing.routed[i]?changed[i]:base[i],lib.family.pairs.at(start+i));
            chunk.seconds[c]=elapsed(begin);
        }
        return chunk;
    },[&](OutputChunk chunk){if(first){chunk.native_seconds=native_seconds;first=false;}consume(std::move(chunk));},"all_cutoffs");
    };
    if(!lib.native.records.empty()){process(lib.native.records,0,0);return;}
    struct NativeChunk {NativeResult result;std::size_t start;double seconds;};
    std::deque<std::future<NativeChunk>> pending;
    std::size_t next=0;
    const auto launch=[&]{
        const auto start=next,n=std::min(native_batch_size,lib.family.pairs.size()-next);next+=n;
        pending.push_back(pool.submit([&,start,n]{const auto begin=Clock::now();auto result=merge_native_range(lib.family.pairs,start,n,lib.native_settings,false,false);return NativeChunk{std::move(result),start,elapsed(begin)};}));
    };
    while(next<lib.family.pairs.size()&&pending.size()<std::size_t(std::max(1,std::min(threads,8))))launch();
    while(!pending.empty()){
        auto chunk=pending.front().get();pending.pop_front();
        if(next<lib.family.pairs.size())launch();
        process(chunk.result.records,chunk.start,chunk.seconds);
    }
}
using OutputFiles=std::array<std::ofstream,5>;
void open_files(OutputFiles& files,const std::string& directory){
    const char* names[]={"merged.fq","r1.fq","r2.fq","orphan.fq","outcomes.jsonl"};
    for(int i=0;i<5;++i){files[i].open(std::filesystem::path(directory)/names[i],std::ios::binary);if(!files[i])throw std::runtime_error("open output");}
}
void write_files(OutputFiles& files,const FinalOutput& out){
    files[0]<<out.merged;files[1]<<out.r1;files[2]<<out.r2;files[3]<<out.orphan;files[4]<<out.outcomes;
}
void close_files(OutputFiles& files){for(auto& file:files){file.close();if(!file)throw std::runtime_error("write output");}}
}
FinalOutput emit_library(const LibraryPipeline& lib,double cutoff){
    FinalOutput result;std::size_t count=0;
    process_records(lib,{cutoff},1,false,nullptr,[&](OutputChunk chunk){
        count+=chunk.count;const auto& out=chunk.outputs[0];
        result.merged+=out.merged;result.r1+=out.r1;result.r2+=out.r2;result.orphan+=out.orphan;result.outcomes+=out.outcomes;
    });
    if(count!=lib.family.pairs.size())throw std::runtime_error("final cardinality");
    return result;
}
void write_library_outputs(const LibraryPipeline& lib,double cutoff,const std::string& directory,int threads,RepairShadow* shadow){
    OutputFiles files;open_files(files,directory);std::size_t count=0;
    process_records(lib,{cutoff},threads,false,shadow,[&](OutputChunk chunk){count+=chunk.count;write_files(files,chunk.outputs[0]);});
    if(count!=lib.family.pairs.size())throw std::runtime_error("final cardinality");
    close_files(files);
}
void write_library_cutoffs(LibraryPipeline& lib,const std::vector<double>& cuts,const std::vector<std::string>& directories,const std::vector<std::string>& labels,int threads,bool debug,bool timings,bool shadow){
    if(cuts.size()!=directories.size()||cuts.size()!=labels.size())throw std::runtime_error("cutoff cardinality");
    std::vector<OutputFiles> files(cuts.size());for(std::size_t c=0;c<cuts.size();++c)open_files(files[c],directories[c]);
    std::vector<RepairShadow> counters(cuts.size());std::vector<double> seconds(cuts.size());
    double router_seconds=0,native_seconds=0;std::size_t count=0;
    process_records(lib,cuts,threads,debug,shadow?counters.data():nullptr,[&](OutputChunk chunk){
        count+=chunk.count;router_seconds+=chunk.router_seconds;native_seconds+=chunk.native_seconds;
        if(debug){lib.routing.routed_tsv+=chunk.routed;lib.routing.audit_tsv+=chunk.audit;lib.native.native_tsv+=chunk.native;}
        for(std::size_t c=0;c<cuts.size();++c){seconds[c]+=chunk.seconds[c];write_files(files[c],chunk.outputs[c]);}
    });
    if(count!=lib.family.pairs.size())throw std::runtime_error("final cardinality");
    for(auto& group:files)close_files(group);
    if(timings||debug)std::cerr<<"MERGE_STAGE streaming_native seconds="<<native_seconds<<" timing_basis=worker_elapsed"<<memory_status()<<'\n';
    if(timings||debug)std::cerr<<"MERGE_STAGE router seconds="<<router_seconds<<" timing_basis=worker_elapsed"<<memory_status()<<'\n';
    for(std::size_t c=0;c<cuts.size();++c){
        if(shadow)std::cerr<<"REPAIR_SHADOW cutoff="<<labels[c]<<" skip_and_unchanged="<<counters[c].skip_and_unchanged.load()<<" skip_but_changed="<<counters[c].skip_but_changed.load()<<" noskip_unchanged="<<counters[c].noskip_unchanged.load()<<" noskip_changed="<<counters[c].noskip_changed.load()<<" clean_and_unchanged="<<counters[c].clean_and_unchanged.load()<<" clean_but_changed="<<counters[c].clean_but_changed.load()<<" notclean="<<counters[c].notclean.load()<<'\n';
        if(timings||debug)std::cerr<<"MERGE_STAGE repair_assemble_write cutoff="<<labels[c]<<" seconds="<<seconds[c]<<" timing_basis=worker_elapsed_excluding_write"<<memory_status()<<'\n';
    }
}
}
