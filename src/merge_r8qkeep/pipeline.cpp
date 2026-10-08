#include "pipeline.hpp"
#include "canonical_pairs.hpp"
#include "prior_sample.hpp"
#include "ordered_chunks.hpp"
#include "memory_status.hpp"
#include "string_input.hpp"
#include <algorithm>
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
std::vector<std::string> fields(const std::string& line){std::vector<std::string> out;std::size_t begin=0;for(;;){const auto end=line.find('\t',begin);out.push_back(line.substr(begin,end==std::string::npos?end:end-begin));if(end==std::string::npos)return out;begin=end+1;}}
std::vector<std::string> rows(const std::string& text){StringInputBuffer buffer(text);std::istream input(&buffer);std::vector<std::string> out;std::string line;while(std::getline(input,line))out.push_back(std::move(line));return out;}
std::string id(const std::string& h){auto space=h.find_first_of(" \t");auto s=h.substr(1,space==std::string::npos?space:space-1);if(s.size()>2&&s[s.size()-2]=='/'&&(s.back()=='1'||s.back()=='2'))s.resize(s.size()-2);return s;}
void json_string(std::string& out,const std::string& s){static const char hex[]="0123456789abcdef";out+='"';for(unsigned char c:s){switch(c){case '"':out+="\\\"";break;case '\\':out+="\\\\";break;case '\b':out+="\\b";break;case '\f':out+="\\f";break;case '\n':out+="\\n";break;case '\r':out+="\\r";break;case '\t':out+="\\t";break;default:if(c<32){out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}else out+=char(c);}}out+='"';}
}
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug,bool timings,bool defer_routing,bool counts){
    auto start=std::chrono::steady_clock::now();
    const auto stage=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(timings||debug)std::cerr<<"MERGE_STAGE "<<name<<" seconds="<<std::chrono::duration<double>(now-start).count()<<memory_status()<<'\n';start=now;};
    if(counts){canonicalize_family_input(r1,r2);stage("canonical_input");}
    LibraryPipeline out;out.family=correct_family(std::move(r1),std::move(r2),threads);stage("family");
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
    const int concurrent=std::max(1,std::min(4,threads)),per_fit=std::max(1,threads/concurrent);
    for(int begin=0;begin<4;begin+=concurrent){std::array<std::future<std::optional<MixtureFit>>,4> futures;
        for(int k=begin;k<std::min(4,begin+concurrent);++k){const int mode=k/2,f=k%2;std::vector<const PairLikelihood*> training;for(const auto& p:pairs[mode])if(p.fold==f)training.push_back(&p);
            if(concurrent>1)futures[k]=std::async(std::launch::async,[training=std::move(training),per_fit,counts](){return counts?fit_mixture_counts(training,per_fit):fit_mixture_views(training,per_fit);});else fits[mode][f]=counts?fit_mixture_counts(training,per_fit):fit_mixture_views(training,per_fit);}
        if(concurrent>1)for(int k=begin;k<std::min(4,begin+concurrent);++k)fits[k/2][k%2]=futures[k].get();
    }
    out.nominal_prior=prior_text(fits[0]);out.full_prior=prior_text(fits[1]);stage("prior_fit");
    if(timings)for(int mode=0;mode<2;++mode)for(int f=0;f<2;++f)if(fits[mode][f]){const auto& fit=*fits[mode][f];std::cerr<<"FIT_CACHE mode="<<mode<<" fold="<<f<<" pairs="<<fit.training_pairs<<" keys="<<fit.unique_likelihoods<<" calls="<<fit.exp_calls<<" hits="<<fit.exp_cache_hits<<'\n';}
    for(auto& observations:pairs)std::vector<PairLikelihood>().swap(observations);
    if(counts){out.native=merge_native(out.family.pairs,settings,timings||debug,debug);stage("native");}
    if(!debug){std::string().swap(out.family.edits_tsv);std::string().swap(out.family.families_tsv);std::string().swap(out.observations.proposed_edits_tsv);}
    out.repair=std::make_unique<RepairPipeline>(out.nominal_prior,out.full_prior);if(!defer_routing){out.routing=out.repair->route(out.native.native_tsv,debug,threads);stage("router");}return out;
}
FinalOutput emit_library(const LibraryPipeline& lib,double cutoff){
    const auto base=rows(lib.repair->repair(lib.native.native_tsv,9,cutoff,0,1,nullptr,&lib.routing.routed,&lib.routing.repair_bounds));const bool rerouted=std::find(lib.routing.routed.begin(),lib.routing.routed.end(),true)!=lib.routing.routed.end();
    std::string sparse_table;if(rerouted&&lib.routing.routed_tsv.empty()){StringInputBuffer buffer(lib.native.native_tsv);std::istream input(&buffer);std::string line;std::size_t index=0;while(std::getline(input,line)){sparse_table+=lib.routing.routed.at(index)?lib.routing.replacements.at(index):line+"\n";++index;}}
    const auto alternative=rerouted?rows(lib.repair->repair(lib.routing.routed_tsv.empty()?sparse_table:lib.routing.routed_tsv,7,cutoff)):std::vector<std::string>{};
    if(base.size()!=lib.family.pairs.size()||lib.routing.routed.size()!=base.size()||(rerouted&&alternative.size()!=base.size()))throw std::runtime_error("stage cardinality");
    FinalOutput out;for(std::size_t i=0;i<base.size();++i){const auto v=fields(lib.routing.routed[i]?alternative[i]:base[i]);if(v.size()!=8||std::stoull(v[0])!=i)throw std::runtime_error("stage order");const auto& pair=lib.family.pairs[i];const auto add=[&](int j,std::string& text){if(id(v[j])!=pair.id)throw std::runtime_error("stage record identity");text+=v[j]+"\n"+v[j+1]+"\n+\n"+v[j+2]+"\n";};if(v[1]=="merged")add(2,out.merged);else if(v[1]=="unmerged"){add(2,out.r1);add(5,out.r2);}else if(v[1]=="orphan_r1"||v[1]=="orphan_r2")add(2,out.orphan);else if(v[1]!="dropped")throw std::runtime_error("stage outcome");out.outcomes+="{\"r1_header\":";json_string(out.outcomes,pair.r[0].h);out.outcomes+=",\"r2_header\":";json_string(out.outcomes,pair.r[1].h);out.outcomes+=",\"outcome\":";json_string(out.outcomes,v[1]);out.outcomes+="}\n";}return out;
}
void write_library_outputs(const LibraryPipeline& lib,double cutoff,const std::string& directory,int threads,RepairShadow* shadow){
    const bool profile=std::getenv("FQDUP_CHUNK_PROFILE")!=nullptr;
    using Clock=std::chrono::steady_clock;double prepare_seconds=0,repair_seconds=0,assemble_seconds=0;std::size_t chunks=0;
    const auto elapsed=[](auto t){return std::chrono::duration<double>(Clock::now()-t).count();};
    const std::filesystem::path path=directory;
    std::array<std::ofstream,5> files;const char* names[]={"merged.fq","r1.fq","r2.fq","orphan.fq","outcomes.jsonl"};
    for(int i=0;i<5;++i){files[i].open(path/names[i],std::ios::binary);if(!files[i])throw std::runtime_error("open output");}
    std::size_t position=0,start=0;
    while(position<lib.native.native_tsv.size()){
        const auto prepare_start=Clock::now();
        std::string native,routed;std::size_t count=0;bool alternative=false;
        while(count<2048&&position<lib.native.native_tsv.size()){
            const auto end=lib.native.native_tsv.find('\n',position);if(end==std::string::npos)throw std::runtime_error("unterminated native row");
            const auto line=lib.native.native_tsv.substr(position,end-position+1);native+=line;
            const auto index=start+count;if(index>=lib.routing.routed.size())throw std::runtime_error("route cardinality");
            if(lib.routing.routed[index]){alternative=true;routed+=lib.routing.replacements.at(index);}else routed+=line;
            position=end+1;++count;
        }
        if(profile)prepare_seconds+=elapsed(prepare_start);const auto repair_start=Clock::now();
        const auto base=rows(lib.repair->repair(native,9,cutoff,start,threads,shadow,&lib.routing.routed,&lib.routing.repair_bounds));
        const auto changed=alternative?rows(lib.repair->repair(routed,7,cutoff,start,threads)):std::vector<std::string>{};
        if(profile)repair_seconds+=elapsed(repair_start);const auto assemble_start=Clock::now();
        if(base.size()!=count||(alternative&&changed.size()!=count))throw std::runtime_error("chunk cardinality");
        for(std::size_t i=0;i<count;++i){const auto index=start+i;const auto v=fields(lib.routing.routed[index]?changed[i]:base[i]);if(v.size()!=8||std::stoull(v[0])!=index)throw std::runtime_error("chunk order");
            const auto& pair=lib.family.pairs.at(index);const auto add=[&](int j,int channel){if(id(v[j])!=pair.id)throw std::runtime_error("record identity");files[channel]<<v[j]<<'\n'<<v[j+1]<<"\n+\n"<<v[j+2]<<'\n';};
            if(v[1]=="merged")add(2,0);else if(v[1]=="unmerged"){add(2,1);add(5,2);}else if(v[1]=="orphan_r1"||v[1]=="orphan_r2")add(2,3);else if(v[1]!="dropped")throw std::runtime_error("chunk outcome");
            std::string outcome="{\"r1_header\":";json_string(outcome,pair.r[0].h);outcome+=",\"r2_header\":";json_string(outcome,pair.r[1].h);outcome+=",\"outcome\":";json_string(outcome,v[1]);outcome+="}\n";files[4]<<outcome;
        }
        start+=count;
        if(profile)assemble_seconds+=elapsed(assemble_start);++chunks;
    }
    if(start!=lib.family.pairs.size())throw std::runtime_error("final cardinality");
    for(auto& file:files){file.close();if(!file)throw std::runtime_error("write output");}
    if(profile)std::fprintf(stderr,"EMIT_PROFILE threads=%d rows=%zu chunks=%zu chunk_rows=2048 prepare=%.6f repair_parse=%.6f ordered_assembly_write=%.6f\n",threads,start,chunks,prepare_seconds,repair_seconds,assemble_seconds);
}
}

namespace fqdup::r8qkeep {
void write_library_cutoffs(LibraryPipeline& lib,const std::vector<double>& cuts,const std::vector<std::string>& directories,const std::vector<std::string>& labels,int threads,bool debug,bool timings,bool shadow) {
    if(cuts.size()!=directories.size()||cuts.size()!=labels.size())throw std::runtime_error("cutoff cardinality");
    using Clock=std::chrono::steady_clock;
    const auto elapsed=[](auto begin){return std::chrono::duration<double>(Clock::now()-begin).count();};
    std::vector<std::array<std::ofstream,5>> files(cuts.size());
    const char* names[]={"merged.fq","r1.fq","r2.fq","orphan.fq","outcomes.jsonl"};
    for(std::size_t c=0;c<cuts.size();++c)for(int j=0;j<5;++j){files[c][j].open(std::filesystem::path(directories[c])/names[j],std::ios::binary);if(!files[c][j])throw std::runtime_error("open output");}
    std::vector<RepairShadow> counters(cuts.size());
    struct Chunk { std::vector<FinalOutput> outputs;std::vector<double> seconds;double router_seconds=0;std::string routed,audit;std::size_t count=0; };
    std::vector<double> seconds(cuts.size());double router_seconds=0;std::size_t count=0;
    ordered_text_chunks(lib.native.native_tsv,threads,0,[&](const std::string& native,std::size_t start){
        Chunk chunk;chunk.count=std::count(native.begin(),native.end(),'\n');chunk.outputs.resize(cuts.size());chunk.seconds.resize(cuts.size());
        // Each worker owns this chunk through every cutoff. Global row identity
        // is restored at assembly; route/repair masks and cache use local rows.
        RepairAlignmentCache cache(0,chunk.count);const auto route_begin=Clock::now();
        const auto routing=lib.repair->route(native,debug,1,0,&cache);chunk.router_seconds=elapsed(route_begin);
        std::string routed;bool alternative=false;const auto native_rows=rows(native);
        for(std::size_t i=0;i<chunk.count;++i){if(routing.routed.at(i)){alternative=true;routed+=routing.replacements.at(i);}else routed+=native_rows[i]+"\n";}
        if(debug){chunk.routed=routing.routed_tsv;for(const auto& line:rows(routing.audit_tsv)){auto tab=line.find('\t');if(tab==std::string::npos)throw std::runtime_error("audit row");chunk.audit+=std::to_string(start+std::stoull(line.substr(0,tab)))+line.substr(tab)+"\n";}}
        for(std::size_t c=0;c<cuts.size();++c){const auto begin=Clock::now();auto& out=chunk.outputs[c];
            const auto base=rows(lib.repair->repair(native,9,cuts[c],0,1,shadow?&counters[c]:nullptr,&routing.routed,&routing.repair_bounds,&cache));
            const auto changed=alternative?rows(lib.repair->repair(routed,7,cuts[c],0,1,nullptr,nullptr,nullptr,&cache)):std::vector<std::string>{};
            if(base.size()!=chunk.count||(alternative&&changed.size()!=chunk.count))throw std::runtime_error("chunk cardinality");
            for(std::size_t i=0;i<chunk.count;++i){const auto v=fields(routing.routed[i]?changed[i]:base[i]);if(v.size()!=8||std::stoull(v[0])!=i)throw std::runtime_error("chunk order");
                const auto& pair=lib.family.pairs.at(start+i);const auto add=[&](int j,std::string& text){if(id(v[j])!=pair.id)throw std::runtime_error("record identity");text+=v[j]+"\n"+v[j+1]+"\n+\n"+v[j+2]+"\n";};
                if(v[1]=="merged")add(2,out.merged);else if(v[1]=="unmerged"){add(2,out.r1);add(5,out.r2);}else if(v[1]=="orphan_r1"||v[1]=="orphan_r2")add(2,out.orphan);else if(v[1]!="dropped")throw std::runtime_error("chunk outcome");
                out.outcomes+="{\"r1_header\":";json_string(out.outcomes,pair.r[0].h);out.outcomes+=",\"r2_header\":";json_string(out.outcomes,pair.r[1].h);out.outcomes+=",\"outcome\":";json_string(out.outcomes,v[1]);out.outcomes+="}\n";
            }
            chunk.seconds[c]=elapsed(begin);
        }
        return chunk;
    },[&](Chunk chunk){count+=chunk.count;router_seconds+=chunk.router_seconds;
        if(debug){lib.routing.routed_tsv+=chunk.routed;lib.routing.audit_tsv+=chunk.audit;}
        for(std::size_t c=0;c<cuts.size();++c){seconds[c]+=chunk.seconds[c];const auto& out=chunk.outputs[c];files[c][0]<<out.merged;files[c][1]<<out.r1;files[c][2]<<out.r2;files[c][3]<<out.orphan;files[c][4]<<out.outcomes;}
    },"all_cutoffs");
    if(count!=lib.family.pairs.size())throw std::runtime_error("final cardinality");
    for(auto& group:files)for(auto& file:group){file.close();if(!file)throw std::runtime_error("write output");}
    if(timings||debug)std::cerr<<"MERGE_STAGE router seconds="<<router_seconds<<" timing_basis=worker_elapsed"<<memory_status()<<'\n';
    for(std::size_t c=0;c<cuts.size();++c){if(shadow)std::cerr<<"REPAIR_SHADOW cutoff="<<labels[c]<<" skip_and_unchanged="<<counters[c].skip_and_unchanged.load()<<" skip_but_changed="<<counters[c].skip_but_changed.load()<<" noskip_unchanged="<<counters[c].noskip_unchanged.load()<<" noskip_changed="<<counters[c].noskip_changed.load()<<'\n';
        if(timings||debug)std::cerr<<"MERGE_STAGE repair_assemble_write cutoff="<<labels[c]<<" seconds="<<seconds[c]<<" timing_basis=worker_elapsed_excluding_write"<<memory_status()<<'\n';}
}
}
