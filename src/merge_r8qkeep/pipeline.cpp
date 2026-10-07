#include "pipeline.hpp"
#include "memory_status.hpp"
#include <algorithm>
#include <array>
#include <future>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <iostream>
namespace fqdup::r8qkeep {
namespace {
std::vector<std::string> fields(const std::string& line){std::vector<std::string> out;std::size_t begin=0;for(;;){const auto end=line.find('\t',begin);out.push_back(line.substr(begin,end==std::string::npos?end:end-begin));if(end==std::string::npos)return out;begin=end+1;}}
std::vector<std::string> rows(const std::string& text){std::istringstream input(text);std::vector<std::string> out;std::string line;while(std::getline(input,line))out.push_back(std::move(line));return out;}
std::string id(const std::string& h){auto space=h.find_first_of(" \t");auto s=h.substr(1,space==std::string::npos?space:space-1);if(s.size()>2&&s[s.size()-2]=='/'&&(s.back()=='1'||s.back()=='2'))s.resize(s.size()-2);return s;}
void json_string(std::string& out,const std::string& s){static const char hex[]="0123456789abcdef";out+='"';for(unsigned char c:s){switch(c){case '"':out+="\\\"";break;case '\\':out+="\\\\";break;case '\b':out+="\\b";break;case '\f':out+="\\f";break;case '\n':out+="\\n";break;case '\r':out+="\\r";break;case '\t':out+="\\t";break;default:if(c<32){out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}else out+=char(c);}}out+='"';}
}
LibraryPipeline prepare_library(std::vector<FamilyRead> r1,std::vector<FamilyRead> r2,const NativeSettings& settings,int threads,bool debug,bool timings){
    auto start=std::chrono::steady_clock::now();
    const auto stage=[&](const char* name){const auto now=std::chrono::steady_clock::now();if(timings||debug)std::cerr<<"MERGE_STAGE "<<name<<" seconds="<<std::chrono::duration<double>(now-start).count()<<memory_status()<<'\n';start=now;};
    LibraryPipeline out;out.family=correct_family(std::move(r1),std::move(r2));stage("family");out.native=merge_native(out.family.pairs,settings,timings||debug);stage("native");out.observations=census(out.family.pairs,out.native,debug);stage("census");
    std::array<std::vector<PairLikelihood>,2> pairs{std::move(out.observations.nominal_pairs),std::move(out.observations.full_pairs)};
    std::array<std::array<std::optional<MixtureFit>,2>,2> fits;
    const int concurrent=std::max(1,std::min(4,threads)),per_fit=std::max(1,threads/concurrent);
    for(int begin=0;begin<4;begin+=concurrent){std::array<std::future<std::optional<MixtureFit>>,4> futures;
        for(int k=begin;k<std::min(4,begin+concurrent);++k){const int mode=k/2,f=k%2;std::vector<PairLikelihood> training;for(const auto& p:pairs[mode])if(p.fold==f)training.push_back(p);
            if(concurrent>1)futures[k]=std::async(std::launch::async,[training=std::move(training),per_fit](){return fit_mixture(training,per_fit);});else fits[mode][f]=fit_mixture(training,per_fit);}
        if(concurrent>1)for(int k=begin;k<std::min(4,begin+concurrent);++k)fits[k/2][k%2]=futures[k].get();
    }
    out.nominal_prior=prior_text(fits[0]);out.full_prior=prior_text(fits[1]);stage("prior_fit");
    if(timings)for(int mode=0;mode<2;++mode)for(int f=0;f<2;++f)if(fits[mode][f]){const auto& fit=*fits[mode][f];std::cerr<<"FIT_CACHE mode="<<mode<<" fold="<<f<<" pairs="<<fit.training_pairs<<" keys="<<fit.unique_likelihoods<<" calls="<<fit.exp_calls<<" hits="<<fit.exp_cache_hits<<'\n';}
    out.repair=std::make_unique<RepairPipeline>(out.nominal_prior,out.full_prior);out.routing=out.repair->route(out.native.native_tsv);stage("router");return out;
}
FinalOutput emit_library(const LibraryPipeline& lib,double cutoff){
    const auto base=rows(lib.repair->repair(lib.native.native_tsv,9,cutoff));const bool rerouted=std::find(lib.routing.routed.begin(),lib.routing.routed.end(),true)!=lib.routing.routed.end();const auto alternative=rerouted?rows(lib.repair->repair(lib.routing.routed_tsv,7,cutoff)):std::vector<std::string>{};
    if(base.size()!=lib.family.pairs.size()||lib.routing.routed.size()!=base.size()||(rerouted&&alternative.size()!=base.size()))throw std::runtime_error("stage cardinality");
    FinalOutput out;for(std::size_t i=0;i<base.size();++i){const auto v=fields(lib.routing.routed[i]?alternative[i]:base[i]);if(v.size()!=8||std::stoull(v[0])!=i)throw std::runtime_error("stage order");const auto& pair=lib.family.pairs[i];const auto add=[&](int j,std::string& text){if(id(v[j])!=pair.id)throw std::runtime_error("stage record identity");text+=v[j]+"\n"+v[j+1]+"\n+\n"+v[j+2]+"\n";};if(v[1]=="merged")add(2,out.merged);else if(v[1]=="unmerged"){add(2,out.r1);add(5,out.r2);}else if(v[1]=="orphan_r1"||v[1]=="orphan_r2")add(2,out.orphan);else if(v[1]!="dropped")throw std::runtime_error("stage outcome");out.outcomes+="{\"r1_header\":";json_string(out.outcomes,pair.r[0].h);out.outcomes+=",\"r2_header\":";json_string(out.outcomes,pair.r[1].h);out.outcomes+=",\"outcome\":";json_string(out.outcomes,v[1]);out.outcomes+="}\n";}return out;
}
}
