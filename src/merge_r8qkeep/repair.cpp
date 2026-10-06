#include "repair.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <numeric>
#include <utility>
#include <vector>
namespace fqdup_r8qkeep_repair_internal {
#include "repair_headers/quality_guard_repair.hpp"
struct Record { std::string header,seq,qual; };
struct MergeRecord {Record merged,unmerged2;bool is_merged=false,is_orphan=false,orphan_r1=false;};
struct Batch {std::vector<MergeRecord> records;int64_t n_merged=0,n_unmerged=0,n_orphan=0,n_orphan_r1=0,n_orphan_r2=0,n_dropped=0;};
std::vector<std::string> split(const std::string& s){std::vector<std::string> v;std::size_t p=0;for(;;){auto e=s.find('\t',p);v.push_back(s.substr(p,e==std::string::npos?e:e-p));if(e==std::string::npos)return v;p=e+1;}}
fqquality::PairMixture model(const std::string& text){std::istringstream input(text);return fqquality::PairMixture(input);}
}
namespace fqdup::r8qkeep {
using namespace ::fqdup_r8qkeep_repair_internal;
struct RepairPipeline::Impl {
    fqquality::PairMixture nominal,full;
    Impl(const std::string& a,const std::string& b):nominal(model(a)),full(model(b)){}
};
RepairPipeline::RepairPipeline(const std::string& nominal,const std::string& full):impl_(std::make_shared<Impl>(nominal,full)){}
RoutingResult RepairPipeline::route(const std::string& native_tsv) const {
    const auto& model=impl_->nominal;fqalign::QualityProcessor processor;std::istringstream input(native_tsv);std::ostringstream out,audit;RoutingResult result;std::string line;std::size_t row=0;
    while(std::getline(input,line)) {
        auto v=split(line);if(v.size()!=13)throw std::runtime_error("route columns");bool routed=false;double worst=0;std::string reason="native_nonmerged";int insert=0;
        if(v[6]=="merged") {
            insert=v[8].size();reason="geometry_abstain";const auto g=fqalign::quality_independent_geometry(v[1],v[4]);
            if(g.eligible&&g.insert==insert){auto p=processor.risks(v[1],v[2],v[4],v[5],insert);if(p.correspondence.eligible){
                fqquality::PairMixture::Audit post;const bool fitted=model.apply(p,v[1],v[2],v[4],v[5],1,&post);bool valid=true;const auto rc=fqalign::reverse_complement(v[4]);
                for(auto [i,j]:p.correspondence.invariant){const int k=int(v[4].size())-1-j;if(i+k!=insert-1){valid=false;break;}if(v[1][i]=='N'||rc[j]=='N'||v[1][i]==rc[j])continue;
                    const char chosen=v[8].at(i);double risk=1;
                    if(chosen==v[1][i]||chosen==rc[j]){const double e1=fqalign::phred_error(v[2][i]),e2=fqalign::phred_error(v[5][k]);auto cp=fqquality::PairMixture::conditional(e1,e2,false);double a=cp.first,b=cp.second;
                        if(fitted){a=b=0;for(int z=0;z<fqquality::PairMixture::K;++z){const double e=fqquality::PairMixture::component(z);cp=fqquality::PairMixture::conditional(z==0?e1:e,z==0?e2:e,false);a+=post.posterior[z]*cp.first;b+=post.posterior[z]*cp.second;}}
                        risk=chosen==v[1][i]?a:b;
                    }else if(chosen!='N'){valid=false;break;}worst=std::max(worst,risk);
                }
                reason=valid?"certain_merge":"correspondence_abstain";
                if(valid&&worst>=.49){routed=true;reason="uncertain_merge";v[6]="unmerged";v[7]=v[0];v[8]=v[1].substr(0,insert);v[9]=v[2].substr(0,insert);v[10]=v[3];v[11]=v[4].substr(0,insert);v[12]=v[5].substr(0,insert);}
            }}
        }
        audit<<row++<<'\t'<<routed<<'\t'<<insert<<'\t'<<worst<<'\t'<<reason<<'\n';for(std::size_t k=0;k<v.size();++k){if(k)out<<'\t';out<<v[k];}out<<'\n';result.routed.push_back(routed);
    }
    result.routed_tsv=out.str();result.audit_tsv=audit.str();return result;
}
std::string RepairPipeline::repair(const std::string& input_tsv,int mode,double cut) const {
    if((mode!=7&&mode!=9)||!std::isfinite(cut)||cut<0||cut>1)throw std::runtime_error("invalid repair arguments");
    const auto& model=impl_->nominal;const auto& full_model=impl_->full;fqalign::QualityProcessor processor;std::istringstream input(input_tsv);std::ostringstream output;std::string line;std::size_t index=0;
    while(std::getline(input,line)){const auto v=split(line);if(v.size()!=13)throw std::runtime_error("repair columns");Record a{v[0],v[1],v[2]},b{v[3],v[4],v[5]};if(a.seq.size()!=a.qual.size()||b.seq.size()!=b.qual.size())throw std::runtime_error("raw lengths");Batch batch;{
        NativeQualityGuard<Batch,Record> guard(batch,a,b,processor,cut,false,0,&model,mode,&full_model,false);
        MergeRecord r;r.merged={v[7],v[8],v[9]};r.unmerged2={v[10],v[11],v[12]};
        if(v[6]=="merged"){r.is_merged=true;batch.n_merged=1;}else if(v[6]=="unmerged")batch.n_unmerged=1;else if(v[6]=="orphan_r1"||v[6]=="orphan_r2"){r.is_orphan=true;r.orphan_r1=v[6]=="orphan_r1";batch.n_orphan=1;if(r.orphan_r1)batch.n_orphan_r1=1;else batch.n_orphan_r2=1;}else if(v[6]!="dropped")throw std::runtime_error("native outcome");
        if(v[6]!="dropped")batch.records.push_back(std::move(r));
    }
    output<<index++<<'\t';if(batch.records.empty()){output<<"dropped\t\t\t\t\t\t\n";continue;}const auto& r=batch.records.back();const char* outcome=r.is_merged?"merged":r.is_orphan?(r.orphan_r1?"orphan_r1":"orphan_r2"):"unmerged";output<<outcome<<'\t'<<r.merged.header<<'\t'<<r.merged.seq<<'\t'<<r.merged.qual<<'\t'<<r.unmerged2.header<<'\t'<<r.unmerged2.seq<<'\t'<<r.unmerged2.qual<<'\n';
    }
    return output.str();
}
}
