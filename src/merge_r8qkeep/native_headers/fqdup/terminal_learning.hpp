#pragma once
#include "fqdup/terminal_construct.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <unordered_map>
namespace fqdup {
// Bounded, experimental learning on ORIGINAL records only. Counts used as
// independent support are group capped; mixture prevalence uses raw records.
// Group proxies and a composition null do not establish biological provenance.
struct TerminalLearningConfig {
    size_t max_records=1000000, sample_bytes=64*1024*1024, max_nodes=100000;
    size_t max_seeds=16, max_candidates=256;
    unsigned max_cut=48, min_cut=12, min_discovery=24, min_validation=8;
    double modeled_false_clips_per_million=10;
    bool apply=false;
};
struct PackedTerminalWindow {
    std::array<uint64_t,2> words{};
    uint64_t valid=0;
    uint8_t length=0;
    static int base(char c) { switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;} }
    void push(char c) {
        if(length==64)return;
        int b=base(c);size_t i=length++;
        if(b>=0){words[i/32]|=static_cast<uint64_t>(b)<<(2*(i%32));valid|=uint64_t{1}<<i;}
    }
    int at(size_t i) const { return i<length && (valid&(uint64_t{1}<<i)) ? static_cast<int>((words[i/32]>>(2*(i%32)))&3) : -1; }
    int code(size_t pos,size_t n=8) const {
        if(pos+n>length)return -1;
        unsigned value=0;for(size_t j=pos;j<pos+n;++j){int b=at(j);if(b<0)return -1;value=(value<<2)|b;}return static_cast<int>(value);
    }
    bool prefix(std::string_view path) const {
        if(path.size()>length)return false;
        for(size_t j=0;j<path.size();++j)if(at(j)!=base(path[j]))return false;return true;
    }
};
struct TerminalLearningSample {
    std::array<PackedTerminalWindow,2> end;
    PackedTerminalWindow middle;
    uint64_t group=0;
    uint32_t length=0;
};
struct LearnedTerminalCandidate {
    TerminalConstruct rule;
    std::string path, reason="insufficient_support";
    std::array<uint64_t,3> groups{};
    uint64_t raw_matches=0;
    std::array<uint64_t,3> context_callable{}, distinct_context{};
    std::array<std::array<uint64_t,4>,3> next_base{};
    double biological_probability=1;
    bool eligible=false;
};
struct FrozenTerminalModel {
    std::vector<LearnedTerminalCandidate> candidates;
    std::array<uint64_t,3> groups{};
    uint64_t sampled_records=0, background_windows=0;
    bool saturated=false;
    std::string group_basis="exact_sequence_or_central_proxy";
    std::vector<TerminalConstruct> rules(bool apply) const {
        std::vector<TerminalConstruct> result;
        for(const auto& c:candidates){auto rule=c.rule;rule.enabled=apply&&c.eligible;result.push_back(std::move(rule));}
        return result;
    }
    void write(const std::string& path,bool apply) const {
        std::ofstream out(path);if(!out)throw std::runtime_error("cannot open terminal model output");
        out<<"#fqdup_terminal_model_v1\n#scores=experimental_not_calibrated_provenance\n#sampled_records="<<sampled_records<<" background_windows="<<background_windows<<" saturated="<<saturated<<" group_basis="<<group_basis<<"\n";
        out<<"id\tend\tsequence\tlog_odds\teligible\tenabled\tdiscovery_groups\tcalibration_groups\tvalidation_groups\traw_matches\tbiological_probability\tdiscovery_contexts\tcalibration_contexts\tvalidation_contexts\treason\n"<<std::setprecision(17);
        for(const auto& c:candidates)out<<c.rule.id<<'\t'<<(c.rule.five_prime?"5":"3")<<'\t'<<c.rule.sequence<<'\t'<<c.rule.log_odds<<'\t'<<c.eligible<<'\t'<<(apply&&c.eligible)<<'\t'<<c.groups[0]<<'\t'<<c.groups[1]<<'\t'<<c.groups[2]<<'\t'<<c.raw_matches<<'\t'<<c.biological_probability<<'\t'<<c.distinct_context[0]<<'\t'<<c.distinct_context[1]<<'\t'<<c.distinct_context[2]<<'\t'<<c.reason<<'\n';
        out.flush();if(!out)throw std::runtime_error("terminal model write failed");
    }
};
inline std::vector<TerminalConstruct> read_terminal_model(const std::string& path) {
    std::ifstream in(path);std::string line;
    if(!in || !std::getline(in,line) || line!="#fqdup_terminal_model_v1")throw std::invalid_argument("unsupported terminal model version");
    std::vector<TerminalConstruct> rules;
    bool header=false;
    while(std::getline(in,line)) {
        if(line.empty()||line[0]=='#')continue;
        if(!header){if(line!="id\tend\tsequence\tlog_odds\teligible\tenabled\tdiscovery_groups\tcalibration_groups\tvalidation_groups\traw_matches\tbiological_probability\tdiscovery_contexts\tcalibration_contexts\tvalidation_contexts\treason")throw std::invalid_argument("invalid terminal model columns");header=true;continue;}
        std::vector<std::string> fields;std::istringstream row(line);std::string field;while(std::getline(row,field,'\t'))fields.push_back(field);
        if(fields.size()!=15 || (fields[1]!="5"&&fields[1]!="3") || (fields[4]!="0"&&fields[4]!="1") || (fields[5]!="0"&&fields[5]!="1") || (fields[5]=="1"&&fields[4]!="1"))throw std::invalid_argument("invalid terminal model row");
        size_t consumed=0;double score=std::stod(fields[3],&consumed);
        if(consumed!=fields[3].size()||!std::isfinite(score))throw std::invalid_argument("invalid terminal model score");
        rules.push_back({fields[0],fields[2],fields[1]=="5",score,fields[5]=="1"});
        if(rules.size()>256)throw std::invalid_argument("terminal model exceeds 256 configurations");
    }
    if(!header||in.bad())throw std::invalid_argument("incomplete terminal model");
    CompiledTerminalCatalogue validate(rules);return rules;
}
class BoundedTerminalLearner {
    TerminalLearningConfig config_;
    std::vector<TerminalLearningSample> sample_;
    bool saturated_=false;
    struct Node { std::array<uint32_t,4> child{};uint32_t support=0; };
    std::vector<Node> nodes_;
    static uint64_t mix(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
    static unsigned partition(uint64_t group){unsigned p=mix(group)&3;return p<2?0:p-1;}
    void extend(std::vector<uint32_t>& indices,size_t begin,size_t end,unsigned side,
                std::string path,FrozenTerminalModel& model) {
        if(end-begin<config_.min_discovery || path.size()>config_.max_cut)return;
        if(nodes_.size()>=config_.max_nodes){model.saturated=true;return;}
        const auto node=nodes_.size();nodes_.push_back({});nodes_[node].support=static_cast<uint32_t>(end-begin);
        std::array<uint64_t,4> count{};uint64_t callable=0;
        for(size_t i=begin;i<end;++i){int b=sample_[indices[i]].end[side].at(path.size());if(b>=0){++count[b];++callable;}}
        auto maximum=*std::max_element(count.begin(),count.end());
        // A diversity change proposes a junction; it never proves technical origin.
        if(path.size()>=config_.min_cut && callable>=config_.min_discovery && maximum<0.90*callable) {
            if(model.candidates.size()==config_.max_candidates){model.saturated=true;return;}
            LearnedTerminalCandidate c;c.path=path;auto sequence=path;if(side==1)std::reverse(sequence.begin(),sequence.end());
            c.rule={"learned_"+std::to_string(model.candidates.size()+1),sequence,side==0,-100,false};model.candidates.push_back(std::move(c));
        }
        if(path.size()==config_.max_cut)return;
        size_t cursor=begin;
        for(int b=0;b<4;++b) {
            auto stop=std::partition(indices.begin()+cursor,indices.begin()+end,[&](uint32_t i){return sample_[i].end[side].at(path.size())==b;});
            size_t next=static_cast<size_t>(stop-indices.begin());
            // Follow a supported conserved continuation, including a major
            // longer-remnant component. Random biological sub-branches stop.
            if(next-cursor>=config_.min_discovery && next-cursor>=0.55*(end-begin)) {
                nodes_[node].child[b]=static_cast<uint32_t>(nodes_.size());
                extend(indices,cursor,next,side,path+"ACGT"[b],model);
            }
            cursor=next;
        }
    }
public:
    explicit BoundedTerminalLearner(TerminalLearningConfig config={}) : config_(config) {
        if(!config_.max_records || config_.max_records>1000000 || config_.sample_bytes<sizeof(TerminalLearningSample) || config_.sample_bytes>128*1024*1024 || !config_.max_nodes || config_.max_nodes>100000 || !config_.max_seeds || config_.max_seeds>16 || !config_.max_candidates || config_.max_candidates>256 || config_.min_cut<12 || config_.min_cut>config_.max_cut || config_.max_cut>48 || config_.min_discovery<8 || config_.min_validation<4 || !std::isfinite(config_.modeled_false_clips_per_million) || config_.modeled_false_clips_per_million<=0)
            throw std::invalid_argument("unsupported bounded terminal learning configuration");
        sample_.reserve(std::min(config_.max_records,config_.sample_bytes/sizeof(TerminalLearningSample)));
        nodes_.reserve(config_.max_nodes);
    }
    bool add(std::string_view seq,uint64_t group=0) {
        if(sample_.size()==config_.max_records || (sample_.size()+1)*sizeof(TerminalLearningSample)>config_.sample_bytes){saturated_=true;return false;}
        TerminalLearningSample r;r.length=static_cast<uint32_t>(std::min(seq.size(),size_t{UINT32_MAX}));
        for(size_t j=0;j<std::min(seq.size(),size_t{64});++j){r.end[0].push(seq[j]);r.end[1].push(seq[seq.size()-1-j]);}
        if(seq.size()>=32){size_t start=(seq.size()-16)/2;for(size_t j=start;j<start+16;++j)r.middle.push(seq[j]);}
        if(group==0) {
            // Long reads use a candidate-independent central proxy. Short
            // reads can only be exact-duplicate grouped; this is disclosed.
            std::string_view key=seq.size()>=112?seq.substr((seq.size()-32)/2,32):seq;
            group=1469598103934665603ULL;for(char c:key){group^=static_cast<unsigned char>(c);group*=1099511628211ULL;}
        }
        r.group=group;sample_.push_back(r);return true;
    }
    size_t size()const{return sample_.size();}
    size_t bytes()const{return sample_.capacity()*sizeof(TerminalLearningSample)+nodes_.capacity()*sizeof(Node);}
    FrozenTerminalModel freeze() {
        FrozenTerminalModel model;model.sampled_records=sample_.size();model.saturated=saturated_;
        std::vector<uint32_t> unique(sample_.size());std::iota(unique.begin(),unique.end(),0);
        std::stable_sort(unique.begin(),unique.end(),[&](uint32_t a,uint32_t b){return sample_[a].group<sample_[b].group;});
        unique.erase(std::unique(unique.begin(),unique.end(),[&](uint32_t a,uint32_t b){return sample_[a].group==sample_[b].group;}),unique.end());
        for(auto i:unique)++model.groups[partition(sample_[i].group)];
        for(unsigned side=0;side<2;++side) {
            std::vector<uint32_t> counts(65536),background(65536);uint64_t background_callable=0,terminal_callable=0;
            for(auto i:unique)if(partition(sample_[i].group)==0) {
                int k=sample_[i].end[side].code(0);if(k>=0){++counts[k];++terminal_callable;}
                int b=sample_[i].middle.code(0);if(b>=0){++background[b];++background_callable;}
            }
            std::vector<uint32_t> seeds;
            for(uint32_t k=0;k<65536;++k)if(counts[k]>=config_.min_discovery) {
                double expected=background_callable ? (background[k]+0.5)/(background_callable+32768.0)*terminal_callable : 0;
                if(counts[k]>=4*(expected+1))seeds.push_back(k);
            }
            std::sort(seeds.begin(),seeds.end(),[&](uint32_t a,uint32_t b){return counts[a]==counts[b]?a<b:counts[a]>counts[b];});
            if(seeds.size()>config_.max_seeds){model.saturated=true;seeds.resize(config_.max_seeds);}
            for(auto k:seeds) {
                std::vector<uint32_t> members;members.reserve(counts[k]);
                for(auto i:unique)if(partition(sample_[i].group)==0 && sample_[i].end[side].code(0)==static_cast<int>(k))members.push_back(i);
                std::string path(8,'A');auto code=k;for(int j=7;j>=0;--j){path[j]="ACGT"[code&3];code>>=2;}
                extend(members,0,members.size(),side,path,model);
            }
        }
        if(model.candidates.empty())return model;
        const size_t n=model.candidates.size();
        std::array<std::vector<int>,2> head;for(auto& h:head)h.assign(65536,-1);std::vector<int> next(n,-1);
        for(size_t j=0;j<n;++j){unsigned side=model.candidates[j].rule.five_prime?0:1;PackedTerminalWindow w;for(char b:model.candidates[j].path)w.push(b);int k=w.code(0);next[j]=head[side][k];head[side][k]=static_cast<int>(j);}
        auto matching=[&](const TerminalLearningSample& r,unsigned side,auto&& fn){int k=r.end[side].code(0);if(k<0)return;for(int j=head[side][k];j>=0;j=next[j])if(r.end[side].prefix(model.candidates[j].path))fn(static_cast<size_t>(j));};
        std::vector<std::array<std::array<uint64_t,1024>,3>> context(n);
        for(auto i:unique) {
            const auto& r=sample_[i];unsigned part=partition(r.group);
            for(unsigned side=0;side<2;++side)matching(r,side,[&](size_t j){auto& c=model.candidates[j];++c.groups[part];int b=r.end[side].at(c.path.size());if(b>=0){++c.context_callable[part];++c.next_base[part][b];}int k=r.end[side].code(c.path.size());if(k>=0)context[j][part][static_cast<unsigned>(k)>>6]|=uint64_t{1}<<(k&63);});
        }
        for(size_t j=0;j<n;++j)for(size_t p=0;p<3;++p){for(auto bits:context[j][p])model.candidates[j].distinct_context[p]+=__builtin_popcountll(bits);}
        for(const auto& r:sample_)for(unsigned side=0;side<2;++side)matching(r,side,[&](size_t j){++model.candidates[j].raw_matches;});
        size_t longest=0;for(const auto& c:model.candidates)longest=std::max(longest,c.path.size());
        std::array<uint64_t,4> middle_base{};std::array<std::array<uint64_t,4>,4> middle_trans{};
        for(auto i:unique)if(partition(sample_[i].group)==1) {
            const auto& r=sample_[i];size_t start=r.length>=16?(r.length-16)/2:0;
            if(start<longest || start+16+longest>r.length || r.middle.length!=16)continue;
            bool valid=true;for(size_t k=0;k<16;++k)if(r.middle.at(k)<0)valid=false;if(!valid)continue;
            ++model.background_windows;
            for(size_t k=0;k<16;++k){int b=r.middle.at(k);++middle_base[b];if(k)++middle_trans[r.middle.at(k-1)][b];}
        }
        uint64_t total_base=std::accumulate(middle_base.begin(),middle_base.end(),uint64_t{0});
        std::array<std::array<uint64_t,4>,2> first_base{};
        std::array<uint64_t,2> terminal_background{};
        for(auto i:unique)if(partition(sample_[i].group)==1)for(unsigned side=0;side<2;++side) {
            bool candidate=false;matching(sample_[i],side,[&](size_t){candidate=true;});
            int b=sample_[i].end[side].at(0);
            if(!candidate && b>=0){++first_base[side][b];++terminal_background[side];}
        }
        std::vector<double> inverse(n),prior(n,1e-12);
        for(size_t j=0;j<n;++j) {
            const auto& path=model.candidates[j].path;
            unsigned side=model.candidates[j].rule.five_prime?0:1;
            int previous=PackedTerminalWindow::base(path[0]);double logp=std::log((first_base[side][previous]+0.5)/(terminal_background[side]+2.0));
            for(size_t k=1;k<path.size();++k){int b=PackedTerminalWindow::base(path[k]);auto denom=std::accumulate(middle_trans[previous].begin(),middle_trans[previous].end(),uint64_t{0});logp+=std::log((middle_trans[previous][b]+0.5)/(denom+2.0));previous=b;}
            model.candidates[j].biological_probability=std::exp(logp);inverse[j]=std::exp(-logp);
        }
        // Bounded aggregated matching signatures; raw records fit prevalence,
        // group-capped observations above are the independent support gates.
        using Mask=std::array<uint64_t,4>;
        struct Hash {size_t operator()(const Mask& m)const{size_t h=0;for(auto v:m)h^=v+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
        std::array<std::unordered_map<Mask,uint64_t,Hash>,2> classes;
        std::array<uint64_t,2> totals{};
        for(const auto& r:sample_)if(partition(r.group)==1)for(unsigned side=0;side<2;++side) {
            Mask mask{};matching(r,side,[&](size_t j){mask[j/64]|=uint64_t{1}<<(j%64);});
            if(classes[side].size()>=4096 && !classes[side].count(mask)){model.saturated=true;continue;}
            ++classes[side][mask];++totals[side];
        }
        for(unsigned side=0;side<2;++side) {
            double none=0.9;size_t count=0;
            for(size_t j=0;j<n;++j)if(model.candidates[j].rule.five_prime==(side==0)){prior[j]=0.1;++count;}
            if(!count||!totals[side])continue;
            for(size_t j=0;j<n;++j)if(model.candidates[j].rule.five_prime==(side==0))prior[j]/=count;
            for(unsigned iteration=0;iteration<5;++iteration) {
                std::vector<double> updated(n,0);double new_none=0;
                for(const auto& item:classes[side]) {
                    double denominator=none;
                    for(size_t j=0;j<n;++j)if(item.first[j/64]&(uint64_t{1}<<(j%64)))denominator+=prior[j]*inverse[j];
                    new_none+=item.second*none/denominator;
                    for(size_t j=0;j<n;++j)if(item.first[j/64]&(uint64_t{1}<<(j%64)))updated[j]+=item.second*prior[j]*inverse[j]/denominator;
                }
                none=std::max(new_none/totals[side],1e-12);
                for(size_t j=0;j<n;++j)if(model.candidates[j].rule.five_prime==(side==0))prior[j]=std::max(updated[j]/totals[side],1e-12);
            }
            for(size_t j=0;j<n;++j)if(model.candidates[j].rule.five_prime==(side==0))model.candidates[j].rule.log_odds=std::max(-100.0,std::min(100.0,std::log(prior[j]*inverse[j]/none)));
        }
        double union_probability=0;
        for(auto& c:model.candidates) {
            bool diverse=true;
            for(size_t p=0;p<3;++p){auto maximum=*std::max_element(c.next_base[p].begin(),c.next_base[p].end());if(c.context_callable[p]<config_.min_validation || maximum>=0.90*c.context_callable[p] || c.distinct_context[p]<config_.min_validation)diverse=false;}
            if(model.saturated)c.reason="discovery_saturated";
            else if(model.background_windows<32 || terminal_background[c.rule.five_prime?0:1]<32)c.reason="missing_background";
            else if([&](){std::array<unsigned,4> counts{};for(char b:c.path)++counts[PackedTerminalWindow::base(b)];double entropy=0;for(auto count:counts)if(count){double p=double(count)/c.path.size();entropy-=p*std::log2(p);}return entropy<1.25;}())c.reason="low_complexity_candidate";
            else if(c.groups[0]<config_.min_discovery || c.groups[1]<config_.min_validation || c.groups[2]<config_.min_validation)c.reason="insufficient_independent_groups";
            else if(!diverse)c.reason="uncertain_junction";
            else {c.eligible=true;c.reason="experimental_sequence_supported";union_probability+=c.biological_probability;}
        }
        if(union_probability*1e6>config_.modeled_false_clips_per_million)for(auto& c:model.candidates)if(c.eligible){c.eligible=false;c.reason="modeled_background_budget_exceeded";}
        for(auto& c:model.candidates)c.rule.enabled=config_.apply&&c.eligible;
        return model;
    }
};
}

