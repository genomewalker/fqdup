#pragma once
#include "native_quality.hpp"
#include "disagreement_risk.hpp"
#include "selective_disagreement.hpp"
#include "pair_mixture_xc.hpp"
#include "joint_pair_intervals.hpp"
#include "exclusive_core_intervals_xc.hpp"
#include "full_trusted_correspondence.hpp"
#include "retained_median_budget_xc.hpp"
#include "merged_span_budget.hpp"
#include "repair_correspondence.hpp"

// Owned by one repair worker; assignment/clear retains capacity between rows.
template<class Record>
struct NativeQualityScratch {
    Record raw1,raw2;
    std::string reverse;
    fqalign::PairRisk pair;
    std::vector<double> risk,proposed,selected,original1,original2,short_risk;
    std::vector<double> a,b,rawa,rawb;
    std::vector<fqquality::JointColumn> columns;
};

// Runs before PairOutcomeGuard on scope exit, including all early continues.
// Batches retain their native bounded queues and original input ordering.
template<class Batch, class Record>
class NativeQualityGuard {
    Batch& batch;
    fqalign::QualityProcessor& processor;
    double cutoff;
    bool raise_disagreements;
    int selective_mode;
    const fqquality::PairMixture* mixture;
    int mixture_mode;
    bool merged_floor;
    const fqquality::PairMixture* full_mixture;
    size_t original_size;
    NativeQualityScratch<Record> owned_scratch;
    NativeQualityScratch<Record>& scratch;
    Record& raw1;
    Record& raw2;
    static void trim(Record& r,fqalign::TrimInterval iv,const char* channel,int mate) {
        const size_t length=r.seq.size();
        if(iv.end<=iv.begin){r.seq.clear();r.qual.clear();return;}
        r.header+=" FQT=v1;ch="+std::string(channel)+";m="+std::to_string(mate)
            +";olen="+std::to_string(length)+";s="+std::to_string(iv.begin)+";e="+std::to_string(iv.end);
        r.seq=r.seq.substr(iv.begin,iv.end-iv.begin);
        r.qual=r.qual.substr(iv.begin,iv.end-iv.begin);
    }
    static void prefix(const Record& emitted,const Record& raw) {
        if(emitted.seq.size()>raw.seq.size()||emitted.qual.size()!=emitted.seq.size()
           ||raw.seq.compare(0,emitted.seq.size(),emitted.seq)!=0
           ||raw.qual.compare(0,emitted.qual.size(),emitted.qual)!=0) {
            std::cerr<<"experimental gapped Mott: native mate is not a raw prefix\n";
            std::abort();
        }
    }
    static void repair(Record& emitted,const Record& recipient,const Record& donor,int mate,const fqalign::Correspondence& corr,size_t parent2_length) {
        if(emitted.seq.empty()||!corr.eligible)return;
        const auto tag=emitted.header.rfind(" FQT=v1;");if(tag==std::string::npos)throw std::runtime_error("repair trace");
        const auto field=emitted.header.find(";s=",tag);if(field==std::string::npos)throw std::runtime_error("repair start");
        const size_t begin=std::stoull(emitted.header.substr(field+3));std::string edits;
        for(auto [i,j]:corr.invariant){const size_t k=parent2_length-1-size_t(j),pos=mate==1?size_t(i):k,other=mate==1?k:size_t(i);if(pos<begin||pos>=begin+emitted.seq.size())continue;
          const char a=recipient.seq.at(pos),b=fqalign::reverse_complement(std::string(1,donor.seq.at(other)))[0];const int q1=(unsigned char)recipient.qual.at(pos)-33,q2=(unsigned char)donor.qual.at(other)-33;
          if(a=='N'||b=='N'||a==b||q2<=q1)continue;
          const int q=std::max(1,q2-q1);emitted.seq.at(pos-begin)=b;emitted.qual.at(pos-begin)=char(q+33);
          if(!edits.empty())edits+=';';edits+=std::to_string(pos)+":"+std::to_string(other)+":"+b+":"+std::to_string(q);
        }
        if(!edits.empty())emitted.header+=" FQC=v1;edits="+edits;
    }
public:
    NativeQualityGuard(Batch& out,const Record& r1,const Record& r2,
                       fqalign::QualityProcessor& p,double c,bool raise=false,int selective=0,const fqquality::PairMixture* model=nullptr,int mode=1,const fqquality::PairMixture* full_model=nullptr,bool floor=false,NativeQualityScratch<Record>* reuse=nullptr)
        :batch(out),processor(p),cutoff(c),raise_disagreements(raise),selective_mode(selective),mixture(model),mixture_mode(mode),full_mixture(full_model),merged_floor(floor),original_size(out.records.size()),scratch(reuse?*reuse:owned_scratch),raw1(scratch.raw1),raw2(scratch.raw2) {
        if(cutoff>=0){raw1=r1;raw2=r2;}
    }
    ~NativeQualityGuard() {
        if(cutoff<0||batch.records.size()==original_size)return;
        if(batch.records.size()!=original_size+1)std::abort();
        auto& r=batch.records.back();
        if(r.is_merged) {
            auto& risk=scratch.risk; risk.clear(); risk.reserve(r.merged.qual.size());
            for(char q:r.merged.qual)risk.push_back(fqalign::phred_error(q));
            if(mixture && mixture_mode>=6) {
                const auto geom=fqalign::quality_independent_geometry(raw1.seq,raw2.seq);
                auto& pr=scratch.pair;processor.risks_into(pr,raw1.seq,raw1.qual,raw2.seq,raw2.qual,geom.eligible?geom.insert:-1);
                fqquality::PairMixture::Audit audit;
                if(mixture->apply(pr,raw1.seq,raw1.qual,raw2.seq,raw2.qual,1,&audit)) {
                    auto& rc=scratch.reverse;fqalign::reverse_complement_into(raw2.seq,rc);auto& proposed=scratch.proposed;proposed=risk;bool valid=true;
                    for(auto [i,j]:pr.correspondence.invariant) {
                        if(raw1.seq.at(i)=='N'||rc.at(j)=='N')continue;
                        const int k=int(raw2.seq.size())-1-j;
                        if(i+k!=int(r.merged.seq.size())-1){valid=false;break;}
                        if(size_t(i)>=proposed.size())continue;
                        const char chosen=r.merged.seq.at(i);int mate=chosen==raw1.seq.at(i)?1:chosen==rc.at(j)?2:chosen=='N'?0:-1;
                        if(mate<0){valid=false;break;}if(mate==0){proposed[i]=1;continue;}
                        const double e1=fqalign::phred_error(raw1.qual.at(i)),e2=fqalign::phred_error(raw2.qual.at(k));double value=0;
                        for(int z=0;z<fqquality::PairMixture::K;++z){double e=fqquality::PairMixture::component(z);auto cp=fqquality::PairMixture::conditional(z==0?e1:e,z==0?e2:e,raw1.seq.at(i)==rc.at(j));value+=audit.posterior[z]*(mate==1?cp.first:cp.second);}
                        proposed[i]=value;
                    }
                    if(valid)risk.swap(proposed);
                }
            }
            auto iv=mixture_mode==9?fqquality::merged_span_budget(risk,cutoff):fqalign::mott_segment(risk,risk.size(),cutoff);
            if(merged_floor && iv.end>iv.begin){auto& selected=scratch.selected;selected.assign(risk.begin()+iv.begin,risk.begin()+iv.end);const auto cb=fqalign::reverse_complement(raw2.seq);for(size_t i=0;i<selected.size();++i){const int c=int(iv.begin+i),k=int(r.merged.seq.size())-1-c;if(c>=int(raw1.seq.size())||k<0||k>=int(raw2.seq.size()))continue;char x=raw1.seq[c],y=cb[raw2.seq.size()-1-k];if(x!='N'&&y!='N'&&x!=y)selected[i]=std::max(selected[i],.5);}auto local=fqalign::mott_segment(selected,selected.size(),cutoff);if(local.end>local.begin)iv={iv.begin+local.begin,iv.begin+local.end};else iv={};}
            trim(r.merged,iv,"merged",0);
            if(r.merged.seq.empty()){--batch.n_merged;++batch.n_dropped;batch.records.pop_back();}
            return;
        }
        const auto g=fqalign::quality_independent_geometry(raw1.seq,raw2.seq);
        auto& p=scratch.pair;processor.risks_into(p,raw1.seq,raw1.qual,raw2.seq,raw2.qual,g.eligible?g.insert:-1);
        bool mixture_applied=false;
        if(mixture && (!(mixture_mode==3||mixture_mode==4)||!r.is_orphan)) {
            auto& original1=scratch.original1;auto& original2=scratch.original2;original1=p.r1;original2=p.r2;
            mixture_applied=mixture->apply(p,raw1.seq,raw1.qual,raw2.seq,raw2.qual,mixture_mode);
            if(mixture_applied && (mixture_mode==3||mixture_mode==4))
                fqquality::guard_pair_mixture(p,original1,original2,r.merged.seq.size(),r.unmerged2.seq.size(),cutoff);
        }
        bool full_applied=false;
        if(mixture_mode>=7 && !mixture_applied){
            for(size_t i=0;i<raw1.qual.size();++i)p.r1[i]=fqalign::phred_error(raw1.qual[i]);
            for(size_t i=0;i<raw2.qual.size();++i)p.r2[i]=fqalign::phred_error(raw2.qual[i]);
            if(full_mixture && !r.is_orphan){fqalign::reverse_complement_into(r.unmerged2.seq,scratch.reverse);auto corr=fqalign::full_trusted_correspondence(r.merged.seq,scratch.reverse);if(corr.eligible){const int shift=int(raw2.seq.size()-r.unmerged2.seq.size());for(auto& c:corr.invariant)c.second+=shift;for(auto& c:corr.certified)c.second+=shift;p.correspondence=std::move(corr);full_applied=full_mixture->apply(p,raw1.seq,raw1.qual,raw2.seq,raw2.qual,6,nullptr,true);mixture_applied=full_applied;}}
        }
        if(mixture && mixture_mode>=6 && !mixture_applied)
            mixture->prior_floor(p,raw1.seq,raw1.qual,raw2.seq,raw2.qual);
        if(raise_disagreements)fqquality::raise_anchored_disagreements(p,raw1.seq,raw2.seq);
        if(selective_mode>0 && (selective_mode!=2 || !r.is_orphan)) {
            auto& original1=scratch.original1;auto& original2=scratch.original2;original1=p.r1;original2=p.r2;
            fqquality::raise_weaker_disagreements(p,raw1.seq,raw2.seq,cutoff);
            if(selective_mode==2)fqquality::protect_unique_coverage(p,original1,original2,r.merged.seq.size(),r.unmerged2.seq.size(),cutoff);
        }
        if(r.is_orphan) {
            const int mate=r.orphan_r1?1:2;
            const auto& raw=mate==1?raw1:raw2;
            const auto& risk=mate==1?p.r1:p.r2;
            prefix(r.merged,raw);
            auto& short_risk=scratch.short_risk;short_risk.assign(risk.begin(),risk.begin()+r.merged.seq.size());const auto rawcut=(mixture_mode==7||mixture_mode==9)&&!mixture_applied?fqquality::retained_median_cutoff(short_risk,raw.qual.substr(0,r.merged.seq.size()),cutoff):cutoff;
            trim(r.merged,fqalign::mott_segment(short_risk,short_risk.size(),rawcut),"orphan",mate);
            if(r.merged.seq.empty()) {
                --batch.n_orphan;
                if(mate==1)--batch.n_orphan_r1;else --batch.n_orphan_r2;
                ++batch.n_dropped;batch.records.pop_back();
            }
            return;
        }
        prefix(r.merged,raw1);prefix(r.unmerged2,raw2);
        const size_t parent2_length=r.unmerged2.seq.size();
        fqalign::Correspondence repair_corr;
        if(mixture_mode==9){fqalign::reverse_complement_into(r.unmerged2.seq,scratch.reverse);repair_corr=fqrepair::full_trusted_correspondence(r.merged.seq,scratch.reverse
#if defined(FQDUP_TRUTH_ABLATION) && FQDUP_TRUTH_ABLATION
,true
#endif
);}
        std::array<fqalign::TrimInterval,2> intervals{fqalign::mott_segment(p.r1,r.merged.seq.size(),cutoff),fqalign::mott_segment(p.r2,r.unmerged2.seq.size(),cutoff)};
        if(mixture_mode==7&&!mixture_applied){auto& a=scratch.a;auto& b=scratch.b;a.assign(p.r1.begin(),p.r1.begin()+r.merged.seq.size());b.assign(p.r2.begin(),p.r2.begin()+r.unmerged2.seq.size());intervals={fqalign::mott_segment(a,a.size(),fqquality::retained_median_cutoff(a,raw1.qual.substr(0,a.size()),cutoff)),fqalign::mott_segment(b,b.size(),fqquality::retained_median_cutoff(b,raw2.qual.substr(0,b.size()),cutoff))};}
        if(mixture_applied && (mixture_mode==5||mixture_mode>=6)) {
            auto& rc=scratch.reverse;fqalign::reverse_complement_into(raw2.seq,rc);auto& columns=scratch.columns;columns.clear();
            for(auto [i,j]:p.correspondence.invariant)if(raw1.seq.at(i)!='N'&&rc.at(j)!='N')columns.push_back({i,int(raw2.seq.size())-1-j,raw1.seq.at(i)==rc.at(j)});
            {size_t hqd=0;for(auto& c:columns)hqd+=!c.agree&&raw1.qual.at(c.i)>=63&&raw2.qual.at(c.k)>=63;if(mixture_mode==9&&hqd>=2){auto& a=scratch.a;auto& b=scratch.b;a.resize(r.merged.seq.size());b.resize(r.unmerged2.seq.size());for(size_t i=0;i<a.size();++i)a[i]=fqalign::phred_error(raw1.qual[i]);for(size_t i=0;i<b.size();++i)b[i]=fqalign::phred_error(raw2.qual[i]);intervals={fqalign::mott_segment(a,a.size(),cutoff),fqalign::mott_segment(b,b.size(),cutoff)};columns.clear();}}bool exclusive=false;if(full_applied){auto& a=scratch.a;auto& b=scratch.b;a.assign(p.r1.begin(),p.r1.begin()+r.merged.seq.size());b.assign(p.r2.begin(),p.r2.begin()+r.unmerged2.seq.size());auto& rawa=scratch.rawa;auto& rawb=scratch.rawb;rawa=a;rawb=b;for(size_t i=0;i<a.size();++i)rawa[i]=fqalign::phred_error(raw1.qual[i]);for(size_t i=0;i<b.size();++i)rawb[i]=fqalign::phred_error(raw2.qual[i]);exclusive=fqquality::exclusive_core_intervals(a,b,rawa,rawb,columns,cutoff,fqquality::PairMixture::fold(raw1.seq,raw2.seq),intervals);}
            if(!exclusive)fqquality::joint_intervals(p.r1,p.r2,r.merged.seq.size(),r.unmerged2.seq.size(),columns,cutoff,intervals);
        }
        trim(r.merged,intervals[0],"r1",1);
        trim(r.unmerged2,intervals[1],"r2",2);
        if(mixture_mode==9){repair(r.merged,raw1,raw2,1,repair_corr,parent2_length);repair(r.unmerged2,raw2,raw1,2,repair_corr,parent2_length);}
        const bool first=!r.merged.seq.empty(),second=!r.unmerged2.seq.empty();
        if(first&&second)return;
        --batch.n_unmerged;
        if(!first&&!second){++batch.n_dropped;batch.records.pop_back();return;}
        r.is_orphan=true;r.orphan_r1=first;
        if(!first)r.merged=std::move(r.unmerged2);
        r.unmerged2=Record{};
        ++batch.n_orphan;
        if(first)++batch.n_orphan_r1;else ++batch.n_orphan_r2;
    }
};
