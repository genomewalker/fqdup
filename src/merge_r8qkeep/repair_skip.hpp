#pragma once
// Included inside the wrapper namespace; uses the frozen repair arithmetic.
struct RepairSkipBound {
    static constexpr int64_t scale=1000000000000LL;
    explicit RepairSkipBound(const fqquality::PairMixture&) {}
    static bool supported(const std::vector<std::string>& v) {
        if(v.size()!=13||v[6]!="merged"||v[8].empty()||v[8].size()>512||v[8].size()!=v[9].size())return false;
        if(v[1].size()!=v[2].size()||v[4].size()!=v[5].size())return false;
        if(!std::numeric_limits<double>::is_iec559||std::fegetround()!=FE_TONEAREST)return false;
        for(int j:{1,4,8})for(char c:v[j])if(std::string_view("ACGTN").find(c)==std::string_view::npos)return false;
        for(int j:{2,5,9})for(unsigned char c:v[j])if(c<33||c>126)return false;
        return true;
    }
    // Called only with the router's fitted nominal posterior and correspondence,
    // when eligible geometry has insert == emitted length. These are exactly the
    // inputs the mode-9 guard recomputes. Cache a sum, never 66 weights per row.
    static int64_t weighted_sum(const std::vector<std::string>& v,
                               const fqalign::PairRisk& pr,
                               const fqquality::PairMixture::Audit& audit) {
        if(!supported(v)||!pr.correspondence.eligible||!audit.supported)return -1;
        if(v[1].size()>512||v[4].size()>512)return -1;
        std::array<double,512> upper{};
        for(size_t i=0;i<v[8].size();++i)upper[i]=fqalign::phred_error(v[9][i]);
        // Retain the original chosen-N fallback even outside invariant columns.
        for(size_t i=0;i<v[8].size();++i){const size_t k=v[8].size()-1-i;
            if(i<v[1].size()&&k<v[4].size()&&v[1][i]!='N'&&v[4][k]!='N'&&v[8][i]=='N')return -1;
        }
        const auto rc=fqalign::reverse_complement(v[4]);
        for(auto [i,j]:pr.correspondence.invariant){
            if(v[1].at(i)=='N'||rc.at(j)=='N')continue;
            const int k=int(v[4].size())-1-j;
            if(i+k!=int(v[8].size())-1)return -1;
            if(size_t(i)>=v[8].size())continue;
            const char chosen=v[8][i];
            const int mate=chosen==v[1][i]?1:chosen==rc[j]?2:0;
            if(!mate)return -1;
            const double e1=fqalign::phred_error(v[2].at(i)),e2=fqalign::phred_error(v[5].at(k));
            double value=0;
            for(int z=0;z<fqquality::PairMixture::K;++z){
                const double e=fqquality::PairMixture::component(z);
                const auto cp=fqquality::PairMixture::conditional(z==0?e1:e,z==0?e2:e,v[1][i]==rc[j]);
                value+=audit.posterior[z]*(mate==1?cp.first:cp.second);
            }
            // Identical normalization and ordered dot product; SKIP.md proves
            // the extra allowance conservatively covers positive arithmetic.
            if(!std::isfinite(value)||value<0||value+1e-12>1)return -1;
            upper[i]=std::max(upper[i],value+1e-12);
        }
        int64_t sum=0;
        for(size_t i=0;i<v[8].size();++i)sum+=int64_t(std::floor(upper[i]*scale+.5));
        return sum;
    }
    bool unchanged(const std::vector<std::string>& v,const fqquality::PairMixture& model,double cut,int64_t cached=-1) const {
        if(v.size()!=13||!std::isfinite(cut)||cut<0||cut>1)return false;
        if(v[6]=="dropped")return true;
        if(!supported(v))return false;
        const int f=1-fqquality::PairMixture::fold(v[1],v[4]);
        const bool update=model.supported[f]&&v[1].size()>=30&&v[1].size()<=512&&v[4].size()>=30&&v[4].size()<=512;
        const int64_t budget=int64_t(std::floor(cut*scale+.5))*int64_t(v[8].size());
        if(update)return cached>=0&&cached<=budget;
        int64_t sum=0;
        for(char q:v[9])sum+=int64_t(std::floor(fqalign::phred_error(q)*scale+.5));
        return sum<=budget;
    }
};
inline std::string serialize_repair(const Batch& batch,size_t index){
    std::ostringstream out;out<<index<<'\t';
    if(batch.records.empty()){out<<"dropped\t\t\t\t\t\t\n";return out.str();}
    const auto& r=batch.records.back();const char* outcome=r.is_merged?"merged":r.is_orphan?(r.orphan_r1?"orphan_r1":"orphan_r2"):"unmerged";
    out<<outcome<<'\t'<<r.merged.header<<'\t'<<r.merged.seq<<'\t'<<r.merged.qual<<'\t'<<r.unmerged2.header<<'\t'<<r.unmerged2.seq<<'\t'<<r.unmerged2.qual<<'\n';return out.str();
}
inline std::string native_full_span(const std::vector<std::string>& v,size_t index){
    Batch batch;if(v[6]=="dropped")return serialize_repair(batch,index);
    MergeRecord r;r.merged={v[7],v[8],v[9]};r.unmerged2={v[10],v[11],v[12]};r.is_merged=v[6]=="merged";r.is_orphan=v[6]=="orphan_r1"||v[6]=="orphan_r2";r.orphan_r1=v[6]=="orphan_r1";
    const auto tag=[](Record& rec,const char* channel,int mate){if(rec.seq.empty())return;rec.header+=" FQT=v1;ch="+std::string(channel)+";m="+std::to_string(mate)+";olen="+std::to_string(rec.seq.size())+";s=0;e="+std::to_string(rec.seq.size());};
    tag(r.merged,r.is_merged?"merged":r.is_orphan?"orphan":"r1",r.is_merged?0:r.is_orphan?(r.orphan_r1?1:2):1);
    if(!r.is_merged&&!r.is_orphan)tag(r.unmerged2,"r2",2);
    batch.records.push_back(std::move(r));return serialize_repair(batch,index);
}
