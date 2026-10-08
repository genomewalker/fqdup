#pragma once
// Included inside the wrapper namespace; uses the frozen repair arithmetic.
struct RepairSkipBound {
    double flat[2][2]{};
    explicit RepairSkipBound(const fqquality::PairMixture& model) {
        for(int f=0;f<2;++f)for(int agree=0;agree<2;++agree)
            for(int z=1;z<fqquality::PairMixture::K;++z)if(model.weights[f][z]>0){
                const double e=fqquality::PairMixture::component(z);
                const auto cp=fqquality::PairMixture::conditional(e,e,agree);
                flat[f][agree]=std::max(flat[f][agree],std::max(cp.first,cp.second));
            }
    }
    bool unchanged(const std::vector<std::string>& v,const fqquality::PairMixture& model,double cut) const {
        if(v[6]=="dropped")return true;
        if(v[6]!="merged"||v[8].empty()||v[8].size()>512||v[8].size()!=v[9].size())return false;
        if(!std::numeric_limits<double>::is_iec559||std::fegetround()!=FE_TONEAREST)return false;
        // Fall back for unsupported input, preserving full-path validation.
        for(int j:{1,4,8})for(char c:v[j])if(std::string_view("ACGTN").find(c)==std::string_view::npos)return false;
        for(int j:{2,5,9})for(unsigned char c:v[j])if(c<33||c>126)return false;
        const int f=1-fqquality::PairMixture::fold(v[1],v[4]);
        const bool update=model.supported[f]&&v[1].size()>=30&&v[1].size()<=512&&v[4].size()>=30&&v[4].size()<=512;
        constexpr int64_t scale=1000000000000LL;
        const int64_t budget=int64_t(std::floor(cut*scale+.5))*int64_t(v[8].size());
        int64_t sum=0;
        for(size_t i=0;i<v[8].size();++i){
            double upper=fqalign::phred_error(v[9][i]);const size_t k=v[8].size()-1-i;
            if(update&&i<v[1].size()&&k<v[4].size()&&v[1][i]!='N'&&v[4][k]!='N'){
                const char a=v[1][i],b=v[4][k]=='A'?'T':v[4][k]=='T'?'A':v[4][k]=='G'?'C':'G';
                const char chosen=v[8][i];
                if(chosen=='N')return false;
                if(chosen==a||chosen==b){
                    double bound=flat[f][a==b];
                    if(model.weights[f][0]>0){const auto cp=fqquality::PairMixture::conditional(fqalign::phred_error(v[2][i]),fqalign::phred_error(v[5][k]),a==b);bound=std::max(bound,chosen==a?cp.first:cp.second);}
                    // 66-term positive normalization and dot product: see SKIP.md.
                    upper=std::max(upper,bound+1e-12);
                    if(upper>1)return false;
                }
            }
            sum+=int64_t(std::floor(upper*scale+.5));if(sum>budget)return false;
        }
        return true;
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
