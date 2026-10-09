#pragma once
// A perfect native path belongs to the best-minus-two path set unless a
// length-L+3 ungapped band edge is also perfect. Invariants must therefore be
// a subset of that native path. See CLEAN.md; this performs no alignment.
struct NativeCleanBound {
    std::array<int,2> component{{-2,-2}};
    explicit NativeCleanBound(const fqquality::PairMixture& model) {
        for(int f=0;f<2;++f){
            if(!model.supported[f]){component[f]=-1;continue;}
            int found=-1;
            for(int z=0;z<fqquality::PairMixture::K;++z)if(model.weights[f][z]>0){
                if(found>=0){found=-2;break;}found=z;
            }
            component[f]=found;
        }
    }
    static char complement(char c) {
        switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';default:return 'N';}
    }
    template<class Row> int64_t operator()(const Row& v) const {
        if(v.size()!=13)return -1;
        if(v[6]=="dropped")return 0;
        if(!RepairSkipBound::supported(v)||v[1].size()>512||v[4].size()>512)return -1;
        const int z=component[1-fqquality::PairMixture::fold(v[1],v[4])];
        if(z==-2)return -1; // Mixed priors retain the existing exact check.
        const int insert=v[8].size(),n=v[1].size(),m=v[4].size();
        const int begin=std::max(0,insert-m),end=std::min(n,insert),length=end-begin;
        if(length<30)return -1;
        for(int i=begin;i<end;++i)
            if(v[1][i]=='N'||v[1][i]!=complement(v[4][insert-1-i])||v[8][i]!=v[1][i])return -1;
        // All banded paths have at most L+3 aligned columns. A gapped path
        // pays at least four; only perfect ungapped L+3 paths can exclude L.
        for(int delta:{-3,3}){
            const int other=insert+delta,lo=std::max(0,other-m),hi=std::min(n,other);
            if(hi-lo!=length+3)continue;
            bool perfect=true;
            for(int i=lo;i<hi;++i)if(v[1][i]=='N'||v[1][i]!=complement(v[4][other-1-i])){perfect=false;break;}
            if(perfect)return -1;
        }
        int64_t sum=0;
        for(int i=0;i<insert;++i){
            double upper=fqalign::phred_error(v[9][i]);
            if(z>=0&&i>=begin&&i<end){
                const double risk=fqquality::PairMixture::lookup(z,v[2][i],v[5][insert-1-i],true).risk.first;
                upper=std::max(upper,risk+1e-12);
            }
            sum+=int64_t(std::floor(upper*RepairSkipBound::scale+.5));
        }
        return sum;
    }
};
// Included beside repair_skip.hpp, using the same frozen model and arithmetic.
// Perfect ungapped overlap forces reference geometry to select this insert
// or abstain: any different winner has this >=30-base, zero-mismatch candidate
// in its second-best scan, violating the required second fraction >= .45.
template<class Row>
int64_t clean_pair_bound(const Row& v,const fqquality::PairMixture& model,
                         fqalign::QualityProcessor& processor) {
    if(v.size()!=13)return -1;
    if(v[6]=="dropped")return 0;
    if(!RepairSkipBound::supported(v)||v[1].size()>512||v[4].size()>512)return -1;
    const int insert=v[8].size(),n=v[1].size(),m=v[4].size();
    const int begin=std::max(0,insert-m),end=std::min(n,insert);
    if(end-begin<30)return -1;
    for(int i=begin;i<end;++i){
        const int j=i-insert+m;
        const char other=fqalign::complement(v[4][m-1-j]);
        if(v[1][i]=='N'||other=='N'||v[1][i]!=other||v[8][i]!=v[1][i])return -1;
    }
    auto p=processor.risks(v[1],v[2],v[4],v[5],insert);
    if(p.correspondence.eligible){
        for(auto [i,j]:p.correspondence.invariant){
            const int k=m-1-j;
            const char other=fqalign::complement(v[4].at(k));
            if(v[1][i]=='N'||other=='N')continue;
            if(i+k!=insert-1||v[1][i]!=other||v[8].at(i)!=v[1][i])return -1;
        }
        fqquality::PairMixture::Audit audit;
        if(model.apply(p,v[1],v[2],v[4],v[5],1,&audit))
            return RepairSkipBound::weighted_sum(v,p,audit);
    }
    int64_t sum=0;
    for(char q:v[9])sum+=int64_t(std::floor(fqalign::phred_error(q)*RepairSkipBound::scale+.5));
    return sum;
}
