#pragma once
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
    const auto rc=fqalign::reverse_complement(v[4]);
    for(int i=begin;i<end;++i){
        const int j=i-insert+m;
        if(v[1][i]=='N'||rc[j]=='N'||v[1][i]!=rc[j]||v[8][i]!=v[1][i])return -1;
    }
    auto p=processor.risks(v[1],v[2],v[4],v[5],insert);
    if(p.correspondence.eligible){
        for(auto [i,j]:p.correspondence.invariant){
            if(v[1][i]=='N'||rc[j]=='N')continue;
            const int k=m-1-j;
            if(i+k!=insert-1||v[1][i]!=rc[j]||v[8].at(i)!=v[1][i])return -1;
        }
        fqquality::PairMixture::Audit audit;
        if(model.apply(p,v[1],v[2],v[4],v[5],1,&audit))
            return RepairSkipBound::weighted_sum(v,p,audit);
    }
    int64_t sum=0;
    for(char q:v[9])sum+=int64_t(std::floor(fqalign::phred_error(q)*RepairSkipBound::scale+.5));
    return sum;
}
