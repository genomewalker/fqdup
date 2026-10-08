class Task24OriginalQualityProcessor {
    BandedCorrespondence verifier;
public:
    // Raw R2 orientation. Caller supplies quality-independent geometry, or -1.
    PairRisk risks(const std::string& s1,const std::string& q1,
                   const std::string& s2,const std::string& q2,
                   int insert,bool enabled=true) {
        if(s1.size()!=q1.size()||s2.size()!=q2.size())
            throw std::invalid_argument("sequence/quality length mismatch");
        PairRisk out;
        out.r1.reserve(q1.size());out.r2.reserve(q2.size());
        for(char c:q1)out.r1.push_back(phred_error(c));
        for(char c:q2)out.r2.push_back(phred_error(c));
        if(!enabled||insert<0){out.correspondence.reason=enabled?"no_geometry":"disabled";return out;}
        out.correspondence=verifier.align(s1,reverse_complement(s2),insert);
        for(auto [i,j]:out.correspondence.certified) {
            const int k=int(s2.size())-1-j;
            const double p=agreement_error(out.r1.at(i),out.r2.at(k));
            if(p<out.r1[i]||p<out.r2[k])out.discounted.emplace_back(i,k);
            out.r1[i]=std::min(out.r1[i],p);out.r2[k]=std::min(out.r2[k],p);
        }
        return out;
    }
};
