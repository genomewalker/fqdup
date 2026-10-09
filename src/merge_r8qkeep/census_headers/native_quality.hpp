#pragma once
#include "native_correspondence.hpp"
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace fqalign {
// Development port of the frozen nominal-risk policy, not calibrated estimates.
inline double phred_error(char c) {
    const unsigned q = static_cast<unsigned char>(c);
    if (q < 33 || q > 126) throw std::invalid_argument("invalid Phred+33 quality");
    static const auto table = [] {
        std::array<double,94> t{};
        for (unsigned i=0;i<t.size();++i)
            t[i]=std::clamp(std::pow(10.0,-double(i)/10.0),1e-10,1.0);
        return t;
    }();
    return table[q-33];
}
inline double agreement_error(double a,double b) {
    a=std::clamp(a,1e-10,1.0);b=std::clamp(b,1e-10,1.0);
    const double num=a*b/3.0, den=(1.0-a)*(1.0-b)+num;
    return den>0.0?num/den:1.0;
}
inline char complement(char c) {
    switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';
        case 'T':return 'A';case 'N':return 'N';default:throw std::invalid_argument("unsupported base");}
}
// Callers supply distinct input and output storage; capacity survives rows.
inline void reverse_complement_into(std::string_view s,std::string& out) {
    out.resize(s.size());
    for(size_t i=0;i<s.size();++i)out[s.size()-1-i]=complement(s[i]);
}
inline std::string reverse_complement(std::string_view s) {
    std::string out;reverse_complement_into(s,out);return out;
}
struct PairRisk {
    std::vector<double> r1,r2;
    std::vector<std::pair<int,int>> discounted;
    Correspondence correspondence;
};
struct Geometry {
    int insert=-1, informative=0, mismatches=0;
    double fraction=1.0, second=9.0;
    bool eligible=false;
};
// Frozen Phase-A search: min overlap 11; selection 30/.30/.45.
// O(read_length) scratch, rather than an offsets-by-positions matrix.
inline Geometry quality_independent_geometry(const std::string& a,const std::string& raw_b) {
    Geometry out;
    if(a.empty()||raw_b.empty()||a.size()>512||raw_b.size()>512)return out;
    static thread_local std::string b;reverse_complement_into(raw_b,b);
    if(!::fqdup_fast_bases::valid(a))
        throw std::invalid_argument("unsupported base");
    auto count=[&](int insert) {
        const int begin=std::max(0,insert-int(b.size()));
        const int end=std::min(int(a.size()),insert);
        if(end<=begin)return std::pair<int,int>{0,0};
        return ::fqdup::r8qkeep::geometry_counts(a.data()+begin,
            b.data()+(begin-insert+int(b.size())),end-begin);
    };
    const int last=int(a.size()+b.size())-11;
    ::fqdup::r8qkeep::GeometryRanking ranking;
    for(int insert=11;insert<=last;++insert) {
        auto [n,mm]=count(insert);
        ranking.observe(insert,n,mm);
    }
    const auto& best=ranking.best;
    if(best.insert<0)return out;
    out.insert=best.insert;out.informative=best.informative;
    out.mismatches=best.mismatches;out.fraction=double(best.mismatches)/best.informative;
    const auto& second=ranking.runner_up();
    if(second.insert>=0)out.second=double(second.mismatches)/second.informative;
    out.eligible=out.informative>=30&&out.fraction<=0.30&&out.second>=0.45;
    return out;
}
class QualityProcessor {
    BandedCorrespondence verifier;
    std::string reverse_buffer;
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
        reverse_complement_into(s2,reverse_buffer);
        out.correspondence=verifier.align(s1,reverse_buffer,insert);
        for(auto [i,j]:out.correspondence.certified) {
            const int k=int(s2.size())-1-j;
            const double p=agreement_error(out.r1.at(i),out.r2.at(k));
            if(p<out.r1[i]||p<out.r2[k])out.discounted.emplace_back(i,k);
            out.r1[i]=std::min(out.r1[i],p);out.r2[k]=std::min(out.r2[k],p);
        }
        return out;
    }
};
struct TrimInterval { size_t begin=0,end=0; };
inline TrimInterval mott_segment(const std::vector<double>& risk,size_t prefix,double cutoff) {
    if(prefix>risk.size()||prefix>size_t(std::numeric_limits<int64_t>::max()/1000000000000LL)
       ||!std::isfinite(cutoff)||cutoff<0.0||cutoff>1.0)
        throw std::invalid_argument("invalid trimming arguments");
    // Integer sums reproduce the intended 12-decimal score grid and avoid drift.
    // Exact historical floating-point parity still requires differential validation.
    constexpr double scale=1e12;
    int64_t sum=0,best=0;size_t start=0,bstart=0,bend=0;
    for(size_t i=0;i<prefix;++i) {
        if(!std::isfinite(risk[i])||risk[i]<0.0||risk[i]>1.0)
            throw std::invalid_argument("invalid risk");
        const int64_t x=static_cast<int64_t>(std::nearbyint((cutoff-risk[i])*scale));
        if(sum<0){sum=x;start=i;}else sum+=x;
        const size_t length=i+1-start,blength=bend-bstart;
        if(sum>best||(sum==best&&sum>0&&(length>blength||(length==blength&&start<bstart))))
            {best=sum;bstart=start;bend=i+1;}
    }
    return best>0?TrimInterval{bstart,bend}:TrimInterval{};
}
}
