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
inline std::string reverse_complement(std::string_view s) {
    std::string out(s.size(),'N');
    for (size_t i=0;i<s.size();++i) {
        char c=s[i];
        switch(c){case 'A':c='T';break;case 'C':c='G';break;
          case 'G':c='C';break;case 'T':c='A';break;case 'N':break;
          default:throw std::invalid_argument("unsupported base");}
        out[s.size()-1-i]=c;
    }
    return out;
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
inline Geometry quality_independent_geometry_uncached(const std::string& a,const std::string& raw_b) {
    Geometry out;
    if(a.empty()||raw_b.empty()||a.size()>512||raw_b.size()>512)return out;
    const auto b=reverse_complement(raw_b);
    for(char c:a)if(std::string_view("ACGTN").find(c)==std::string_view::npos)
        throw std::invalid_argument("unsupported base");
    auto count=[&](int insert) {
        int informative=0,mismatches=0;
        for(int i=std::max(0,insert-int(b.size()));i<std::min(int(a.size()),insert);++i) {
            int j=i-insert+int(b.size());
            if(a[i]!='N'&&b[j]!='N'){++informative;mismatches+=a[i]!=b[j];}
        }
        return std::pair<int,int>{informative,mismatches};
    };
    const int last=int(a.size()+b.size())-11;
    for(int insert=11;insert<=last;++insert) {
        auto [n,mm]=count(insert);
        if(n<11)continue;
        const double frac=double(mm)/n;
        if(out.insert<0||frac<out.fraction)
            {out.insert=insert;out.informative=n;out.mismatches=mm;out.fraction=frac;}
    }
    if(out.insert<0)return out;
    for(int insert=11;insert<=last;++insert) {
        if(insert==out.insert)continue;
        auto [n,mm]=count(insert);
        if(n>=30)out.second=std::min(out.second,double(mm)/n);
    }
    out.eligible=out.informative>=30&&out.fraction<=0.30&&out.second>=0.45;
    return out;
}
inline Geometry quality_independent_geometry(const std::string& a,const std::string& b) {
    return alignment_cache::get<Geometry>(a,b,0,[&]{return quality_independent_geometry_uncached(a,b);});
}
class QualityProcessor {
    BandedCorrespondence verifier;
public:
    // Raw R2 orientation. Caller supplies quality-independent geometry, or -1.
    PairRisk risks(const std::string& s1,const std::string& q1,
                   const std::string& s2,const std::string& q2,
                   int insert,bool enabled=true) {
        PairRisk out;
        risks_into(out,s1,q1,s2,q2,insert,enabled);
        return out;
    }
    void risks_into(PairRisk& out,const std::string& s1,const std::string& q1,
                    const std::string& s2,const std::string& q2,int insert,bool enabled=true) {
        if(s1.size()!=q1.size()||s2.size()!=q2.size())
            throw std::invalid_argument("sequence/quality length mismatch");
        out.r1.clear();out.r2.clear();out.discounted.clear();
        out.correspondence.best=-1000000000;out.correspondence.paths=0;
        out.correspondence.eligible=false;out.correspondence.reason="no_alignment";
        out.correspondence.invariant.clear();out.correspondence.certified.clear();
        out.r1.reserve(q1.size());out.r2.reserve(q2.size());
        for(char c:q1)out.r1.push_back(phred_error(c));
        for(char c:q2)out.r2.push_back(phred_error(c));
        if(!enabled||insert<0){out.correspondence.reason=enabled?"no_geometry":"disabled";return;}
        const auto corr=verifier.align(s1,reverse_complement(s2),insert);
        out.correspondence=corr;
        for(auto [i,j]:out.correspondence.certified) {
            const int k=int(s2.size())-1-j;
            const double p=agreement_error(out.r1.at(i),out.r2.at(k));
            if(p<out.r1[i]||p<out.r2[k])out.discounted.emplace_back(i,k);
            out.r1[i]=std::min(out.r1[i],p);out.r2[k]=std::min(out.r2[k],p);
        }
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

