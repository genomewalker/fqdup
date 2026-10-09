#pragma once
// Frozen from 1a0ba92; function renamed only.
namespace frozen_geometry {
using fqalign::Geometry;
using fqalign::reverse_complement;
inline Geometry geometry(const std::string& a,const std::string& raw_b) {
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
}
