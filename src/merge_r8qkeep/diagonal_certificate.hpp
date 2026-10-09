#pragma once
#include <algorithm>
#include <cstdlib>
#include <string>
namespace fqdup_fast_diagonal {
// No matrix or heap scratch. False means unknown, never rejection of a pair.
// The optimistic suffix awards +1 to every possible remaining aligned column
// and charges no future gaps. Exhaustively excluding all near-best alternatives
// proves the perfect diagonal is the complete best-minus-two path set.
inline bool unique(const std::string& a,const std::string& b,int offset,int radius) {
    const int n=a.size(),m=b.size(),ai=std::max(0,-offset),bj=std::max(0,offset);
    const int length=std::min(n-ai,m-bj);
    if(length<30)return false;
    for(int k=0;k<length;++k)if(a[ai+k]=='N'||a[ai+k]!=b[bj+k])return false;
    const int threshold=length-2,limit=4*(length+7);int visited=0;
    const auto inside=[&](int i,int j){return i>=1&&i<=n&&j>=1&&j<=m&&std::abs(j-i-offset)<=radius;};
    const auto sub=[&](int i,int j){return a[i]==b[j]&&a[i]!='N'?1:-1;};
    const auto visit=[&](auto&& self,int i,int j,int state,int score,bool diagonal)->bool {
        if(score+std::min(n-i,m-j)<threshold)return true;
        if(++visited>limit)return false;
        if(state==0&&(i==n||j==m)&&score>=threshold&&!diagonal)return false;
        if(i==n||j==m)return true; // No later aligned terminal is possible.
        if(inside(i+1,j+1)&&!self(self,i+1,j+1,0,score+sub(i,j),diagonal))return false;
        if(state!=2&&inside(i+1,j)&&!self(self,i+1,j,1,score-(state==1?1:4),false))return false;
        if(state!=1&&inside(i,j+1)&&!self(self,i,j+1,2,score-(state==2?1:4),false))return false;
        return true;
    };
    for(int j=std::max(1,1+offset-radius);j<=std::min(m,1+offset+radius);++j)
        if(!visit(visit,1,j,0,sub(0,j-1),j-1==offset))return false;
    for(int i=std::max(2,1-offset-radius);i<=std::min(n,1-offset+radius);++i)
        if(!visit(visit,i,1,0,sub(i-1,0),1-i==offset))return false;
    return true;
}
inline bool full_offset(const std::string& a,const std::string& b,int& offset) {
    if(std::min(a.size(),b.size())<30)return false;
    if(a.size()<=b.size()){
        const auto p=b.find(a);if(p==std::string::npos)return false;offset=int(p);
    }else{
        const auto p=a.find(b);if(p==std::string::npos)return false;offset=-int(p);
    }
    return unique(a,b,offset,1024);
}
template<class Correspondence>
inline void fill(Correspondence& out,int n,int m,int offset) {
    const int ai=std::max(0,-offset),bj=std::max(0,offset),length=std::min(n-ai,m-bj);
    out.best=length;out.paths=1;out.eligible=true;out.reason="ok";
    out.invariant.reserve(length);out.certified.reserve(length-10);
    for(int k=0;k<length;++k){
        out.invariant.emplace_back(ai+k,bj+k);
        if(k>=5&&k<length-5)out.certified.emplace_back(ai+k,bj+k);
    }
}
}
