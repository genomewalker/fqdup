#include "one_mate_adapter.hpp"
#include "deferred_adapter_clip.hpp"
#include <iostream>
using namespace fqdup::r8qkeep;
void check(bool x) { if (!x) throw std::runtime_error("one-mate regression"); }
int main() {
 const std::string a="AGATCGGAAGAGCACACGTCTGAACTCCAG",b="AGATCGGAAGAGCGTCGTGTAGGGAAAGAG";
 FamilyPair p; p.r[0].s=std::string(25,'C')+a; p.r[1].s=std::string(55,'T');
 for(auto& r:p.r)r.q=std::string(r.s.size(),'I');
 auto clip=[&]{return one_mate_adapter_boundary(p,a,b,1);};
 check(clip().boundary==25&&clip().mate==2);
 p.r[0].s[27]='T';check(clip().boundary==-1);
 p.r[0].q[27]='4';check(clip().boundary==25); // Q19 mismatch allowed
 p.r[0].q[27]='5';check(clip().boundary==-1); // Q20 rejected
 p.r[0].s[27]='N';p.r[0].q[27]='!';check(clip().boundary==-1);
 p.r[0].s=std::string(25,'C')+a;
 p.r[0].s+=std::string(10,'T')+a;p.r[0].q=std::string(p.r[0].s.size(),'I');check(clip().boundary==-1);
 p.r[0].s.back()='T';p.r[0].s[p.r[0].s.size()-2]='T';check(clip().boundary==-1); // second near placement
 p.r[0].s=std::string(25,'C')+a;p.r[0].q=std::string(55,'I');
 p.r[1].s=std::string(24,'T')+b; p.r[1].q=std::string(54,'I');check(clip().boundary==-1);
 p.r[1].s[24]='T';p.r[1].s[25]='T';check(clip().boundary==-1); // conflicting near
 p.r[1].s=std::string(25,'T');p.r[1].q=std::string(25,'I');check(clip().boundary==-1);
 p.r[1].s=std::string(55,'T');p.r[1].q=std::string(55,'I');
 check(one_mate_adapter_boundary(p,a,b,26).boundary==-1);
 p.r[0].s.resize(36);p.r[0].q.resize(36);check(clip().boundary==-1); // 11-base stub
 p.r[0].s+=a[11];p.r[0].q+='I';check(clip().boundary==25); // exactly 12
 std::swap(p.r[0],p.r[1]);check(one_mate_adapter_boundary(p,b,a,1).mate==1);
 RepairedRecord r;r.outcome="unmerged";
 r.first={"@p FQT=v1;ch=r1;m=1;s=0;e=30",std::string(30,'A'),std::string(30,'I')};
 r.second={"@p FQT=v1;ch=r2;m=2;s=4;e=30 FQC=v1;edits=10:1:A:40;28:1:C:40",std::string(26,'C'),std::string(26,'H')};
 auto c=final_adapter_clip(r,25,2);check(c.first.seq==r.first.seq&&c.second.seq.size()==21&&c.second.qual==std::string(21,'H')&&c.second.header.find("28:")==std::string::npos);
 r.outcome="merged";check(final_adapter_clip(r,25,2).first.seq==r.first.seq);
 r.outcome="orphan_r1";check(final_adapter_clip(r,25,2).first.seq==r.first.seq);
 r.outcome="unmerged";check(final_adapter_clip(r,4,2).outcome=="unmerged");
 std::cout<<"ONE_MATE_ADAPTER_PASS\n";
}
