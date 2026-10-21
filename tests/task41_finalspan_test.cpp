#include "deferred_adapter_clip.hpp"
#include <cassert>
using namespace fqdup::r8qkeep;
int main(){
    RepairedRecord input;input.outcome="unmerged";
    input.first={"@pair/1 FQT=v1;ch=r1;m=1;olen=12;s=2;e=10 FQC=v1;edits=3:4:A:30;8:9:C:31","ACGTACGT","IIIIIIII"};
    input.second={"@pair/2 FQT=v1;ch=r2;m=2;olen=12;s=7;e=12","ACGTA","IIIII"};
    auto clipped=final_adapter_clip(input,6);
    assert(clipped.outcome=="orphan_r1");
    assert(clipped.first.seq=="ACGT"&&clipped.first.qual=="IIII");
    assert(clipped.first.header=="@pair/1 FQT=v1;ch=orphan;m=1;olen=12;s=2;e=6 FQC=v1;edits=3:4:A:30");
    assert(input.first.seq=="ACGTACGT");
    assert(final_adapter_clip(input,12).first.header==input.first.header);
    assert(final_adapter_clip(input,1).outcome=="dropped");
    input.first=input.second;input.second={"@pair/2 FQT=v1;ch=r2;m=2;olen=12;s=0;e=5","ACGTA","IIIII"};
    clipped=final_adapter_clip(input,6);assert(clipped.outcome=="orphan_r2");
    assert(clipped.first.header.find("ch=orphan;m=2")!=std::string::npos);
    input.outcome="merged";assert(final_adapter_clip(input,1).first.seq==input.first.seq);
    FamilyPair pair;NativeRecord native;native.raw=&pair;native.outcome="unmerged";
    native.first={"@pair/1","ACGTACGT","IIIIIIII"};native.second={"@pair/2","ACGTACGT","IIIIIIII"};
    RepairedRecord direct;direct.native=&native;clipped=final_adapter_clip(direct,6);
    assert(!clipped.native&&clipped.first.seq=="ACGTAC"&&clipped.second.seq=="ACGTAC");
    assert(clipped.first.header.find("olen=8;s=0;e=6")!=std::string::npos);
}
