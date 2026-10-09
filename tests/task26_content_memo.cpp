#include "../src/merge_r8qkeep/content_memo.hpp"
#include <iostream>
#include <random>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("content memo mismatch");}
int main(){
    FamilyPair pair;pair.id="x";
    for(int m=0;m<2;++m){pair.r[m].h="@x/"+std::to_string(m+1);pair.r[m].id="x";pair.r[m].s="ACGTNry";pair.r[m].q="!~I9-+I";}
    std::string key;memo_content(pair,key);ContentMemo cache;cache.context(1);
    MemoValue value;value.record.outcome="merged";value.record.first={"@task26 FQT=v1;ch=merged;m=0;olen=7;s=0;e=7","ACGTNry","!~I9-+I"};value.native_clean=true;
    const auto h=ContentMemo::hash(key);cache.put(key,h,value);require(cache.find(key,h)!=nullptr);
    require(cache.find(key+"different",h)==nullptr); // Deliberate full-hash collision.
    FamilyPair renamed=pair;renamed.id="y";renamed.r[0].h="@y/1 comment FQT=v1;old";renamed.r[1].h="@y/2 other";
    std::string other;memo_content(renamed,other);require(key==other);
    auto output=value;memo_restore_header(output.record.first.header,renamed);
    require(output.record.first.header==renamed.r[0].h+" FQT=v1;ch=merged;m=0;olen=7;s=0;e=7");
    for(int m=0;m<2;++m){std::string header="@task26/"+std::to_string(m+1)+" FQC=verbatim";memo_restore_header(header,renamed);require(header==renamed.r[m].h+" FQC=verbatim");}
    renamed.r[0].h="@y/1";std::string merged="@task26 tag";memo_restore_header(merged,renamed);require(merged=="@y tag");
    cache.discard();require(cache.find(key,h)==nullptr);cache.put(key,h,value);require(cache.find(key,h)!=nullptr);
    cache.context(2);require(cache.find(key,h)==nullptr);
    std::mt19937 rng(2603);
    for(unsigned length=0;length<=512;++length){
        for(auto& r:pair.r){r.s.clear();r.q.clear();for(unsigned j=0;j<length;++j){r.s+="ACGTNry"[rng()%7];r.q+=char(rng()%256);}}
        memo_content(pair,key);
        if(length){auto changed=pair;changed.r[0].q[length/2]^=1;memo_content(changed,other);require(key!=other);
            changed=pair;changed.r[1].s[length/2]=changed.r[1].s[length/2]=='N'?'A':'N';memo_content(changed,other);require(key!=other);}
    }
    NativeRecord native;native.outcome="unmerged";native.first={"@task26/1","AC","!!"};native.second={"@task26/2","GT","~~"};
    RepairedRecord borrowed;borrowed.native=&native;const auto owned=memo_owned(borrowed,true,true);
    require(!owned.record.native&&owned.record.first.header=="@task26/1 FQT=v1;ch=r1;m=1;olen=2;s=0;e=2");
    require(owned.record.second.header=="@task26/2 FQT=v1;ch=r2;m=2;olen=2;s=0;e=2");
    std::cout<<"content keys, collision checks, quality edges, header provenance and owned results PASS\n";
}
