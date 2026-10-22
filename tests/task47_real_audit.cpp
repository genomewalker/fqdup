#define TASK47_SPAN_LIBRARY
#include "task47_stacked_clip_check.cpp"
#undef TASK47_SPAN_LIBRARY
#include "one_mate_adapter.hpp"
using namespace fqdup::r8qkeep;
int main(int argc,char**argv){try{
 require(argc==6,"real-audit parent candidate adapter2 prefix library");
 auto p=records(argv[1]),c=records(argv[2]);ifstream in(string(argv[2])+"/debug/native.tsv");require(bool(in),"raw native open");
 ofstream out(string(argv[4])+"-fires.tsv"),dist(string(argv[4])+"-n.tsv");require(out&&dist,"output open");
 out<<"id\tn\tsilent_mate\told_end\tnew_end\tsilent_adapter\tmismatches\n";
 long long pairs=0,eligible=0,fires=0,quality=0;map<int,long long> hist;string line;
 while(getline(in,line)){
  ++pairs;vector<string> f;istringstream row(line);string v;while(getline(row,v,'\t'))f.push_back(v);require(f.size()>=7,"native columns");
  FamilyPair raw;raw.id=idof(f[0]);raw.r[0].s=f[1];raw.r[0].q=f[2];raw.r[1].s=f[4];raw.r[1].q=f[5];
  auto rule=one_mate_adapter_boundary(raw,"AGATCGGAAGAGCACACGTCTGAACTCCAG",argv[3],1);eligible+=rule.boundary>=0;
  for(string mate:{"merged","1","2"}){
   const string key=raw.id+"\t"+mate;auto pi=p.find(key),ci=c.find(key);if(pi==p.end()&&ci==c.end())continue;
   require(pi!=p.end()&&ci!=c.end(),"channel transition "+key);
   const auto& a=pi->second;const auto& b=ci->second;
   if(mate=="merged")require(a.h==b.h&&a.s==b.s&&a.q==b.q,"merged changed");
   else if(a.h!=b.h||a.s!=b.s||a.q!=b.q){
    auto os=span(a),ns=span(b);require(a.ch==b.ch&&a.ch!="orphan","nonpair changed");
    require(rule.boundary==int(ns.second)&&rule.mate==stoi(mate)&&ns.first==os.first&&ns.second<os.second,"wrong clip "+key);
    require(a.s.substr(0,b.s.size())==b.s,"retained sequence changed");quality+=a.q.substr(0,b.q.size())!=b.q;
    ++fires;++hist[rule.boundary];const string adapter=rule.mate==1?"AGATCGGAAGAGCACACGTCTGAACTCCAG":argv[3];
    const string tail=raw.r[rule.mate-1].s.substr(rule.boundary,adapter.size());int mm=0;for(size_t k=0;k<tail.size();++k)mm+=tail[k]!=adapter[k];
    out<<raw.id<<'\t'<<rule.boundary<<'\t'<<rule.mate<<'\t'<<os.second<<'\t'<<ns.second<<'\t'<<tail<<'\t'<<mm<<'\n';
   }
   p.erase(pi);c.erase(ci);
  }
 }
 require(p.empty()&&c.empty(),"unknown IDs");dist<<"n\tfires\n";for(auto v:hist)dist<<v.first<<'\t'<<v.second<<'\n';
 cout<<argv[5]<<" pairs="<<pairs<<" eligible_raw="<<eligible<<" changed_mates="<<fires<<" fire_rate="<<double(fires)/pairs<<" retained_quality_changes="<<quality<<" verdict="<<(quality?"FAIL":"PASS")<<'\n';return quality?1:0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 2;}}
