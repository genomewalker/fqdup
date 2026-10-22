#define TASK47_SPAN_LIBRARY
#include "task47_stacked_clip_check.cpp"
#undef TASK47_SPAN_LIBRARY
#include "one_mate_adapter.hpp"
using namespace fqdup::r8qkeep;
int main(int argc,char**argv) {try {
 require(argc==7,"audit truth parent-root candidate-root cutoffs report-prefix library");
 ifstream cutsin(argv[4]);require(bool(cutsin),"cutoffs open");vector<string> cuts;string cut;
 vector<map<string,Fastq>> before,after;
 while(getline(cutsin,cut)) {cuts.push_back(cut);before.push_back(records(string(argv[2])+"/"+cut));after.push_back(records(string(argv[3])+"/"+cut));}
 require(!cuts.empty(),"empty cutoffs");
 ofstream out(string(argv[5])+"-audit.tsv"),changes(string(argv[5])+"-changes.tsv"),noise(string(argv[5])+"-noise.tsv");
 require(out&&changes&&noise,"audit output");
 out<<"library\tcutoff\tpairs\trule_candidates\tchanged_mates\tadapter_removed\tmolecule_removed\tn_lt_frag\tlong_insert\tdiscordant\tquality_changed\tfirst10_loss\tlast10_loss\tC5_first10_loss\tG3_last10_loss\tverdict\n";
 changes<<"library\tcutoff\tid\tmate\tn\tfrag_len\told_end\tnew_end\tadapter_removed\n";
 struct Counts {long long changed=0,adapter=0,molecule=0,shortclip=0,longclip=0,discordant=0,quality=0,first=0,last=0,c5=0,g3=0;};
 vector<Counts> count(cuts.size());long long pairs=0,candidates=0;map<int,long long> hist;vector<map<char,long long>> consensus(30);
 long long noise_pairs=0;map<int,long long> mismatches;
 const string a="AGATCGGAAGAGCACACGTCTGAACTCCAG",b="AGATCGGAAGAGCGTCGTGTAGGGAAAGAG";
 ifstream in(argv[1]);require(bool(in),"truth open");string line;
 while(getline(in,line)) {
  Parser parser{line};auto j=parser.value();++pairs;const string id=j.at("id").s;const int frag=stoi(j.at("frag_len").s);
  const auto& reads=j.at("reads").a;FamilyPair raw;raw.id=id;
  map<int,char> truebase;
  for(int m=0;m<2;++m){raw.r[m].s=reads[m].at("sequence").s;raw.r[m].q=reads[m].at("quality").s;
   for(const auto& o:reads[m].at("origins").a)if(o.at("kind").s=="molecule"){int pos=stoi(o.at("coordinate").s);char base=o.at("truth_base").s[0];if(m==1){const string bases="ACGT",rc="TGCA";auto at=bases.find(base);if(at!=string::npos)base=rc[at];}if(!truebase.count(pos))truebase[pos]=base;}}
  const bool c5=truebase.count(0)&&truebase[0]=='C',g3=truebase.count(frag-1)&&truebase[frag-1]=='G';
  auto rule=one_mate_adapter_boundary(raw,a,b,1);
  if(rule.boundary>=0) {++candidates;++hist[rule.boundary];
   if(rule.mate==2) {int mm=0,n=min(int(b.size()),int(raw.r[1].s.size())-rule.boundary);
    for(int k=0;k<n;++k){char base=raw.r[1].s[rule.boundary+k];mm+=base!=b[k];++consensus[k][base];}++mismatches[mm];++noise_pairs;}}
  for(size_t ci=0;ci<cuts.size();++ci) {
   auto& p=before[ci];auto& c=after[ci];auto& s=count[ci];
   for(string mate:{"merged","1","2"}) {
    const string key=id+"\t"+mate;auto pi=p.find(key),qi=c.find(key);
    if(pi==p.end()&&qi==c.end())continue;
    require(pi!=p.end()&&qi!=c.end(),"channel added/dropped "+id+" "+cuts[ci]);
    const auto& old=pi->second;const auto& now=qi->second;
    if(mate=="merged")require(old.h==now.h&&old.s==now.s&&old.q==now.q,"merged changed "+id);
    else if(old.h!=now.h||old.s!=now.s||old.q!=now.q||old.ch!=now.ch) {
     ++s.changed;auto os=span(old),ns=span(now);const int m=stoi(mate)-1;
     require(ns.first>=os.first&&ns.second<=os.second,"span expansion "+id);
     s.quality+=old.q.substr(ns.first-os.first,ns.second-ns.first)!=now.q;
     require(old.s.substr(ns.first-os.first,ns.second-ns.first)==now.s,"retained bases changed "+id);
     require(old.ch==now.ch&&p.count(id+"\t1")&&p.count(id+"\t2"),"not a retained proper pair "+id);
     s.shortclip+=int(ns.second)<frag;s.longclip+=frag>int(raw.r[m].s.size());
     s.discordant+=rule.boundary!=int(ns.second)||rule.mate!=stoi(mate);
     long long ad=0;
     const auto& orig=reads[m].at("origins").a;
     for(size_t k=os.first;k<os.second;++k)if(k<ns.first||k>=ns.second){const auto& o=orig.at(k);
      if(o.at("kind").s=="adapter"){++ad;++s.adapter;}else{++s.molecule;int pos=stoi(o.at("coordinate").s);s.first+=pos<10;s.last+=pos>=frag-10;s.c5+=c5&&pos<10;s.g3+=g3&&pos>=frag-10;}}
     changes<<argv[6]<<'\t'<<cuts[ci]<<'\t'<<id<<'\t'<<mate<<'\t'<<rule.boundary<<'\t'<<frag<<'\t'<<os.second<<'\t'<<ns.second<<'\t'<<ad<<'\n';
    }
    // Keep both mates in the map until all per-pair checks have completed.
   }
   for(string mate:{"merged","1","2"}){p.erase(id+"\t"+mate);c.erase(id+"\t"+mate);}
  }
 }
 bool bad=false;
 for(size_t ci=0;ci<cuts.size();++ci){require(before[ci].empty()&&after[ci].empty(),"unknown output IDs");const auto& s=count[ci];bool ok=!(s.molecule||s.shortclip||s.longclip||s.discordant||s.quality||s.first||s.last||s.c5||s.g3);bad|=!ok;
  out<<argv[6]<<'\t'<<cuts[ci]<<'\t'<<pairs<<'\t'<<candidates<<'\t'<<s.changed<<'\t'<<s.adapter<<'\t'<<s.molecule<<'\t'<<s.shortclip<<'\t'<<s.longclip<<'\t'<<s.discordant<<'\t'<<s.quality<<'\t'<<s.first<<'\t'<<s.last<<'\t'<<s.c5<<'\t'<<s.g3<<'\t'<<(ok?"PASS":"FAIL")<<'\n';}
 noise<<"R1_anchor_R2_silent_pairs\t"<<noise_pairs<<"\nposition\tA\tC\tG\tT\tN\tconsensus\tdeclared\n";
 for(size_t k=0;k<consensus.size();++k){char best='N';long long top=0;noise<<k;for(char base:string("ACGTN")){noise<<'\t'<<consensus[k][base];if(consensus[k][base]>top){top=consensus[k][base];best=base;}}noise<<'\t'<<best<<'\t'<<b[k]<<'\n';}
 noise<<"mismatches\tpairs\n";for(auto v:mismatches)noise<<v.first<<'\t'<<v.second<<'\n';
 noise<<"n\trule_candidates\n";for(auto v:hist)noise<<v.first<<'\t'<<v.second<<'\n';
 cout<<argv[6]<<" pairs="<<pairs<<" candidates="<<candidates<<" cutoffs="<<cuts.size()<<" verdict="<<(bad?"FAIL":"PASS")<<'\n';return bad?1:0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 2;}}
