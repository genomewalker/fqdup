// Task 38: per-truth-pair, all-channel attribution of fqdup and three AR3 arms (read-only on existing outputs).
#include <array>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
void require(bool b,const string& m){if(!b)throw runtime_error(m);}
struct J {string s;map<string,J> o;vector<J> a;const J& at(const string& k)const{return o.at(k);}};
struct Parser {
 const string& t;size_t p=0;
 void ws(){while(p<t.size()&&isspace((unsigned char)t[p]))++p;}
 char take(){require(p<t.size(),"JSON end");return t[p++];}
 string str(){require(take()=='"',"JSON string");string s;for(char c=take();c!='"';c=take()){if(c=='\\'){c=take();require(c=='"'||c=='\\'||c=='/',"unsupported JSON escape");}s+=c;}return s;}
 J value(){ws();J v;char c=take();if(c=='"'){--p;v.s=str();}else if(c=='{'||c=='['){ws();char end=c=='{'?'}':']';if(p<t.size()&&t[p]==end){++p;return v;}do{ws();string k;if(c=='{'){k=str();ws();require(take()==':',"JSON colon");}auto x=value();if(c=='{')require(v.o.emplace(k,move(x)).second,"duplicate JSON key");else v.a.push_back(move(x));ws();char sep=take();if(sep==end)break;require(sep==',',"JSON separator");}while(true);}else{--p;while(p<t.size()&&t[p]!=','&&t[p]!=']'&&t[p]!='}'&&!isspace((unsigned char)t[p]))v.s+=t[p++];require(!v.s.empty(),"JSON scalar");}return v;}
};
string idof(string h){h=h.substr(0,h.find_first_of(" \t"));if(!h.empty()&&h[0]=='@')h.erase(0,1);if(h.size()>2&&h[h.size()-2]=='/'&&(h.back()=='1'||h.back()=='2'))h.resize(h.size()-2);return h;}
map<string,string> fields(const string& s){map<string,string> m;istringstream in(s);string x;while(getline(in,x,';')){auto p=x.find('=');if(p!=string::npos)m[x.substr(0,p)]=x.substr(p+1);}return m;}
template<class F> void readfq(const string& path,const F& fn){ifstream in(path);require(bool(in),"open "+path);string h,s,p,q;while(getline(in,h)){require(bool(getline(in,s))&&bool(getline(in,p))&&bool(getline(in,q)),"truncated FASTQ "+path);require(!h.empty()&&h[0]=='@'&&s.size()==q.size(),"invalid FASTQ "+path);fn(h,s,q);}}
size_t occurrences(const string& hay,const string& needle){size_t n=0;for(size_t p=hay.find(needle);p!=string::npos;p=hay.find(needle,p+1))++n;return n;}

struct Rec{string h,s,q;int mate;string channel;};
struct Edit{char before,after;};
struct Arm{
 string name;bool tagged;
 unordered_map<string,string> merged;
 unordered_map<string,vector<Rec>> mates;
 unordered_map<string,array<map<size_t,Edit>,3>> family;
 unordered_map<string,string> route;
};
struct MateStats{long long n=0,k=0,em=0,ea=0,t=0,at=0,raw=0,bad=0;
 void add(const MateStats& o){n+=o.n;k+=o.k;em+=o.em;ea+=o.ea;t+=o.t;at+=o.at;raw+=o.raw;bad+=o.bad;}};
struct Channel{long long opportunity=0,cov_ok=0,cov_bad=0,cov_conflict=0,omitted=0,copies_ok=0,copies_bad=0,adapter=0,synthetic=0,ambiguous_pairs=0;};
struct Amb{long long spans=0,prefix_also_later=0,orphans=0;};

void load(Arm& a,const string& dir,const string& ledger,const string& debug){
 readfq(dir+"/merged.fq",[&](const string& h,const string& s,const string&){require(a.merged.emplace(idof(h),s).second,"duplicate merged "+a.name);});
 for(string ch:{"r1","r2","orphan"})readfq(dir+"/"+ch+".fq",[&](const string& h,const string& s,const string& q){
  int mate=ch=="r1"?1:ch=="r2"?2:0;auto pos=h.rfind(" FQT=v1;");
  if(pos!=string::npos){auto f=fields(h.substr(pos+1,h.find(' ',pos+1)-pos-1));int m=stoi(f.at("m"));require(!mate||m==mate,"FQT mate mismatch");mate=m;}
  if(!mate){string token=h.substr(0,h.find(' '));if(token.size()>2&&token[token.size()-2]=='/')mate=token.back()-'0';}
  a.mates[idof(h)].push_back({h,s,q,mate,ch});});
 if(!ledger.empty()){ifstream in(ledger);require(bool(in),"ledger open");string row;getline(in,row);
  while(getline(in,row)){istringstream r(row);string id;int mate,oq,qq;size_t p;char b,af;require(bool(r>>id>>mate>>p>>b>>oq>>af>>qq),"ledger syntax");
   require((mate==1||mate==2)&&a.family[id][mate].emplace(p,Edit{b,af}).second,"duplicate family edit");}}
 if(!debug.empty()){ifstream nat(debug+"/native.tsv"),aud(debug+"/route-audit.tsv");require(nat&&aud,"debug open");string n,r;
  while(getline(nat,n)){require(bool(getline(aud,r)),"route-audit short");string cls=r.substr(r.rfind('\t')+1);a.route[idof(n.substr(0,n.find('\t')))]=cls;}}
}

string ovbin(long long f,long long rl){long long ov=f<=rl?f:2*rl-f;return ov<=0?"none":ov<11?"lt11":ov<=22?"11-22":ov<=50?"23-50":"gt50";}

int main(int argc,char** argv){try{
 require(argc==8,"usage: attribution lib truth.jsonl fqdup_dir AR3-speed_dir AR3-t29_dir AR3-default_dir outdir");
 string lib=argv[1],out=argv[7];
 vector<Arm> arms(4);arms[0].name="fqdup";arms[0].tagged=true;
 const char* names[]={"parent","fixed-main","fixed-main-copy"};
 load(arms[0],string(argv[3])+"/0.316",string(argv[3])+"/debug/family-edits.tsv",string(argv[3])+"/debug");
 for(int i=1;i<4;++i){arms[i].name=names[i-1];arms[i].tagged=true; string root=argv[3+i]; load(arms[i],root+"/0.316",root+"/debug/family-edits.tsv",root+"/debug");}

 map<string,long long> transitions,cohorts_pairs;map<string,MateStats> cohorts,attribution,totals_mates;map<string,long long> totals_merged;
 vector<Channel> chan(4);vector<Amb> amb(4);long long target_mismatch=0,truth_pairs=0;
 ofstream pairs_out(out+"/pairs.tsv");require(bool(pairs_out),"pairs open");
 ifstream truth(argv[2]);require(bool(truth),"truth open");string line;
 while(getline(truth,line)){
  Parser parser{line};J j=parser.value();++truth_pairs;
  string id=j.at("id").s,target=j.at("merged_target").s;bool should=target!="null"&&target!="-"&&!target.empty();
  long long frag=stoll(j.at("frag_len").s);
  const auto& reads=j.at("reads").a;require(reads.size()==2,"truth mates");
  long long rl=reads[0].at("sequence").s.size();
  map<long long,char> opportunity;
  for(const auto& r:reads)for(const auto& o:r.at("origins").a)if(o.at("kind").s=="molecule"){long long c=stoll(o.at("coordinate").s);opportunity[c]=1;}
  if(should)for(const auto& o:reads[0].at("origins").a)if(o.at("kind").s=="molecule"){long long c=stoll(o.at("coordinate").s);if(c>=(long long)target.size()||target[c]!=o.at("truth_base").s[0])++target_mismatch;}
  vector<string> state(4);vector<array<MateStats,3>> pm(4);vector<bool> ambiguous(4,false);
  for(int ai=0;ai<4;++ai){
   Arm& a=arms[ai];Channel& ch=chan[ai];
   unordered_map<long long,array<long long,2>> copies;long long adapter=0,synthetic=0;
   auto mit=a.merged.find(id);auto rit=a.mates.find(id);
   if(mit!=a.merged.end()){
    require(rit==a.mates.end(),"merged/mate collision "+a.name+" "+id);
    const string& m=mit->second;
    if(!should){state[ai]="FP";synthetic+=m.size();}
    else if(m.size()!=target.size()){state[ai]="wrong_len";synthetic+=m.size();}
    else{size_t e=0;for(size_t i=0;i<m.size();++i){bool ok=m[i]==target[i];e+=!ok;++copies[i][ok?0:1];}state[ai]=e?"base_err":"exact";}
    a.merged.erase(mit);
   }else if(rit==a.mates.end())state[ai]="dropped";
   else{
    state[ai]=rit->second.size()==2?"unmerged":"orphan";
    for(auto& r:rit->second){
     if(!r.mate){size_t in1=occurrences(reads[0].at("sequence").s,r.s),in2=occurrences(reads[1].at("sequence").s,r.s);
      require(in1||in2,"orphan matches neither mate "+a.name+" "+id);if(in1&&in2){++amb[ai].orphans;ambiguous[ai]=true;continue;}r.mate=in1?1:2;}
     const J& raw=reads.at(r.mate-1);string seq=raw.at("sequence").s;const auto& orig=raw.at("origins").a;
     auto fit=a.family.find(id);if(fit!=a.family.end())for(const auto& e:fit->second[r.mate]){require(e.first<seq.size()&&seq[e.first]==e.second.before,"family raw base mismatch");seq[e.first]=e.second.after;}
     size_t begin=0,end=r.s.size();auto tag=r.h.rfind(" FQT=v1;");
     if(tag!=string::npos){auto f=fields(r.h.substr(tag+1,r.h.find(' ',tag+1)-tag-1));begin=stoull(f.at("s"));end=stoull(f.at("e"));}
     else{size_t n=occurrences(seq,r.s);require(n>0,"untagged span not found "+a.name+" "+id);
      if(n>1){++amb[ai].spans;if(seq.compare(0,r.s.size(),r.s)==0)++amb[ai].prefix_also_later;ambiguous[ai]=true;continue;}
      begin=seq.find(r.s);end=begin+r.s.size();}
     require(end>=begin&&end-begin==r.s.size()&&end<=seq.size(),"retained span "+a.name+" "+id);
     map<size_t,char> edits;auto qc=r.h.rfind(" FQC=v1;edits=");
     if(qc!=string::npos){istringstream in(r.h.substr(qc+14));string ed;while(getline(in,ed,';')){for(char& c:ed)if(c==':')c=' ';istringstream e(ed);size_t pos,donor;char base;int qual;require(bool(e>>pos>>donor>>base>>qual),"FQC syntax");edits.emplace(pos,base);}}
     int slot=r.channel=="r1"?0:r.channel=="r2"?1:2;MateStats& s=pm[ai][slot];++s.n;s.raw+=seq.size();bool bad=false;
     for(size_t i=0;i<seq.size();++i){const auto& o=orig[i];bool mol=o.at("kind").s=="molecule";
      if(i>=begin&&i<end){char b=r.s[i-begin];auto e=edits.find(i);require(b==(e==edits.end()?seq[i]:e->second),"untraced base change "+a.name+" "+id);
       if(mol){bool ok=b==o.at("truth_base").s[0];ok?++s.k:++s.em;bad|=!ok;++copies[stoll(o.at("coordinate").s)][ok?0:1];}else{++s.ea;++adapter;bad=true;}}
      else if(mol)++s.t;else ++s.at;}
     require(s.k+s.em+s.ea+s.t+s.at==s.raw,"conservation "+a.name+" "+id);s.bad+=bad;
    }
    a.mates.erase(rit);
   }
   if(ambiguous[ai]){++ch.ambiguous_pairs;continue;}
   ch.opportunity+=opportunity.size();ch.adapter+=adapter;ch.synthetic+=synthetic;
   for(auto& c:copies)require(opportunity.count(c.first),"coordinate outside opportunity "+a.name+" "+id);
   for(auto& c:opportunity){auto it=copies.find(c.first);if(it==copies.end()){++ch.omitted;continue;}long long ok=it->second[0],bad=it->second[1];
    ch.copies_ok+=ok;ch.copies_bad+=bad;if(ok&&bad)++ch.cov_conflict;else if(ok)++ch.cov_ok;else ++ch.cov_bad;}
  }
  string sm=should?"1":"0";
  {MateStats s0;for(int k=0;k<3;++k)s0.add(pm[0][k]);auto rt=arms[0].route.find(id);pairs_out<<lib<<'\t'<<id<<'\t'<<sm<<'\t'<<frag<<'\t'<<ovbin(frag,rl)<<'\t'<<(rt==arms[0].route.end()?"missing":rt->second)<<'\t'<<state[0]<<'\t'<<state[1]<<'\t'<<state[2]<<'\t'<<state[3]<<'\t'<<s0.em<<'\t'<<s0.ea<<'\t'<<s0.t<<'\n';}
  for(int ai=0;ai<4;++ai){totals_merged[arms[ai].name+"\t"+state[ai]]++;if(!ambiguous[ai])for(int k=0;k<3;++k)totals_mates[arms[ai].name+"\t"+to_string(k)].add(pm[ai][k]);}
  for(int ai=1;ai<4;++ai)transitions[lib+"\t"+sm+"\t"+state[0]+"\t"+arms[ai].name+"\t"+state[ai]]++;
  bool shared=!ambiguous[0]&&!ambiguous[1]&&!ambiguous[2]&&!ambiguous[3];
  for(int ai=0;ai<4;++ai)shared=shared&&state[ai]=="unmerged";
  auto addc=[&](const string& cohort,int ai){MateStats s;for(int k=0;k<3;++k)s.add(pm[ai][k]);string key=lib+"\t"+cohort+"\t"+arms[ai].name;cohorts[key].add(s);cohorts_pairs[key]++;};
  if(shared)for(int ai=0;ai<4;++ai)addc("shared_unmerged",ai);
  for(int ai=1;ai<=2;++ai)if(!ambiguous[0]&&!ambiguous[ai]&&state[0]=="unmerged"&&state[ai]=="unmerged"){string cohort="shared_"+arms[ai].name;addc(cohort,0);addc(cohort,ai);}
  bool fq_un=state[0]=="unmerged"||state[0]=="orphan";
  for(int ai=1;ai<4;++ai){bool ar_merged=state[ai]=="exact"||state[ai]=="base_err"||state[ai]=="wrong_len"||state[ai]=="FP";
   bool ar_un=state[ai]=="unmerged"||state[ai]=="orphan";bool fq_merged=state[0]=="exact"||state[0]=="base_err"||state[0]=="wrong_len"||state[0]=="FP";
   if(fq_un&&ar_merged&&!ambiguous[0])addc(string("fqdup_unmerged_")+arms[ai].name+"_merged",0);
   if(ar_un&&fq_merged&&!ambiguous[ai])addc(arms[ai].name+"_unmerged_fqdup_merged",ai);}
  if(fq_un&&!ambiguous[0]){MateStats s;for(int k=0;k<3;++k)s.add(pm[0][k]);
   auto rt=arms[0].route.find(id);string route=rt==arms[0].route.end()?"missing":rt->second;
   attribution[lib+"\t"+sm+"\t"+route+"\t"+ovbin(frag,rl)+"\t"+(frag<rl?"adapter_in_read":"no_adapter")+"\t"+state[2]].add(s);}
 }
 for(auto& a:arms){require(a.merged.empty(),"merged IDs missing from truth "+a.name);require(a.mates.empty(),"mate IDs missing from truth "+a.name);}

 auto open=[&](const string& name){ofstream f(out+"/"+name,ios::app);require(bool(f),"write "+name);return f;};
 {auto f=open("transitions.tsv");for(auto& t:transitions)f<<t.first<<'\t'<<t.second<<'\n';}
 {auto f=open("allchannel.tsv");for(int ai=0;ai<4;++ai){auto& c=chan[ai];f<<lib<<'\t'<<arms[ai].name<<'\t'<<c.opportunity<<'\t'<<c.cov_ok<<'\t'<<c.cov_bad<<'\t'<<c.cov_conflict<<'\t'<<c.omitted<<'\t'<<c.copies_ok<<'\t'<<c.copies_bad<<'\t'<<c.adapter<<'\t'<<c.synthetic<<'\t'<<c.ambiguous_pairs<<'\n';}}
 auto ms=[](ostream& f,const MateStats& s){f<<s.n<<'\t'<<s.k<<'\t'<<s.em<<'\t'<<s.ea<<'\t'<<s.t<<'\t'<<s.at<<'\t'<<s.bad;};
 {auto f=open("cohorts.tsv");for(auto& c:cohorts){f<<c.first<<'\t'<<cohorts_pairs[c.first]<<'\t';ms(f,c.second);f<<'\n';}}
 {auto f=open("fqdup_attribution.tsv");for(auto& c:attribution){f<<c.first<<'\t';ms(f,c.second);f<<'\n';}}
 {auto f=open("ambiguity.tsv");for(int ai=0;ai<4;++ai)f<<lib<<'\t'<<arms[ai].name<<'\t'<<amb[ai].spans<<'\t'<<amb[ai].prefix_also_later<<'\t'<<amb[ai].orphans<<'\t'<<chan[ai].ambiguous_pairs<<'\n';}
 {auto f=open("totals.tsv");for(auto& t:totals_merged)f<<lib<<"\tmerged_state\t"<<t.first<<'\t'<<t.second<<'\n';
  const char* slot[]={"r1","r2","orphan"};for(auto& t:totals_mates){auto tab=t.first.find('\t');f<<lib<<"\tmates\t"<<t.first.substr(0,tab)<<'\t'<<slot[stoi(t.first.substr(tab+1))]<<'\t';ms(f,t.second);f<<'\n';}
  f<<lib<<"\tcheck\ttruth_pairs\t"<<truth_pairs<<"\n"<<lib<<"\tcheck\tr1_truth_base_vs_target_mismatch\t"<<target_mismatch<<'\n';}
 cerr<<"ATTRIBUTION_DONE "<<lib<<" pairs="<<truth_pairs<<'\n';
}catch(const exception& e){cerr<<"ERROR "<<e.what()<<'\n';return 1;}}
