// Targeted runtime-evidence diagnostic; truth is used only to score clipping.
#define main task41_attribution_main
#include "task41_attribution.cpp"
#undef main
#include <algorithm>
#include <set>
vector<string> cols(string s) {
 for(size_t p=0;(p=s.find("\\t",p))!=string::npos;)s.replace(p,2,"\t");
 vector<string> v;istringstream in(s);string x;while(getline(in,x,'\t'))v.push_back(x);return v;
}
int start(const string&s,const string&q,const string&a,int boundary,int mode) {
 int result=-1;
 for(int p=max(0,boundary);p<int(s.size());++p) {
  int n=min(int(a.size()),int(s.size())-p);if(n<11||(mode==0&&n<int(a.size())))continue;
  bool ok=true;for(int j=0;j<n;++j)if(s[p+j]!=a[j]) {
   if(mode!=2||s[p+j]=='N'||int((unsigned char)q[p+j])-33>=20){ok=false;break;}
  }
  if(ok){if(result>=0)return -1;result=p;}
 }
 return result;
}
int main(int argc,char**argv){try{
 require(argc==5,"usage: probe q0-dir truth-dir retained-output-dir report-dir");
 const string a[2]={"AGATCGGAAGAGCACACGTCTGAACTCCAG","AGATCGGAAGAGCGTCGTGTAGGGAAAGAG"};
 const char* names[]={"complete_exact","visible_exact_min11","visible_lowQ20_min11"};
 ofstream detail(string(argv[4])+"/pairs.tsv"),summary(string(argv[4])+"/summary.tsv");
 require(detail&&summary,"report open");
 detail<<"library\tid\tmate\tanchored_boundary\ttruth_adapter_start\tvisible\tmismatches\tN\tlowQ\tcomplete_exact\tvisible_exact\tvisible_lowQ\n";
 summary<<"rule\tpairs\tchanged_mates\tadapter_removed\tmolecule_removed\tinterior_molecule_removed\n";
 long long changed[3]={},adapter[3]={},molecule[3]={},interior[3]={},pairs=0;
 for(string lib:{"L1","L2","L3","L4"}) {
  map<string,int> targets;ifstream q0(string(argv[1])+"/"+lib+".fullpass.tsv");require(bool(q0),"q0 open");string line;
  while(getline(q0,line)){auto v=cols(line);if(v.size()>6&&v[3]=="unweighted_sequence_floor")require(targets.emplace(v[1],stoi(v[6])).second,"duplicate target");}
  // The 148 set-A pairs are floor-rejected mates retaining adapter, not every floor rejection.
  map<string,array<pair<int,int>,2>> spans;
  string base=string(argv[3])+"/outputs/"+lib+"-parent/0.316";
  for(int m=0;m<2;++m)readfq(base+"/r"+to_string(m+1)+".fq",[&](const string&h,const string&,const string&){string id=idof(h);if(!targets.count(id))return;auto f=fields(h.substr(h.rfind(" FQT=v1;")+1));spans[id][m]={stoi(f.at("s")),stoi(f.at("e"))};});
  map<string,array<string,4>> runtime;
  ifstream nat(string(argv[3])+"/outputs/"+lib+"-parent/debug/native.tsv");require(bool(nat),"native open");
  while(getline(nat,line)){auto v=cols(line);string id=idof(v.at(0));if(targets.count(id))runtime[id]={v.at(1),v.at(2),v.at(4),v.at(5)};}
  ifstream truth(string(argv[2])+"/"+lib+".jsonl");require(bool(truth),"truth open");
  while(getline(truth,line)) {
   // Avoid parsing unrelated, large per-base truth records.
   auto endid=line.find('"',7);string id=line.substr(7,endid-7);if(!targets.count(id)||!spans.count(id))continue;
   Parser parser{line};auto j=parser.value();long long frag=stoll(j.at("frag_len").s);
   bool retains=false;for(int m=0;m<2;++m){auto&orig=j.at("reads").a[m].at("origins").a;for(int i=spans[id][m].first;i<spans[id][m].second;++i)retains|=orig.at(i).at("kind").s=="adapter";}
   if(!retains)continue;++pairs;int clips[3][2];
   for(int m=0;m<2;++m){auto&r=runtime.at(id);auto&s=r[m*2];auto&q=r[m*2+1];auto&orig=j.at("reads").a[m].at("origins").a;
    int p=0;while(p<int(orig.size())&&orig[p].at("kind").s!="adapter")++p;
    int n=min(int(a[m].size()),int(s.size())-p),mm=0,ns=0,low=0;
    for(int k=0;k<n;++k){mm+=s[p+k]!=a[m][k];ns+=s[p+k]=='N';low+=int((unsigned char)q[p+k])-33<20;}
    detail<<lib<<'\t'<<id<<'\t'<<m+1<<'\t'<<targets[id]<<'\t'<<p<<'\t'<<n<<'\t'<<mm<<'\t'<<ns<<'\t'<<low;
    for(int rule=0;rule<3;++rule){clips[rule][m]=start(s,q,a[m],targets[id],rule);detail<<'\t'<<clips[rule][m];}detail<<'\n';
   }
   for(int rule=0;rule<3;++rule){int clip=clips[rule][0];if(clip<0||clip!=clips[rule][1])continue;
    for(int m=0;m<2;++m){auto&orig=j.at("reads").a[m].at("origins").a;auto span=spans[id][m];changed[rule]+=clip<span.second;
     for(int p=max(clip,span.first);p<span.second;++p){auto&o=orig.at(p);if(o.at("kind").s=="adapter")++adapter[rule];else{++molecule[rule];auto c=stoll(o.at("coordinate").s);interior[rule]+=c>=0&&c<frag;}}
    }
   }
  }
 }
 require(pairs==148,"target set must contain exactly 148 pairs, got "+to_string(pairs));
 for(int r=0;r<3;++r)summary<<names[r]<<'\t'<<pairs<<'\t'<<changed[r]<<'\t'<<adapter[r]<<'\t'<<molecule[r]<<'\t'<<interior[r]<<'\n';

 return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
