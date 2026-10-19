// Task35: strict raw-coordinate scoring; no alignment-based remapping.
#include <array>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
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
struct Rec{string h,s,q;int mate;string channel;};
struct Counts{long long n=0,k=0,em=0,ea=0,t=0,bad=0,at=0,raw=0;};
map<string,string> fields(string s,char sep=';'){map<string,string> m;istringstream in(s);string x;while(getline(in,x,sep)){auto p=x.find('=');if(p!=string::npos)m[x.substr(0,p)]=x.substr(p+1);}return m;}
template<class F> void readfq(const string& path,const F& fn){ifstream in(path);require(bool(in),"open "+path);string h,s,p,q;while(getline(in,h)){require(bool(getline(in,s))&&bool(getline(in,p))&&bool(getline(in,q)),"truncated FASTQ");require(!h.empty()&&h[0]=='@'&&!p.empty()&&p[0]=='+'&&s.size()==q.size(),"invalid FASTQ");fn(h,s,q);}}
int main(int argc,char** argv){try{
 require(argc==6,"usage: scorer truth.jsonl output-dir arm library family-edits.tsv");string dir=argv[2],arm=argv[3],lib=argv[4];
 struct Edit {char before,after;int quality;};
 map<string,map<int,map<size_t,Edit>>> family;
 ifstream ledger(argv[5]);require(bool(ledger),"family ledger open");string row;getline(ledger,row);
 while(getline(ledger,row)){istringstream in(row);string id;int mate,oq,q;size_t pos;char before,after;
 require(bool(in>>id>>mate>>pos>>before>>oq>>after>>q),"family ledger syntax");
 require((mate==1||mate==2)&&family[id][mate].emplace(pos,Edit{before,after,q}).second,"duplicate family edit");}
 set<string> truth_ids;

 unordered_map<string,vector<Rec>> outputs;set<string> merged;
 readfq(dir+"/merged.fq",[&](auto h,auto,auto){require(merged.insert(idof(h)).second,"duplicate merged ID");});
 for(string ch:{"r1","r2","orphan"})readfq(dir+"/"+ch+".fq",[&](auto h,auto s,auto q){int mate=ch=="r1"?1:ch=="r2"?2:0;auto pos=h.rfind(" FQT=v1;");if(pos!=string::npos){auto f=fields(h.substr(pos+1,h.find(' ',pos+1)-pos-1));int m=stoi(f.at("m"));require(!mate||m==mate,"FQT mate mismatch");mate=m;}if(!mate){string token=h.substr(0,h.find(' '));if(token.size()>2&&token[token.size()-2]=='/')mate=token.back()-'0';}require(mate>=0&&mate<=2,"mate index");auto id=idof(h);require(!merged.count(id),"merged/mate collision");for(auto& r:outputs[id])require(r.mate!=mate,"duplicate output mate");outputs[id].push_back({h,s,q,mate,ch});});
 ifstream truth(argv[1]);require(bool(truth),"truth open");string line;map<string,Counts> counts;long long dropped=0,total=0,ambiguous_orphans=0,ambiguous_spans=0;bool example=false;
 while(getline(truth,line)){Parser p{line};J j=p.value();p.ws();require(p.p==line.size(),"JSON trailing data");string id=j.at("id").s;require(truth_ids.insert(id).second,"duplicate truth ID");++total;const auto& reads=j.at("reads").a;require(reads.size()==2,"truth mate count");auto it=outputs.find(id);if(it==outputs.end()){if(!merged.erase(id))++dropped;continue;}
 for(auto& r:it->second){if(!r.mate){bool in1=reads.at(0).at("sequence").s.find(r.s)!=string::npos,in2=reads.at(1).at("sequence").s.find(r.s)!=string::npos;require(in1||in2,"orphan matches neither mate: "+id);r.mate=in1?1:2;if(in1&&in2)++ambiguous_orphans;}const J& raw=reads.at(r.mate-1);string seq=raw.at("sequence").s;const auto& orig=raw.at("origins").a;require(seq.size()==orig.size(),"truth origin length");for(const auto& item:family[id][r.mate]){require(item.first<seq.size()&&seq[item.first]==item.second.before,"family raw base mismatch");seq[item.first]=item.second.after;}
 size_t begin=0,end=r.s.size();auto tag=r.h.rfind(" FQT=v1;");if(tag!=string::npos){auto f=fields(r.h.substr(tag+1,r.h.find(' ',tag+1)-tag-1));begin=stoull(f.at("s"));end=stoull(f.at("e"));require(stoull(f.at("olen"))<=seq.size(),"FQT original length");}else if(seq.compare(0,r.s.size(),r.s)!=0){auto found=seq.find(r.s);require(found!=string::npos,"untagged span not found: "+id);if(seq.find(r.s,found+1)!=string::npos)++ambiguous_spans;begin=found;end=found+r.s.size();}require(end>=begin&&end-begin==r.s.size()&&end<=seq.size(),"retained span");map<size_t,char> edits;auto qc=r.h.rfind(" FQC=v1;edits=");if(qc!=string::npos){istringstream in(r.h.substr(qc+14));string edit;while(getline(in,edit,';')){for(char& c:edit)if(c==':')c=' ';istringstream e(edit);size_t pos,donor;char base;int qual;require(bool(e>>pos>>donor>>base>>qual),"FQC edit syntax");require(pos>=begin&&pos<end&&edits.emplace(pos,base).second,"FQC edit range");require(r.q.at(pos-begin)==char(qual+33),"FQC quality");}}
 auto& c=counts[r.channel+"/"+to_string(r.mate)];++c.n;c.raw+=seq.size();bool bad=false;long long rawmol=0,rawad=0,keepmol=0,trim=0,keepad=0,trimad=0;
 for(size_t i=0;i<seq.size();++i){const auto& o=orig[i];string kind=o.at("kind").s;bool mol=kind=="molecule";require(mol||kind=="adapter","unknown origin");mol?++rawmol:++rawad;if(i>=begin&&i<end){char b=r.s[i-begin];auto e=edits.find(i);require(b==(e==edits.end()?seq[i]:e->second),"untraced base change: "+id+"/"+to_string(r.mate)+" at "+to_string(i));if(mol){++keepmol;if(string(1,b)==o.at("truth_base").s)++c.k;else{++c.em;bad=true;}}else{++keepad;++c.ea;bad=true;}}else if(mol){++c.t;++trim;}else{++c.at;++trimad;}}
 require(keepmol+trim==rawmol&&keepad+trimad==rawad,"raw conservation");c.bad+=bad;
 if(!example&&(begin||end<seq.size())){cerr<<"EXAMPLE\t"<<lib<<'\t'<<arm<<'\t'<<id<<'\t'<<r.mate<<"\traw="<<seq.size()<<"\tspan=["<<begin<<','<<end<<")\tmolecule="<<rawmol<<"\tadapter="<<rawad<<"\tkept_molecule="<<keepmol<<"\ttrimmed_molecule="<<trim<<"\tkept_adapter="<<keepad<<"\ttrimmed_adapter="<<trimad<<'\n';example=true;}
 }outputs.erase(it);
 }
 require(outputs.empty()&&merged.empty(),"output IDs missing from truth");for(string key:{"r1/1","r2/2","orphan/1","orphan/2"}){auto c=counts[key];require(c.k+c.em+c.ea+c.t+c.at==c.raw,"aggregate conservation");cout<<arm<<'\t'<<lib<<'\t'<<key<<'\t'<<c.n<<'\t'<<c.k<<'\t'<<c.em<<'\t'<<c.ea<<'\t'<<c.t<<'\t'<<c.bad<<'\t'<<c.at<<'\t'<<c.raw<<'\n';}cerr<<"PAIRS\t"<<lib<<'\t'<<arm<<"\ttotal="<<total<<"\tdropped="<<dropped<<"\tambiguous_orphans="<<ambiguous_orphans<<"\tambiguous_spans="<<ambiguous_spans<<'\n';
 }catch(const exception& e){cerr<<e.what()<<'\n';return 1;}}
