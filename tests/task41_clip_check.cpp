// T41-1: exhaustively check changed records against raw truth origins.
#define main task41_attribution_main
#include "task41_attribution.cpp"
#undef main
struct Fastq { string h,s,q; };
map<string,Fastq> records(const string& dir) {
 map<string,Fastq> result;
 for(string ch:{"merged","r1","r2","orphan"})
  readfq(dir+"/"+ch+".fq",[&](const string& h,const string& s,const string& q){
   string mate=ch;
   if(ch=="orphan") { auto p=h.rfind(" FQT=v1;");require(p!=string::npos,"orphan provenance");mate+=fields(h.substr(p+1)).at("m"); }
   require(result.emplace(idof(h)+"\t"+mate,Fastq{h,s,q}).second,"duplicate record");
  });
 return result;
}
int main(int argc,char**argv) {try {
 require(argc==5,"usage: clip_check truth parent-dir candidate-dir report");
 auto parent=records(argv[2]),candidate=records(argv[3]);
 ofstream out(argv[4]); require(bool(out),"report open");
 out<<"id\tchannel\tsequence_changed\tquality_changed\tlength_changed\tadapter_removed\tmolecule_removed\tinterior_molecule_removed\n";
 long long changed=0,adapter=0,molecule=0,interior=0,failed=0;
 ifstream in(argv[1]);require(bool(in),"truth open");string line;
 while(getline(in,line)) {
  Parser parser{line};auto j=parser.value();string id=j.at("id").s;auto frag=stoll(j.at("frag_len").s);
  for(string ch:{"merged","r1","r2","orphan1","orphan2"}) {
   auto key=id+"\t"+ch;auto p=parent.find(key),c=candidate.find(key);
   if(p==parent.end()&&c==candidate.end())continue;
   require(p!=parent.end()&&c!=candidate.end(),"channel transition "+key);
   auto a=p->second,b=c->second;
   if(a.s!=b.s||a.q!=b.q) {
    ++changed;require(ch!="merged","merged changed");
    // Compare retained raw coordinates; downstream repair may restore a prefix.
    auto pa=a.h.rfind(" FQT=v1;"),pb=b.h.rfind(" FQT=v1;");require(pa!=string::npos&&pb!=string::npos,"FQT required");
    auto fa=fields(a.h.substr(pa+1)),fb=fields(b.h.substr(pb+1));
    require(fa.at("m")==fb.at("m"),"changed mate");
    size_t start=stoull(fa.at("s")),end=stoull(fa.at("e")),newstart=stoull(fb.at("s")),newend=stoull(fb.at("e"));require(end>=start&&newend>=newstart&&end-start==a.s.size()&&newend-newstart==b.s.size(),"span lengths");int mate=stoi(fa.at("m"));
    const auto& orig=j.at("reads").a.at(mate-1).at("origins").a;
    long long ad=0,mo=0,inside=0;
    for(size_t i=start;i<end;++i){if(i>=newstart&&i<newend)continue;auto& o=orig.at(i);if(o.at("kind").s=="adapter")++ad;else{++mo;auto coord=stoll(o.at("coordinate").s);inside+=coord>=0&&coord<frag;}}
    adapter+=ad;molecule+=mo;interior+=inside;failed+=inside>0;
    out<<id<<'\t'<<ch<<'\t'<<(a.s!=b.s)<<'\t'<<(a.q!=b.q)<<'\t'<<(a.s.size()!=b.s.size())<<'\t'<<ad<<'\t'<<mo<<'\t'<<inside<<'\n';
   }
   parent.erase(p);candidate.erase(c);
  }
 }
 require(parent.empty()&&candidate.empty(),"output IDs missing from truth");
 cout<<"changed="<<changed<<" adapter_removed="<<adapter<<" molecule_removed="<<molecule<<" interior_molecule_removed="<<interior<<" verdict="<<(failed?"FAIL":"PASS")<<'\n';
 return failed?1:0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 2;}}
