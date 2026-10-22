// T41-1: exhaustively check changed records against raw truth origins.
#define main task47_attribution_main
#include "task47_attribution.cpp"
#undef main
struct Fastq { string h,s,q,ch; };
map<string,Fastq> records(const string& dir) {
 map<string,Fastq> result;
 for(string ch:{"merged","r1","r2","orphan"})
  readfq(dir+"/"+ch+".fq",[&](const string& h,const string& s,const string& q){
   string mate="merged";
   if(ch!="merged") {
    auto p=h.rfind(" FQT=v1;");require(p!=string::npos,"mate provenance");
    auto f=fields(h.substr(p+1,h.find(' ',p+1)-p-1));mate=f.at("m");
    require(mate=="1"||mate=="2","invalid mate");
    require(ch=="orphan"||ch=="r"+mate,"channel/mate mismatch");
   }
   require(result.emplace(idof(h)+"\t"+mate,Fastq{h,s,q,ch}).second,"duplicate record");
  });
 return result;
}
pair<size_t,size_t> span(const Fastq& r) {
 auto p=r.h.rfind(" FQT=v1;");require(p!=string::npos,"FQT required");
 auto f=fields(r.h.substr(p+1,r.h.find(' ',p+1)-p-1));
 size_t s=stoull(f.at("s")),e=stoull(f.at("e"));
 require(e>=s&&e-s==r.s.size(),"span lengths");return {s,e};
}
#ifndef TASK47_SPAN_LIBRARY
int main(int argc,char**argv) {try {
 require(argc==5||argc==6,"usage: clip_check truth parent-dir candidate-dir report [--channel-changes]");
 bool channel_changes=argc==6;require(!channel_changes||string(argv[5])=="--channel-changes","invalid mode");
 auto parent=records(argv[2]),candidate=records(argv[3]);
 ofstream out(argv[4]); require(bool(out),"report open");
 out<<"id\tchannel\tsequence_changed\tquality_changed\tlength_changed\tadapter_removed\tmolecule_removed\tinterior_molecule_removed\tnew_channel\tmate\n";
 long long changed=0,adapter=0,molecule=0,interior=0,failed=0;
 ifstream in(argv[1]);require(bool(in),"truth open");string line;
 while(getline(in,line)) {
  Parser parser{line};auto j=parser.value();string id=j.at("id").s;auto frag=stoll(j.at("frag_len").s);
  if(channel_changes&&(parent.count(id+"\tmerged")||candidate.count(id+"\tmerged"))){for(string mate:{"merged","1","2"}){parent.erase(id+"\t"+mate);candidate.erase(id+"\t"+mate);}continue;}
  for(string mate:{"merged","1","2"}) {
   auto key=id+"\t"+mate;auto p=parent.find(key),c=candidate.find(key);
   if(p==parent.end()&&c==candidate.end())continue;
   if(channel_changes&&p==parent.end()){candidate.erase(c);continue;}
   require(p!=parent.end(),"unexpected addition "+key);
   const auto& a=p->second;
   if(mate=="merged") {
    require(c!=candidate.end()&&a.h==c->second.h&&a.s==c->second.s&&a.q==c->second.q,"merged changed "+key);
   } else {
    auto oldspan=span(a);pair<size_t,size_t> newspan{0,0};
    bool missing=c==candidate.end();if(!missing)newspan=span(c->second);
    const auto& orig=j.at("reads").a.at(stoi(mate)-1).at("origins").a;
    require(oldspan.second<=orig.size()&&newspan.second<=orig.size(),"span outside raw read");
    if(missing||a.s!=c->second.s||a.q!=c->second.q||a.ch!=c->second.ch||oldspan!=newspan) {
     ++changed;long long ad=0,mo=0,inside=0;
     for(size_t i=oldspan.first;i<oldspan.second;++i){
      if(!missing&&i>=newspan.first&&i<newspan.second)continue;
      const auto& o=orig.at(i);
      if(o.at("kind").s=="adapter")++ad;
      else {require(o.at("kind").s=="molecule","unknown origin");++mo;auto coord=stoll(o.at("coordinate").s);inside+=coord>=0&&coord<frag;}
     }
     adapter+=ad;molecule+=mo;interior+=inside;failed+=inside>0;
     out<<id<<'\t'<<a.ch<<'\t'<<(missing||a.s!=c->second.s)<<'\t'<<(missing||a.q!=c->second.q)<<'\t'<<(missing||a.s.size()!=c->second.s.size())<<'\t'<<ad<<'\t'<<mo<<'\t'<<inside<<'\t'<<(missing?"dropped":c->second.ch)<<'\t'<<mate<<'\n';
    }
   }
   parent.erase(p);if(c!=candidate.end())candidate.erase(c);
  }
 }
 require(parent.empty()&&candidate.empty(),"output IDs missing from truth");
 cout<<"changed="<<changed<<" adapter_removed="<<adapter<<" molecule_removed="<<molecule<<" interior_molecule_removed="<<interior<<" verdict="<<(failed?"FAIL":"PASS")<<'\n';
 return failed?1:0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 2;}}

#endif
