#include "../src/merge_r8qkeep/compact_store.hpp"
#include <iostream>
#include <random>
using namespace fqdup::r8qkeep;
static void check(bool ok){if(!ok)throw std::runtime_error("absolute edit parity");}
struct Edit {std::uint32_t pair,position;unsigned char mate;char base,quality;bool changed;};
int main(){WorkerPool pool(8);std::mt19937 rng(260108);
 for(bool canonical:{false,true})for(int workers:{1,8}){
  CompactPairStore scalar;CompactPairStore fast;std::vector<FamilyPair> expected(4101);
  for(unsigned p=0;p<expected.size();++p){auto& row=expected[p];const auto id=std::to_string(expected.size()-p);
   for(int m=0;m<2;++m){row.r[m].h="@"+id+"/"+std::to_string(m+1)+" comment";row.r[m].s.resize(1+rng()%512);row.r[m].q.resize(row.r[m].s.size());
    for(std::size_t i=0;i<row.r[m].s.size();++i){row.r[m].s[i]="ACGTN"[rng()%5];row.r[m].q[i]=char(rng()%256);}}
   scalar.append(row.r[0],row.r[1]);fast.append(row.r[0],row.r[1]);
  }
  if(canonical){scalar.canonicalize();fast.canonicalize();}
  CompactPairStore::Decoder cursor;for(std::size_t p=0;p<expected.size();++p)scalar.read(p,expected[p],cursor);
  for(int round=0;round<3;++round){std::vector<Edit> edits;std::vector<QualityBlocks::Edit> qualities,bases;
   for(unsigned i=0;i<20000;++i){const auto p=rng()%expected.size();const unsigned m=rng()%2;auto& r=expected[p].r[m];const auto pos=rng()%r.s.size();
    const char base=i%4?r.s[pos]:char(rng()%256),quality=char(rng()%256);const bool changed=base!=r.s[pos];r.s[pos]=base;r.q[pos]=quality;
    edits.push_back({std::uint32_t(p),std::uint32_t(pos),static_cast<unsigned char>(m),base,quality,changed});
    const auto at=fast.edit_address(p,m,pos);qualities.push_back({at,static_cast<unsigned char>(quality)});if(changed)bases.push_back({at,static_cast<unsigned char>(base)});
   }
   scalar.apply_edits(std::move(edits),&pool,workers,[](const Edit& e){return e.changed;});fast.apply_absolute_edits(std::move(qualities),std::move(bases),&pool,workers);
   CompactPairStore::Decoder a;CompactPairStore::Decoder b;FamilyPair x,y;
   for(std::size_t p=0;p<expected.size();++p){scalar.read(p,x,a);fast.read(p,y,b);check(x.id==y.id);for(int m=0;m<2;++m)check(x.r[m].h==y.r[m].h&&x.r[m].s==y.r[m].s&&x.r[m].q==y.r[m].q&&x.r[m].s==expected[p].r[m].s&&x.r[m].q==expected[p].r[m].q);}
  }
  fast.apply_absolute_edits({}, {},&pool,workers);
  for(int bad=0;bad<3;++bad){bool caught=false;try{fast.edit_address(bad==0?fast.size():0,bad==1?2:0,bad==2?100000:0);}catch(const std::out_of_range&){caught=true;}check(caught);}
  for(bool base:{false,true}){bool caught=false;try{std::vector<QualityBlocks::Edit> bad{{fast.symbols(),0}};if(base)fast.apply_absolute_edits({},bad,&pool,workers);else fast.apply_absolute_edits(bad,{},&pool,workers);}catch(const std::out_of_range&){caught=true;}check(caught);}
 }
 std::cout<<"ABSOLUTE_EDIT_PASS scalar equality canonical order duplicate edits arbitrary bytes unchanged bases p1/p8 bounds empty\n";
}
