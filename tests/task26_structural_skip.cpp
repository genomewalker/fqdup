#include "../src/merge_r8qkeep/repair.cpp"
#include <random>
#include <iostream>
using namespace fqdup::r8qkeep;
static void require(bool b){if(!b)throw std::runtime_error("structural skip parity");}
static std::string prior(int component){std::ostringstream o;o.precision(17);o<<"FQDUP_PAIR_MIXTURE_V1 64\n";for(int f=0;f<2;++f){o<<f<<' '<<(component!=-1);for(int z=0;z<66;++z)o<<' '<<(component==-2?(z==0?.99:z==65?.01:0):(z==component?1:0));o<<'\n';}return o.str();}
static std::size_t checked=0,direct_total=0;static bool all_clean=false,no_clean=false,mixed=false;
static void check(RepairPipeline& pipeline,NativeRecords rows,const std::vector<double>& cuts,std::size_t start){
 const auto minimum=*std::min_element(cuts.begin(),cuts.end());const auto budget=int64_t(std::floor(minimum*1000000000000LL+.5));
 const auto bounds=pipeline.native_clean_bounds(rows,minimum);std::vector<std::size_t> unresolved;std::vector<int64_t> pending;
 for(std::size_t i=0;i<rows.size();++i)if(bounds[i]<0||bounds[i]>budget*int64_t(rows[i][8].size())){unresolved.push_back(i);pending.push_back(bounds[i]);}
 const auto work=rows.select(unresolved);RepairAlignmentCache reference_cache(start,rows.size()),cache(0,work.size());
 const auto reference_route=pipeline.route_records(rows,true,start,&reference_cache,true,minimum);
 const auto route=pipeline.route_records(work,false,0,&cache,false,minimum,&pending);
 all_clean|=rows.size()&&work.size()==0;no_clean|=rows.size()&&work.size()==rows.size();mixed|=work.size()>0&&work.size()<rows.size();
 direct_total+=rows.size()-work.size();
 for(double cut:cuts){RepairShadow shadow;
  const auto reference=pipeline.repair_records(rows,9,cut,start,&shadow,&reference_route.routed,&reference_route.repair_bounds,&reference_cache,&reference_route.clean_bounds);
  const auto reference_changed=pipeline.repair_records(rows,7,cut,start,nullptr,nullptr,nullptr,&reference_cache,nullptr,&reference_route.routed,&reference_route.record_replacements);
  const auto result=pipeline.repair_records(work,9,cut,0,nullptr,&route.routed,&route.repair_bounds,&cache,&route.clean_bounds);
  const auto changed=pipeline.repair_records(work,7,cut,0,nullptr,nullptr,nullptr,&cache,nullptr,&route.routed,&route.record_replacements);
  require(shadow.skip_but_changed==0&&shadow.clean_but_changed==0);std::size_t slot=0;RepairedRecord direct;
  for(std::size_t i=0;i<rows.size();++i){const RepairedRecord* actual;
   if(slot==unresolved.size()||unresolved[slot]!=i){require(!reference_route.routed[i]);direct.native=&rows[i];actual=&direct;}
   else{require(route.routed[slot]==reference_route.routed[i]);actual=route.routed[slot]?&changed[slot]:&result[slot];++slot;}
   const auto& expected=reference_route.routed[i]?reference_changed[i]:reference[i];
   require(repaired_record_tsv(*actual,start+i)==repaired_record_tsv(expected,start+i));++checked;
  }
 }
}
int main(){std::mt19937 rng(261045);
 for(int component:{-2,-1,0,1,32,65}){const auto p=prior(component);RepairPipeline pipeline(p,p);std::string input;
  for(int row=0;row<48;++row){std::string a(40,'A'),q(40,'I');for(auto&c:a)c="ACGT"[rng()%4];auto b=fqdup_r8qkeep_repair_internal::fqalign::reverse_complement(a);auto seq=a,qual=q;
   if(row%3==0)b[19]=b[19]=='A'?'C':'A';if(row%5==0)qual[0]='!';if(row%7==0)seq[20]='N';if(row%11==0){a[3]='N';b[5]='N';}
   const auto id="@pair"+std::to_string(row);const auto h1=id+"/1 comment FQT=user",h2=id+"/2 comment";
   input+=h1+'\t'+a+'\t'+q+'\t'+h2+'\t'+b+'\t'+q+"\tmerged\t"+h1+'\t'+seq+'\t'+qual+"\t\t\t\n";
  }
  input+="@drop/1\tA\tI\t@drop/2\tT\tI\tdropped\t\t\t\t\t\t\n";
  input+="@orphan/1\tACGT\tIIII\t@orphan/2\tACGT\tIIII\torphan_r1\t@orphan/1\tACGT\tIIII\t\t\t\n";
  input+="@unmerged/1\tACGT\tIIII\t@unmerged/2\tACGT\tIIII\tunmerged\t@unmerged/1\tACGT\tIIII\t@unmerged/2\tACGT\tIIII\n";
  auto parsed=parse_native_records(input);NativeRecords rows(parsed.records);
  for(double cut:{0.,.0001,.01,.316,1.})check(pipeline,rows,{cut},37);
  check(pipeline,rows,{.316,.01,1.,.1},413);check(pipeline,rows.slice(0,0),{.316},999);
  check(pipeline,rows.slice(48,1),{.316},85);check(pipeline,rows.slice(49,2),{.316},86);
  auto bad=parse_native_records("@bad/1\tAA\tI\t@bad/2\tT\tI\tdropped\t\t\t\t\t\t\n");require(pipeline.native_clean_bounds(bad.records,.316)[0]<0);bool failed=false;try{pipeline.repair_records(bad.records,9,.316);}catch(const std::runtime_error&e){failed=std::string(e.what())=="raw lengths";}require(failed);
 }
 require(all_clean&&no_clean&&mixed&&direct_total>0);std::cout<<"STRUCTURAL_SKIP_PASS compared="<<checked<<" classified_direct="<<direct_total<<" all_clean="<<all_clean<<" no_clean="<<no_clean<<" mixed="<<mixed<<'\n';
}
