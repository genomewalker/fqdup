#include "../src/merge_r8qkeep/stage_records.hpp"
#include <random>
#include <iostream>
using namespace fqdup::r8qkeep;
int main(){
 std::mt19937 rng(261044);std::vector<FamilyPair> parents(2048);std::vector<NativeRecord> rows(2048);
 for(std::size_t i=0;i<rows.size();++i){parents[i].id=std::to_string(i);rows[i].raw=&parents[i];rows[i].outcome=i%2?"merged":"dropped";rows[i].first.seq.assign(i%513,'A');}
 NativeRecords empty;if(empty.begin()!=empty.end()||empty.slice(0,0).size())throw std::runtime_error("empty view");
 auto transient=NativeRecords(rows).begin();if(&*transient!=&rows[0])throw std::runtime_error("iterator ownership");
 for(int trial=0;trial<2000;++trial){const std::size_t offset=rng()%2049,n=rng()%(2049-offset);auto view=NativeRecords(rows).slice(offset,n);std::vector<std::size_t> selected;for(std::size_t i=0;i<n;++i)if(rng()%3==0)selected.push_back(i);if(trial%2)std::shuffle(selected.begin(),selected.end(),rng);
  auto indexed=view.select(selected);if(indexed.size()!=selected.size()||std::size_t(std::distance(indexed.begin(),indexed.end()))!=selected.size())throw std::runtime_error("indexed count");
  std::size_t p=0;for(const auto&row:indexed){if(&row!=&rows[offset+selected[p++]]||row.raw!=&parents[offset+selected[p-1]])throw std::runtime_error("indexed identity");}
  for(std::size_t start=0;start<=selected.size();++start){auto tail=indexed.slice(start,selected.size()-start);if(tail.size()&&&tail[0]!=&rows[offset+selected[start]])throw std::runtime_error("indexed slice");}
 }
 bool bad=false;try{NativeRecords(rows).select(std::vector<std::size_t>{rows.size()});}catch(const std::out_of_range&){bad=true;}if(!bad)throw std::runtime_error("index bounds");
 std::cout<<"INDEXED_NATIVE_RECORDS_PASS trials=2000\n";
}
