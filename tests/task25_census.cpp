#include "../src/merge_r8qkeep/census.hpp"
#include "../src/merge_r8qkeep/ordered_chunks.hpp"
#include <cstring>
#include <iostream>
#include <stdexcept>
using namespace fqdup::r8qkeep;
static void require(bool value){if(!value)throw std::runtime_error("census range parity");}
static void compare_pairs(const std::vector<PairLikelihood>& a,const std::vector<PairLikelihood>& b){
    require(a.size()==b.size());
    for(std::size_t i=0;i<a.size();++i){
        require(a[i].id==b[i].id&&a[i].fold==b[i].fold&&a[i].agreements==b[i].agreements&&a[i].disagreements==b[i].disagreements&&a[i].sufficient_statistics==b[i].sufficient_statistics);
        require(std::memcmp(a[i].log_likelihood.data(),b[i].log_likelihood.data(),sizeof(a[i].log_likelihood))==0);
    }
}
static void compare(const CensusResult& a,const CensusResult& b){
    require(a.nominal_tsv==b.nominal_tsv&&a.full_tsv==b.full_tsv&&a.proposed_edits_tsv==b.proposed_edits_tsv);
    compare_pairs(a.nominal_pairs,b.nominal_pairs);compare_pairs(a.full_pairs,b.full_pairs);
}
int main(){
    // Native preparation requires fixed batch boundaries at every worker count.
    for(int workers:{1,3,8}){
        std::vector<std::pair<std::size_t,std::size_t>> ranges;
        ordered_index_chunks(4097,workers,11,[](std::size_t begin,std::size_t n){return std::pair{begin,n};},
            [&](auto range){ranges.push_back(range);},"batch_test",2048);
        require((ranges==std::vector<std::pair<std::size_t,std::size_t>>{{11,2048},{2059,2048},{4107,1}}));
    }
    std::vector<FamilyPair> pairs(1031); // More than one 1024-row queue window.
    const std::string seq="ACGTTGCAAGCTAGGCTAACGTAGCTAGTCGA";
    std::string rc;for(auto it=seq.rbegin();it!=seq.rend();++it)rc+="TGCA"[std::string("ACGT").find(*it)];
    for(std::size_t i=0;i<pairs.size();++i){
        auto& p=pairs[i];p.id="pair"+std::to_string(i);
        for(int mate=0;mate<2;++mate)p.r[mate]={"@"+p.id+"/"+std::to_string(mate+1),p.id,mate?rc:seq,std::string(seq.size(),char('I'+i%3))};
    }
    NativeResult native;
    for(std::size_t i=0;i<pairs.size();++i){
        const auto& p=pairs[i];native.records.push_back({&p,i%8?"merged":"unmerged",{p.r[0].h,p.r[0].s,p.r[0].q},{p.r[1].h,p.r[1].s,p.r[1].q}});
    }
    for(bool counts:{false,true}){
        const auto serial=census(pairs,native,true,1,137,counts);
        require(!serial.nominal_pairs.empty());
        for(int workers:{3,8})compare(serial,census(pairs,native,true,workers,137,counts));
        NativeResult text;for(const auto& row:native.records)text.native_tsv+=native_record_tsv(row);
        compare(serial,census(pairs,text,true,8,137,counts));
        auto empty=census({},NativeResult{},true,1,137,counts);
        compare(empty,census({},NativeResult{},true,8,137,counts));
        // Duplicate IDs crossing chunk boundaries retain serial accumulation.
        const auto saved=pairs.back().id;pairs.back().id=pairs.front().id;
        compare(census(pairs,native,true,1,137,counts),census(pairs,native,true,8,137,counts));pairs.back().id=saved;
    }
    native.records[700].raw=&pairs[699];
    bool failed=false;try{census(pairs,native,false,8);}catch(const std::runtime_error&){failed=true;}require(failed);
    native.records[700].raw=&pairs[700];native.records.pop_back();
    failed=false;try{census(pairs,native,false,8);}catch(const std::runtime_error&){failed=true;}require(failed);
    std::cout<<"TASK25_CENSUS_PASS\n";
}
