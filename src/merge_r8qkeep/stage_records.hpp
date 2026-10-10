#pragma once
#include "family.hpp"
#include <algorithm>
#include <iterator>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
struct StageRead {std::string header,seq,qual;};
struct NativeRecord {
    const FamilyPair* raw=nullptr;
    std::string outcome;
    StageRead first,second;
    std::size_t size()const{return 13;}
    const std::string& operator[](std::size_t i)const{
        switch(i){
        case 0:return raw->r[0].h;case 1:return raw->r[0].s;case 2:return raw->r[0].q;
        case 3:return raw->r[1].h;case 4:return raw->r[1].s;case 5:return raw->r[1].q;
        case 6:return outcome;case 7:return first.header;case 8:return first.seq;case 9:return first.qual;
        case 10:return second.header;case 11:return second.seq;case 12:return second.qual;
        default:throw std::out_of_range("native field");
        }
    }
};
class NativeRecords {
 const NativeRecord* data_=nullptr;std::size_t size_=0;const std::size_t* indices_=nullptr;
 NativeRecords(const NativeRecord* data,std::size_t n,const std::size_t* indices):data_(data),size_(n),indices_(indices){}
public:
 struct Iterator {
  using iterator_category=std::forward_iterator_tag;using value_type=NativeRecord;using difference_type=std::ptrdiff_t;using pointer=const NativeRecord*;using reference=const NativeRecord&;
  const NativeRecord* data=nullptr;const std::size_t* indices=nullptr;std::size_t position=0;
  reference operator*()const{return data[indices?indices[position]:position];}pointer operator->()const{return &**this;}
  Iterator& operator++(){++position;return *this;}Iterator operator++(int){auto old=*this;++*this;return old;}
  bool operator==(const Iterator& b)const{return data==b.data&&indices==b.indices&&position==b.position;}bool operator!=(const Iterator& b)const{return !(*this==b);}
 };
 NativeRecords()=default;
 NativeRecords(const std::vector<NativeRecord>& rows):data_(rows.data()),size_(rows.size()){}
 NativeRecords(const NativeRecord* data,std::size_t n):data_(data),size_(n){}
 std::size_t size()const{return size_;}
 const NativeRecord& operator[](std::size_t i)const{return data_[indices_?indices_[i]:i];}
 Iterator begin()const{return {data_,indices_,0};}Iterator end()const{return {data_,indices_,size_};}
 NativeRecords slice(std::size_t first,std::size_t n)const{
  if(first>size_||n>size_-first)throw std::out_of_range("native slice");
  if(indices_)return {data_,n,first?indices_+first:indices_};
  return {first?data_+first:data_,n};
 }
 // Both record storage and index storage must outlive the borrowed view.
 NativeRecords select(const std::vector<std::size_t>& indices)const{
  if(indices_)throw std::logic_error("nested native selection");
  for(auto i:indices)if(i>=size_)throw std::out_of_range("native selection");
  return {data_,indices.size(),indices.data()};
 }
};
// Full-span fast results borrow native sequence/quality instead of allocating
// one set of strings per row/cutoff. Other results own the guard's moved output.
struct RepairedRecord {
    const NativeRecord* native=nullptr;
    std::string outcome="dropped";
    StageRead first,second;
};
inline std::string native_record_tsv(const NativeRecord& row){
    std::string out;for(std::size_t i=0;i<row.size();++i){if(i)out+='\t';out+=row[i];}out+='\n';return out;
}
// Legacy TSV adapters are used by regression tests and explicit diagnostics;
// the production native/census/router/repair path passes records directly.
struct ParsedNativeRecords {
    std::vector<FamilyPair> parents;
    std::vector<NativeRecord> records;
    ParsedNativeRecords()=default;
    ParsedNativeRecords(ParsedNativeRecords&&)=default;
    ParsedNativeRecords& operator=(ParsedNativeRecords&&)=default;
    ParsedNativeRecords(const ParsedNativeRecords&)=delete;
};
inline ParsedNativeRecords parse_native_records(const std::string& text){
    ParsedNativeRecords out;
    const auto n=std::count(text.begin(),text.end(),'\n')+(!text.empty()&&text.back()!='\n');
    out.parents.reserve(n);out.records.reserve(n);
    std::istringstream input(text);std::string line;
    while(std::getline(input,line)){
        std::vector<std::string> v;std::size_t begin=0;
        for(;;){const auto end=line.find('\t',begin);v.push_back(line.substr(begin,end==std::string::npos?end:end-begin));if(end==std::string::npos)break;begin=end+1;}
        if(v.size()!=13)throw std::runtime_error("native column count");
        FamilyPair p;p.r[0].h=std::move(v[0]);p.r[0].s=std::move(v[1]);p.r[0].q=std::move(v[2]);
        p.r[1].h=std::move(v[3]);p.r[1].s=std::move(v[4]);p.r[1].q=std::move(v[5]);out.parents.push_back(std::move(p));
        NativeRecord row;row.raw=&out.parents.back();row.outcome=std::move(v[6]);
        row.first={std::move(v[7]),std::move(v[8]),std::move(v[9])};row.second={std::move(v[10]),std::move(v[11]),std::move(v[12])};
        out.records.push_back(std::move(row));
    }
    return out;
}
}
