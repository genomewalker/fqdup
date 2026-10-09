#pragma once
#include "packed_kernels.hpp"
#include "stage_records.hpp"
#include <functional>
#include <vector>
namespace fqdup::r8qkeep {
// Cache keys exclude IDs. Length framing and exact exception bytes make this
// encoding injective for every sequence/quality byte, including non-ACGT.
inline void memo_content(const FamilyPair& pair,std::string& key){
    key.clear();
    const auto integer=[&](std::uint64_t n){for(unsigned i=0;i<8;++i)key+=char(n>>(8*i));};
    for(const auto& r:pair.r){
        integer(r.s.size());const auto begin=key.size();key.resize(begin+(r.s.size()+3)/4,0);
        for(std::size_t i=0;i<r.s.size();++i)key[begin+i/4]|=char((packed::code_table[static_cast<unsigned char>(r.s[i])]&3)<<(2*(i%4)));
        std::size_t n=0;for(unsigned char c:r.s)n+=packed::code_table[c]==4;
        integer(n);for(std::size_t i=0;i<r.s.size();++i)if(packed::code_table[static_cast<unsigned char>(r.s[i])]==4){integer(i);key+=r.s[i];}
        integer(r.q.size());key+=r.q;
    }
}
struct MemoValue {
    RepairedRecord record;
    bool native_clean=false,clean_pass=false;
    std::size_t bytes()const{return sizeof(*this)+record.outcome.capacity()+record.first.header.capacity()+record.first.seq.capacity()+record.first.qual.capacity()+record.second.header.capacity()+record.second.seq.capacity()+record.second.qual.capacity();}
};
// Results are generated with these distinct placeholders, so header provenance
// is unambiguous even when the caller's two headers happen to be identical.
inline void memo_headers(FamilyPair& pair){pair.id="task26";pair.r[0].id=pair.r[1].id=pair.id;pair.r[0].h="@task26/1";pair.r[1].h="@task26/2";}
inline void memo_restore_header(std::string& header,const FamilyPair& pair){
    if(header.empty())return;
    std::string prefix;std::size_t length=0;
    if(header.compare(0,9,"@task26/1")==0){prefix=pair.r[0].h;length=9;}
    else if(header.compare(0,9,"@task26/2")==0){prefix=pair.r[1].h;length=9;}
    else if(header.compare(0,7,"@task26")==0){prefix=pair.r[0].h;length=7;
        if(prefix.size()>=2&&prefix[prefix.size()-2]=='/'&&(prefix.back()=='1'||prefix.back()=='2'))prefix.resize(prefix.size()-2);}
    else throw std::runtime_error("memo header provenance");
    header.replace(0,length,prefix);
}
inline MemoValue memo_owned(const RepairedRecord& input,bool native_clean,bool clean_pass){
    MemoValue value;value.native_clean=native_clean;value.clean_pass=clean_pass;
    if(!input.native){value.record=input;return value;}
    const auto& n=*input.native;auto& r=value.record;r.outcome=n.outcome;r.first=n.first;r.second=n.second;
    const auto tag=[](StageRead& read,const char* channel,int mate){if(read.seq.empty())return;
        const auto length=std::to_string(read.seq.size());read.header+=" FQT=v1;ch="+std::string(channel)+";m="+std::to_string(mate)+";olen="+length+";s=0;e="+length;};
    if(r.outcome=="merged")tag(r.first,"merged",0);
    else if(r.outcome=="unmerged"){tag(r.first,"r1",1);tag(r.second,"r2",2);}
    else if(r.outcome=="orphan_r1"||r.outcome=="orphan_r2")tag(r.first,"orphan",r.outcome=="orphan_r1"?1:2);
    return value;
}
class ContentMemo {
    struct Entry {std::size_t hash=0;std::string key;MemoValue value;};
    static constexpr std::size_t slots=2048,max_entry=8192;
    std::vector<Entry> entries_{slots};
    std::uint64_t context_=0;
public:
    void context(std::uint64_t id){if(context_!=id){entries_.clear();entries_.resize(slots);context_=id;}}
    void discard(){entries_.clear();entries_.resize(slots);}
    static std::size_t hash(std::string_view key){return std::hash<std::string_view>{}(key);}
    const MemoValue* find(std::string_view key,std::size_t h)const{
        const auto& e=entries_[h%slots];return e.hash==h&&e.key==key?&e.value:nullptr;
    }
    void put(const std::string& key,std::size_t h,const MemoValue& value){
        if(key.capacity()+value.bytes()>max_entry)return;
        auto& e=entries_[h%slots];e.hash=h;e.key=key;e.value=value;
        if(e.key.capacity()+e.value.bytes()>2*max_entry)e=Entry{};
    }
    static constexpr std::size_t capacity_bound(){return slots*(sizeof(Entry)+2*max_entry);}
};
}
