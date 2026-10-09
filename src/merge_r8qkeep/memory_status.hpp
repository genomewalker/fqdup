#pragma once
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>
namespace fqdup::r8qkeep {
inline bool memory_enabled(){return std::getenv("FQDUP_MEMORY")!=nullptr;}
// Count owned capacity, including the terminator, but not inline string storage
// (already included in the containing object). No allocator-specific SSO size.
inline std::size_t string_heap(const std::string& s){
    const auto data=reinterpret_cast<std::uintptr_t>(s.data());
    const auto object=reinterpret_cast<std::uintptr_t>(&s);
    return data>=object&&data<object+sizeof(s)?0:s.capacity()+1;
}
template<class T> inline std::size_t vector_bytes(const std::vector<T>& v){return v.capacity()*sizeof(T);}
inline std::size_t vector_bytes(const std::vector<bool>& v){return (v.capacity()+7)/8;}
inline std::size_t strings_bytes(const std::vector<std::string>& v){
    std::size_t bytes=vector_bytes(v);for(const auto& s:v)bytes+=string_heap(s);return bytes;
}
inline std::string memory_status(){
    std::ifstream input("/proc/self/status");std::string line;long rss=-1,hwm=-1;
    while(std::getline(input,line)){std::istringstream fields(line);std::string key;long value;
        if(fields>>key>>value){if(key=="VmRSS:")rss=value;else if(key=="VmHWM:")hwm=value;}}
    return " rss_kib="+std::to_string(rss)+" hwm_kib="+std::to_string(hwm);
}
struct MemoryUsage {
    std::size_t raw=0,ids=0,family=0,native=0,census=0,routing=0,output=0,diagnostics=0;
    template<class Read> void read(const Read& r){
        raw+=sizeof(Read)-sizeof(r.id)+string_heap(r.h)+string_heap(r.s)+string_heap(r.q);
        ids+=sizeof(r.id)+string_heap(r.id);
    }
    template<class Read> void reads(const std::vector<Read>& rows){
        raw+=(rows.capacity()-rows.size())*(sizeof(Read)-sizeof(rows[0].id));
        ids+=(rows.capacity()-rows.size())*sizeof(rows[0].id);
        for(const auto& r:rows)read(r);
    }
    template<class Pair> void pairs(const std::vector<Pair>& rows){
        raw+=(rows.capacity()-rows.size())*(sizeof(Pair)-3*sizeof(std::string));
        ids+=(rows.capacity()-rows.size())*3*sizeof(std::string);
        for(const auto& p:rows){ids+=sizeof(p.id)+string_heap(p.id);for(const auto& r:p.r)read(r);}
    }
    void report(const char* stage,std::size_t pairs,const char* scope="owned_capacity")const{
        std::cerr<<"MERGE_MEMORY stage="<<stage<<" scope="<<scope<<" pairs="<<pairs
            <<" raw_bytes="<<raw<<" id_bytes="<<ids<<" family_bytes="<<family
            <<" native_bytes="<<native<<" census_bytes="<<census<<" routing_bytes="<<routing
            <<" output_bytes="<<output<<" diagnostic_bytes="<<diagnostics
            <<" total_bytes="<<raw+ids+family+native+census+routing+output+diagnostics
            <<memory_status()<<'\n';
    }
};
}
