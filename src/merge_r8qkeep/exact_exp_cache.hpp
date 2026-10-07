#pragma once
#include "fdlibm.hpp"
#include <array>
#include <cstdint>
#include <cstring>
namespace fqdup::r8qkeep {
// Bounded per-fit cache: exact argument bits, exact fdlibm result, collision eviction.
class ExactExpCache {
    struct Entry { std::uint64_t bits=0; double value=0; bool valid=false; };
    std::array<Entry,16384> entries_{};
    std::uint64_t calls_=0,hits_=0;
public:
    std::uint64_t calls()const{return calls_;}
    std::uint64_t hits()const{return hits_;}
    double operator()(double x) {
        ++calls_;
        std::uint64_t bits;std::memcpy(&bits,&x,sizeof(bits));
        std::uint64_t hash=bits;hash^=hash>>30;hash*=0xbf58476d1ce4e5b9ULL;
        hash^=hash>>27;hash*=0x94d049bb133111ebULL;hash^=hash>>31;
        auto& entry=entries_[hash&(entries_.size()-1)];
        if(entry.valid&&entry.bits==bits){++hits_;return entry.value;}
        const double result=fqdup::fdlibm::exp(x);
        entry={bits,result,true};return result;
    }
};
}
