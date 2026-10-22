#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
namespace fqdup::r8qkeep {
// T47-B3: frozen bins and support bound, declared before measurement.
inline constexpr std::size_t calibration_min_bin=1000;
struct HoldoutCell {std::uint64_t n=0,disagree=0;};
struct HoldoutBins {
    std::array<HoldoutCell,2*3*4*94> cells{};
    static std::size_t index(int mate,int window,int base,int q){return ((mate*3+window)*4+base)*94+q;}
    void add(const HoldoutBins& other){for(std::size_t i=0;i<cells.size();++i){cells[i].n+=other.cells[i].n;cells[i].disagree+=other.cells[i].disagree;}}
    std::array<HoldoutCell,94> qualities()const{
        std::array<HoldoutCell,94> bins{};
        for(std::size_t i=0;i<cells.size();++i){bins[i%94].n+=cells[i].n;bins[i%94].disagree+=cells[i].disagree;}
        return bins;
    }
    std::string tsv()const{
        std::ostringstream out;out.precision(17);out<<"mate\twindow\tbase\tQ\tn\tdisagree\tphred_expected\n";
        const char* windows[]={"first10","mid","last10"};
        for(int m=0;m<2;++m)for(int w=0;w<3;++w)for(int b=0;b<4;++b)for(int q=0;q<94;++q){
            const auto& c=cells[index(m,w,b,q)];if(c.n)out<<m+1<<'\t'<<windows[w]<<'\t'<<"ACGT"[b]<<'\t'<<q<<'\t'<<c.n<<'\t'<<c.disagree<<'\t'<<c.n*std::pow(10.,-q/10.)<<'\n';
        }
        return out.str();
    }
};
}
