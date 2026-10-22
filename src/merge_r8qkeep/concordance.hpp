#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
namespace fqdup::r8qkeep {
inline constexpr std::uint64_t calibration_min_columns=10000;
inline constexpr double calibration_ratio_low=.5,calibration_ratio_high=2.;
inline constexpr double correction_ratio_low=.67,correction_ratio_high=1.5;
inline constexpr double correction_min_expected=25.;
struct ConcordanceCell {
    std::uint64_t n=0,disagree=0;double expected=0;
    void add(const ConcordanceCell& b){n+=b.n;disagree+=b.disagree;expected+=b.expected;}
    double ratio()const{return expected>0?disagree/expected:std::numeric_limits<double>::infinity();}
};
struct ConcordanceBins {
    std::array<ConcordanceCell,9*9*4> cells{};
    void add(const ConcordanceBins& b){for(std::size_t i=0;i<cells.size();++i)cells[i].add(b.cells[i]);}
    ConcordanceCell total()const{ConcordanceCell all;for(const auto& c:cells)all.add(c);return all;}
    bool cutoff_ok()const{const auto all=total();return all.n>=calibration_min_columns&&all.ratio()>=calibration_ratio_low&&all.ratio()<=calibration_ratio_high;}
    bool correction_ok()const{
        const auto all=total();if(all.n<calibration_min_columns||all.ratio()<correction_ratio_low||all.ratio()>correction_ratio_high)return false;
        for(const auto& c:cells)if(c.expected>=correction_min_expected&&(c.ratio()<calibration_ratio_low||c.ratio()>calibration_ratio_high))return false;
        return true;
    }
    std::string tsv()const{
        std::ostringstream out;out.precision(17);out<<"q1_bin\tq2_bin\tbase\tn\tdisagree\texpected\tratio\tstatus\n";
        for(int q1=0;q1<9;++q1)for(int q2=0;q2<9;++q2)for(int b=0;b<4;++b){const auto& c=cells[(q1*9+q2)*4+b];if(!c.n)continue;
            out<<q1<<'\t'<<q2<<'\t'<<"ACGT"[b]<<'\t'<<c.n<<'\t'<<c.disagree<<'\t'<<c.expected<<'\t'<<c.ratio()<<'\t'<<(c.expected<correction_min_expected?"SPARSE":c.ratio()>=.5&&c.ratio()<=2?"PASS":"FAIL")<<'\n';}
        const auto all=total();out<<"ALL\tALL\tALL\t"<<all.n<<'\t'<<all.disagree<<'\t'<<all.expected<<'\t'<<all.ratio()<<'\t'<<(cutoff_ok()?"ok":"bad")<<'\n';return out.str();
    }
};
struct RawConcordance {ConcordanceBins bins;std::string nominal_tsv;};
}
