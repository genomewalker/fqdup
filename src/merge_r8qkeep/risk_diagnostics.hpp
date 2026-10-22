#pragma once
#include "holdout_bins.hpp"
#include <algorithm>
#include <map>
#include <optional>
#include <string_view>
#include <stdexcept>
#include <vector>
namespace fqdup::r8qkeep {
struct OutputRisk {
    std::array<std::array<std::uint64_t,94>,2> qualities{};
    std::array<std::uint64_t,2> records{};
    std::map<std::size_t,std::uint64_t> merged_lengths;
    std::uint64_t adapter_bases=0,adapter_mates=0;
    void scan(std::string_view text,int channel){
        std::size_t line=0,begin=0;
        while(begin<text.size()){
            const auto end=text.find('\n',begin);if(end==std::string_view::npos)throw std::runtime_error("risk FASTQ line");
            if(line%4==3){++records[channel];if(channel==0)++merged_lengths[end-begin];
                for(std::size_t i=begin;i<end;++i){const int q=int(static_cast<unsigned char>(text[i]))-33;if(q<0||q>=94)throw std::runtime_error("risk quality");++qualities[channel][q];}}
            ++line;begin=end+1;
        }
        if(line%4)throw std::runtime_error("risk FASTQ cardinality");
    }
};
inline std::string risk_report(const std::vector<OutputRisk>& risks,const std::vector<std::string>& labels,const HoldoutBins& holdout){
    const auto bins=holdout.qualities();
    const std::optional<double> copy_floor=bins[60].n>=calibration_min_bin?std::optional<double>(double(bins[60].disagree)/bins[60].n):std::nullopt;
    std::ostringstream out;out.precision(17);
    out<<"cutoff\tchannel\trecords\tbases\traw_risk\tcalibrated_risk\tuncalibrated_bases\tcopy_error_floor\tadapter_residual_bases\tadapter_residual_mates\n";
    for(std::size_t c=0;c<risks.size();++c)for(int channel=0;channel<2;++channel){
        std::uint64_t bases=0,uncalibrated=0;double raw=0,calibrated=0;
        for(int q=0;q<94;++q){const auto n=risks[c].qualities[channel][q];bases+=n;const double phred=std::pow(10.,-q/10.);raw+=n*phred;
            if(channel==1){const bool supported=bins[q].n>=calibration_min_bin;calibrated+=n*(supported?double(bins[q].disagree)/bins[q].n:phred);if(!supported)uncalibrated+=n;}
            else if(copy_floor)calibrated+=n*std::max(phred,*copy_floor);
        }
        out<<labels.at(c)<<'\t'<<(channel==0?"merged":"mates")<<'\t'<<risks[c].records[channel]<<'\t'<<bases<<'\t'<<(bases?raw/bases:0)<<'\t';
        if(channel==1||copy_floor)out<<(bases?calibrated/bases:0);else out<<"NA";
        out<<'\t'<<uncalibrated<<'\t';if(copy_floor)out<<*copy_floor;else out<<"NA";out<<'\t'<<(channel==1?risks[c].adapter_bases:0)<<'\t'<<(channel==1?risks[c].adapter_mates:0)<<'\n';
    }
    return out.str();
}
}
