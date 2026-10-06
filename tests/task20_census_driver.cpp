#include "../src/merge_r8qkeep/census.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        if(argc!=5)throw std::runtime_error("usage: census_driver CORRECTED_R1 CORRECTED_R2 NATIVE.tsv NEW_OUTPUT_DIR");
        if(!std::filesystem::create_directory(argv[4]))throw std::runtime_error("output directory exists");
        const auto r1=fqdup::r8qkeep::read_family_fastq(argv[1]),r2=fqdup::r8qkeep::read_family_fastq(argv[2]);
        if(r1.size()!=r2.size())throw std::runtime_error("mate count");std::vector<fqdup::r8qkeep::FamilyPair> pairs;
        for(std::size_t i=0;i<r1.size();++i){if(r1[i].id!=r2[i].id)throw std::runtime_error("mate identity");pairs.push_back({r1[i].id,{r1[i],r2[i]}});}
        std::ifstream input(argv[3]);if(!input)throw std::runtime_error("native input open");fqdup::r8qkeep::NativeResult native;native.native_tsv.assign(std::istreambuf_iterator<char>(input),{});
        const auto result=fqdup::r8qkeep::census(pairs,native);const std::string out=argv[4];
        const auto write=[&](const std::string& file,const std::string& text){std::ofstream f(out+"/"+file);f<<text;if(!f)throw std::runtime_error("write failed");};
        write("nominal.tsv",result.nominal_tsv);write("full.tsv",result.full_tsv);write("proposed-edits.tsv",result.proposed_edits_tsv);return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
