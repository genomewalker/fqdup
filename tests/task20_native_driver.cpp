#include "../src/merge_r8qkeep/native.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        if(argc!=4)throw std::runtime_error("usage: native_driver RAW_R1 RAW_R2 NEW_OUTPUT_DIR");
        if(!std::filesystem::create_directory(argv[3]))throw std::runtime_error("output directory exists");
        auto family=fqdup::r8qkeep::correct_family(fqdup::r8qkeep::read_family_fastq(argv[1]),fqdup::r8qkeep::read_family_fastq(argv[2]));
        auto result=fqdup::r8qkeep::merge_native(family.pairs);
        const std::string out=argv[3];
        const auto write=[&](const std::string& file,const std::string& text){std::ofstream f(out+"/"+file);f<<text;if(!f)throw std::runtime_error("write failed: "+file);};
        write("merged.fq",result.merged_fq);write("r1.fq",result.r1_fq);write("r2.fq",result.r2_fq);write("orphan.fq",result.orphan_fq);
        write("outcomes.jsonl",result.outcomes_jsonl);write("native.tsv",result.native_tsv);write("native.outcomes.json",result.native_outcomes_json);
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
