#include "../src/merge_r8qkeep/family.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: family_driver R1 R2 NEW_OUTPUT_DIR");
        if (!std::filesystem::create_directory(argv[3])) throw std::runtime_error("output directory exists");
        auto result = fqdup::r8qkeep::correct_family(fqdup::r8qkeep::read_family_fastq(argv[1]), fqdup::r8qkeep::read_family_fastq(argv[2]));
        const std::string out = argv[3];
        for (int mate = 0; mate < 2; ++mate) {
            std::ofstream f(out + (mate ? "/R2.fastq" : "/R1.fastq"));
            for (const auto& pair : result.pairs) {
                const auto& read = pair.r[mate];
                f << read.h << '\n' << read.s << "\n+\n" << read.q << '\n';
            }
            if (!f) throw std::runtime_error("FASTQ write failed");
        }
        std::ofstream(out + "/family-edits.tsv") << result.edits_tsv;
        std::ofstream(out + "/families.tsv") << result.families_tsv;
        std::ofstream(out + "/stats.json") << result.stats_json;
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
