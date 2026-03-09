// fqdup - FASTQ utilities for paired-end deduplication
// Subcommands:
//   fqdup sort        - External merge sort by read ID
//   fqdup extend      - De Bruijn graph extension of merged reads (replaces Tadpole)
//   fqdup derep_pairs - Two-pass paired deduplication (no damage/error modeling)
//   fqdup derep       - Single-file deduplication with damage-aware hashing + error correction

#include <iostream>
#include <string>

int sort_main(int argc, char** argv);
int extend_main(int argc, char** argv);
int derep_pairs_main(int argc, char** argv);
int derep_main(int argc, char** argv);
int damage_main(int argc, char** argv);

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <subcommand> [options]\n"
              << "\nSubcommands:\n"
              << "  sort        Sort FASTQ by read ID (external merge sort)\n"
              << "  extend      Extend merged reads via self-referential k-mer graph\n"
              << "  derep_pairs Deduplicate sorted paired-end FASTQ (structural dedup)\n"
              << "  derep       Deduplicate sorted single-file FASTQ (damage + error correction)\n"
              << "  damage      Profile deamination damage and compute position mask\n"
              << "\nRun '" << prog << " <subcommand> --help' for subcommand options.\n"
              << "\nTypical workflow:\n"
              << "  " << prog << " sort        -i merged.fq.gz  -o merged.sorted.fq.gz  --max-memory 64G\n"
              << "  " << prog << " extend      -i merged.sorted.fq.gz -o extended.fq.gz\n"
              << "  " << prog << " sort        -i extended.fq.gz -o extended.sorted.fq.gz --max-memory 64G\n"
              << "  " << prog << " derep_pairs -n merged.sorted.fq.gz -e extended.sorted.fq.gz \\\n"
              << "                            -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz\n"
              << "  " << prog << " derep       -i nonext.deduped.fq.gz -o nonext.final.fq.gz \\\n"
              << "                            --damage-auto --error-correct\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    std::string sub = argv[1];

    if (sub == "sort")        return sort_main(argc - 1, argv + 1);
    if (sub == "extend")      return extend_main(argc - 1, argv + 1);
    if (sub == "derep_pairs") return derep_pairs_main(argc - 1, argv + 1);
    if (sub == "derep")       return derep_main(argc - 1, argv + 1);
    if (sub == "damage")      return damage_main(argc - 1, argv + 1);

    if (sub == "-h" || sub == "--help") {
        usage(argv[0]);
        return 0;
    }

    std::cerr << "Error: unknown subcommand '" << sub << "'\n\n";
    usage(argv[0]);
    return 1;
}
