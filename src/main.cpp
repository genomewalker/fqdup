// fqdup - FASTQ utilities for paired-end deduplication
// Subcommands:
//   fqdup sort   - External merge sort by read ID
//   fqdup derep  - Memory-efficient two-pass deduplication (requires sorted input)

#include <iostream>
#include <string>

int sort_main(int argc, char** argv);
int derep_main(int argc, char** argv);

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <subcommand> [options]\n"
              << "\nSubcommands:\n"
              << "  sort   Sort FASTQ by read ID (external merge sort)\n"
              << "  derep  Deduplicate sorted paired-end FASTQ (two-pass, ~16B/read)\n"
              << "\nRun '" << prog << " <subcommand> --help' for subcommand options.\n"
              << "\nTypical workflow:\n"
              << "  " << prog << " sort   -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G\n"
              << "  " << prog << " sort   -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G\n"
              << "  " << prog << " derep  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \\\n"
              << "                       -o-non out.non.fq.gz -o-ext out.ext.fq.gz\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    std::string sub = argv[1];

    if (sub == "sort")  return sort_main(argc - 1, argv + 1);
    if (sub == "derep") return derep_main(argc - 1, argv + 1);

    if (sub == "-h" || sub == "--help") {
        usage(argv[0]);
        return 0;
    }

    std::cerr << "Error: unknown subcommand '" << sub << "'\n\n";
    usage(argv[0]);
    return 1;
}
