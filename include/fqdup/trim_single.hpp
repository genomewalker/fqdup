#pragma once
#include "fqdup/fastq_types.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
namespace fqdup::raw_se {
struct Options {
 std::string input,output,report,adapter5,adapter3;
 std::size_t min_length=0,max_length=0,clip5=0,clip3=0,overlap=0;
 unsigned threads=1,min_phred=0;bool quality_ends=false,help=false;
};
struct Counts {
 std::uint64_t input=0,retained=0,empty=0,short_reads=0,long_reads=0,input_bases=0,retained_bases=0,discarded_bases=0;
 std::uint64_t hard5_reads=0,hard3_reads=0,adapter5_reads=0,adapter3_reads=0,quality_reads=0;
 std::uint64_t hard5_bases=0,hard3_bases=0,adapter5_bases=0,adapter3_bases=0,quality_bases=0;
 void add(const Counts&);void check()const;
};
struct Result {FastqRecord record;Counts counts;std::size_t begin=0,end=0;bool keep=false;};
Options parse_options(int argc,const char* const* argv);
void validate_options(const Options&);
void validate_record(const FastqRecord&);
// Options must already satisfy validate_options. Moves the raw record; new
// sequence/quality values are aligned contiguous slices no longer than input.
Result transform(FastqRecord record,const Options&);
// Reuses the qualified FASTQ factory/writer. Regular immutable input required
// for a fixed-buffer preflight NUL scan (not a second FASTQ parser).
Counts run(const Options&);
std::string report(const Options&,const Counts&);
}
int trim_single_main(int argc,char** argv);
