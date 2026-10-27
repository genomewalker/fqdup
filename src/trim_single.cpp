#include "fqdup/trim_single.hpp"
#include "fqdup/fastq_common.hpp"
#include "merge_r8qkeep/worker_pool.hpp"
#include <algorithm>
#include <array>
#include <charconv>
#include <deque>
#include <filesystem>
#include <fstream>
#include <future>
#include <limits>
#include <set>
#include <sstream>
#include <string_view>
#include <sys/stat.h>
#ifndef RAW_SE_MUTANT
#define RAW_SE_MUTANT 0
#endif
namespace fqdup::raw_se {
namespace {
[[noreturn]] void fail(const std::string& why){throw std::runtime_error("single-end: "+why);}
std::uint64_t add(std::uint64_t a,std::uint64_t b){
#if RAW_SE_MUTANT != 18
 if(a>std::numeric_limits<std::uint64_t>::max()-b)fail("counter overflow");
#endif
 return a+b;
}
char upper(char c){return c>='a'&&c<='z'?char(c-'a'+'A'):c;}
bool canonical(char c){c=upper(c);return c=='A'||c=='C'||c=='G'||c=='T';}
bool base(char c){switch(upper(c)){case 'A':case 'C':case 'G':case 'T':case 'N':case 'R':case 'Y':case 'S':case 'W':case 'K':case 'M':case 'B':case 'D':case 'H':case 'V':return true;default:return false;}}
std::size_t number(std::string_view text){std::size_t value=0;auto p=std::from_chars(text.data(),text.data()+text.size(),value);
#if RAW_SE_MUTANT == 16
 if(p.ec!=std::errc{})fail("invalid integer");
#else
 if(text.empty()||p.ec!=std::errc{}||p.ptr!=text.data()+text.size())fail("invalid integer: "+std::string(text));
#endif
 return value;
}
bool match(char read,char adapter){
#if RAW_SE_MUTANT == 7
 if(upper(read)=='N')return true;
#endif
 return canonical(read)&&upper(read)==upper(adapter);
}
bool exact(std::string_view read,std::string_view adapter){if(read.size()!=adapter.size())return false;for(std::size_t i=0;i<read.size();++i)if(!match(read[i],adapter[i]))return false;return true;}
bool traced(const std::string& s){return s.find("FQT=")!=std::string::npos||s.find("FQC=")!=std::string::npos;}
void printable(const std::string& text){for(unsigned char c:text)if(c<32||c>126)fail("header/plus contains non-printable byte");}
bool same_file(const std::string& a,const std::string& b){std::error_code ec;
 if(std::filesystem::weakly_canonical(a)==std::filesystem::weakly_canonical(b))return true;
 const bool result=std::filesystem::equivalent(a,b,ec);return !ec&&result;
}
struct Stamp {struct stat value{};explicit Stamp(const std::string& path){if(::stat(path.c_str(),&value)!=0||!S_ISREG(value.st_mode))fail("regular input file required for NUL preflight");}
 bool same(const Stamp& other)const{return value.st_dev==other.value.st_dev&&value.st_ino==other.value.st_ino&&value.st_size==other.value.st_size&&value.st_mtim.tv_sec==other.value.st_mtim.tv_sec&&value.st_mtim.tv_nsec==other.value.st_mtim.tv_nsec&&value.st_ctim.tv_sec==other.value.st_ctim.tv_sec&&value.st_ctim.tv_nsec==other.value.st_ctim.tv_nsec;}};
void reject_nul(const std::string& path){
 gzFile file=gzopen(path.c_str(),"rb");if(!file)fail("cannot open input preflight");
 try{std::array<char,65536> buffer{};for(;;){const int n=gzread(file,buffer.data(),static_cast<unsigned>(buffer.size()));if(n<0)fail("gzip read/trailer error in NUL preflight");if(n==0){int error=Z_OK;(void)gzerror(file,&error);if(error!=Z_OK&&error!=Z_STREAM_END)fail("truncated gzip preflight");break;}
#if RAW_SE_MUTANT != 21
  if(std::memchr(buffer.data(),0,static_cast<std::size_t>(n)))fail("NUL byte in raw FASTQ");
#endif
 }}catch(...){gzclose(file);throw;}
 if(gzclose(file)!=Z_OK)fail("gzip close/trailer error in preflight");
}
struct Batch {std::vector<FastqRecord> records;Counts counts;};
}
Options parse_options(int argc,const char* const* argv){
 if(argc<1||argc>64)fail("argument count exceeds bounded CLI (64)");Options out;bool single=false,min=false,overlap=false,phred=false;std::set<std::string> seen;
#if RAW_SE_MUTANT == 22
 out.quality_ends=true;out.min_phred=20;
#endif
 for(int i=1;i<argc;++i){std::string key=argv[i];if(key.size()>1048576)fail("CLI token exceeds 1MiB");
  if(key=="--input")key="-i";if(key=="--output")key="-o";if(key=="--threads")key="-p";
  if(!seen.insert(key).second)fail("duplicate option: "+key);
  if(key=="--single-end"){single=true;continue;}if(key=="--help"||key=="-h"){out.help=true;continue;}
  if(key=="--stub5"||key=="--stub3"||key=="--scan-reads")fail("legacy stub/scan options conflict with --single-end");
  const auto value=[&](){if(i+1>=argc)fail("missing value: "+key);std::string_view text(argv[++i]);if(text.size()>1048576)fail("CLI token exceeds 1MiB");return text;};
  if(key=="-i")out.input=value();else if(key=="-o")out.output=value();else if(key=="--report")out.report=value();
  else if(key=="--adapter5"){out.adapter5=value();if(out.adapter5.empty())fail("empty adapter");}else if(key=="--adapter3"){out.adapter3=value();if(out.adapter3.empty())fail("empty adapter");}
  else if(key=="--min-length"){out.min_length=number(value());min=true;}else if(key=="--max-length")out.max_length=number(value());
  else if(key=="--clip5")out.clip5=number(value());else if(key=="--clip3")out.clip3=number(value());
  else if(key=="--adapter-min-overlap"){out.overlap=number(value());overlap=true;}
  else if(key=="--min-phred"){const auto n=number(value());if(n>93)fail("min-phred must be in 0..93");out.min_phred=static_cast<unsigned>(n);phred=true;}
  else if(key=="-p"){const auto n=number(value());if(n<1||n>64)fail("threads must be in 1..64");out.threads=static_cast<unsigned>(n);}
  else if(key=="--quality-trim"){const auto mode=value();if(mode!="none"&&mode!="ends")fail("quality-trim must be none or ends");out.quality_ends=mode=="ends";}
  else fail("unknown option: "+key);
 }
 if(out.help)return out;
 if(!single||!min||out.input.empty()||out.output.empty())fail("--single-end, -i, -o and explicit --min-length are required");
 if((!out.adapter5.empty()||!out.adapter3.empty())&&!overlap)fail("adapter requires explicit --adapter-min-overlap");
 if(overlap&&out.adapter5.empty()&&out.adapter3.empty())fail("overlap requires an adapter");
 if(out.quality_ends&&!phred)fail("quality ends requires explicit --min-phred");if(phred&&!out.quality_ends)fail("min-phred requires quality ends");
 validate_options(out);return out;
}
void validate_options(const Options& o){if(o.threads<1||o.threads>64||o.min_phred>93)fail("invalid thread/quality bound");if(o.max_length&&o.max_length<o.min_length)fail("max-length must be zero or >=min-length");
 for(const auto* adapter:{&o.adapter5,&o.adapter3})if(!adapter->empty()){if(o.overlap==0||o.overlap>adapter->size())fail("adapter overlap outside supplied adapter");for(char c:*adapter)if(!canonical(c))fail("adapters require literal A/C/G/T DNA");}}
void validate_record(const FastqRecord& r){
#if RAW_SE_MUTANT != 11
 if(r.header.size()<2||r.header.front()!='@'||r.plus.empty()||r.plus.front()!='+')fail("invalid FASTQ header/separator");printable(r.header);printable(r.plus);
#endif
#if RAW_SE_MUTANT != 12
 if(r.seq.size()!=r.qual.size())fail("FASTQ sequence/quality lengths differ");
#endif
#if RAW_SE_MUTANT != 13
 for(char c:r.seq)if(!base(c))fail("unsupported literal DNA byte");for(unsigned char q:r.qual)if(q<33||q>126)fail("illegal Phred+33 byte");
#endif
#if RAW_SE_MUTANT != 14
 if(traced(r.header)||traced(r.plus))fail("raw-only mode rejects FQT/FQC traced input; supply original raw reads");
#endif
}
Result transform(FastqRecord record,const Options& o){validate_record(record);Result out;auto& c=out.counts;c.input=1;c.input_bases=record.seq.size();std::size_t begin=0,end=record.seq.size();
 auto cut5=[&](std::size_t n,std::uint64_t& reads,std::uint64_t& bases){if(n){begin+=n;++reads;bases=n;}};
 auto cut3=[&](std::size_t n,std::uint64_t& reads,std::uint64_t& bases){if(n){end-=n;++reads;bases=n;}};
 cut5(std::min(o.clip5,end-begin),c.hard5_reads,c.hard5_bases);
#if RAW_SE_MUTANT == 1
 cut3(std::min(o.clip3,record.seq.size()),c.hard3_reads,c.hard3_bases);
#else
 cut3(std::min(o.clip3,end-begin),c.hard3_reads,c.hard3_bases);
#endif
 if(!o.adapter5.empty()){std::size_t length=0;for(std::size_t n=o.overlap;n<=std::min(end-begin,o.adapter5.size());++n){
#if RAW_SE_MUTANT == 8
  const auto adapter=std::string_view(o.adapter5).substr(0,n);
#else
  const auto adapter=std::string_view(o.adapter5).substr(o.adapter5.size()-n);
#endif
  if(exact(std::string_view(record.seq).substr(begin,n),adapter)){
#if RAW_SE_MUTANT == 3
   length=n;break;
#else
   length=n;
#endif
  }}cut5(length,c.adapter5_reads,c.adapter5_bases);}
 if(!o.adapter3.empty())for(std::size_t i=begin;i<end;++i){auto n=std::min(o.adapter3.size(),end-i);
#if RAW_SE_MUTANT == 5
  if(n<=o.overlap)continue;
#else
  if(n<o.overlap)continue;
#endif
#if RAW_SE_MUTANT == 6
  n=o.overlap;
#endif
  if(exact(std::string_view(record.seq).substr(i,n),std::string_view(o.adapter3).substr(0,n))){
#if RAW_SE_MUTANT == 4
   cut3(end-i-1,c.adapter3_reads,c.adapter3_bases);
#else
   cut3(end-i,c.adapter3_reads,c.adapter3_bases);
#endif
   break;
  }}
 if(o.quality_ends){const auto first=begin,last=end;
#if RAW_SE_MUTANT == 9
  while(begin<end&&unsigned(static_cast<unsigned char>(record.qual[begin])-33)<=o.min_phred)++begin;
  while(end>begin&&unsigned(static_cast<unsigned char>(record.qual[end-1])-33)<=o.min_phred)--end;
#else
  while(begin<end&&unsigned(static_cast<unsigned char>(record.qual[begin])-33)<o.min_phred)++begin;
  while(end>begin&&unsigned(static_cast<unsigned char>(record.qual[end-1])-33)<o.min_phred)--end;
#endif
  c.quality_bases=(begin-first)+(last-end);c.quality_reads=c.quality_bases!=0;
 }
 out.begin=begin;out.end=end;
 if(end==begin)c.empty=1;else if(end-begin<o.min_length)c.short_reads=1;
#if RAW_SE_MUTANT == 10
 else if(o.max_length&&end-begin>=o.max_length)c.long_reads=1;
#else
 else if(o.max_length&&end-begin>o.max_length)c.long_reads=1;
#endif
 else{out.keep=true;c.retained=1;c.retained_bases=end-begin;}
 if(!out.keep)c.discarded_bases=end-begin;
#if RAW_SE_MUTANT == 23
 record.seq.reserve(record.seq.capacity()+1024);
#endif
 record.seq=record.seq.substr(begin,end-begin);
#if RAW_SE_MUTANT != 2
 record.qual=record.qual.substr(begin,end-begin);
#endif
#if RAW_SE_MUTANT == 15
 record.plus="+";
#endif
 out.record=std::move(record);
#if RAW_SE_MUTANT == 19
 ++c.retained;
#endif
 c.check();return out;
}
void Counts::add(const Counts& b){
#define FIELD(f) f=raw_se::add(f,b.f)
 FIELD(input);FIELD(retained);FIELD(empty);FIELD(short_reads);FIELD(long_reads);FIELD(input_bases);FIELD(retained_bases);FIELD(discarded_bases);FIELD(hard5_reads);FIELD(hard3_reads);FIELD(adapter5_reads);FIELD(adapter3_reads);FIELD(quality_reads);FIELD(hard5_bases);FIELD(hard3_bases);FIELD(adapter5_bases);FIELD(adapter3_bases);FIELD(quality_bases);
#undef FIELD
 check();}
void Counts::check()const{const auto dropped=raw_se::add(raw_se::add(empty,short_reads),long_reads);
#if RAW_SE_MUTANT != 19
 if(input!=raw_se::add(retained,dropped))fail("read conservation");
#endif
 auto clipped=raw_se::add(raw_se::add(hard5_bases,hard3_bases),raw_se::add(adapter5_bases,adapter3_bases));clipped=raw_se::add(clipped,quality_bases);if(input_bases!=raw_se::add(raw_se::add(clipped,retained_bases),discarded_bases))fail("base conservation");
 for(auto count:{hard5_reads,hard3_reads,adapter5_reads,adapter3_reads,quality_reads})if(count>input)fail("clipped-read count exceeds input");}
std::string report(const Options& o,const Counts& c){c.check();std::ostringstream s;s<<"raw_single_end_v1\nquality_policy\t"<<(o.quality_ends?"ends":"none")<<"\nmin_phred\t"<<(o.quality_ends?std::to_string(o.min_phred):"unused")<<"\nmin_length\t"<<o.min_length<<"\nmax_length\t"<<o.max_length<<"\nclip5\t"<<o.clip5<<"\nclip3\t"<<o.clip3<<"\nadapter5\t"<<o.adapter5<<"\nadapter3\t"<<o.adapter3<<"\nadapter_min_overlap\t"<<o.overlap<<"\nread_clip_categories\toverlap_allowed\n";
#define FIELD(f) s<<#f<<'\t'<<c.f<<'\n'
 FIELD(input);FIELD(retained);FIELD(empty);FIELD(short_reads);FIELD(long_reads);FIELD(input_bases);FIELD(retained_bases);FIELD(discarded_bases);FIELD(hard5_reads);FIELD(hard3_reads);FIELD(adapter5_reads);FIELD(adapter3_reads);FIELD(quality_reads);FIELD(hard5_bases);FIELD(hard3_bases);FIELD(adapter5_bases);FIELD(adapter3_bases);FIELD(quality_bases);
#undef FIELD
 s<<"dropped\t"<<raw_se::add(raw_se::add(c.empty,c.short_reads),c.long_reads)<<'\n';
 if(!s)fail("report formatting");return s.str();}
Counts run(const Options& o){validate_options(o);
#if RAW_SE_MUTANT != 26
 if(o.input=="-"||o.input=="/dev/stdin")fail("stdin aliases '-' and '/dev/stdin' are unsupported in raw mode");
#endif
#if RAW_SE_MUTANT != 17
 if(same_file(o.input,o.output)||(!o.report.empty()&&(same_file(o.input,o.report)||same_file(o.output,o.report))))fail("input/output/report identify the same file");
#endif
 const Stamp stamp(o.input);reject_nul(o.input);if(!stamp.same(Stamp(o.input)))fail("input changed during preflight");
 auto reader=make_fastq_reader(o.input,1);FastqWriter writer(o.output,output_wants_gzip(o.output),1);fqdup::r8qkeep::WorkerPool pool(static_cast<int>(o.threads));std::deque<std::future<Batch>> pending;Counts totals;
 const std::size_t bound=2*o.threads;const auto drain=[&](){auto batch=pending.front().get();pending.pop_front();
#if RAW_SE_MUTANT == 20
  std::reverse(batch.records.begin(),batch.records.end());
#endif
  for(const auto& r:batch.records)writer.write(r);totals.add(batch.counts);};
 const auto read_one=[&](FastqRecord& r){
#if RAW_SE_MUTANT == 24
  try{return reader->read(r);}catch(...){return false;}
#else
  return reader->read(r);
#endif
 };
 try{bool eof=false;while(!eof){Batch batch;batch.records.reserve(64);FastqRecord r;for(unsigned n=0;n<64;++n){if(!read_one(r)){eof=true;break;}batch.records.push_back(std::move(r));}
  if(!batch.records.empty()){pending.push_back(pool.submit([batch=std::move(batch),&o]()mutable{Batch result;result.records.reserve(batch.records.size());for(auto& r:batch.records){
#if RAW_SE_MUTANT == 25
   try{
#endif
   auto value=transform(std::move(r),o);result.counts.add(value.counts);if(value.keep)result.records.push_back(std::move(value.record));
#if RAW_SE_MUTANT == 25
   }catch(const std::exception&){}
#endif
  }return result;}));}
  if(pending.size()>=bound)drain();}
  while(!pending.empty())drain();pool.stop();if(!stamp.same(Stamp(o.input)))fail("input changed during processing");writer.close();
 }catch(...){pool.stop();throw;}
 totals.check();return totals;
}
} // namespace
int trim_single_main(int argc,char** argv){try{const auto o=fqdup::raw_se::parse_options(argc,argv);if(o.help){std::cerr<<"Usage: fqdup trim --single-end -i FILE -o FILE --min-length N [--clip5 N] [--clip3 N] [--max-length N]\n  [--adapter5 DNA] [--adapter3 DNA] [--adapter-min-overlap K]\n  [--quality-trim none|ends] [--min-phred 0..93] [-p 1..64] [--report FILE]\nRaw regular FASTQ only; known FQT/FQC traces rejected. Quality defaults to none.\nExact matches can clip an indistinguishable biological adapter motif.\n";return 0;}
 const auto totals=fqdup::raw_se::run(o);const auto text=fqdup::raw_se::report(o,totals);if(!o.report.empty()){std::ofstream file(o.report,std::ios::binary);file<<text;file.close();if(!file)throw std::runtime_error("single-end: report write failed");}std::cerr<<text;return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
