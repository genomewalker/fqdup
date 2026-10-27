#include "fqdup/trim_single.hpp"
#include "raw_se_new_counter.hpp"
#include <zlib.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <sys/wait.h>
#include <unistd.h>
using namespace fqdup::raw_se;
static void need(bool yes,const char* reason){if(!yes)throw std::runtime_error(reason);}
template<class F>static bool rejects(F f){try{f();return false;}catch(const std::exception&){return true;}}
static FastqRecord record(std::string seq,std::string qual=""){if(qual.empty())qual=std::string(seq.size(),'I');return {"@raw metadata",std::move(seq),"+raw metadata",std::move(qual)};}
static Result apply(FastqRecord r,const Options& o,const char* reason){try{return transform(std::move(r),o);}catch(const std::exception&){throw std::runtime_error(reason);}}
static std::vector<std::string> args(){return {"trim","--single-end","-i","input.fq","-o","output.fq","--min-length","0"};}
static Options parse(const std::vector<std::string>& a){std::vector<const char*> p;for(const auto& s:a)p.push_back(s.c_str());return parse_options(static_cast<int>(p.size()),p.data());}
struct Temp{std::string path;Temp(){char name[]="/tmp/fqdup-raw-se-XXXXXX";const auto p=::mkdtemp(name);if(!p)throw std::runtime_error("SETUP_TEMP");path=p;}~Temp(){std::error_code ec;std::filesystem::remove_all(path,ec);}};
static void put(const std::string& path,const std::string& text){std::ofstream f(path,std::ios::binary);f.write(text.data(),static_cast<std::streamsize>(text.size()));f.close();if(!f)throw std::runtime_error("SETUP_WRITE");}
static std::string get(const std::string& path){std::ifstream f(path,std::ios::binary);std::ostringstream s;s<<f.rdbuf();if(!f.good()&&!f.eof())throw std::runtime_error("SETUP_READ");return s.str();}
static std::string encode(const FastqRecord& r){return r.header+'\n'+r.seq+'\n'+r.plus+'\n'+r.qual+'\n';}
// Real production entry, isolated so inherited stdin/cwd cannot affect other
// checks. Called only after earlier worker pools have stopped and joined.
static int stdin_entry(const std::string& cwd,const std::string& input,const std::string& diagnostic,std::vector<std::string> arguments){
 const pid_t child=::fork();need(child>=0,"SETUP_STDIN_ALIAS");
 if(child==0){
  const int in=::open(input.c_str(),O_RDONLY|O_CLOEXEC);
  const int err=::open(diagnostic.c_str(),O_WRONLY|O_CREAT|O_TRUNC|O_CLOEXEC,0600);
  if(in<0||err<0||::chdir(cwd.c_str())!=0||::dup2(in,STDIN_FILENO)<0||::dup2(err,STDERR_FILENO)<0)::_exit(98);
  if(in!=STDIN_FILENO)::close(in);if(err!=STDERR_FILENO)::close(err);
  std::vector<char*> values;for(auto& argument:arguments)values.push_back(argument.data());
  const int status=trim_single_main(static_cast<int>(values.size()),values.data());::_exit(status);
 }
 int status=0;pid_t waited;do{waited=::waitpid(child,&status,0);}while(waited<0&&errno==EINTR);
 need(waited==child&&WIFEXITED(status),"STDIN_ALIAS_CHILD_STATUS");need(WEXITSTATUS(status)!=98,"SETUP_STDIN_ALIAS");return WEXITSTATUS(status);
}
static void stdin_alias_regression(){
 Temp temporary;const std::string named=temporary.path+"/-",redirected=temporary.path+"/actual.fq",separate=temporary.path+"/output.fq",summary=temporary.path+"/report.tsv",diagnostic=temporary.path+"/stderr.txt";
 auto a=record("ACGT");a.header="@named-dash metadata";const auto named_bytes=encode(a);
 auto b=record("NNNN");b.header="@redirected-stdin metadata";const auto redirected_bytes=encode(b);
 const std::string output_before="output sentinel\n",report_before="report sentinel\n";
 for(const std::string alias:{"-","/dev/stdin"}){
  put(named,named_bytes);put(redirected,redirected_bytes);put(separate,output_before);put(summary,report_before);
  // '-' uses the actual redirected input itself as output, reproducing the
  // destructive ownership counterexample. /dev/stdin uses a separate output.
  const auto output=alias=="-"?redirected:separate;
  const int status=stdin_entry(temporary.path,redirected,diagnostic,{"trim","--single-end","-i",alias,"-o",output,"--min-length","1","--report",summary});
  need(status==1&&get(diagnostic)=="single-end: stdin aliases '-' and '/dev/stdin' are unsupported in raw mode\n","STDIN_ALIAS");
  need(get(named)==named_bytes&&get(redirected)==redirected_bytes&&get(separate)==output_before&&get(summary)==report_before,"STDIN_ALIAS");
 }
 // './-' is an ordinary filename in the actual factory; it must read A, not
 // redirected B, and must not be rejected by an overbroad spelling check.
 put(named,named_bytes);put(redirected,redirected_bytes);put(separate,output_before);put(summary,report_before);
 const int status=stdin_entry(temporary.path,redirected,diagnostic,{"trim","--single-end","-i","./-","-o",separate,"--min-length","1","--report",summary});
 need(status==0&&get(separate)==named_bytes&&get(named)==named_bytes&&get(redirected)==redirected_bytes,"ORDINARY_DASH_PATH");
 need(get(summary).find("raw_single_end_v1\n")==0,"ORDINARY_DASH_PATH");
}
static bool legal_number(const std::string& text){if(text.empty())return false;std::size_t n=0;for(unsigned char c:text){if(c<'0'||c>'9'||n>(std::numeric_limits<std::size_t>::max()-(c-'0'))/10)return false;n=n*10+(c-'0');}return true;}
// Independent interval reference for random bounded valid records: enumerate all
// legal end matches before selecting, then use a quality mask on surviving ends.
static std::pair<std::size_t,std::size_t> naive(const FastqRecord& r,const Options& o){
 std::size_t left=std::min(o.clip5,r.seq.size()),right=r.seq.size();right-=std::min(o.clip3,right-left);
 auto up=[](char c){return c>='a'&&c<='z'?char(c-'a'+'A'):c;};
 auto eq=[&](std::size_t start,const std::string& a,std::size_t offset,std::size_t n){for(std::size_t j=0;j<n;++j){const char c=up(r.seq[start+j]);if((c!='A'&&c!='C'&&c!='G'&&c!='T')||c!=up(a[offset+j]))return false;}return true;};
 std::size_t longest=0;if(!o.adapter5.empty())for(std::size_t n=1;n<=right-left&&n<=o.adapter5.size();++n)if(n>=o.overlap&&eq(left,o.adapter5,o.adapter5.size()-n,n))longest=n;left+=longest;
 if(!o.adapter3.empty()){std::vector<std::size_t> matches;for(std::size_t start=left;start<right;++start){const auto n=std::min(o.adapter3.size(),right-start);if(n>=o.overlap&&eq(start,o.adapter3,0,n))matches.push_back(start);}if(!matches.empty())right=*std::min_element(matches.begin(),matches.end());}
 if(o.quality_ends){std::vector<bool> low(right-left);for(std::size_t i=left;i<right;++i)low[i-left]=unsigned(static_cast<unsigned char>(r.qual[i])-33)<o.min_phred;const auto base=left;while(left<right&&low[left-base])++left;while(right>left&&low[right-1-base])--right;}
 return {left,right};
}
static void truth_panel(const std::string& path){std::ofstream out(path);if(!out)throw std::runtime_error("SETUP_TRUTH");out<<"case\tdamage5\tdamage3\tlow_quality\tquality_ends\tfirst10_retained\tlast10_retained\tendogenous_wrongly_clipped\tadapter_bases_removed\tretained_length\tfirst10_mask_5prime_inward\tlast10_mask_3prime_inward\n";std::map<std::size_t,unsigned> histogram;unsigned id=0;
 for(bool d5:{false,true})for(bool d3:{false,true})for(bool low:{false,true})for(bool quality:{false,true})for(bool a5:{false,true})for(bool a3:{false,true}){
  std::string molecule="CCACGTCGACGTACGTACGTACGTACGTACGTACGTGGGG";if(d5)molecule.front()='T';if(d3)molecule.back()='A';const std::string left=a5?"AGTCGA":"",right=a3?"TTGACA":"";
  auto r=record(left+molecule+right);if(low)for(unsigned j=0;j<5;++j){r.qual[left.size()+j]='!';r.qual[left.size()+molecule.size()-1-j]='!';}
  Options o;o.min_length=1;o.adapter5="AGTCGA";o.adapter3="TTGACA";o.overlap=4;o.quality_ends=quality;o.min_phred=20;const auto reference=naive(r,o);const auto v=apply(r,o,"TRUTH");need(v.begin==reference.first&&v.end==reference.second,"TRUTH");
  auto count=[&](std::size_t a,std::size_t b){return std::min(b,v.end)>std::max(a,v.begin)?std::min(b,v.end)-std::max(a,v.begin):0;};const auto kept=count(left.size(),left.size()+molecule.size());const auto first=count(left.size(),left.size()+10),last=count(left.size()+molecule.size()-10,left.size()+molecule.size());const auto adapter_kept=count(0,left.size())+count(left.size()+molecule.size(),r.seq.size());
  std::string first_mask,last_mask;for(unsigned position=0;position<10;++position){const auto f=left.size()+position,l=left.size()+molecule.size()-1-position;first_mask+=f>=v.begin&&f<v.end?'1':'0';last_mask+=l>=v.begin&&l<v.end?'1':'0';}need(std::count(first_mask.begin(),first_mask.end(),'1')==static_cast<long>(first)&&std::count(last_mask.begin(),last_mask.end(),'1')==static_cast<long>(last),"TRUTH");need(quality||kept==molecule.size(),"TRUTH");out<<id++<<'\t'<<d5<<'\t'<<d3<<'\t'<<low<<'\t'<<quality<<'\t'<<first<<'\t'<<last<<'\t'<<molecule.size()-kept<<'\t'<<left.size()+right.size()-adapter_kept<<'\t'<<(v.keep?v.record.seq.size():0)<<'\t'<<first_mask<<'\t'<<last_mask<<'\n';++histogram[v.keep?v.record.seq.size():0];
 }
 // Deliberately indistinguishable biological adapter motif: report wrong clipping.
 {Options o;o.min_length=1;o.adapter3="TTGACA";o.overlap=4;auto v=apply(record("CCCCCTTGACA"),o,"TRUTH");need(v.record.seq=="CCCCC","TRUTH");out<<"indistinguishable_biological_motif\t0\t0\t0\t0\t5\t4\t6\t0\t5\t1111100000\t0000001111\n";++histogram[5];}
 for(auto x:histogram)out<<"length_histogram\t"<<x.first<<'\t'<<x.second<<'\n';out.close();if(!out)throw std::runtime_error("TRUTH_WRITE");}
int main(int argc,char** argv){try{
 Options o;o.min_length=1;
 {auto v=apply(record("ACGT"),o,"COUNTS");need(v.counts.input==1&&v.counts.retained==1&&v.counts.empty==0&&v.counts.retained_bases==4,"COUNTS");}
 {auto p=o;p.clip5=1;auto r=record("ACGT","!5I#");auto v=apply(r,p,"ALIGNMENT");need(v.record.seq=="CGT"&&v.record.qual=="5I#","ALIGNMENT");}
 {auto p=o;p.clip5=100;p.clip3=100;auto v=apply(record("ACGT"),p,"HARD_CLIP");need(!v.keep&&v.counts.empty==1&&v.counts.hard5_bases==4&&v.counts.hard3_bases==0,"HARD_CLIP");}
 {auto p=o;p.clip5=2;p.clip3=2;auto v=apply(record("ACGTN","12345"),p,"HARD_CLIP");need(v.record.seq=="G"&&v.record.qual=="3"&&v.counts.hard5_bases==2&&v.counts.hard3_bases==2,"HARD_CLIP");}
 {auto p=o;p.adapter5="ACAC";p.overlap=2;auto v=apply(record("ACACG"),p,"LONGEST5");need(v.record.seq=="G","LONGEST5");}
 {auto p=o;p.adapter3="ACGT";p.overlap=2;auto v=apply(record("GACGTGG"),p,"OFF_BY_ONE");need(v.record.seq=="G","OFF_BY_ONE");
  need(apply(record("GAC"),p,"MIN_OVERLAP").record.seq=="G","MIN_OVERLAP");need(apply(record("GA"),p,"MIN_OVERLAP").record.seq=="GA","MIN_OVERLAP");
  need(apply(record("GACCTG"),p,"INTERNAL_PARTIAL").record.seq=="GACCTG","INTERNAL_PARTIAL");need(apply(record("GANG"),p,"N_WILDCARD").record.seq=="GANG","N_WILDCARD");}
 {auto p=o;p.adapter5="ACGT";p.overlap=2;need(apply(record("GTAA"),p,"ORIENTATION").record.seq=="AA","ORIENTATION");need(apply(record("ACAA"),p,"ORIENTATION").record.seq=="ACAA","ORIENTATION");}
 {auto p=o;p.adapter3="ACGT";p.overlap=2;need(apply(record("GGT"),p,"ORIENTATION").record.seq=="GGT","ORIENTATION");}
 {auto p=o;p.quality_ends=true;p.min_phred=20;need(apply(record("AAACAAA","!!!5!!!"),p,"QUALITY_EQUAL").record.seq=="C","QUALITY_EQUAL");auto v=apply(record("ACGT","!!!!"),p,"QUALITY_EQUAL");need(!v.keep&&v.counts.empty==1,"QUALITY_EQUAL");}
 {auto p=o;p.max_length=4;need(apply(record("ACGT"),p,"MAX_EQUAL").keep,"MAX_EQUAL");need(!apply(record("ACGTA"),p,"MAX_EQUAL").keep,"MAX_EQUAL");}
 for(std::size_t n=0;n<=6;++n)for(std::size_t min=0;min<=6;++min)for(std::size_t max=0;max<=7;++max){if(max&&max<min)continue;auto p=o;p.min_length=min;p.max_length=max;auto v=apply(record(std::string(n,'A')),p,"DROP_REASONS");need(v.counts.empty==(n==0)&&v.counts.short_reads==(n>0&&n<min)&&v.counts.long_reads==(n>0&&n>=min&&max>0&&n>max),"DROP_REASONS");}
 {auto r=record("ACGT");r.header="bad";need(rejects([&]{validate_record(r);}),"HEADER");r=record("ACGT");r.plus="bad";need(rejects([&]{validate_record(r);}),"HEADER");r=record("ACGT");r.header+='\t';need(rejects([&]{validate_record(r);}),"HEADER");}
 {auto r=record("ACGT");r.qual.pop_back();need(rejects([&]{validate_record(r);}),"LENGTH");}
 {auto r=record("ACGT");r.qual[0]=' ';need(rejects([&]{validate_record(r);}),"ALPHABET");r=record("ACGU");need(rejects([&]{validate_record(r);}),"ALPHABET");}
 {auto r=record("ACGT");r.header+=" FQT=v1;ch=r1";need(rejects([&]{validate_record(r);}),"TRACE");r.header="@raw";r.plus+=" FQC=v1;edits=1:A>G";need(rejects([&]{validate_record(r);}),"TRACE");}
 {auto r=record("AcgNRY");auto v=apply(r,o,"METADATA");need(v.record.header==r.header&&v.record.plus==r.plus&&v.record.seq==r.seq,"METADATA");}
 {auto a=args();Options p;try{p=parse(a);}catch(...){throw std::runtime_error("DEFAULT_NONE");}need(!p.quality_ends,"DEFAULT_NONE");}
 for(const std::string bad:{"2x","-1","+1","NaN","nan","inf","infinity","1e3","1.0","184467440737095516160",""}){auto a=args();a.back()=bad;need(rejects([&]{(void)parse(a);}),"CLI_FULL_TOKEN");}
 {auto a=args();a.back()=std::to_string(std::numeric_limits<std::size_t>::max());need(!rejects([&]{(void)parse(a);}),"CLI_FULL_TOKEN");}
 for(const std::string bad:{"NaN","inf","1.5","-1","94"}){auto a=args();a.insert(a.end(),{"--quality-trim","ends","--min-phred",bad});need(rejects([&]{(void)parse(a);}),"CLI_CONTRACT");}
 {auto a=args();a.back()="10";a.insert(a.end(),{"--max-length","5"});need(rejects([&]{(void)parse(a);}),"CLI_CONTRACT");}
 for(const auto extra:std::vector<std::vector<std::string>>{{"--quality-trim","ends"},{"--min-phred","20"},{"--adapter3","ACGT"},{"--adapter3","ACNT","--adapter-min-overlap","2"},{"--adapter3","ACGT","--adapter-min-overlap","5"},{"--stub5","ACGT"},{"--scan-reads","2"},{"-p","0"},{"-p","65"},{"--quality-trim","ends","--min-phred","94"},{"--clip5","-1"},{"--max-length","-1"}}){auto a=args();a.insert(a.end(),extra.begin(),extra.end());need(rejects([&]{(void)parse(a);}),"CLI_CONTRACT");}
 {Counts a,b;a.input=a.retained=std::numeric_limits<std::uint64_t>::max();b.input=b.retained=1;need(rejects([&]{a.add(b);}),"OVERFLOW");}
 {auto r=record(std::string(4096,'A'));raw_se_alloc::start();auto v=transform(std::move(r),o);raw_se_alloc::stop();need(v.record.seq.size()==4096&&raw_se_alloc::calls<=2&&raw_se_alloc::bytes<=8194,"ALLOCATION");}
 std::mt19937_64 rng(0x5241575345323032ULL);const char alphabet[]="ACGTNacgtRY";
 for(unsigned trial=0;trial<8192;++trial){std::string seq(trial%257,'A'),qual(seq.size(),'!');for(char& c:seq)c=alphabet[rng()%(sizeof alphabet-1)];for(char& q:qual)q=char(33+rng()%94);auto r=record(seq,qual);Options p;p.min_length=trial%12;p.max_length=trial%3?0:300;p.clip5=rng()%20;p.clip3=rng()%20;p.adapter5="ACGT";p.adapter3="TGCA";p.overlap=1+rng()%4;p.quality_ends=trial%2;p.min_phred=rng()%94;const auto expected=naive(r,p);auto v=apply(r,p,"FUZZ_TRANSFORM");need(v.begin==expected.first&&v.end==expected.second&&v.record.seq==r.seq.substr(v.begin,v.end-v.begin)&&v.record.qual==r.qual.substr(v.begin,v.end-v.begin),"FUZZ_TRANSFORM");need(v.record.header==r.header&&v.record.plus==r.plus,"FUZZ_TRANSFORM");v.counts.check();}
 // Every quality byte and sequence byte gets a real validation expectation.
 for(unsigned byte=0;byte<256;++byte){auto r=record("A");r.qual[0]=char(byte);need(rejects([&]{validate_record(r);})==(byte<33||byte>126),"BYTE_VALIDATION");}
 for(unsigned byte=0;byte<256;++byte){auto r=record("A");r.seq[0]=char(byte);const bool legal=std::string("ACGTRYSWKMBDHVNacgtryswkmbdhvn").find(char(byte))!=std::string::npos;need(rejects([&]{validate_record(r);})!=legal,"BYTE_VALIDATION");r=record("A");r.header+=char(byte);need(rejects([&]{validate_record(r);})==(byte<32||byte>126),"BYTE_VALIDATION");}
 for(unsigned trial=0;trial<8192;++trial){auto a=args();std::string text(trial%65,' ');for(char& c:text)c=char(1+rng()%255);a.back()=text;const bool expected=legal_number(text);std::size_t size=0;for(const auto& s:a)size+=s.size();raw_se_alloc::start();const bool rejected=rejects([&]{(void)parse(a);});raw_se_alloc::stop();need(rejected!=expected,"CLI_FUZZ");need(raw_se_alloc::calls<=10*a.size()+32&&raw_se_alloc::bytes<=32*(size+a.size()*128)+4096,"CLI_ALLOCATION");}
 Temp tmp;Options io=o;io.input=tmp.path+"/in.fq";io.output=tmp.path+"/out.fq";
 {const std::string text="@a\nAC\n+metadata\nII\n";put(io.input,text);auto p=io;p.output=p.input;(void)rejects([&]{(void)run(p);});need(get(io.input)==text,"SAME_FILE");need(rejects([&]{(void)run(p);}),"SAME_FILE");}
 {const std::string text=std::string("@x",2)+char(0)+"junk\npadding\nAC\n+\nII\n";put(io.input,text);need(rejects([&]{(void)run(io);}),"NUL_PREFLIGHT");}
 for(const std::string truncated:{"@a\n","@a\nAC\n","@a\nAC\n+\n"}){put(io.input,truncated);need(rejects([&]{(void)run(io);}),"TRUNCATED");}
 {put(io.input,"@bad\nA\n+\n \n");need(rejects([&]{(void)run(io);}),"DROPPED_VALIDATION");auto p=io;p.min_length=999;need(rejects([&]{(void)run(p);}),"DROPPED_VALIDATION");}
 for(unsigned length=1;length<=8;++length){const auto text=encode(record(std::string(length,'A')));for(std::size_t cut=1;cut+1<text.size();++cut){put(io.input,text.substr(0,cut));need(rejects([&]{(void)run(io);}),"TRUNCATED_FUZZ");}}
 {auto p=io;p.input=tmp.path+"/input.fq.gz";p.output=tmp.path+"/output.fq.gz";const auto text=encode(record("ACGTRYNN"));gzFile f=gzopen(p.input.c_str(),"wb");need(f!=nullptr,"SETUP_GZIP");need(gzwrite(f,text.data(),static_cast<unsigned>(text.size()))==static_cast<int>(text.size()),"SETUP_GZIP");need(gzclose(f)==Z_OK,"SETUP_GZIP");const auto compressed=get(p.input);(void)run(p);f=gzopen(p.output.c_str(),"rb");need(f!=nullptr,"SETUP_GZIP");std::array<char,256> buf{};const int n=gzread(f,buf.data(),buf.size());need(n>=0&&std::string(buf.data(),static_cast<std::size_t>(n))==text,"GZIP_IDENTITY");need(gzclose(f)==Z_OK,"GZIP_IDENTITY");
  for(unsigned removed=1;removed<=8;++removed){put(p.input,compressed.substr(0,compressed.size()-removed));need(rejects([&]{(void)run(p);}),"GZIP_TRUNCATED");}}
 {std::string input,expected;for(unsigned i=0;i<513;++i){auto r=record("ACGTACGT");r.header="@r"+std::to_string(i)+" metadata";r.plus="+keep "+std::to_string(i);input+=encode(r);auto p=o;p.clip5=1;p.clip3=2;expected+=encode(apply(r,p,"ORDER").record);}put(io.input,input);std::string totals;for(unsigned threads:{1U,8U,16U}){auto p=io;p.threads=threads;p.clip5=1;p.clip3=2;auto counts=run(p);need(get(p.output)==expected,"ORDER");const auto text=report(p,counts);if(totals.empty())totals=text;need(text==totals,"THREAD_COUNTS");}}
 {std::string input;for(unsigned i=0;i<70;++i)input+=encode(record("ACGT"));input+="@bad\nACGT\n+\n!!! \n";put(io.input,input);auto p=io;p.threads=8;need(rejects([&]{(void)run(p);}),"WORKER_EXCEPTION");put(io.input,input+"@truncated\nACGT\n");need(rejects([&]{(void)run(p);}),"READER_EXCEPTION");}
 {put(io.input,encode(record("ACGT")));auto p=io;p.output="/dev/full";need(std::filesystem::exists(p.output),"SETUP_DEV_FULL");need(rejects([&]{(void)run(p);}),"WRITER_EXCEPTION");}
 truth_panel(argc==2?argv[1]:tmp.path+"/truth.tsv");
 stdin_alias_regression();
 std::cout<<"RAW_SE_PROPERTIES_PASS seed=0x5241575345323032 transform_fuzz=8192 cli_fuzz=8192 qualities=256 p=1,8,16 truth=64+motif\n";return 0;
 }catch(const std::exception& e){raw_se_alloc::stop();std::cerr<<e.what()<<'\n';return 1;}}
