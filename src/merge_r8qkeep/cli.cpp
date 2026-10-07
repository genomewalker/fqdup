#include "pipeline.hpp"
#include "fqdup/fastq_types.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
int merge_main(int argc,char** argv);
namespace {
using namespace fqdup::r8qkeep;
std::vector<std::string> paths(const std::string& text){std::vector<std::string> out;std::istringstream in(text);std::string p;while(std::getline(in,p,',')){if(p.empty())throw std::runtime_error("empty input lane");out.push_back(p);}return out;}
std::vector<FamilyRead> input(const std::vector<std::string>& paths){auto reader=make_chained_fastq_reader(paths,1);std::vector<FamilyRead> out;FastqRecord r;while(reader->read(r)){auto end=r.header.find_first_of(" \t");auto id=r.header.substr(1,end==std::string::npos?end:end-1);if(id.size()>2&&id[id.size()-2]=='/'&&(id.back()=='1'||id.back()=='2'))id.resize(id.size()-2);if(r.seq.size()!=r.qual.size())throw std::runtime_error("FASTQ lengths");out.push_back({r.header,id,r.seq,r.qual});}return out;}
void write(const std::filesystem::path& path,const std::string& text){std::ofstream out(path,std::ios::binary);out<<text;out.close();if(!out)throw std::runtime_error("write failed: "+path.string());}
double number(const std::string& text){std::size_t n;const double x=std::stod(text,&n);if(n!=text.size()||!std::isfinite(x))throw std::runtime_error("invalid numeric argument");return x;}
}
int merge_r8qkeep_main(int argc,char** argv){
    for(int i=1;i<argc;++i)if(std::string(argv[i])=="--legacy-merge"){std::vector<char*> args;for(int j=0;j<argc;++j)if(j!=i)args.push_back(argv[j]);return merge_main(args.size(),args.data());}
    try{std::vector<std::string> r1,r2;std::string outdir;double cut=.316;int threads=1;bool debug=false;NativeSettings settings;
        for(int i=1;i<argc;++i){const std::string a=argv[i];const auto value=[&](){if(i+1>=argc)throw std::runtime_error("missing value: "+a);return std::string(argv[++i]);};
            if(a=="--help"||a=="-h"){std::cout<<"Usage: fqdup merge -1 R1.fq[,LANE...] -2 R2.fq[,LANE...] -o NEW_OUTPUT_DIR [--cutoff 0.316] [-p THREADS] [--debug-stages]\nOutputs: merged.fq r1.fq r2.fq orphan.fq outcomes.jsonl\nNative options: --adapter1 --adapter2 --library-type --min-length --min-overlap --max-mm-rate\nThe r8qkeep recipe enables short-tail/gapped clipping and disables poly-G/internal-panel filtering.\nUse --legacy-merge for the previous repository merge command and its options.\n";return 0;}
            if(a=="-1"){auto p=paths(value());r1.insert(r1.end(),p.begin(),p.end());}else if(a=="-2"){auto p=paths(value());r2.insert(r2.end(),p.begin(),p.end());}else if(a=="-o")outdir=value();else if(a=="--cutoff")cut=number(value());else if(a=="-p"||a=="--threads")threads=std::stoi(value());else if(a=="--adapter1")settings.adapter1=value();else if(a=="--adapter2")settings.adapter2=value();else if(a=="--library-type")settings.library_type=value();else if(a=="--min-length")settings.min_length=std::stoi(value());else if(a=="--min-overlap")settings.min_overlap=std::stoi(value());else if(a=="--max-mm-rate")settings.max_mismatch_rate=number(value());else if(a=="--debug-stages")debug=true;else if(a=="--no-poly-g"||a=="--no-internal-panel"||a=="--experimental-short-tail"||a=="--experimental-gapped-clipping"){}else throw std::runtime_error("unsupported r8qkeep option: "+a+"; use --legacy-merge for the previous command");
        }
        if(r1.empty()||r2.empty()||r1.size()!=r2.size()||outdir.empty()||cut<0||cut>1||threads<1||settings.min_length<1||settings.min_overlap<1||settings.max_mismatch_rate<0||settings.max_mismatch_rate>1||(settings.library_type!="ds"&&settings.library_type!="ss"))throw std::runtime_error("invalid merge arguments");
        if(!std::filesystem::create_directory(outdir))throw std::runtime_error("output directory exists");auto library=prepare_library(input(r1),input(r2),settings,threads);auto final=emit_library(library,cut);const std::filesystem::path out=outdir;
        write(out/"merged.fq",final.merged);write(out/"r1.fq",final.r1);write(out/"r2.fq",final.r2);write(out/"orphan.fq",final.orphan);write(out/"outcomes.jsonl",final.outcomes);
        if(debug){std::filesystem::create_directory(out/"debug");write(out/"debug"/"family-edits.tsv",library.family.edits_tsv);write(out/"debug"/"families.tsv",library.family.families_tsv);write(out/"debug"/"stats.json",library.family.stats_json);write(out/"debug"/"native.tsv",library.native.native_tsv);write(out/"debug"/"nominal.tsv",library.observations.nominal_tsv);write(out/"debug"/"full.tsv",library.observations.full_tsv);write(out/"debug"/"proposed-edits.tsv",library.observations.proposed_edits_tsv);write(out/"debug"/"nominal.txt",library.nominal_prior);write(out/"debug"/"full.txt",library.full_prior);write(out/"debug"/"routed.tsv",library.routing.routed_tsv);write(out/"debug"/"route-audit.tsv",library.routing.audit_tsv);}
        return 0;
    }catch(const std::exception& e){std::cerr<<"fqdup merge: "<<e.what()<<'\n';return 1;}
}
