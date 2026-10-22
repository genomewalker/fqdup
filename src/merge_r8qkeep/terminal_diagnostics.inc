// Raw-base terminal census, sharing one input pass across the complete grid.
// Input sequence is read again so family edits cannot change the opportunity bases.
struct TerminalProducts {std::string bins,summary;};
TerminalProducts output_geometry_diagnostics(const std::vector<std::string>& first,const std::vector<std::string>& second,const std::vector<std::string>& dirs,const std::vector<std::string>& labels,const AdapterVocabulary& vocabulary,std::vector<OutputRisk>& risks){
    struct Span {std::uint32_t begin=0,end=0,length=0;bool present=false;};
    struct Cell {std::uint64_t total=0,trimmed=0;};
    const auto stride=dirs.size()*2;
    std::unordered_map<std::string,std::size_t> ids;std::vector<Span> spans;
    std::vector<std::array<std::array<std::array<Cell,256>,2>,2>> bins(dirs.size());
    const auto field=[](const std::string& h,const char* name){const auto tag=h.rfind(" FQT=v1;");if(tag==std::string::npos)throw std::runtime_error("terminal missing FQT");const auto at=h.find(std::string(";")+name+"=",tag);if(at==std::string::npos)throw std::runtime_error("terminal FQT field");return std::stoul(h.substr(at+std::char_traits<char>::length(name)+2));};
    for(std::size_t cut=0;cut<dirs.size();++cut)for(int channel=0;channel<3;++channel){
        const char* names[]={"r1.fq","r2.fq","orphan.fq"};std::ifstream in(std::filesystem::path(dirs[cut])/names[channel]);if(!in)throw std::runtime_error("terminal output open");
        std::string h,s,plus,q;
        while(std::getline(in,h)){
            if(!std::getline(in,s)||!std::getline(in,plus)||!std::getline(in,q))throw std::runtime_error("terminal truncated FASTQ");
            const int mate=channel<2?channel+1:int(field(h,"m"));if(mate<1||mate>2)throw std::runtime_error("terminal mate");
            const auto id=raw_id(h);auto [it,inserted]=ids.emplace(id,ids.size());if(inserted)spans.resize(spans.size()+stride);
            auto& span=spans[it->second*stride+cut*2+mate-1];if(span.present)throw std::runtime_error("terminal duplicate mate");
            span={std::uint32_t(field(h,"s")),std::uint32_t(field(h,"e")),std::uint32_t(field(h,"olen")),true};
            if(span.end<span.begin||span.end-span.begin!=s.size())throw std::runtime_error("terminal span");
            const auto residual=vocabulary.residual(s,mate);risks.at(cut).adapter_bases+=residual;risks.at(cut).adapter_mates+=residual>0;
        }
    }
    for(int mate=0;mate<2;++mate){
        auto reader=make_chained_fastq_reader(mate==0?first:second,1);FastqRecord r;
        while(reader->read(r)){
            const auto it=ids.find(raw_id(r.header));if(it==ids.end())continue;
            for(std::size_t cut=0;cut<dirs.size();++cut){const auto& span=spans[it->second*stride+cut*2+mate];if(!span.present)continue;
                if(span.length>r.seq.size())throw std::runtime_error("terminal raw length");
                for(int window=0;window<2;++window){
                    const auto start=window==0?0u:(span.length>10?span.length-10:0u),end=window==0?std::min(10u,span.length):span.length;
                    for(auto p=start;p<end;++p){auto& cell=bins[cut][mate][window][static_cast<unsigned char>(r.seq[p])];++cell.total;cell.trimmed+=window==0?p<span.begin:p>=span.end;}
                }
            }
        }
    }
    std::ostringstream out,summary;out.precision(17);summary.precision(17);
    out<<"cutoff\tmate\twindow\tbase\ttotal\ttrimmed\tfraction\n";
    summary<<"cutoff\twindow\ttarget\ttarget_total\ttarget_trimmed\ttarget_fraction\tother_total\tother_trimmed\tother_fraction\n";
    for(std::size_t cut=0;cut<dirs.size();++cut)for(int window=0;window<2;++window){
        const char target=window==0?'C':'G';Cell wanted,other;
        for(int mate=0;mate<2;++mate)for(int base=0;base<256;++base){const auto& c=bins[cut][mate][window][base];if(!c.total)continue;
            out<<labels[cut]<<'\t'<<mate+1<<'\t'<<(window==0?"first10":"last10")<<'\t'<<char(base)<<'\t'<<c.total<<'\t'<<c.trimmed<<'\t'<<double(c.trimmed)/c.total<<'\n';
            auto& aggregate=base==target?wanted:other;aggregate.total+=c.total;aggregate.trimmed+=c.trimmed;
        }
        summary<<labels[cut]<<'\t'<<(window==0?"first10":"last10")<<'\t'<<target<<'\t'<<wanted.total<<'\t'<<wanted.trimmed<<'\t'<<(wanted.total?double(wanted.trimmed)/wanted.total:0)<<'\t'<<other.total<<'\t'<<other.trimmed<<'\t'<<(other.total?double(other.trimmed)/other.total:0)<<'\n';
    }
    return {out.str(),summary.str()};
}
