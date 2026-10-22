// Included in cli.cpp's anonymous namespace. Stream raw input again instead
// of retaining a second complete library through family correction.
std::string raw_id(const std::string& header){
    auto id=header.substr(1,header.find_first_of(" \t")-1);
    if(id.size()>2&&id[id.size()-2]=='/'&&(id.back()=='1'||id.back()=='2'))id.resize(id.size()-2);
    return id;
}
void calibrate_raw_inputs(LibraryPipeline& library,const std::vector<std::string>& first,const std::vector<std::string>& second,const std::filesystem::path& out,int threads,bool counts){
    std::unordered_set<std::string> excluded;
    CompactPairStore::Decoder decoder;std::string id;
    for(int index:library.family.calibration_excluded){
        if(library.compact)library.compact->id(index,id,decoder);else id=library.family.pairs.at(index).id;
        excluded.insert(id);
    }
    const auto count=library.pair_count();
    const auto threshold=count?(static_cast<unsigned __int128>(std::min<std::size_t>(25000,count))<<64)/count:0;
    std::filesystem::create_directories(out/"debug");
    std::ofstream nominal(out/"debug"/"nominal.tsv",std::ios::binary);bool header=false;
    auto a=make_chained_fastq_reader(first,1),b=make_chained_fastq_reader(second,1);FastqRecord x,y;
    std::vector<FamilyPair> batch;batch.reserve(8192);
    const auto flush=[&](){
        if(batch.empty())return;
        const auto native=merge_native(batch,library.native_settings,false,false,threads,library.pool.get());
        const auto data=raw_concordance(batch,native);library.concordance.add(data.bins);
        nominal<<(header?data.nominal_tsv.substr(data.nominal_tsv.find('\n')+1):data.nominal_tsv);header=true;batch.clear();
    };
    for(;;){
        const bool one=a->read(x),two=b->read(y);if(one!=two)throw std::runtime_error("calibration raw pair count");if(!one)break;
        id=raw_id(x.header);if(raw_id(y.header)!=id)throw std::runtime_error("calibration raw pair order");
        if(excluded.count(id)||!counts||static_cast<unsigned __int128>(XXH64(id.data(),id.size(),prior_sample_seed))<threshold)continue;
        batch.push_back(FamilyPair{id,{FamilyRead{x.header,id,x.seq,x.qual},FamilyRead{y.header,id,y.seq,y.qual}}});
        if(batch.size()==8192)flush();
    }
    flush();if(!header)nominal<<"id\tfold\ti\tk\tq1_bin\tq2_bin\te1\te2\tbase1\trc_base2\tcertified\n";
    nominal.close();if(!nominal)throw std::runtime_error("write heldout nominal");
}
