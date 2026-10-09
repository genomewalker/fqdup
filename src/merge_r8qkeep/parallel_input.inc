// Included in cli.cpp's anonymous namespace. Each mate reader is accessed by
// exactly one task at a time. Raw arenas and encoded futures are bounded.
struct CompactRawBatch {
    struct Row {std::size_t offset,header,sequence,quality;};
    std::string bytes;
    std::vector<Row> rows;
    void clear(){bytes.clear();rows.clear();}
    void append(const FastqRecord& r){rows.push_back({bytes.size(),r.header.size(),r.seq.size(),r.qual.size()});bytes+=r.header;bytes+=r.seq;bytes+=r.qual;}
    std::array<std::string_view,3> read(std::size_t i)const{const auto&r=rows[i];const auto*p=bytes.data()+r.offset;return {std::string_view(p,r.header),std::string_view(p+r.header,r.sequence),std::string_view(p+r.header+r.sequence,r.quality)};}
};
class CompactRawPool:public std::enable_shared_from_this<CompactRawPool>{
    std::mutex mutex_;
    std::vector<std::unique_ptr<CompactRawBatch>> free_;
public:
    std::shared_ptr<CompactRawBatch> acquire(){
        std::unique_ptr<CompactRawBatch> value;
        {std::lock_guard<std::mutex> lock(mutex_);if(!free_.empty()){value=std::move(free_.back());free_.pop_back();}}
        if(!value)value=std::make_unique<CompactRawBatch>();value->clear();
        return {value.release(),[self=shared_from_this()](CompactRawBatch* raw)noexcept{
            std::unique_ptr<CompactRawBatch> owned(raw);
            try{std::lock_guard<std::mutex> lock(self->mutex_);self->free_.push_back(std::move(owned));}catch(...){}
        }};
    }
};
std::shared_ptr<CompactPairStore> parallel_compact_input(const std::vector<std::string>& first,const std::vector<std::string>& second,int threads){
    constexpr std::size_t batch=4096;
    const bool diagnostics=std::getenv("FQDUP_FAMILY_DIAGNOSTICS")!=nullptr;
    const auto cpu_ns=[&]()->std::uint64_t {if(!diagnostics)return 0;struct timespec value{};if(::clock_gettime(CLOCK_THREAD_CPUTIME_ID,&value))throw std::runtime_error("thread CPU clock");return std::uint64_t(value.tv_sec)*1000000000ULL+value.tv_nsec;};
    std::atomic<std::uint64_t> parse_cpu{0},encode_cpu{0};std::uint64_t join_cpu=0;
    struct Reader {std::unique_ptr<FastqReaderBase> input;FastqRecord record;std::unique_ptr<BufferedPlainInput> plain;};
    const auto source=[](const std::vector<std::string>& paths){Reader out;
        if(paths.size()==1)out.plain=BufferedPlainInput::open(paths[0]);
        if(!out.plain)out.input=make_chained_fastq_reader(paths,1);return out;};
    std::array<Reader,2> reader{{source(first),source(second)}};
    auto raw_pool=std::make_shared<CompactRawPool>();
    WorkerPool readers(2),encoders(std::max(1,std::min(threads,8)-2));
    auto result=std::make_shared<CompactPairStore>();
    std::deque<std::future<std::unique_ptr<CompactPairStore>>> pending;
    const auto launch=[&]{
        std::array<std::future<std::shared_ptr<CompactRawBatch>>,2> next;
        for(unsigned side=0;side<2;++side)next[side]=readers.submit([&,side]{
            const auto started=cpu_ns();auto raw=raw_pool->acquire();auto& source=reader[side];
            if(source.plain)source.plain->read(*raw,batch);
            else while(raw->rows.size()<batch&&source.input->read(source.record))raw->append(source.record);
            parse_cpu.fetch_add(cpu_ns()-started,std::memory_order_relaxed);return raw;
        });return next;
    };
    const auto consume=[&]{auto part=pending.front().get();pending.pop_front();const auto started=cpu_ns();result->append_store(*part);join_cpu+=cpu_ns()-started;};
    try{
        auto next=launch();
        for(;;){
            auto a=next[0].get(),b=next[1].get();
            if(a->rows.size()!=b->rows.size())throw std::runtime_error("pair count");
            if(a->rows.empty())break;
            next=launch();
            pending.push_back(encoders.submit([&,a=std::move(a),b=std::move(b)]{
                const auto started=cpu_ns();
                auto store=std::make_unique<CompactPairStore>();
                for(std::size_t i=0;i<a->rows.size();++i){const auto x=a->read(i),y=b->read(i);store->append(x[0],x[1],x[2],y[0],y[1],y[2]);}
                encode_cpu.fetch_add(cpu_ns()-started,std::memory_order_relaxed);return store;
            }));
            if(pending.size()>=std::size_t(std::max(2,std::min(threads,8))))consume();
        }
        while(!pending.empty())consume();
    }catch(...){readers.stop();encoders.stop();throw;}
    readers.stop();encoders.stop();
    if(diagnostics)std::cerr<<"COMPACT_INPUT pairs="<<result->size()<<" parse_cpu_seconds="<<parse_cpu.load()/1e9<<" encode_cpu_seconds="<<encode_cpu.load()/1e9<<" join_cpu_seconds="<<join_cpu/1e9<<'\n';
    return result;
}
