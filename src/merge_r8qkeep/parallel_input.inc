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
    struct Reader {std::unique_ptr<FastqReaderBase> input;FastqRecord record;};
    std::array<Reader,2> reader{{{make_chained_fastq_reader(first,1),{}},{make_chained_fastq_reader(second,1),{}}}};
    auto raw_pool=std::make_shared<CompactRawPool>();
    WorkerPool readers(2),encoders(std::max(1,std::min(threads,8)-2));
    auto result=std::make_shared<CompactPairStore>();
    std::deque<std::future<std::unique_ptr<CompactPairStore>>> pending;
    const auto launch=[&]{
        std::array<std::future<std::shared_ptr<CompactRawBatch>>,2> next;
        for(unsigned side=0;side<2;++side)next[side]=readers.submit([&,side]{
            auto raw=raw_pool->acquire();auto& source=reader[side];
            while(raw->rows.size()<batch&&source.input->read(source.record))raw->append(source.record);
            return raw;
        });return next;
    };
    const auto consume=[&]{auto part=pending.front().get();pending.pop_front();result->append_store(*part);};
    try{
        auto next=launch();
        for(;;){
            auto a=next[0].get(),b=next[1].get();
            if(a->rows.size()!=b->rows.size())throw std::runtime_error("pair count");
            if(a->rows.empty())break;
            next=launch();
            pending.push_back(encoders.submit([a=std::move(a),b=std::move(b)]{
                auto store=std::make_unique<CompactPairStore>();
                for(std::size_t i=0;i<a->rows.size();++i){const auto x=a->read(i),y=b->read(i);store->append(x[0],x[1],x[2],y[0],y[1],y[2]);}
                return store;
            }));
            if(pending.size()>=std::size_t(std::max(2,std::min(threads,8))))consume();
        }
        while(!pending.empty())consume();
    }catch(...){readers.stop();encoders.stop();throw;}
    readers.stop();encoders.stop();return result;
}
