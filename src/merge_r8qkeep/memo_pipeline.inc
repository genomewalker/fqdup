// Included inside stream_pipeline after the standard process callback. The
// cache context is unique for this immutable library/prior/settings/cutoff.
static std::atomic<std::uint64_t> memo_context_counter{0};
const auto memo_context=++memo_context_counter;
// One fixed 512-row chunk per worker; fewer than 1% hits disables probing.
constexpr std::size_t memo_pilot_rows=512;
std::atomic<std::uint64_t> memo_hits{0},memo_rows{0},memo_bypassed{0};
const auto memo_process=[&](std::size_t start,std::size_t n){
    constexpr std::size_t batch=512,table_size=1024;
    struct Scratch {
        ContentMemo cache;
        std::uint64_t context=0,seen=0,hits=0;bool bypass=false;
        std::vector<FamilyPair> raw=std::vector<FamilyPair>(batch),misses=std::vector<FamilyPair>(batch);
        std::vector<std::string> keys=std::vector<std::string>(batch);
        std::vector<MemoValue> values=std::vector<MemoValue>(batch),computed;
        std::array<std::size_t,batch> hash{},which{},representative{};
        std::array<int,table_size> pending{};
    };
    thread_local Scratch scratch;
    scratch.cache.context(memo_context);scratch.pending.fill(-1);
    if(scratch.context!=memo_context){scratch.context=memo_context;scratch.seen=scratch.hits=0;scratch.bypass=false;}
    const auto account=[&](const NativeResult& native,const OutputChunk& out){
        if(!memory)return;
        MemoryUsage raw;raw.pairs(scratch.raw);raw.pairs(scratch.misses);
        auto bytes=raw.raw+raw.ids+native_bytes(native)+strings_bytes(scratch.keys)+vector_bytes(scratch.values)+vector_bytes(scratch.computed)+sizeof(scratch.hash)+sizeof(scratch.which)+sizeof(scratch.representative)+sizeof(scratch.pending);
        for(const auto* rows:{&scratch.values,&scratch.computed})for(const auto& value:*rows)bytes+=string_heap(value.record.outcome)+stage_read_bytes(value.record.first)+stage_read_bytes(value.record.second);
        maximum(max_native,bytes);
        bytes=vector_bytes(out.outputs)+vector_bytes(out.seconds)+vector_bytes(out.native_clean_pass)+string_heap(out.routed)+string_heap(out.audit)+string_heap(out.native);
        for(const auto& row:out.outputs)bytes+=string_heap(row.merged)+string_heap(row.r1)+string_heap(row.r2)+string_heap(row.orphan)+string_heap(row.outcomes);
        maximum(max_output,bytes);
    };
    CompactPairStore::Decoder decoder;std::size_t unique=0,hits=0;
    const auto native_begin=Clock::now();
    if(scratch.bypass){
        for(std::size_t i=0;i<n;++i)lib.compact->read(start+i,scratch.raw[i],decoder);
        auto native=merge_native_range(scratch.raw,0,n,lib.native_settings,false,false);
        auto out=process(native.records,start,elapsed(native_begin),nullptr,&native.alignments);account(native,out);
        memo_rows.fetch_add(n,std::memory_order_relaxed);memo_bypassed.fetch_add(n,std::memory_order_relaxed);return out;
    };
    for(std::size_t i=0;i<n;++i){
        lib.compact->read(start+i,scratch.raw[i],decoder);
        memo_content(scratch.raw[i],scratch.keys[i]);const auto h=ContentMemo::hash(scratch.keys[i]);scratch.hash[i]=h;
        if(const auto* cached=scratch.cache.find(scratch.keys[i],h)){
            scratch.values[i]=*cached;scratch.which[i]=batch;++hits;continue;
        }
        auto bucket=h%table_size;
        while(scratch.pending[bucket]>=0&&scratch.keys[scratch.pending[bucket]]!=scratch.keys[i])bucket=(bucket+1)%table_size;
        if(scratch.pending[bucket]>=0){scratch.which[i]=scratch.which[scratch.pending[bucket]];++hits;continue;}
        scratch.pending[bucket]=int(i);scratch.which[i]=unique;scratch.representative[unique]=i;
        scratch.misses[unique]=scratch.raw[i];memo_headers(scratch.misses[unique]);++unique;
    }
    NativeResult native;
    if(unique)native=merge_native_range(scratch.misses,0,unique,lib.native_settings,false,false);
    const auto native_seconds=elapsed(native_begin);
    OutputChunk out;
    if(unique){
        out=process(native.records,0,native_seconds,&scratch.computed,&native.alignments);
        for(std::size_t j=0;j<unique;++j){const auto i=scratch.representative[j];scratch.cache.put(scratch.keys[i],scratch.hash[i],scratch.computed[j]);}
    }else{out.outputs.resize(1);out.seconds.resize(1);out.native_seconds=native_seconds;}
    out.count=n;out.native_clean=0;out.native_clean_pass.assign(1,0);
    const auto emit_begin=Clock::now();
    for(std::size_t i=0;i<n;++i){
        auto& value=scratch.values[i];if(scratch.which[i]!=batch)value=scratch.computed[scratch.which[i]];
        out.native_clean+=value.native_clean;out.native_clean_pass[0]+=value.clean_pass;
        memo_restore_header(value.record.first.header,scratch.raw[i]);memo_restore_header(value.record.second.header,scratch.raw[i]);
        append_record(out.outputs[0],value.record,scratch.raw[i]);
    }
    out.seconds[0]+=elapsed(emit_begin);memo_hits.fetch_add(hits,std::memory_order_relaxed);memo_rows.fetch_add(n,std::memory_order_relaxed);
    account(native,out);scratch.seen+=n;scratch.hits+=hits;
    if(!force_memo&&scratch.seen==memo_pilot_rows&&scratch.hits*100<scratch.seen){scratch.bypass=true;scratch.cache.release();}
    return out;
};
