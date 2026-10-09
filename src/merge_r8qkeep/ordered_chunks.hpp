#pragma once
#include <algorithm>
#include "worker_pool.hpp"
#include <future>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <ctime>
namespace fqdup::r8qkeep {
// Index ranges borrow records directly; jobs finish out of order, consumption
// remains ordered. At most workers futures (and outputs) are retained.
template<class Work,class Consume>
void ordered_index_chunks(WorkerPool& pool,std::size_t count,int threads,std::size_t row_offset,Work work,Consume consume,const char* label="records",std::size_t rows_per_chunk=0) {
    const int workers=std::max(1,std::min(threads,8));
    const std::size_t chunk_rows=rows_per_chunk?rows_per_chunk:std::max<std::size_t>(1,1024/workers);
    using Result=decltype(work(row_offset,std::size_t{}));
    const bool profile=std::getenv("FQDUP_CHUNK_PROFILE")!=nullptr;
    const auto start=std::chrono::steady_clock::now();std::size_t chunks=0,pos=0;
    std::deque<std::future<Result>> pending;
    const auto launch=[&]{
        const auto begin=pos,n=std::min(chunk_rows,count-pos);pos+=n;++chunks;
        pending.push_back(pool.submit([work,begin,n,row_offset]{return work(row_offset+begin,n);}));
    };
    try {
        while(pos<count&&pending.size()<std::size_t(workers))launch();
        while(!pending.empty()){
            auto result=pending.front().get();pending.pop_front();
            if(pos<count)launch();
            consume(std::move(result));
        }
    }catch(...){
        // All references captured by jobs remain alive until they have drained.
        // Also drain a submitted job if storing its future failed to allocate.
        pool.stop();
        throw;
    }
    if(profile)std::fprintf(stderr,"CHUNK_PROFILE stage=%s workers=%d chunk_rows=%zu max_inflight=%d rows=%zu chunks=%zu wall=%.6f\n",label,workers,chunk_rows,workers,count,chunks,std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
}
template<class Work,class Consume>
void ordered_index_chunks(std::size_t count,int threads,std::size_t row_offset,Work work,Consume consume,const char* label="records",std::size_t rows_per_chunk=0) {
    WorkerPool pool(std::max(1,std::min(threads,8)));
    ordered_index_chunks(pool,count,threads,row_offset,std::move(work),std::move(consume),label,rows_per_chunk);
}

// Bound outstanding row chunks and consume results in original input order.
template<class Work,class Consume>
void ordered_text_chunks(const std::string& text,int threads,std::size_t row_offset,Work work,Consume consume,const char* label="chunks") {
    const int workers=std::max(1,std::min(threads,8));
    const std::size_t chunk_rows=std::max<std::size_t>(1,1024/workers);
    using Result=decltype(work(std::string(),row_offset));
    struct TimedResult { Result value;double wall,cpu; };
    const bool profile=std::getenv("FQDUP_CHUNK_PROFILE")!=nullptr;
    using Clock=std::chrono::steady_clock;
    const auto begin_time=Clock::now();
    const auto seconds=[](auto begin){return std::chrono::duration<double>(Clock::now()-begin).count();};
    const auto cpu=[](){timespec t{};clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);return double(t.tv_sec)+double(t.tv_nsec)*1e-9;};
    double busy=0,worker_cpu=0,wait=0,consume_time=0,launch_time=0;std::size_t chunks=0;
    std::size_t pos=0,row=row_offset;
    while(pos<text.size()) {
        std::vector<std::future<TimedResult>> pending;
        const auto launch_begin=Clock::now();
        for(int w=0;w<workers&&pos<text.size();++w) {
            const auto begin=pos,start=row;
            for(std::size_t n=0;n<chunk_rows&&pos<text.size();++n,++row) {
                const auto end=text.find('\n',pos);pos=end==std::string::npos?text.size():end+1;
            }
            pending.push_back(std::async(std::launch::async,[&,chunk=text.substr(begin,pos-begin),start](){const auto t=Clock::now();const auto c=profile?cpu():0;auto value=work(chunk,start);return TimedResult{std::move(value),profile?seconds(t):0,profile?cpu()-c:0};}));++chunks;
        }
        if(profile)launch_time+=seconds(launch_begin);
        for(auto& result:pending){const auto t=Clock::now();auto part=result.get();if(profile){wait+=seconds(t);busy+=part.wall;worker_cpu+=part.cpu;}const auto c=Clock::now();consume(std::move(part.value));if(profile)consume_time+=seconds(c);}
    }
    if(profile)std::fprintf(stderr,"CHUNK_PROFILE stage=%s workers=%d chunk_rows=%zu max_inflight=%d rows=%zu chunks=%zu wall=%.6f worker_wall=%.6f worker_cpu=%.6f wait=%.6f consume=%.6f launch=%.6f\n",label,workers,chunk_rows,workers,row-row_offset,chunks,seconds(begin_time),busy,worker_cpu,wait,consume_time,launch_time);
}
}
