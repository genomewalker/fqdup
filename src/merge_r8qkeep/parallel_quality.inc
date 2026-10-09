// Included in CompactPairStore's private section. A wave writes only above all
// still-unread source bits. Barriers between waves make in-place widening safe.
void widen_qualities(unsigned old_width,WorkerPool* pool,int workers){
    const auto convert=[&](std::uint64_t begin,std::uint64_t end){
        for(auto i=end;i>begin;){--i;const auto old=get(qualities_,i,old_width);const auto c=old_width==8?old:alphabet_[old];put(qualities_,i,width_,width_==8?c:codes_[c]);}
    };
    auto remaining=starts_.back();
    if(pool&&workers>1){
        while(remaining>65536){
            const auto begin=((remaining*old_width+width_-1)/width_+7)/8*8;
            if(begin>=remaining)break;
            const auto stride=((remaining-begin+workers-1)/workers+7)/8*8;
            std::vector<std::future<void>> jobs;jobs.reserve(workers);
            try{for(auto lo=begin;lo<remaining;lo+=stride){const auto hi=std::min(remaining,lo+stride);jobs.push_back(pool->submit([&,lo,hi]{convert(lo,hi);}));}
                for(auto& job:jobs)job.get();
            }catch(...){for(auto& job:jobs)if(job.valid())job.wait();throw;}
            remaining=begin;
        }
    }
    convert(0,remaining);
}
