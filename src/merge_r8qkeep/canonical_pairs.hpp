#pragma once
#include "family.hpp"
#include "worker_pool.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace fqdup::r8qkeep {
// Counts-mode prerequisite. IDs must be normalized by read_family_fastq.
// Validate before sorting so that sorting cannot conceal mispaired inputs.
inline void canonicalize_family_input(std::vector<FamilyRead>& r1,
                                      std::vector<FamilyRead>& r2,
                                      WorkerPool* pool=nullptr,int threads=1) {
    if (r1.size() != r2.size()) throw std::runtime_error("pair count");
    for (std::size_t i = 0; i < r1.size(); ++i)
        if (r1[i].id != r2[i].id) throw std::runtime_error("pair order " + r1[i].id);

    std::vector<std::size_t> order(r1.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    const auto key = [&](std::size_t i) {
        return std::tie(r1[i].id, r1[i].h, r1[i].s, r1[i].q,
                        r2[i].h, r2[i].s, r2[i].q);
    };
    const auto less=[&](std::size_t a,std::size_t b){return key(a)<key(b);};
    const int workers=std::max(1,std::min(threads,8));
    if(!pool||workers==1||order.size()<4096)std::sort(order.begin(),order.end(),less);
    else{
        // Sort only indices. Equal keys identify byte-identical mate records,
        // so a different stable ordering within a tie is unobservable.
        std::vector<std::size_t> scratch(order.size());
        const auto bound=[&](int part){return order.size()*part/workers;};
        const auto parallel=[&](int count,auto work){
            std::vector<std::future<void>> pending;pending.reserve(count);
            try{
                for(int task=0;task<count;++task)pending.push_back(pool->submit([&,task]{work(task);}));
                for(auto& result:pending)result.get();
            }catch(...){pool->stop();throw;}
        };
        parallel(workers,[&](int part){std::sort(order.begin()+bound(part),order.begin()+bound(part+1),less);});
        for(int width=1;width<workers;width*=2){
            parallel((workers+2*width-1)/(2*width),[&](int part){
                const auto begin=bound(part*2*width),middle=bound(std::min((part*2+1)*width,workers)),
                    end=bound(std::min((part+1)*2*width,workers));
                std::merge(order.begin()+begin,order.begin()+middle,order.begin()+middle,
                           order.begin()+end,scratch.begin()+begin,less);
            });
            order.swap(scratch);
        }
    }
    // order[destination] = source. Move cycles with two saved reads, retaining
    // only one size_t per pair instead of duplicating the sequence storage.
    for (std::size_t start = 0; start < order.size(); ++start) {
        if (order[start] == start) continue;
        auto saved1 = std::move(r1[start]);
        auto saved2 = std::move(r2[start]);
        auto at = start;
        while (order[at] != start) {
            const auto from = order[at];
            r1[at] = std::move(r1[from]);
            r2[at] = std::move(r2[from]);
            order[at] = at;
            at = from;
        }
        r1[at] = std::move(saved1);
        r2[at] = std::move(saved2);
        order[at] = at;
    }
}
}
