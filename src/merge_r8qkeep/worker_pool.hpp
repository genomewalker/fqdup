#pragma once
#include "worker_placement.hpp"
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>
namespace fqdup::r8qkeep {
// Callers bound outstanding futures. Nested waits require reserved capacity:
// counts fitting reserves at most two coordinators and uses the remaining
// workers for their jobs; p1/p2 coordinators submit no nested work.
class WorkerPool {
    std::mutex mutex_;
    std::condition_variable ready_;
    std::deque<std::function<void()>> queue_;
    std::vector<std::thread> workers_;
    bool stopping_=false,placement_reported_=false;
    std::vector<int> placement_;
    std::atomic<int> placed_{0};
public:
    explicit WorkerPool(int count):placement_(worker_cpu_plan(count)) {
        try {
            for(int i=0;i<count;++i)workers_.emplace_back([this,i]{
                if(!placement_.empty()&&place_worker(placement_[i]))placed_.fetch_add(1,std::memory_order_relaxed);
                for(;;){
                    std::function<void()> job;
                    {std::unique_lock<std::mutex> lock(mutex_);
                     ready_.wait(lock,[this]{return stopping_||!queue_.empty();});
                     if(queue_.empty())return;
                     job=std::move(queue_.front());queue_.pop_front();}
                    job();
                }
            });
        }catch(...){stop();throw;}
    }
    WorkerPool(const WorkerPool&)=delete;
    WorkerPool& operator=(const WorkerPool&)=delete;
    ~WorkerPool(){stop();}
    void stop(){
        {std::lock_guard<std::mutex> lock(mutex_);stopping_=true;}
        ready_.notify_all();
        for(auto& worker:workers_)if(worker.joinable())worker.join();
        if(!placement_reported_&&!placement_.empty()&&std::getenv("FQDUP_FAMILY_DIAGNOSTICS")){
            std::cerr<<"WORKER_PLACEMENT workers="<<workers_.size()<<" pinned="<<placed_.load()<<" cpus=";
            for(auto cpu:placement_)std::cerr<<cpu<<',';std::cerr<<'\n';
        }
        placement_reported_=true;
    }
    template<class Work>
    auto submit(Work work)->std::future<decltype(work())>{
        using Result=decltype(work());
        auto task=std::make_shared<std::packaged_task<Result()>>(std::move(work));
        auto future=task->get_future();
        {std::lock_guard<std::mutex> lock(mutex_);
         if(stopping_)throw std::runtime_error("submit to stopped worker pool");
         queue_.emplace_back([task]{(*task)();});}
        ready_.notify_one();return future;
    }
};
}
