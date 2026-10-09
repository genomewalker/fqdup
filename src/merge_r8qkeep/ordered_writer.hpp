#pragma once
#include "worker_pool.hpp"
namespace fqdup::r8qkeep {
// One writer owns consumption. Credits include the active write as well as
// queued futures, bounding retained records even behind a slow first chunk.
template<class Result> class OrderedWriter {
    std::mutex mutex_;
    std::condition_variable ready_;
    std::deque<std::future<Result>> pending_;
    std::function<void(Result)> consume_;
    std::size_t limit_;
    bool active_=false,done_=false;
    std::exception_ptr error_;
    std::thread thread_;
    void run() noexcept {
        try {
            for(;;){
                std::future<Result> next;
                {std::unique_lock<std::mutex> lock(mutex_);
                 ready_.wait(lock,[&]{return done_||!pending_.empty();});
                 if(pending_.empty())return;
                 next=std::move(pending_.front());pending_.pop_front();active_=true;}
                consume_(next.get());
                {std::lock_guard<std::mutex> lock(mutex_);active_=false;}
                ready_.notify_all();
            }
        }catch(...){
            {std::lock_guard<std::mutex> lock(mutex_);error_=std::current_exception();done_=true;active_=false;}
            ready_.notify_all();
        }
    }
    void join() noexcept {
        {std::lock_guard<std::mutex> lock(mutex_);done_=true;}
        ready_.notify_all();
        if(thread_.joinable())thread_.join();
    }
public:
    OrderedWriter(std::size_t limit,std::function<void(Result)> consume):consume_(std::move(consume)),limit_(limit){
        if(!limit_)throw std::invalid_argument("zero writer capacity");
        thread_=std::thread([this]{run();});
    }
    ~OrderedWriter(){join();}
    void submit(std::future<Result> next){
        std::unique_lock<std::mutex> lock(mutex_);
        ready_.wait(lock,[&]{return done_||pending_.size()+std::size_t(active_)<limit_;});
        if(error_)std::rethrow_exception(error_);
        if(done_)throw std::runtime_error("submit to closed writer");
        pending_.push_back(std::move(next));lock.unlock();ready_.notify_all();
    }
    void finish(){join();if(error_)std::rethrow_exception(error_);}
};
}
