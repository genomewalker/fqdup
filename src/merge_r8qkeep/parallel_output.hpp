#pragma once
#include "streaming_output.hpp"
#include "ordered_writer.hpp"
#include <array>
#include <chrono>
#include <numeric>
namespace fqdup::r8qkeep {
// Files are independent; each actor preserves its own byte order. The caller
// owns streams and must retain them until this object has joined every actor.
class ParallelOutputFiles {
    std::array<StreamingOutput,5>& files_;
    bool measure_,finished_=false;
    std::size_t limit_;
    std::array<double,5> seconds_{};
    std::array<std::size_t,5> max_text_{};
    std::array<std::unique_ptr<OrderedWriter<std::string>>,5> writers_;
public:
    ParallelOutputFiles(std::array<StreamingOutput,5>& files,bool measure=false,std::size_t limit=16):files_(files),measure_(measure),limit_(limit){
        for(std::size_t i=0;i<writers_.size();++i)writers_[i]=std::make_unique<OrderedWriter<std::string>>(limit_,[this,i](std::string text){
            const auto start=measure_?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
            files_[i].append(text);
            if(measure_)seconds_[i]+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        });
    }
    ParallelOutputFiles(const ParallelOutputFiles&)=delete;
    ParallelOutputFiles& operator=(const ParallelOutputFiles&)=delete;
    void append(std::array<std::string*,5> text){
        if(finished_)throw std::logic_error("finished parallel outputs");
        for(std::size_t i=0;i<text.size();++i){
            if(!text[i])throw std::invalid_argument("null output text");
            if(text[i]->empty())continue;
            max_text_[i]=std::max(max_text_[i],text[i]->capacity());
            std::promise<std::string> ready;auto future=ready.get_future();ready.set_value(std::move(*text[i]));writers_[i]->submit(std::move(future));
        }
    }
    void finish(){
        std::exception_ptr error;
        for(auto& writer:writers_)try{writer->finish();}catch(...){if(!error)error=std::current_exception();}
        finished_=true;if(error)std::rethrow_exception(error);
    }
    double seconds()const{if(!finished_)throw std::logic_error("unfinished writer timing");return std::accumulate(seconds_.begin(),seconds_.end(),0.0);}
    std::size_t text_capacity_bound()const{return (limit_+1)*std::accumulate(max_text_.begin(),max_text_.end(),std::size_t{0});}
};
}
