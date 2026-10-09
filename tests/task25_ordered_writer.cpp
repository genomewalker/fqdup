#include "../src/merge_r8qkeep/ordered_writer.hpp"
#include <iostream>
using namespace fqdup::r8qkeep;
static void require(bool value){if(!value)throw std::runtime_error("ordered writer failure");}
int main(){
    // The later result is ready first; the writer must still consume 0,1.
    std::vector<int> result;
    {OrderedWriter<int> writer(2,[&](int value){result.push_back(value);});
     std::promise<int> first,second;writer.submit(first.get_future());writer.submit(second.get_future());
     second.set_value(1);first.set_value(0);writer.finish();}
    require(result==std::vector<int>({0,1}));
    for(bool consumer_error:{false,true}){
        bool caught=false;OrderedWriter<int> writer(1,[&](int){if(consumer_error)throw std::runtime_error("write error");});
        std::promise<int> value;writer.submit(value.get_future());
        if(consumer_error)value.set_value(1);else value.set_exception(std::make_exception_ptr(std::runtime_error("work error")));
        try{writer.finish();}catch(const std::runtime_error&){caught=true;}require(caught);
    }
    {OrderedWriter<int> empty(1,[](int){throw std::runtime_error("unexpected output");});empty.finish();}
    // Producer failure/destruction drains the outstanding write.
    result.clear();
    {OrderedWriter<int> writer(1,[&](int value){result.push_back(value);});
     std::promise<int> value;value.set_value(7);writer.submit(value.get_future());}
    require(result==std::vector<int>({7}));
    // Fit coordinators may wait only when another pool worker can run jobs.
    for(int threads:{1,2,3,8}){
        WorkerPool pool(threads);std::vector<std::future<int>> fits;
        for(int f=0;f<std::min(2,threads);++f)fits.push_back(pool.submit([&,f]{
            if(threads<=2)return f;
            std::vector<std::future<int>> parts;
            for(int j=0;j<threads-2;++j)parts.push_back(pool.submit([j]{return j;}));
            int sum=f;for(auto& part:parts)sum+=part.get();return sum;
        }));
        for(int f=0;f<int(fits.size());++f)require(fits[f].get()==f+(threads>2?(threads-3)*(threads-2)/2:0));
    }
    std::cout<<"TASK25_ORDERED_WRITER_PASS\n";
}
