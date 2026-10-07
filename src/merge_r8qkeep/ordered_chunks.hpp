#pragma once
#include <algorithm>
#include <future>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
// Bound outstanding row chunks and consume results in original input order.
template<class Work,class Consume>
void ordered_text_chunks(const std::string& text,int threads,std::size_t row_offset,Work work,Consume consume) {
    const int workers=std::max(1,std::min(threads,8));
    const std::size_t chunk_rows=std::max<std::size_t>(1,1024/workers);
    using Result=decltype(work(std::string(),row_offset));
    std::size_t pos=0,row=row_offset;
    while(pos<text.size()) {
        std::vector<std::future<Result>> pending;
        for(int w=0;w<workers&&pos<text.size();++w) {
            const auto begin=pos,start=row;
            for(std::size_t n=0;n<chunk_rows&&pos<text.size();++n,++row) {
                const auto end=text.find('\n',pos);pos=end==std::string::npos?text.size():end+1;
            }
            pending.push_back(std::async(std::launch::async,[&,chunk=text.substr(begin,pos-begin),start](){return work(chunk,start);}));
        }
        for(auto& result:pending)consume(result.get());
    }
}
}
