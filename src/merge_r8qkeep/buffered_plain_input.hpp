#pragma once
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
namespace fqdup::r8qkeep {
// Optional regular-file reader for the compact arena pipeline. Compressed,
// chained, and streaming inputs retain the existing backend. No record strings
// are constructed: header/sequence/quality lines append directly to the arena.
class BufferedPlainInput {
    using File=std::unique_ptr<std::FILE,int(*)(std::FILE*)>;
    File file_;
    std::vector<char> buffer_;
    std::size_t at_=0,end_=0;
    std::string ignored_;
    explicit BufferedPlainInput(File file):file_(std::move(file)),buffer_(1<<20){}
    bool fill(){
        if(at_!=end_)return true;
        at_=0;end_=std::fread(buffer_.data(),1,buffer_.size(),file_.get());
        if(std::ferror(file_.get()))throw std::runtime_error("plain FASTQ read failure");
        return end_!=0;
    }
    bool line(std::string& out){
        const auto initial=out.size();std::size_t segment=0;bool visible=true;
        // Match the existing gzgets(8192)/strlen loop even for embedded NULs:
        // a physical segment is consumed in full, but only its visible prefix
        // is appended. CR bytes and unterminated final lines are preserved.
        for(;;){
            if(!fill())return out.size()!=initial;
            const auto* begin=buffer_.data()+at_;const auto available=std::min<std::size_t>(end_-at_,8191-segment);
            const auto* newline=static_cast<const char*>(std::memchr(begin,'\n',available));
            const auto count=newline?std::size_t(newline-begin)+1:available;
            if(visible){const auto* zero=static_cast<const char*>(std::memchr(begin,0,count));
                out.append(begin,zero?std::size_t(zero-begin):count);if(zero)visible=false;}
            at_+=count;segment+=count;
            if(newline&&visible){out.pop_back();return true;}
            if(newline||segment==8191){segment=0;visible=true;}
        }
    }
public:
    static std::unique_ptr<BufferedPlainInput> open(const std::string& path){
        if(path=="-"||path=="/dev/stdin")return {};
        struct stat before{};if(::stat(path.c_str(),&before)||!S_ISREG(before.st_mode))return {};
        File file(std::fopen(path.c_str(),"rb"),std::fclose);if(!file)return {};
        struct stat status{};if(::fstat(::fileno(file.get()),&status)||!S_ISREG(status.st_mode))return {};
        unsigned char magic[2]{};const auto n=std::fread(magic,1,2,file.get());
        if((n==2&&magic[0]==0x1f&&magic[1]==0x8b)||std::ferror(file.get()))return {};
        if(std::fseek(file.get(),0,SEEK_SET))return {};return std::unique_ptr<BufferedPlainInput>(new BufferedPlainInput(std::move(file)));
    }
    template<class Batch> void read(Batch& batch,std::size_t limit){
        while(batch.rows.size()<limit){
            const auto start=batch.bytes.size();if(!line(batch.bytes))break;const auto h=batch.bytes.size()-start;
            if(!line(batch.bytes))throw std::runtime_error("truncated plain FASTQ sequence");const auto s=batch.bytes.size()-start-h;
            ignored_.clear();if(!line(ignored_))throw std::runtime_error("truncated plain FASTQ plus line");
            if(!line(batch.bytes))throw std::runtime_error("truncated plain FASTQ quality");const auto q=batch.bytes.size()-start-h-s;
            batch.rows.push_back({start,h,s,q});
        }
    }
};
}
