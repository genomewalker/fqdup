#pragma once
#include <cerrno>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#if defined(__linux__)
#include <fcntl.h>
#endif
namespace fqdup::r8qkeep {
// A single writer owns each stream. Schedule completed byte ranges while
// workers compute the next chunks; ordinary close still completes and reports
// errors. Range writeback is a throughput hint, not a durability primitive.
class StreamingOutput {
public:
    using Writeback=int(*)(int,std::uint64_t,std::uint64_t);
private:
    static int platform_writeback(int fd,std::uint64_t offset,std::uint64_t count){
#if defined(__linux__)
        if(offset>std::uint64_t(std::numeric_limits<off64_t>::max())||count>std::uint64_t(std::numeric_limits<off64_t>::max())-offset){errno=EINVAL;return -1;}
        return ::sync_file_range(fd,static_cast<off64_t>(offset),static_cast<off64_t>(count),SYNC_FILE_RANGE_WRITE);
#else
        (void)fd;(void)offset;(void)count;errno=ENOSYS;return -1;
#endif
    }
    std::FILE* file_=nullptr;
    Writeback writeback_;
    std::uint64_t written_=0,scheduled_=0,submitted_=0;
    std::size_t window_;
    bool good_=false,enabled_=true;
    void failed(){good_=false;throw std::runtime_error("write output");}
public:
    explicit StreamingOutput(std::size_t window=8U<<20,Writeback writeback=platform_writeback):writeback_(writeback),window_(window){
        if(!window_||!writeback_)throw std::invalid_argument("output writeback settings");
    }
    StreamingOutput(const StreamingOutput&)=delete;
    StreamingOutput& operator=(const StreamingOutput&)=delete;
    StreamingOutput(StreamingOutput&& other)noexcept:file_(std::exchange(other.file_,nullptr)),writeback_(other.writeback_),written_(other.written_),scheduled_(other.scheduled_),submitted_(other.submitted_),window_(other.window_),good_(other.good_),enabled_(other.enabled_){other.good_=false;}
    ~StreamingOutput(){if(file_)std::fclose(file_);}
    explicit operator bool()const{return good_;}
    void open(const std::filesystem::path& path){
        if(file_)throw std::logic_error("output already open");
        written_=scheduled_=submitted_=0;enabled_=true;file_=std::fopen(path.c_str(),"wb");good_=file_!=nullptr;
    }
    void append(std::string_view text){
        if(!file_||!good_)failed();
        if(text.empty())return;
        if(std::fwrite(text.data(),1,text.size(),file_)!=text.size())failed();
        written_+=text.size();
        if(enabled_&&written_-scheduled_>=window_){
            if(std::fflush(file_))failed();
            if(writeback_(::fileno(file_),scheduled_,written_-scheduled_)){
                if(errno==ENOSYS||errno==EINVAL||errno==EOPNOTSUPP||errno==ESPIPE)enabled_=false;
                else failed();
            }
            else submitted_+=written_-scheduled_;
            scheduled_=written_;
        }
    }
    std::uint64_t writeback_bytes()const{return submitted_;}
    bool writeback_enabled()const{return enabled_;}
    void close(){if(file_){const auto owned=std::exchange(file_,nullptr);if(std::fclose(owned))good_=false;}}
};
}
