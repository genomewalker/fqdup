#pragma once
#include "worker_pool.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <string_view>
#include <vector>
namespace fqdup::r8qkeep {
// Secondary lossless compression: bounded independent blocks use local packed
// codes or run tokens. Every original quality byte remains reconstructible.
class QualityBlocks {
public:
    static constexpr std::size_t block_size=4096;
    struct Cursor {const QualityBlocks* owner=nullptr;std::uint64_t generation=0;std::size_t block=std::numeric_limits<std::size_t>::max();bool decoded=false;std::vector<char> text;};
    struct Edit {
        std::uint64_t position:56,value:8;
        Edit(std::uint64_t p=0,unsigned char v=0):position(p),value(v){if(p>=(std::uint64_t(1)<<56))throw std::out_of_range("quality edit address");}
    };
    static_assert(sizeof(Edit)==8);
private:
    struct Block {std::uint64_t end,offset;std::uint32_t bytes;};
    std::vector<unsigned char> data_;
    std::vector<Block> blocks_;
    std::vector<std::uint32_t> pages_;
    std::string pending_;
    std::uint64_t generation_=1;
    std::size_t live_bytes_=0;
    static void encode(std::string_view text,std::vector<unsigned char>& output){
        if(text.empty()||text.size()>block_size)throw std::logic_error("quality block size");
        std::array<int,256> codes;codes.fill(-1);std::array<unsigned char,16> alphabet{};unsigned count=0;
        for(unsigned char c:text)if(codes[c]<0){codes[c]=count;if(count<16)alphabet[count]=c;if(++count==17)break;}
        unsigned width=0;while((1U<<width)<count)++width;if(width>4)width=8;
        std::array<unsigned char,2*block_size+16> runs{};std::size_t run_bytes=0;
        if(width)for(std::size_t i=0;i<text.size();){auto end=i+1;while(end<text.size()&&text[end]==text[i])++end;
            auto token=std::uint32_t((end-i-1)<<width)|(width==8?static_cast<unsigned char>(text[i]):unsigned(codes[static_cast<unsigned char>(text[i])]));
            do{runs[run_bytes++]=static_cast<unsigned char>((token&127)|(token>=128?128:0));token>>=7;}while(token);i=end;
        }
        const auto raw_bytes=(text.size()*width+7)/8;const bool rle=width&&run_bytes<raw_bytes;
        output.push_back(static_cast<unsigned char>(width|(rle?128:0)));
        if(width<8){output.push_back(static_cast<unsigned char>(count));output.insert(output.end(),alphabet.begin(),alphabet.begin()+count);}
        if(rle){output.insert(output.end(),runs.begin(),runs.begin()+run_bytes);return;}
        if(width==8){output.insert(output.end(),text.begin(),text.end());return;}
        const auto start=output.size();output.resize(start+raw_bytes,0);
        if(width)for(std::size_t i=0;i<text.size();++i){const auto code=unsigned(codes[static_cast<unsigned char>(text[i])]);const auto bit=i*width,index=start+bit/8;const unsigned shift=bit%8;
            output[index]|=static_cast<unsigned char>(code<<shift);if(shift+width>8)output[index+1]|=static_cast<unsigned char>(code>>(8-shift));}
    }
    std::uint64_t committed()const{return blocks_.empty()?0:blocks_.back().end;}
    void index_block(Block block){
        if(blocks_.size()>=std::numeric_limits<std::uint32_t>::max())throw std::length_error("quality block count");
        const auto index=static_cast<std::uint32_t>(blocks_.size());blocks_.push_back(block);live_bytes_+=block.bytes;
        while(std::uint64_t(pages_.size())*block_size<block.end)pages_.push_back(index);
    }
    std::size_t find(std::uint64_t position)const{
        auto i=std::size_t(pages_.at(position/block_size));while(blocks_[i].end<=position)++i;return i;
    }
    void decode(std::size_t index,Cursor& cursor)const{
        const auto& block=blocks_[index];const auto start=index?blocks_[index-1].end:0;const auto n=block.end-start;
        const auto* at=data_.data()+block.offset;const auto* end=at+block.bytes;const unsigned mode=*at++,width=mode&15;
        std::array<unsigned char,16> alphabet{};unsigned count=256;
        if(width<8){count=*at++;if(!count||count>16||std::size_t(end-at)<count)throw std::logic_error("quality alphabet");std::copy_n(at,count,alphabet.begin());at+=count;}
        cursor.text.resize(n);
        if(!width)std::fill(cursor.text.begin(),cursor.text.end(),char(alphabet[0]));
        else if(mode&128){
            std::size_t i=0;const unsigned mask=(1U<<width)-1;
            while(i<n){std::uint32_t token=0;unsigned shift=0,byte;
                do{if(at==end||shift>28)throw std::logic_error("quality run");byte=*at++;token|=(byte&127)<<shift;shift+=7;}while(byte&128);
                const auto length=std::size_t(token>>width)+1;const auto code=token&mask;
                if(length>n-i||code>=count)throw std::logic_error("quality run range");
                std::fill_n(cursor.text.begin()+i,length,char(width==8?code:alphabet[code]));i+=length;
            }
        }else if(width==8){if(std::size_t(end-at)<n)throw std::logic_error("quality raw range");std::memcpy(cursor.text.data(),at,n);}
        else{
            if(std::size_t(end-at)<(n*width+7)/8)throw std::logic_error("quality packed range");
            for(std::size_t i=0;i<n;++i){const auto bit=i*width;const unsigned shift=bit%8;unsigned code=at[bit/8]>>shift;
                if(shift+width>8)code|=unsigned(at[bit/8+1])<<(8-shift);code&=(1U<<width)-1;
                if(code>=count)throw std::logic_error("quality code");cursor.text[i]=char(alphabet[code]);}
        }
        cursor.owner=this;cursor.generation=generation_;cursor.block=index;cursor.decoded=true;
    }
    void compact(){
        std::vector<unsigned char> next;next.reserve(live_bytes_);
        for(auto& block:blocks_){const auto offset=next.size();next.insert(next.end(),data_.begin()+block.offset,data_.begin()+block.offset+block.bytes);block.offset=offset;}
        data_=std::move(next);
    }
public:
    std::uint64_t size()const{return committed()+pending_.size();}
    std::size_t owned_bytes()const{return data_.capacity()+blocks_.capacity()*sizeof(Block)+pages_.capacity()*sizeof(std::uint32_t)+pending_.capacity();}
    void append(std::string_view text){
        while(!text.empty()){const auto n=std::min(block_size-pending_.size(),text.size());pending_.append(text.data(),n);text.remove_prefix(n);if(pending_.size()==block_size)finish();}
    }
    void finish(){
        if(pending_.empty())return;const auto end=committed()+pending_.size(),offset=data_.size();encode(pending_,data_);index_block({end,offset,static_cast<std::uint32_t>(data_.size()-offset)});pending_.clear();
    }
    void join(const QualityBlocks& part){
        if(this==&part)throw std::logic_error("quality self join");finish();const auto start=committed();
        if(part.live_bytes_==part.data_.size()){
            const auto offset=data_.size();data_.insert(data_.end(),part.data_.begin(),part.data_.end());
            for(const auto& b:part.blocks_)index_block({start+b.end,offset+b.offset,b.bytes});
        }else for(const auto& b:part.blocks_){const auto offset=data_.size();data_.insert(data_.end(),part.data_.begin()+b.offset,part.data_.begin()+b.offset+b.bytes);index_block({start+b.end,offset,b.bytes});}
        append(part.pending_);finish();
    }
    void read(std::uint64_t at,std::size_t n,std::string& out,Cursor& cursor)const{
        if(at>size()||n>size()-at)throw std::out_of_range("quality read");out.resize(n);std::size_t written=0;
        while(written<n){
            if(at>=committed()){const auto count=n-written;std::memcpy(out.data()+written,pending_.data()+at-committed(),count);return;}
            std::size_t index=cursor.block;
            const bool cached=cursor.owner==this&&cursor.generation==generation_&&index<blocks_.size()&&at<blocks_[index].end&&(!index||at>=blocks_[index-1].end);
            if(!cached)index=find(at);
            const auto start=index?blocks_[index-1].end:0;const auto count=std::min<std::size_t>(n-written,blocks_[index].end-at);
            const auto* bytes=data_.data()+blocks_[index].offset;const unsigned mode=*bytes++,width=mode&15;
            if(mode&128){
                if(!cached||!cursor.decoded)decode(index,cursor);
                std::memcpy(out.data()+written,cursor.text.data()+at-start,count);
            }else{
                // Raw blocks support direct slices; random record order must
                // never force 4096 decoded bytes for a short quality string.
                const unsigned char* alphabet=nullptr;unsigned symbols=256;
                if(width<8){symbols=*bytes++;alphabet=bytes;bytes+=symbols;}
                if(!width)std::fill_n(out.data()+written,count,char(alphabet[0]));
                else if(width==8)std::memcpy(out.data()+written,bytes+at-start,count);
                else for(std::size_t i=0;i<count;++i){const auto bit=(at-start+i)*width;const unsigned shift=bit%8;unsigned code=bytes[bit/8]>>shift;
                    if(shift+width>8)code|=unsigned(bytes[bit/8+1])<<(8-shift);code&=(1U<<width)-1;
                    if(code>=symbols)throw std::logic_error("quality slice code");out[written+i]=char(alphabet[code]);}
                cursor.owner=this;cursor.generation=generation_;cursor.block=index;cursor.decoded=false;
            }
            written+=count;at+=count;
        }
    }
    void edit(std::uint64_t position,unsigned char value){
        if(position>=size())throw std::out_of_range("quality edit");finish();const auto i=find(position);Cursor cursor;decode(i,cursor);
        const auto start=i?blocks_[i-1].end:0;cursor.text[position-start]=char(value);
        const auto offset=data_.size();encode(std::string_view(cursor.text.data(),cursor.text.size()),data_);
        live_bytes_-=blocks_[i].bytes;blocks_[i].offset=offset;blocks_[i].bytes=data_.size()-offset;live_bytes_+=blocks_[i].bytes;
        ++generation_;if(data_.size()>2*live_bytes_+block_size)compact();
    }
    void rewrite(std::vector<Edit> edits,WorkerPool* pool=nullptr,int workers=1){
        if(edits.empty())return;finish();
        // Blocks are independent. Stable bucketing preserves the journal order
        // within each block (including repeated positions) in linear time;
        // position ordering inside a block is unnecessary for exact replay.
        std::vector<std::size_t> offsets(blocks_.size()+1);
        for(const auto& edit:edits){if(edit.position>=size())throw std::out_of_range("quality bulk edit");++offsets[find(edit.position)+1];}
        for(std::size_t i=1;i<offsets.size();++i)offsets[i]+=offsets[i-1];
        {
            auto next=offsets;std::vector<Edit> grouped(edits.size());
            for(const auto& edit:edits)grouped[next[find(edit.position)]++]=edit;
            edits.swap(grouped);
        }
        struct Chunk {std::vector<Block> blocks;std::vector<unsigned char> data;};
        const auto produce=[&](std::size_t begin,std::size_t end){
            Chunk out;out.blocks.reserve(end-begin);Cursor cursor;
            for(auto i=begin;i<end;++i){const auto offset=out.data.size();const auto& block=blocks_[i];
                if(offsets[i]==offsets[i+1])out.data.insert(out.data.end(),data_.begin()+block.offset,data_.begin()+block.offset+block.bytes);
                else{decode(i,cursor);const auto start=i?blocks_[i-1].end:0;for(auto j=offsets[i];j<offsets[i+1];++j)cursor.text[edits[j].position-start]=char(edits[j].value);encode(std::string_view(cursor.text.data(),cursor.text.size()),out.data);}
                out.blocks.push_back({block.end,offset,static_cast<std::uint32_t>(out.data.size()-offset)});
            }return out;
        };
        std::vector<Block> next_blocks;next_blocks.reserve(blocks_.size());std::vector<unsigned char> next_data;next_data.reserve(live_bytes_);
        const auto consume=[&](Chunk chunk){const auto offset=next_data.size();next_data.insert(next_data.end(),chunk.data.begin(),chunk.data.end());for(auto b:chunk.blocks){b.offset+=offset;next_blocks.push_back(b);}};
        constexpr std::size_t batch=32;
        if(pool&&workers>1){
            std::deque<std::future<Chunk>> pending;std::size_t next=0;
            try{while(next<blocks_.size()||!pending.empty()){
                while(next<blocks_.size()&&pending.size()<std::size_t(workers*2)){const auto begin=next,end=std::min(blocks_.size(),begin+batch);next=end;pending.push_back(pool->submit([&,begin,end]{return produce(begin,end);}));}
                auto future=std::move(pending.front());pending.pop_front();consume(future.get());
            }}catch(...){for(auto& future:pending)if(future.valid())future.wait();throw;}
        }else for(std::size_t begin=0;begin<blocks_.size();begin+=batch)consume(produce(begin,std::min(blocks_.size(),begin+batch)));
        blocks_=std::move(next_blocks);data_=std::move(next_data);live_bytes_=data_.size();++generation_;
    }
};
}
