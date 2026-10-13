#pragma once
#include "family.hpp"
#include "worker_pool.hpp"
#include "packed_kernels.hpp"
#include "quality_blocks.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <tuple>

namespace fqdup::r8qkeep {
// Columnar input storage. Sequence offsets also address independently compressed quality
// blocks. Only exceptional (non-ACGT) symbols carry a position and byte.
// Header checkpoints bound random access; sequential decoders reuse their prefix.
class CompactPairStore {
public:
    static constexpr std::size_t header_interval=32;
    struct Decoder {
        std::size_t next=std::numeric_limits<std::size_t>::max(),offset=0;
        std::string header;
        QualityBlocks::Cursor quality;
    };
private:
    std::vector<std::uint8_t> bases_,headers_;
    QualityBlocks qualities_;
    std::vector<std::uint64_t> starts_{0},checkpoints_,exception_positions_;
    std::vector<unsigned char> exception_values_;
    std::array<unsigned char,256> alphabet_{};
    std::array<int,256> codes_{};
    unsigned alphabet_size_=0,width_=0;
    std::string previous_header_,id_common_;
    std::vector<std::uint32_t> order_;

    static unsigned get(const std::vector<std::uint8_t>& bytes,std::uint64_t at,unsigned width) {
        if(!width)return 0;
        const auto bit=at*width,index=bit/8;const unsigned shift=bit%8;
        unsigned value=bytes[index]>>shift;
        if(shift+width>8)value|=unsigned(bytes[index+1])<<(8-shift);
        return value&((1U<<width)-1);
    }
    static void put(std::vector<std::uint8_t>& bytes,std::uint64_t at,unsigned width,unsigned value) {
        if(!width)return;
        const auto bit=at*width,index=bit/8;const unsigned shift=bit%8,mask=(1U<<width)-1;
        bytes[index]=std::uint8_t((bytes[index]&~(mask<<shift))|(value<<shift));
        if(shift+width>8)bytes[index+1]=std::uint8_t((bytes[index+1]&~(mask>>(8-shift)))|(value>>(8-shift)));
    }
    static unsigned base_code(unsigned char c) {
        switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return 4;}
    }
    static std::string_view id_view(std::string_view header) {
        if(header.empty())throw std::runtime_error("empty FASTQ header");
        header.remove_prefix(1);header=header.substr(0,header.find_first_of(" \t"));
        if(header.size()>2&&header[header.size()-2]=='/'&&(header.back()=='1'||header.back()=='2'))header.remove_suffix(2);
        return header;
    }
    void common_id(std::string_view id){
        if(empty())id_common_.assign(id.data(),id.size());
        else{std::size_t n=0;while(n<std::min(id_common_.size(),id.size())&&id_common_[n]==id[n])++n;id_common_.resize(n);}
    }
    void varint(std::size_t value) {
        while(value>=128){headers_.push_back(std::uint8_t(value)|128);value>>=7;}
        headers_.push_back(std::uint8_t(value));
    }
    std::size_t varint_at(std::size_t& offset)const {
        std::size_t value=0;unsigned shift=0;
        for(;;){const auto byte=headers_.at(offset++);value|=std::size_t(byte&127)<<shift;if(!(byte&128))return value;shift+=7;}
    }
    void add_alphabet(std::string_view text,WorkerPool* =nullptr,int =1) {
        bool changed=false;
        for(unsigned char c:text)if(codes_[c]<0){codes_[c]=int(alphabet_size_);alphabet_[alphabet_size_++]=c;changed=true;}
        if(!changed)return;
        width_=0;
        while((1U<<width_)<alphabet_size_)++width_;
        if(width_>4)width_=8;
    }
    template<class Read> void append_read(const Read& read) {
        if(read.s.size()!=read.q.size())throw std::runtime_error("FASTQ lengths");
        const auto row=starts_.size()-1;
        std::size_t prefix=0;
        if(row%header_interval==0)checkpoints_.push_back(headers_.size());
        else while(prefix<std::min(previous_header_.size(),read.h.size())&&previous_header_[prefix]==read.h[prefix])++prefix;
        varint(prefix);varint(read.h.size()-prefix);
        headers_.insert(headers_.end(),read.h.begin()+prefix,read.h.end());previous_header_=read.h;
        add_alphabet(read.q);
        const auto start=starts_.back(),end=start+read.s.size();
        bases_.resize((end*2+7)/8);
        std::size_t i=0;
        for(;i+32<=read.s.size();i+=32){
            const auto block=packed::block(read.s.data()+i);
            packed::store32(bases_.data()+(start+i)/4,2*((start+i)%4),block.codes);
            auto exceptions=block.exceptions;
            while(exceptions){const auto bit=__builtin_ctz(exceptions);exceptions&=exceptions-1;
                exception_positions_.push_back(start+i+bit);exception_values_.push_back(read.s[i+bit]);}
        }
        for(;i<read.s.size();++i){
            const unsigned code=base_code(read.s[i]);put(bases_,start+i,2,code&3);
            if(code==4){exception_positions_.push_back(start+i);exception_values_.push_back(read.s[i]);}
        }
        qualities_.append(read.q);
        starts_.push_back(end);
    }
    std::size_t physical(std::size_t index)const {return order_.empty()?index:order_.at(index);}
    void decode_sequence(std::size_t row,std::string& sequence)const {
        const auto start=starts_.at(row),end=starts_.at(row+1),length=end-start;
        sequence.resize(length);
        std::size_t i=0;
        for(;i+32<=length;i+=32)packed::decode32(bases_.data()+(start+i)/4,2*((start+i)%4),sequence.data()+i);
        for(;i<length;++i)sequence[i]="ACGT"[get(bases_,start+i,2)];
        auto e=std::lower_bound(exception_positions_.begin(),exception_positions_.end(),start);
        for(;e!=exception_positions_.end()&&*e<end;++e)sequence[*e-start]=char(exception_values_[e-exception_positions_.begin()]);
    }
    void decode_read(std::size_t row,FamilyRead& read,Decoder& decoder)const {
        header(row,read.h,decoder);const auto id=id_view(read.h);read.id.assign(id.data(),id.size());
        decode_sequence(row,read.s);qualities_.read(starts_[row],starts_[row+1]-starts_[row],read.q,decoder.quality);
    }
    void set_base(std::uint64_t at,char base) {
        const auto code=base_code(base);put(bases_,at,2,code&3);
        const auto found=std::lower_bound(exception_positions_.begin(),exception_positions_.end(),at);
        const auto slot=found-exception_positions_.begin();
        if(found!=exception_positions_.end()&&*found==at)exception_values_[slot]=base;
        else if(code==4){exception_positions_.insert(found,at);exception_values_.insert(exception_values_.begin()+slot,base);}
    }
public:
    #include "compact_join.inc"
    #include "parallel_quality.inc"
#include "absolute_edits.inc"
    CompactPairStore(){codes_.fill(-1);}
    std::size_t size()const{return (starts_.size()-1)/2;}
    bool empty()const{return size()==0;}
    void prepare_qualities(std::string_view symbols,WorkerPool* pool=nullptr,int workers=1){add_alphabet(symbols,pool,workers);}
    std::array<std::size_t,5> capacity_bytes()const{return {bases_.capacity(),qualities_.owned_bytes(),headers_.capacity()+id_common_.capacity()+8*checkpoints_.capacity(),8*starts_.capacity()+4*order_.capacity(),8*exception_positions_.capacity()+exception_values_.capacity()};}
    unsigned quality_width()const{return width_;}
    unsigned quality_alphabet_size()const{return alphabet_size_;}
    std::uint64_t symbols()const{return starts_.back();}
    void append(const FamilyRead& first,const FamilyRead& second){append(first.h,first.s,first.q,second.h,second.s,second.q);}
    void append(std::string_view h1,std::string_view s1,std::string_view q1,std::string_view h2,std::string_view s2,std::string_view q2) {
        if(!order_.empty())throw std::logic_error("append after canonical order");
        if(size()>=std::numeric_limits<std::uint32_t>::max())throw std::length_error("compact pair count");
        if(s1.size()!=q1.size()||s2.size()!=q2.size())throw std::runtime_error("FASTQ lengths");
        if(id_view(h1)!=id_view(h2))throw std::runtime_error("pair order "+std::string(id_view(h1)));
        common_id(id_view(h1));
        struct View {std::string_view h,s,q;};
        append_read(View{h1,s1,q1});append_read(View{h2,s2,q2});
    }
    void header(std::size_t row,std::string& out,Decoder& decoder)const {
        if(row>=starts_.size()-1)throw std::out_of_range("compact header");
        if(decoder.next>row||row-decoder.next>=header_interval){decoder.next=row/header_interval*header_interval;decoder.offset=checkpoints_.at(row/header_interval);decoder.header.clear();}
        while(decoder.next<=row){
            const auto prefix=varint_at(decoder.offset),suffix=varint_at(decoder.offset);
            if(prefix>decoder.header.size()||suffix>headers_.size()-decoder.offset)throw std::logic_error("compact header bounds");
            decoder.header.resize(prefix);decoder.header.append(reinterpret_cast<const char*>(headers_.data()+decoder.offset),suffix);
            decoder.offset+=suffix;++decoder.next;
        }
        out=decoder.header;
    }
    void id(std::size_t index,std::string& out,Decoder& decoder)const {
        header(2*physical(index),out,decoder);const auto v=id_view(out);out.assign(v.data(),v.size());
    }
    void read(std::size_t index,FamilyPair& pair,Decoder& decoder)const {
        const auto row=2*physical(index);decode_read(row,pair.r[0],decoder);decode_read(row+1,pair.r[1],decoder);pair.id=pair.r[0].id;
    }
    void sequences(std::size_t index,FamilyPair& pair)const {
        const auto row=2*physical(index);decode_sequence(row,pair.r[0].s);decode_sequence(row+1,pair.r[1].s);
    }
    void content(std::size_t index,FamilyPair& pair,Decoder& decoder)const {
        sequences(index,pair);const auto row=2*physical(index);
        for(unsigned m=0;m<2;++m)qualities_.read(starts_[row+m],starts_[row+m+1]-starts_[row+m],pair.r[m].q,decoder.quality);
    }
    void read(std::size_t index,FamilyPair& pair)const {Decoder decoder;read(index,pair,decoder);}
    // All workers finish reading before these sparse corrections are applied.
    // This keeps quality alphabet expansion and cross-word writes race-free.
    void edit(std::size_t index,unsigned mate,std::size_t position,char base,char quality) {
        if(mate>1)throw std::out_of_range("compact mate");
        const auto row=2*physical(index)+mate,start=starts_.at(row);
        if(position>=starts_.at(row+1)-start)throw std::out_of_range("compact edit");
        add_alphabet(std::string_view(&quality,1));const auto at=start+position;
        qualities_.edit(at,static_cast<unsigned char>(quality));set_base(at,base);
    }
    bool canonicalize(WorkerPool* pool=nullptr,int workers=1,
                      std::vector<std::pair<std::string,std::size_t>>* sample=nullptr,
                      std::uint64_t (*sample_hash)(std::string_view)=nullptr) {
        auto measured=std::chrono::steady_clock::now();
        const bool diagnostics=std::getenv("FQDUP_CANONICAL_DIAGNOSTICS")!=nullptr;
        const auto phase=[&](const char* name){const auto now=std::chrono::steady_clock::now();
            if(diagnostics)std::cerr<<"CANONICAL_PHASE "<<name<<" seconds="<<std::chrono::duration<double>(now-measured).count()<<'\n';
            measured=now;};
        if(!order_.empty())throw std::logic_error("already canonical");
        // A library often shares a long instrument/run prefix. It contributes
        // nothing to lexicographic order; index the first varying bytes.
        workers=pool?bounded_worker_count(workers):1;
        const auto skip=id_common_.size();
        // Three big-endian words cover the varying tile/x/y suffix of lane
        // headers. Short generator IDs retain the original single-word passes.
        std::vector<std::array<std::uint64_t,3>> prefix(size());
        std::vector<std::uint64_t> selected(sample?(size()+63)/64:0);
        const auto threshold=size()?(static_cast<unsigned __int128>(std::min<std::size_t>(25000,size()))<<64)/size():0;
        const auto collect=[&](std::size_t begin,std::size_t end){
            Decoder decoder;std::string h;
            for(auto i=begin;i<end;++i){header(2*i,h,decoder);const auto id=id_view(h);
                if(sample&&static_cast<unsigned __int128>(sample_hash(id))<threshold)selected[i/64]|=std::uint64_t{1}<<(i%64);
                for(std::size_t j=0;j<24;++j){auto& word=prefix[i][j/8];
                    word=(word<<8)|(skip+j<id.size()?static_cast<unsigned char>(id[skip+j]):0);}}
        };
        if(pool&&workers>1&&size()>=4096){
            std::vector<std::future<void>> pending;pending.reserve(workers);
            // Keep sample-bit words private to each collector.
            const auto span=((size()+workers*64-1)/(workers*64))*64;
            try{for(std::size_t begin=0;begin<size();begin+=span){const auto end=std::min(size(),begin+span);pending.push_back(pool->submit([&,begin,end]{collect(begin,end);}));}
                for(auto& task:pending)task.get();
            }catch(...){for(auto& task:pending)if(task.valid())task.wait();throw;}
        }else collect(0,size());
        phase("collect");
        std::vector<std::uint32_t> order(size());std::iota(order.begin(),order.end(),std::uint32_t{0});
        FamilyPair a,b;Decoder da,db;std::string ha,hb;
        const auto compare_ids=[&](auto i,auto j){header(2*i,ha,da);header(2*j,hb,db);return id_view(ha).compare(id_view(hb));};
        const auto key=[](const FamilyPair& p){return std::tie(p.id,p.r[0].h,p.r[0].s,p.r[0].q,p.r[1].h,p.r[1].s,p.r[1].q);};
        const auto less=[&](auto i,auto j){
            if(prefix[i]!=prefix[j])return prefix[i]<prefix[j];
            const auto ids=compare_ids(i,j);if(ids)return ids<0;
            read(i,a,da);read(j,b,db);return key(a)<key(b);
        };
        if(size()<4096)std::sort(order.begin(),order.end(),less);
        else{
            // Stable radix passes order only the fixed-size prefix. Full IDs
            // and the frozen seven-field tuple still resolve every collision.
            constexpr unsigned digit_bits=11,buckets=1U<<digit_bits;
            std::array<std::uint64_t,3> varying{};
            for(const auto& value:prefix)for(unsigned word=0;word<3;++word)varying[word]|=value[word]^prefix[0][word];
            std::vector<std::uint32_t> scratch(order.size());
            std::vector<std::array<std::size_t,buckets>> positions(workers);
            const auto ranges=[&](auto work){
                if(workers==1){work(0,0,order.size());return;}
                std::vector<std::future<void>> pending;pending.reserve(workers);
                try{for(int worker=0;worker<workers;++worker){const auto begin=order.size()*worker/workers,end=order.size()*(worker+1)/workers;
                    pending.push_back(pool->submit([&,worker,begin,end]{work(worker,begin,end);}));}
                    for(auto& task:pending)task.get();
                }catch(...){for(auto& task:pending)if(task.valid())task.wait();throw;}
            };
            for(int word=2;word>=0;--word)for(unsigned shift=0;shift<64;shift+=digit_bits){
                if(!((varying[word]>>shift)&(buckets-1)))continue;
                ranges([&](int worker,std::size_t begin,std::size_t end){auto& counts=positions[worker];counts.fill(0);
                    for(auto i=begin;i<end;++i)++counts[(prefix[order[i]][word]>>shift)&(buckets-1)];});
                std::size_t offset=0;
                for(unsigned bucket=0;bucket<buckets;++bucket)for(auto& counts:positions){const auto next=offset+counts[bucket];counts[bucket]=offset;offset=next;}
                ranges([&](int worker,std::size_t begin,std::size_t end){auto& at=positions[worker];
                    for(auto i=begin;i<end;++i){const auto id=order[i];scratch[at[(prefix[id][word]>>shift)&(buckets-1)]++]=id;}});
                order.swap(scratch);
            }
            for(std::size_t begin=0;begin<order.size();){auto end=begin+1;while(end<order.size()&&prefix[order[begin]]==prefix[order[end]])++end;
                if(end-begin>1)std::sort(order.begin()+begin,order.begin()+end,less);begin=end;}
        }

        phase("ordering");
        bool unique=true;
        for(std::size_t i=1;i<order.size();++i)
            if(prefix[order[i-1]]==prefix[order[i]]&&!compare_ids(order[i-1],order[i]))unique=false;
        if(sample){
            sample->clear();Decoder decoder;std::string id;
            for(std::size_t i=0;i<order.size();++i)if(selected[order[i]/64]&(std::uint64_t{1}<<(order[i]%64))){
                header(2*order[i],id,decoder);sample->emplace_back(id_view(id),i);
            }
        }
        order_=std::move(order);phase("duplicates");return unique;
    }
    void input_order(){std::vector<std::uint32_t>().swap(order_);}
    std::size_t owned_bytes()const {
        return bases_.capacity()+qualities_.owned_bytes()+headers_.capacity()+id_common_.capacity()+exception_values_.capacity()
            +8*(starts_.capacity()+checkpoints_.capacity()+exception_positions_.capacity())+4*order_.capacity();
    }
};
}
