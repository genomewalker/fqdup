#pragma once
#include "family.hpp"
#include "packed_kernels.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <tuple>

namespace fqdup::r8qkeep {
// Columnar input storage. Sequence offsets also address the library-wide quality
// code stream. Only exceptional (non-ACGT) symbols carry a position and byte.
// Header checkpoints bound random access; sequential decoders reuse their prefix.
class CompactPairStore {
public:
    static constexpr std::size_t header_interval=32;
    struct Decoder {
        std::size_t next=std::numeric_limits<std::size_t>::max(),offset=0;
        std::string header;
    };
private:
    std::vector<std::uint8_t> bases_,qualities_,headers_;
    std::vector<std::uint64_t> starts_{0},checkpoints_,exception_positions_;
    std::vector<unsigned char> exception_values_;
    std::array<unsigned char,256> alphabet_{};
    std::array<int,256> codes_{};
    unsigned alphabet_size_=0,width_=0;
    std::string previous_header_;
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
    void varint(std::size_t value) {
        while(value>=128){headers_.push_back(std::uint8_t(value)|128);value>>=7;}
        headers_.push_back(std::uint8_t(value));
    }
    std::size_t varint_at(std::size_t& offset)const {
        std::size_t value=0;unsigned shift=0;
        for(;;){const auto byte=headers_.at(offset++);value|=std::size_t(byte&127)<<shift;if(!(byte&128))return value;shift+=7;}
    }
    void add_alphabet(std::string_view text) {
        const unsigned old_width=width_;bool changed=false;
        for(unsigned char c:text)if(codes_[c]<0){codes_[c]=int(alphabet_size_);alphabet_[alphabet_size_++]=c;changed=true;}
        if(!changed)return;
        width_=0;
        while((1U<<width_)<alphabet_size_)++width_;
        if(width_>4)width_=8;
        if(width_==old_width)return;
        const auto n=starts_.back();
        qualities_.resize((n*width_+7)/8);
        // Widen backwards: destination bits cannot overwrite an unread code.
        for(std::uint64_t i=n;i>0;){--i;const auto old=get(qualities_,i,old_width);const auto c=old_width==8?old:alphabet_[old];put(qualities_,i,width_,width_==8?c:codes_[c]);}
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
        bases_.resize((end*2+7)/8);qualities_.resize((end*width_+7)/8);
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
        if(width_==8){if(!read.q.empty())std::memcpy(qualities_.data()+start,read.q.data(),read.q.size());}
        else if(width_)for(std::size_t i=0;i<read.q.size();++i)put(qualities_,start+i,width_,codes_[static_cast<unsigned char>(read.q[i])]);
        starts_.push_back(end);
    }
    std::size_t physical(std::size_t index)const {return order_.empty()?index:order_.at(index);}
    void decode_read(std::size_t row,FamilyRead& read,Decoder& decoder)const {
        header(row,read.h,decoder);const auto id=id_view(read.h);read.id.assign(id.data(),id.size());
        const auto start=starts_.at(row),end=starts_.at(row+1),length=end-start;
        read.s.resize(length);read.q.resize(length);
        std::size_t i=0;
        for(;i+32<=length;i+=32)packed::decode32(bases_.data()+(start+i)/4,2*((start+i)%4),read.s.data()+i);
        for(;i<length;++i)read.s[i]="ACGT"[get(bases_,start+i,2)];
        if(width_==8){if(length)std::memcpy(read.q.data(),qualities_.data()+start,length);}
        else if(!width_)std::fill(read.q.begin(),read.q.end(),char(alphabet_[0]));
        else for(std::size_t i=0;i<length;++i)read.q[i]=char(alphabet_[get(qualities_,start+i,width_)]);
        auto e=std::lower_bound(exception_positions_.begin(),exception_positions_.end(),start);
        for(;e!=exception_positions_.end()&&*e<end;++e)read.s[*e-start]=char(exception_values_[e-exception_positions_.begin()]);
    }
public:
    CompactPairStore(){codes_.fill(-1);}
    std::size_t size()const{return (starts_.size()-1)/2;}
    bool empty()const{return size()==0;}
    void prepare_qualities(std::string_view symbols){add_alphabet(symbols);}
    std::array<std::size_t,5> capacity_bytes()const{return {bases_.capacity(),qualities_.capacity(),headers_.capacity()+8*checkpoints_.capacity(),8*starts_.capacity()+4*order_.capacity(),8*exception_positions_.capacity()+exception_values_.capacity()};}
    unsigned quality_width()const{return width_;}
    unsigned quality_alphabet_size()const{return alphabet_size_;}
    std::uint64_t symbols()const{return starts_.back();}
    void append(const FamilyRead& first,const FamilyRead& second){append(first.h,first.s,first.q,second.h,second.s,second.q);}
    void append(std::string_view h1,std::string_view s1,std::string_view q1,std::string_view h2,std::string_view s2,std::string_view q2) {
        if(!order_.empty())throw std::logic_error("append after canonical order");
        if(size()>=std::numeric_limits<std::uint32_t>::max())throw std::length_error("compact pair count");
        if(s1.size()!=q1.size()||s2.size()!=q2.size())throw std::runtime_error("FASTQ lengths");
        if(id_view(h1)!=id_view(h2))throw std::runtime_error("pair order "+std::string(id_view(h1)));
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
    void read(std::size_t index,FamilyPair& pair)const {Decoder decoder;read(index,pair,decoder);}
    // All workers finish reading before these sparse corrections are applied.
    // This keeps quality alphabet expansion and cross-word writes race-free.
    void edit(std::size_t index,unsigned mate,std::size_t position,char base,char quality) {
        if(mate>1)throw std::out_of_range("compact mate");
        const auto row=2*physical(index)+mate,start=starts_.at(row);
        if(position>=starts_.at(row+1)-start)throw std::out_of_range("compact edit");
        add_alphabet(std::string_view(&quality,1));const auto at=start+position;
        put(qualities_,at,width_,width_==8?static_cast<unsigned char>(quality):codes_[static_cast<unsigned char>(quality)]);
        const auto code=base_code(base);put(bases_,at,2,code&3);
        const auto found=std::lower_bound(exception_positions_.begin(),exception_positions_.end(),at);
        const auto slot=found-exception_positions_.begin();
        if(found!=exception_positions_.end()&&*found==at)exception_values_[slot]=base;
        else if(code==4){exception_positions_.insert(found,at);exception_values_.insert(exception_values_.begin()+slot,base);}
    }
    bool canonicalize() {
        if(!order_.empty())throw std::logic_error("already canonical");
        // A library often shares a long instrument/run prefix. It contributes
        // nothing to lexicographic order; index the first varying bytes.
        Decoder decoder;std::string h,common;
        for(std::size_t i=0;i<size();++i){
            header(2*i,h,decoder);const auto id=id_view(h);
            if(!i)common.assign(id.data(),id.size());
            else{std::size_t n=0;while(n<std::min(common.size(),id.size())&&common[n]==id[n])++n;common.resize(n);}
            if(common.empty())break;
        }
        const auto skip=common.size();decoder={};
        std::vector<std::uint64_t> prefix(size());
        for(std::size_t i=0;i<size();++i){
            header(2*i,h,decoder);const auto id=id_view(h);
            for(std::size_t j=0;j<8;++j)prefix[i]=(prefix[i]<<8)|(skip+j<id.size()?static_cast<unsigned char>(id[skip+j]):0);
        }
        std::vector<std::uint32_t> order(size());std::iota(order.begin(),order.end(),std::uint32_t{0});
        FamilyPair a,b;Decoder da,db;std::string ha,hb;
        const auto compare_ids=[&](auto i,auto j){header(2*i,ha,da);header(2*j,hb,db);return id_view(ha).compare(id_view(hb));};
        const auto key=[](const FamilyPair& p){return std::tie(p.id,p.r[0].h,p.r[0].s,p.r[0].q,p.r[1].h,p.r[1].s,p.r[1].q);};
        std::sort(order.begin(),order.end(),[&](auto i,auto j){
            if(prefix[i]!=prefix[j])return prefix[i]<prefix[j];
            const auto ids=compare_ids(i,j);if(ids)return ids<0;
            // Only genuinely repeated normalized IDs need sequence/quality
            // tie-breaks. Decoding ordinary reads during sort is unnecessary.
            read(i,a,da);read(j,b,db);return key(a)<key(b);
        });
        bool unique=true;
        for(std::size_t i=1;i<order.size();++i)
            if(prefix[order[i-1]]==prefix[order[i]]&&!compare_ids(order[i-1],order[i]))unique=false;
        order_=std::move(order);return unique;
    }
    void input_order(){std::vector<std::uint32_t>().swap(order_);}
    std::size_t owned_bytes()const {
        return bases_.capacity()+qualities_.capacity()+headers_.capacity()+exception_values_.capacity()
            +8*(starts_.capacity()+checkpoints_.capacity()+exception_positions_.capacity())+4*order_.capacity();
    }
};
}
