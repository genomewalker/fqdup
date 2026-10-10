#pragma once
// Byte-identical local quality codec: bounded run scans and eight-symbol packing.
#include <array>
#include <vector>
#include <string_view>
#include <cstdint>
#include <stdexcept>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
namespace fqdup::r8qkeep::quality_encoding {
constexpr std::size_t block_size=4096;
static std::size_t run_end(std::string_view text,std::size_t begin){
    auto end=begin+1;const char c=text[begin];
    if(end==text.size()||text[end]!=c)return end;
    ++end;
#if defined(__AVX2__)
    const auto value=_mm256_set1_epi8(c);
    for(;end+32<=text.size();end+=32){
        const auto bytes=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(text.data()+end));
        const auto equal=std::uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(bytes,value)));
        if(equal!=UINT32_MAX)return end+__builtin_ctz(~equal);
    }
#endif
    while(end<text.size()&&text[end]==c)++end;
    return end;
}

    template<unsigned Width>
    static void pack_codes(std::string_view text,const std::array<int,256>& codes,unsigned char* out){
        std::size_t i=0;
        const auto code=[&](std::size_t n){return std::uint32_t(codes[static_cast<unsigned char>(text[n])]);};
        for(;i+8<=text.size();i+=8){
            const std::uint32_t bits=code(i)|(code(i+1)<<Width)|(code(i+2)<<(2*Width))|(code(i+3)<<(3*Width))|
                (code(i+4)<<(4*Width))|(code(i+5)<<(5*Width))|(code(i+6)<<(6*Width))|(code(i+7)<<(7*Width));
            for(unsigned b=0;b<Width;++b)out[i*Width/8+b]=static_cast<unsigned char>(bits>>(8*b));
        }
        for(;i<text.size();++i){const auto bit=i*Width,index=bit/8;const unsigned shift=bit%8;
            out[index]|=static_cast<unsigned char>(code(i)<<shift);if(shift+Width>8)out[index+1]|=static_cast<unsigned char>(code(i)>>(8-shift));}
    }
    static void encode(std::string_view text,std::vector<unsigned char>& output){
        if(text.empty()||text.size()>block_size)throw std::logic_error("quality block size");
        std::array<int,256> codes;codes.fill(-1);std::array<unsigned char,16> alphabet{};unsigned count=0;
        // Discover the alphabet at run boundaries while retaining run ends.
        // Dense prefixes fall back to a byte alphabet scan; the first32 saved
        // runs are still reused. This emits exactly the same block bytes.
        std::array<std::uint16_t,block_size> ends;std::size_t stored=0;
        for(std::size_t i=0;i<text.size();){
            const auto c=static_cast<unsigned char>(text[i]);
            if(count<17&&codes[c]<0){codes[c]=count;if(count<16)alphabet[count]=c;++count;}
            auto end=run_end(text,i);
            ends[stored++]=end;i=end;
            if(stored==32&&i<128){
                for(;i<text.size()&&count<17;++i){const auto c=static_cast<unsigned char>(text[i]);if(codes[c]<0){codes[c]=count;if(count<16)alphabet[count]=c;++count;}}
                break;
            }
        }

        unsigned width=0;while((1U<<width)<count)++width;if(width>4)width=8;
        std::array<unsigned char,2*block_size+16> runs;std::size_t run_bytes=0;
        std::array<std::uint16_t,16> seek_position,seek_offset;std::size_t seek_count=0,run=0;
        if(width)for(std::size_t i=0;i<text.size();){
            auto end=i+1;if(run<stored)end=ends[run++];else end=run_end(text,i);
            while(seek_count*256<end){seek_position[seek_count]=i;seek_offset[seek_count++]=run_bytes;}
            auto token=std::uint32_t((end-i-1)<<width)|(width==8?static_cast<unsigned char>(text[i]):unsigned(codes[static_cast<unsigned char>(text[i])]));
            do{runs[run_bytes++]=static_cast<unsigned char>((token&127)|(token>=128?128:0));token>>=7;}while(token);i=end;
        }
        const auto raw_bytes=(text.size()*width+7)/8;const bool rle=width&&run_bytes<raw_bytes;
        const bool seek=rle&&run_bytes>64;
        output.push_back(static_cast<unsigned char>(width|(rle?128:0)|(seek?64:0)));
        if(width<8){output.push_back(static_cast<unsigned char>(count));output.insert(output.end(),alphabet.begin(),alphabet.begin()+count);}
        if(seek)for(std::size_t i=0;i<seek_count;++i){output.push_back(seek_position[i]&255);output.push_back(seek_position[i]>>8);output.push_back(seek_offset[i]&255);output.push_back(seek_offset[i]>>8);}
        if(rle){output.insert(output.end(),runs.begin(),runs.begin()+run_bytes);return;}
        if(width==8){output.insert(output.end(),text.begin(),text.end());return;}
        const auto start=output.size();output.resize(start+raw_bytes,0);
        switch(width){
        case 1:pack_codes<1>(text,codes,output.data()+start);break;
        case 2:pack_codes<2>(text,codes,output.data()+start);break;
        case 3:pack_codes<3>(text,codes,output.data()+start);break;
        case 4:pack_codes<4>(text,codes,output.data()+start);break;
        }
    }
}
