#pragma once
// Frozen scalar block encoder from12d9083.
namespace frozen_quality_encoding {
constexpr std::size_t block_size=4096;
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
            auto end=i+1;while(end<text.size()&&text[end]==text[i])++end;
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
            auto end=i+1;if(run<stored)end=ends[run++];else while(end<text.size()&&text[end]==text[i])++end;
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
        if(width)for(std::size_t i=0;i<text.size();++i){const auto code=unsigned(codes[static_cast<unsigned char>(text[i])]);const auto bit=i*width,index=start+bit/8;const unsigned shift=bit%8;
            output[index]|=static_cast<unsigned char>(code<<shift);if(shift+width>8)output[index+1]|=static_cast<unsigned char>(code>>(8-shift));}
    }
}
