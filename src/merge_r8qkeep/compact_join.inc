// Included in CompactPairStore's public section. Worker chunks end at header
// checkpoint boundaries; their encoded arenas can be joined without reads.
void append_store(const CompactPairStore& part){
    if(part.empty())return;
    if(!order_.empty()||!part.order_.empty()||(starts_.size()-1)%header_interval)throw std::logic_error("compact join boundary");
    if(size()+part.size()>std::numeric_limits<std::uint32_t>::max())throw std::length_error("compact pair count");
    const auto start=starts_.back(),count=part.symbols();
    add_alphabet(std::string_view(reinterpret_cast<const char*>(part.alphabet_.data()),part.alphabet_size_));
    const auto append_bits=[](std::vector<std::uint8_t>& dst,std::uint64_t at,const std::vector<std::uint8_t>& src,std::uint64_t bits){
        if(!bits)return;
        const auto index=at/8;const unsigned shift=at%8;const auto bytes=(bits+7)/8;
        dst.resize((at+bits+7)/8);
        if(!shift)std::memcpy(dst.data()+index,src.data(),bytes);
        else{
            dst[index]=std::uint8_t((dst[index]&((1U<<shift)-1))|(src[0]<<shift));
            for(std::size_t i=1;i<bytes;++i)dst[index+i]=std::uint8_t((src[i-1]>>(8-shift))|(src[i]<<shift));
            if(index+bytes<dst.size())dst[index+bytes]=src[bytes-1]>>(8-shift);
        }
        const unsigned last=(at+bits)%8;if(last)dst.back()&=(1U<<last)-1;
    };
    append_bits(bases_,start*2,part.bases_,count*2);
    bool direct=width_==part.width_;
    if(direct&&width_!=8)for(unsigned i=0;i<part.alphabet_size_;++i)direct&=alphabet_[i]==part.alphabet_[i];
    if(direct)append_bits(qualities_,start*width_,part.qualities_,count*width_);
    else{
        std::vector<std::uint8_t> recoded((count*width_+7)/8);
        if(width_==part.width_&&(width_==1||width_==2||width_==4)){
            std::array<unsigned char,256> table{};const unsigned mask=(1U<<width_)-1;
            for(unsigned b=0;b<256;++b)for(unsigned shift=0;shift<8;shift+=width_){
                const auto old=(b>>shift)&mask;table[b]|=static_cast<unsigned char>((old<part.alphabet_size_?codes_[part.alphabet_[old]]:0)<<shift);}
            for(std::size_t i=0;i<recoded.size();++i)recoded[i]=table[part.qualities_[i]];
        }else for(std::uint64_t i=0;i<count;++i){const auto code=get(part.qualities_,i,part.width_);const auto c=part.width_==8?code:part.alphabet_[code];put(recoded,i,width_,width_==8?c:codes_[c]);}
        append_bits(qualities_,start*width_,recoded,count*width_);
    }
    const auto header_start=headers_.size();
    headers_.insert(headers_.end(),part.headers_.begin(),part.headers_.end());
    for(auto offset:part.checkpoints_)checkpoints_.push_back(header_start+offset);
    for(std::size_t i=1;i<part.starts_.size();++i)starts_.push_back(start+part.starts_[i]);
    for(auto pos:part.exception_positions_)exception_positions_.push_back(start+pos);
    exception_values_.insert(exception_values_.end(),part.exception_values_.begin(),part.exception_values_.end());
    previous_header_=part.previous_header_;
}
