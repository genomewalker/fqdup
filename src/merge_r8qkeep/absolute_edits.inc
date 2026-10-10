// Logical positions are resolved while the compact store is immutable. Base
// and quality changes commute; each separate journal retains its original order.
std::uint64_t edit_address(std::size_t pair,unsigned mate,std::size_t position)const{
    if(mate>1)throw std::out_of_range("compact mate");
    if(pair>=size())throw std::out_of_range("compact pair");
    const auto row=2*physical(pair)+mate,start=starts_.at(row);
    if(position>=starts_.at(row+1)-start)throw std::out_of_range("compact edit");
    return start+position;
}
void apply_absolute_edits(std::vector<QualityBlocks::Edit> qualities,std::vector<QualityBlocks::Edit> bases,WorkerPool* pool=nullptr,int workers=1){
    for(const auto& edit:bases){if(edit.position>=symbols())throw std::out_of_range("compact base edit");set_base(edit.position,char(edit.value));}
    std::vector<QualityBlocks::Edit>().swap(bases);
    std::array<bool,256> observed{};std::string alphabet;
    for(const auto& edit:qualities){
        if(edit.position>=symbols())throw std::out_of_range("compact quality edit");
        if(!observed[edit.value]){observed[edit.value]=true;alphabet+=char(edit.value);}
    }
    add_alphabet(alphabet);qualities_.rewrite(std::move(qualities),pool,workers);
}
