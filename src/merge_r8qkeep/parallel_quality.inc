// Included in CompactPairStore's public section. Family reads finish before
// edits: required base changes are applied in journal order, then independent quality
// blocks are rebuilt in bounded parallel chunks. Duplicate positions retain
// the last journal value, including after canonical ordering of pairs.
// The predicate may omit only a base operation that is a no-op at that point
// in the ordered journal. Quality operations are never omitted.
template<class Edits,class ChangeBase> void apply_edits(Edits edits,WorkerPool* pool,int workers,ChangeBase change_base){
    std::vector<QualityBlocks::Edit> qualities;qualities.reserve(edits.size());
    std::array<bool,256> observed{};std::string alphabet;
    for(const auto& edit:edits){
        if(edit.mate>1)throw std::out_of_range("compact mate");
        const auto row=2*physical(edit.pair)+edit.mate,start=starts_.at(row);
        if(edit.position>=starts_.at(row+1)-start)throw std::out_of_range("compact edit");
        const auto at=start+edit.position;if(change_base(edit))set_base(at,edit.base);
        const auto q=static_cast<unsigned char>(edit.quality);
        if(!observed[q]){observed[q]=true;alphabet+=edit.quality;}
        qualities.push_back({at,q});
    }
    add_alphabet(alphabet);Edits().swap(edits);
    qualities_.rewrite(std::move(qualities),pool,workers);
}
template<class Edits> void apply_edits(Edits edits,WorkerPool* pool=nullptr,int workers=1){
    apply_edits(std::move(edits),pool,workers,[](const auto&){return true;});
}
