#pragma once
#include "content_memo.hpp"
namespace fqdup::r8qkeep {
// Called after parent routing, interval selection and repair, before caching.
// The confirmed adapter boundary is in raw mate coordinates, like FQT s/e.
inline RepairedRecord final_adapter_clip(const RepairedRecord& selected,int boundary){
    auto record=memo_owned(selected,false,false).record;
    if(record.outcome=="merged"||record.outcome=="dropped"||boundary<0)return record;
    const auto clip=[&](StageRead& read){
        if(read.seq.empty())return;
        const auto tag=read.header.rfind(" FQT=v1;");
        if(tag==std::string::npos)throw std::runtime_error("adapter clip missing FQT");
        const auto s=read.header.find(";s=",tag),e=read.header.find(";e=",tag);
        if(s==std::string::npos||e==std::string::npos)throw std::runtime_error("adapter clip span");
        const auto begin=std::stoull(read.header.substr(s+3));
        const auto end=std::stoull(read.header.substr(e+3));
        if(end<begin||end-begin!=read.seq.size())throw std::runtime_error("adapter clip invalid span");
        const auto stop=std::max(begin,std::min(end,static_cast<unsigned long long>(boundary)));
        if(stop==end)return;
        read.seq.resize(stop-begin);read.qual.resize(stop-begin);
        const auto end_field=read.header.find_first_not_of("0123456789",e+3);
        read.header.replace(e+3,end_field==std::string::npos?std::string::npos:end_field-e-3,std::to_string(stop));
        // Repair edits use raw coordinates and must lie inside the new FQT span.
        const auto qc=read.header.rfind(" FQC=v1;edits=");
        if(qc!=std::string::npos){
            std::istringstream in(read.header.substr(qc+14));std::string entry,edits;
            while(std::getline(in,entry,';'))if(std::stoull(entry)<stop){
                if(!edits.empty())edits+=';';edits+=entry;
            }
            read.header.resize(qc);if(!edits.empty())read.header+=" FQC=v1;edits="+edits;
        }
    };
    clip(record.first);if(record.outcome=="unmerged")clip(record.second);
    if(record.outcome=="unmerged"){
        if(record.first.seq.empty()&&record.second.seq.empty())record.outcome="dropped";
        else if(record.first.seq.empty()){
            record.first=std::move(record.second);record.outcome="orphan_r2";
        }else if(record.second.seq.empty())record.outcome="orphan_r1";
        if(record.outcome=="orphan_r1"||record.outcome=="orphan_r2"){
            const auto ch=record.first.header.rfind(" FQT=v1;ch=");
            if(ch==std::string::npos)throw std::runtime_error("adapter orphan trace");
            const auto start=ch+11,stop=record.first.header.find(';',start);
            record.first.header.replace(start,stop-start,"orphan");
        }
    }else if(record.first.seq.empty())record.outcome="dropped";
    return record;
}
}
