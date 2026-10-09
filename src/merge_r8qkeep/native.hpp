#pragma once
#include "family.hpp"
#include "stage_records.hpp"
namespace fqdup::r8qkeep {
struct NativeSettings {
    std::string adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAG";
    std::string adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAG";
    std::string library_type = "ds";
    int min_overlap = 11;
    int min_length = 1;
    float max_mismatch_rate = 0.08f;
};
inline constexpr std::size_t native_batch_size=2048;
struct NativeResult {
    std::vector<NativeRecord> records;
    std::string merged_fq, r1_fq, r2_fq, orphan_fq, outcomes_jsonl;
    std::string native_tsv, native_outcomes_json;
};
NativeResult merge_native_range(const std::vector<FamilyPair>& pairs,std::size_t begin,std::size_t count,const NativeSettings& settings,bool trace=false,bool audit=false);
NativeResult merge_native(const std::vector<FamilyPair>& pairs, const NativeSettings& settings = {},bool trace=false,bool audit=true,int threads=1,WorkerPool* shared_pool=nullptr);
}
