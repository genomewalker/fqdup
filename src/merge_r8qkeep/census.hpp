#pragma once
#include "native.hpp"
#include "mixture_fit.hpp"
#include "concordance.hpp"
namespace fqdup::r8qkeep {
struct CensusResult {
    std::string nominal_tsv, full_tsv, proposed_edits_tsv;
    std::vector<PairLikelihood> nominal_pairs,full_pairs;
};
CensusResult census(const std::vector<FamilyPair>& pairs,const NativeResult& native,bool debug=true,int threads=1,std::size_t start=0,bool counts=false,WorkerPool* shared_pool=nullptr);
RawConcordance raw_concordance(const std::vector<FamilyPair>& pairs,const NativeResult& native);
}
