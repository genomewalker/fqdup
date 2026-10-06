#pragma once
#include "native.hpp"
namespace fqdup::r8qkeep {
struct CensusResult {
    std::string nominal_tsv, full_tsv, proposed_edits_tsv;
};
CensusResult census(const std::vector<FamilyPair>& pairs,const NativeResult& native);
}
