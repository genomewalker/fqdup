#pragma once
#include "compact_store.hpp"
namespace fqdup::r8qkeep {
struct UnsupportedCompactAlphabet {};
FamilyResult correct_compact_family(CompactPairStore& store,int threads,bool debug,WorkerPool* pool,bool calibration=false);
}
