#pragma once
#include <cstddef>
#include <cstdint>
#include <utility>
#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace fqdup::r8qkeep {
// Only exact integer overlap counts; callers retain offset order and tie rules.
// Include globally before the contextually included alignment headers.
inline std::pair<int,int> geometry_counts(const char* a,const char* b,std::size_t n) {
    int informative=0,mismatches=0;
    std::size_t i=0;
#if defined(__AVX2__)
    const auto unknown=_mm256_set1_epi8('N');
    for(;i+32<=n;i+=32) {
        const auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a+i));
        const auto y=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b+i));
        const auto unknowns=_mm256_or_si256(_mm256_cmpeq_epi8(x,unknown),
                                          _mm256_cmpeq_epi8(y,unknown));
        const std::uint32_t known=~std::uint32_t(_mm256_movemask_epi8(unknowns));
        const std::uint32_t equal=std::uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(x,y)));
        informative+=__builtin_popcount(known);
        mismatches+=__builtin_popcount(known & ~equal);
    }
#endif
    for(;i<n;++i) {
        const int known=(a[i]!='N') & (b[i]!='N');
        informative+=known;
        mismatches+=known & (a[i]!=b[i]);
    }
    return {informative,mismatches};
}
}
