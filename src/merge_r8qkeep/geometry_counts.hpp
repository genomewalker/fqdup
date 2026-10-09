#pragma once
#include <cstddef>
#include <cstdint>
#include <utility>
#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace fqdup::r8qkeep {
struct GeometryCandidate {
    int insert=-1,informative=0,mismatches=0;
    bool precedes(const GeometryCandidate& other) const {
        return other.insert<0 || mismatches*other.informative<other.mismatches*informative;
    }
};
// Denominators <=512: cross-products fit int, and distinct ratios differ by
// >=1/512^2, far beyond binary64 rounding. Strict ties retain the first offset.
struct GeometryRanking {
    GeometryCandidate best,first_long,second_long;
    void observe(int insert,int n,int mm){
        if(n<11)return;
        const GeometryCandidate candidate{insert,n,mm};
        if(candidate.precedes(best))best=candidate;
        if(n<30)return;
        if(candidate.precedes(first_long)){second_long=first_long;first_long=candidate;}
        else if(candidate.precedes(second_long))second_long=candidate;
    }
    const GeometryCandidate& runner_up() const {
        return best.insert==first_long.insert?second_long:first_long;
    }
};
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
    // The last block ends at n; mask positions already counted above. Loads
    // remain wholly inside both ranges, including unaligned and short inputs.
    if(i<n&&n>=32){
        const auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a+n-32));
        const auto y=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b+n-32));
        const auto unknowns=_mm256_or_si256(_mm256_cmpeq_epi8(x,unknown),_mm256_cmpeq_epi8(y,unknown));
        const std::uint32_t mask=~std::uint32_t(0)<<(32-(n-i));
        const std::uint32_t known=~std::uint32_t(_mm256_movemask_epi8(unknowns))&mask;
        const std::uint32_t equal=std::uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(x,y)));
        informative+=__builtin_popcount(known);mismatches+=__builtin_popcount(known&~equal);i=n;
    }
    if(i<n&&n>=16){
        const auto unknown16=_mm_set1_epi8('N');
        const auto add=[&](std::size_t pos,unsigned mask){
            const auto x=_mm_loadu_si128(reinterpret_cast<const __m128i*>(a+pos));
            const auto y=_mm_loadu_si128(reinterpret_cast<const __m128i*>(b+pos));
            const auto unknowns=_mm_or_si128(_mm_cmpeq_epi8(x,unknown16),_mm_cmpeq_epi8(y,unknown16));
            const unsigned known=~unsigned(_mm_movemask_epi8(unknowns))&mask;
            const unsigned equal=unsigned(_mm_movemask_epi8(_mm_cmpeq_epi8(x,y)));
            informative+=__builtin_popcount(known);mismatches+=__builtin_popcount(known&~equal);
        };
        add(0,65535U);if(n>16)add(n-16,(65535U<<(32-n))&65535U);i=n;
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
