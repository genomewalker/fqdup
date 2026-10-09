#pragma once
#include <cstddef>
#include <cstdint>
#include <string_view>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
namespace fqdup_fast_bases {
// Exact uppercase A/C/G/T/N membership, including arbitrary input bytes.
// Loads stay within the supplied view; the scalar tail uses the same alphabet.
inline bool valid(std::string_view sequence){
 std::size_t i=0;
#if defined(__AVX2__)
 const auto a=_mm256_set1_epi8('A'),c=_mm256_set1_epi8('C');
 const auto g=_mm256_set1_epi8('G'),t=_mm256_set1_epi8('T'),n=_mm256_set1_epi8('N');
 for(;i+32<=sequence.size();i+=32){
  const auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(sequence.data()+i));
  const auto ac=_mm256_or_si256(_mm256_cmpeq_epi8(x,a),_mm256_cmpeq_epi8(x,c));
  const auto gt=_mm256_or_si256(_mm256_cmpeq_epi8(x,g),_mm256_cmpeq_epi8(x,t));
  const auto accepted=_mm256_or_si256(_mm256_or_si256(ac,gt),_mm256_cmpeq_epi8(x,n));
  if(std::uint32_t(_mm256_movemask_epi8(accepted))!=UINT32_MAX)return false;
 }
#endif
 constexpr std::uint32_t alphabet=(1u<<('A'-'A'))|(1u<<('C'-'A'))|
     (1u<<('G'-'A'))|(1u<<('T'-'A'))|(1u<<('N'-'A'));
 for(;i<sequence.size();++i){
  const unsigned position=static_cast<unsigned char>(sequence[i])-'A';
  if(position>unsigned('T'-'A')||((alphabet>>position)&1u)==0)return false;
 }
 return true;
}
}
