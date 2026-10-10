#pragma once
#include <string>
#include <string_view>
#include <stdexcept>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
namespace fqdup_fast_bases {
inline char complement(char c){switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';case 'N':return 'N';default:throw std::invalid_argument("unsupported base");}}
// Input and output storage must be distinct. Invalid blocks use the scalar
// loop so both the first exception and the partial output remain identical.
inline void reverse_complement_into(std::string_view s,std::string& out){
 out.resize(s.size());std::size_t i=0;
#if defined(__AVX2__)
 const auto letters=_mm256_setr_epi8(0,'A',0,'C','T',0,0,'G',0,0,0,0,0,0,'N',0,0,'A',0,'C','T',0,0,'G',0,0,0,0,0,0,'N',0);
 const auto complements=_mm256_setr_epi8(0,'T',0,'G','A',0,0,'C',0,0,0,0,0,0,'N',0,0,'T',0,'G','A',0,0,'C',0,0,0,0,0,0,'N',0);
 const auto order=_mm256_setr_epi8(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
 for(;i+32<=s.size();i+=32){
  const auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(s.data()+i)),index=_mm256_and_si256(x,_mm256_set1_epi8(15));
  const auto expected=_mm256_shuffle_epi8(letters,index);
  const auto valid=_mm256_and_si256(_mm256_cmpeq_epi8(x,expected),_mm256_cmpgt_epi8(expected,_mm256_setzero_si256()));
  if(static_cast<unsigned>(_mm256_movemask_epi8(valid))!=0xffffffffu)break;
  const auto mapped=_mm256_shuffle_epi8(complements,index);
  const auto reversed=_mm256_permute2x128_si256(_mm256_shuffle_epi8(mapped,order),_mm256_shuffle_epi8(mapped,order),1);
  _mm256_storeu_si256(reinterpret_cast<__m256i*>(out.data()+s.size()-i-32),reversed);
 }
#endif
 for(;i<s.size();++i)out[s.size()-1-i]=complement(s[i]);
}
}
