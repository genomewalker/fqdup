#pragma once
#include <cstddef>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
namespace fqdup_fast_dp {
#if defined(__AVX2__)
// Eight full-matrix or seven banded interior cells of the +1/-1,
// gap-open 4 / gap-extend 1 suffix DP.
// Horizontal(j)=max(diag(j), Horizontal(j+1)-1). Three integer max scans
// enumerate within-tile extensions; full tiles include the right carry.
// Banded tiles mask off-band loads and never extend beyond seven cells.
// No choice/tie metadata is stored here. The original scalar traceback runs
// afterwards against identical scores. Valid lengths <=512 keep all arithmetic
// far from overflow. Boundary cells and partial tiles stay scalar.
template<bool banded,class Cell>
inline void suffix_tile(Cell* out,const Cell* below,char a,const char* b) {
 constexpr int neg=-1000000000;
 static_assert(sizeof(Cell)==16&&offsetof(Cell,suffix)==0&&offsetof(Cell,through)==12);
 const auto offsets=_mm256_setr_epi32(0,4,8,12,16,20,24,28);
 const auto zero=_mm256_setzero_si256();
 const auto absent=_mm256_set1_epi32(-2147482624);
 const auto equal=a=='N'?_mm_setzero_si128():_mm_cmpeq_epi8(_mm_set1_epi8(a),_mm_loadl_epi64(reinterpret_cast<const __m128i*>(b)));
 const auto substitution=_mm256_sub_epi32(_mm256_set1_epi32(-1),_mm256_slli_epi32(_mm256_cvtepi8_epi32(equal),1));
 __m256i diag,down;
 if constexpr(banded){
  // Below starts at (i+1,j+1), the next row's band slot zero. The
  // masked lanes never load the off-band vertical neighbor or slot seven.
  const auto valid=_mm256_setr_epi32(-1,-1,-1,-1,-1,-1,-1,0);
  diag=_mm256_add_epi32(_mm256_mask_i32gather_epi32(absent,&below[0].suffix[0],offsets,valid,4),substitution);
  diag=_mm256_blend_epi32(diag,absent,0x80);
  down=_mm256_mask_i32gather_epi32(_mm256_set1_epi32(neg),&below[0].suffix[1],
       _mm256_setr_epi32(-4,0,4,8,12,16,20,24),_mm256_setr_epi32(0,-1,-1,-1,-1,-1,-1,0),4);
 }else{
  diag=_mm256_add_epi32(_mm256_i32gather_epi32(&below[1].suffix[0],offsets,4),substitution);
  down=_mm256_i32gather_epi32(&below[0].suffix[1],offsets,4);
 }
 const auto one=_mm256_set1_epi32(1),four=_mm256_set1_epi32(4);
 const auto vertical=_mm256_max_epi32(diag,_mm256_sub_epi32(down,one));
 auto horizontal=diag;
 auto shifted=_mm256_sub_epi32(_mm256_permutevar8x32_epi32(horizontal,_mm256_setr_epi32(1,2,3,4,5,6,7,7)),one);
 horizontal=_mm256_max_epi32(horizontal,_mm256_blend_epi32(shifted,absent,0x80));
 shifted=_mm256_sub_epi32(_mm256_permutevar8x32_epi32(horizontal,_mm256_setr_epi32(2,3,4,5,6,7,7,7)),_mm256_set1_epi32(2));
 horizontal=_mm256_max_epi32(horizontal,_mm256_blend_epi32(shifted,absent,0xc0));
 shifted=_mm256_sub_epi32(_mm256_permutevar8x32_epi32(horizontal,_mm256_setr_epi32(4,5,6,7,7,7,7,7)),four);
 horizontal=_mm256_max_epi32(horizontal,_mm256_blend_epi32(shifted,absent,0xf0));
 __m256i right;
 if constexpr(banded){
  // No horizontal extension may leave the original seven-cell band.
  right=_mm256_permutevar8x32_epi32(horizontal,_mm256_setr_epi32(1,2,3,4,5,6,7,7));
 }else{
  const int carry=out[8].suffix[2];
  horizontal=_mm256_max_epi32(horizontal,_mm256_sub_epi32(_mm256_set1_epi32(carry),_mm256_setr_epi32(8,7,6,5,4,3,2,1)));
  right=_mm256_blend_epi32(_mm256_permutevar8x32_epi32(horizontal,_mm256_setr_epi32(1,2,3,4,5,6,7,7)),_mm256_set1_epi32(carry),0x80);
 }
 auto aligned=_mm256_max_epi32(diag,_mm256_set1_epi32(neg));
 aligned=_mm256_max_epi32(aligned,_mm256_sub_epi32(down,four));
 aligned=_mm256_max_epi32(aligned,_mm256_sub_epi32(right,four));
 // Transpose the three score vectors and a zero through vector into AoS cells.
 const auto ab0=_mm256_unpacklo_epi32(aligned,vertical),ab1=_mm256_unpackhi_epi32(aligned,vertical);
 const auto cz0=_mm256_unpacklo_epi32(horizontal,zero),cz1=_mm256_unpackhi_epi32(horizontal,zero);
 const auto r0=_mm256_unpacklo_epi64(ab0,cz0),r1=_mm256_unpackhi_epi64(ab0,cz0);
 const auto r2=_mm256_unpacklo_epi64(ab1,cz1),r3=_mm256_unpackhi_epi64(ab1,cz1);
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+0),_mm256_castsi256_si128(r0));
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+1),_mm256_castsi256_si128(r1));
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+2),_mm256_castsi256_si128(r2));
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+3),_mm256_castsi256_si128(r3));
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+4),_mm256_extracti128_si256(r0,1));
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+5),_mm256_extracti128_si256(r1,1));
 _mm_storeu_si128(reinterpret_cast<__m128i*>(out+6),_mm256_extracti128_si256(r2,1));
 if constexpr(!banded)_mm_storeu_si128(reinterpret_cast<__m128i*>(out+7),_mm256_extracti128_si256(r3,1));
}
template<class Cell>
inline void suffix_eight(Cell* out,const Cell* below,char a,const char* b){suffix_tile<false>(out,below,a,b);}
template<class Cell>
inline void suffix_seven(Cell* out,const Cell* below,char a,const char* b){suffix_tile<true>(out,below,a,b);}

#endif
}
