#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
namespace fqdup_fast_dp {
// Full suffix M scores are within [-512,512]. Unreachable I/D states use
// -16384; no reachable prefix can bring them near the best-minus-two cutoff.
struct ShortPlanarCells {
 static constexpr int neg=-16384;
 struct Scores {std::int16_t* values;std::size_t stride;std::int16_t& operator[](std::size_t s)const{return values[s*stride];}};
 struct Cell {Scores suffix;};
 std::int16_t* values=nullptr;std::size_t plane=0;
 void resize(std::size_t required){
  constexpr std::size_t maximum=513*513;
  if(required>maximum)throw std::length_error("full alignment scratch length");
  struct Buffer {std::unique_ptr<std::int16_t[]> values;std::size_t capacity=0;};
  static thread_local Buffer buffer;
  if(required>buffer.capacity){
   const auto capacity=std::min(maximum,std::max(required,std::max<std::size_t>(64,buffer.capacity*2)));
   std::unique_ptr<std::int16_t[]> next(new std::int16_t[3*capacity]);
   buffer.values=std::move(next);buffer.capacity=capacity;
  }
  values=buffer.values.get();plane=required;
 }
 std::int16_t* score(std::size_t state,std::size_t index){return values+state*plane+index;}
 Cell at(std::size_t index){return {{values+index,plane}};}
};
#if defined(__AVX2__)
template<int width>
inline void suffix_short_tail(ShortPlanarCells& cells,std::size_t out,std::size_t below,char a,const char* b){
 static_assert(width==4||width==8);
 const auto load=[](const std::int16_t* p){
  if constexpr(width==8)return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
  else return _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p));
 };
 __m128i letters;
 if constexpr(width==8)letters=_mm_loadl_epi64(reinterpret_cast<const __m128i*>(b));
 else{std::int32_t bytes;std::memcpy(&bytes,b,sizeof(bytes));letters=_mm_cvtsi32_si128(bytes);}
 const auto equal=a=='N'?_mm_setzero_si128():_mm_cmpeq_epi8(_mm_set1_epi8(a),letters);
 const auto substitution=_mm_sub_epi16(_mm_set1_epi16(-1),_mm_slli_epi16(_mm_cvtepi8_epi16(equal),1));
 const auto diag=_mm_add_epi16(load(cells.score(0,below+1)),substitution),down=load(cells.score(1,below));
 const auto one=_mm_set1_epi16(1),four=_mm_set1_epi16(4),absent=_mm_set1_epi16(-30000);
 const auto vertical=_mm_max_epi16(diag,_mm_sub_epi16(down,one));
 auto horizontal=diag;
 if constexpr(width==4)horizontal=_mm_blend_epi16(horizontal,absent,0xf0);
 horizontal=_mm_max_epi16(horizontal,_mm_sub_epi16(_mm_alignr_epi8(absent,horizontal,2),one));
 horizontal=_mm_max_epi16(horizontal,_mm_sub_epi16(_mm_alignr_epi8(absent,horizontal,4),_mm_set1_epi16(2)));
 if constexpr(width==8)horizontal=_mm_max_epi16(horizontal,_mm_sub_epi16(_mm_alignr_epi8(absent,horizontal,8),four));
 const auto carry=_mm_set1_epi16(*cells.score(2,out+width));
 horizontal=_mm_max_epi16(horizontal,_mm_sub_epi16(carry,_mm_setr_epi16(width,width-1,width-2,width-3,width-4,width-5,width-6,width-7)));
 auto next=horizontal;if constexpr(width==4)next=_mm_blend_epi16(next,carry,0xf0);
 const auto right=_mm_alignr_epi8(carry,next,2);
 auto aligned=_mm_max_epi16(diag,_mm_set1_epi16(ShortPlanarCells::neg));
 aligned=_mm_max_epi16(aligned,_mm_sub_epi16(down,four));
 aligned=_mm_max_epi16(aligned,_mm_sub_epi16(right,four));
 const auto store=[](std::int16_t* p,__m128i value){
  if constexpr(width==8)_mm_storeu_si128(reinterpret_cast<__m128i*>(p),value);
  else _mm_storel_epi64(reinterpret_cast<__m128i*>(p),value);
 };
 store(cells.score(0,out),aligned);store(cells.score(1,out),vertical);store(cells.score(2,out),horizontal);
}
template<int bytes>
inline __m256i shift_right_short(__m256i value,__m256i padding){
 const auto next=_mm256_permute2x128_si256(value,padding,0x21);
 if constexpr(bytes==16)return next;
 else return _mm256_alignr_epi8(next,value,bytes);
}
inline void suffix_sixteen_planar(ShortPlanarCells& cells,std::size_t out,std::size_t below,char a,const char* b){
 const auto equal=a=='N'?_mm_setzero_si128():_mm_cmpeq_epi8(_mm_set1_epi8(a),_mm_loadu_si128(reinterpret_cast<const __m128i*>(b)));
 const auto substitution=_mm256_sub_epi16(_mm256_set1_epi16(-1),_mm256_slli_epi16(_mm256_cvtepi8_epi16(equal),1));
 const auto diag=_mm256_add_epi16(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(cells.score(0,below+1))),substitution);
 const auto down=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(cells.score(1,below)));
 const auto one=_mm256_set1_epi16(1),four=_mm256_set1_epi16(4),absent=_mm256_set1_epi16(-30000);
 const auto vertical=_mm256_max_epi16(diag,_mm256_sub_epi16(down,one));
 auto horizontal=diag;
 horizontal=_mm256_max_epi16(horizontal,_mm256_sub_epi16(shift_right_short<2>(horizontal,absent),one));
 horizontal=_mm256_max_epi16(horizontal,_mm256_sub_epi16(shift_right_short<4>(horizontal,absent),_mm256_set1_epi16(2)));
 horizontal=_mm256_max_epi16(horizontal,_mm256_sub_epi16(shift_right_short<8>(horizontal,absent),four));
 horizontal=_mm256_max_epi16(horizontal,_mm256_sub_epi16(shift_right_short<16>(horizontal,absent),_mm256_set1_epi16(8)));
 const auto carry=_mm256_set1_epi16(*cells.score(2,out+16));
 horizontal=_mm256_max_epi16(horizontal,_mm256_sub_epi16(carry,_mm256_setr_epi16(16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)));
 const auto right=shift_right_short<2>(horizontal,carry);
 auto aligned=_mm256_max_epi16(diag,_mm256_set1_epi16(ShortPlanarCells::neg));
 aligned=_mm256_max_epi16(aligned,_mm256_sub_epi16(down,four));
 aligned=_mm256_max_epi16(aligned,_mm256_sub_epi16(right,four));
 _mm256_storeu_si256(reinterpret_cast<__m256i*>(cells.score(0,out)),aligned);
 _mm256_storeu_si256(reinterpret_cast<__m256i*>(cells.score(1,out)),vertical);
 _mm256_storeu_si256(reinterpret_cast<__m256i*>(cells.score(2,out)),horizontal);
}
#endif
}
