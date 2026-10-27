#include "raw_se_new_counter.hpp"
#include <new>
#include <cstdlib>
#include <limits>
void* operator new(std::size_t n){if(raw_se_alloc::active){if(raw_se_alloc::calls==std::numeric_limits<std::size_t>::max()||n>std::numeric_limits<std::size_t>::max()-raw_se_alloc::bytes)throw std::bad_alloc();++raw_se_alloc::calls;raw_se_alloc::bytes+=n;}if(void* p=std::malloc(n?n:1))return p;throw std::bad_alloc();}
void* operator new[](std::size_t n){return ::operator new(n);}void operator delete(void* p)noexcept{std::free(p);}void operator delete[](void* p)noexcept{std::free(p);}void operator delete(void* p,std::size_t)noexcept{std::free(p);}void operator delete[](void* p,std::size_t)noexcept{std::free(p);}
