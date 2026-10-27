#pragma once
#include <cstddef>
namespace raw_se_alloc {inline thread_local bool active=false;inline thread_local std::size_t calls=0,bytes=0;inline void start(){calls=bytes=0;active=true;}inline void stop(){active=false;}}
