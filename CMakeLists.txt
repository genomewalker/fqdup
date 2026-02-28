cmake_minimum_required(VERSION 3.10)
project(fqdup)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

# Required dependencies
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# xxHash (XXH3) - required
find_path(XXHASH_INCLUDE_DIR xxhash.h)
find_library(XXHASH_LIB NAMES xxhash)
if(NOT XXHASH_INCLUDE_DIR OR NOT XXHASH_LIB)
    message(FATAL_ERROR "xxhash not found. Install libxxhash-dev, xxhash-devel, or conda package 'xxhash'")
endif()
message(STATUS "Found xxhash: include=${XXHASH_INCLUDE_DIR} lib=${XXHASH_LIB}")

# ISA-L (optional - hardware-accelerated decompression)
find_path(ISAL_INCLUDE_DIR isa-l.h PATH_SUFFIXES include)
find_library(ISAL_LIB NAMES isal)
if(ISAL_INCLUDE_DIR AND ISAL_LIB)
    message(STATUS "Found ISA-L: ${ISAL_LIB} (hardware-accelerated decompression enabled)")
    add_definitions(-DHAVE_ISAL)
else()
    message(STATUS "ISA-L not found - using zlib (slower decompression)")
endif()

# jemalloc (optional - better memory release to OS)
find_library(JEMALLOC_LIB NAMES jemalloc)
if(JEMALLOC_LIB)
    message(STATUS "Found jemalloc: ${JEMALLOC_LIB}")
else()
    message(STATUS "jemalloc not found - using system allocator")
endif()

# ska::flat_hash_map (header-only, fetched automatically)
include(FetchContent)
FetchContent_Declare(
    flat_hash_map
    GIT_REPOSITORY https://github.com/skarupke/flat_hash_map.git
    GIT_TAG master
)
FetchContent_MakeAvailable(flat_hash_map)

# rapidgzip (via indexed_bzip2) — in-process parallel gzip decompression.
# Replaces the pigz subprocess; no PATH dependency, true multi-threaded inflate.
set(LIBRAPIDARCHIVE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LIBRAPIDARCHIVE_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(LIBRAPIDARCHIVE_WITH_ISAL   OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    librapidarchive
    GIT_REPOSITORY https://github.com/mxmlnkn/indexed_bzip2.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
    UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(librapidarchive)

# libdart-damage — DART's full damage estimation pipeline as a shared library.
# Provides dart::FrameSelector::{update,finalize}_sample_profile and
# SampleDamageProfile with d_max_combined (mixture model), all channels, GC bins.
FetchContent_Declare(
    libdart_damage
    GIT_REPOSITORY https://github.com/genomewalker/libdart-damage.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
    UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(libdart_damage)

# Build fqdup
add_executable(fqdup
    src/main.cpp
    src/sort.cpp
    src/derep_pairs.cpp
    src/derep.cpp
    src/fastq_io_backend.cpp
    src/logger.cpp
)

target_include_directories(fqdup PRIVATE
    ${XXHASH_INCLUDE_DIR}
    ${flat_hash_map_SOURCE_DIR}
)

target_link_libraries(fqdup PRIVATE
    ZLIB::ZLIB
    Threads::Threads
    ${XXHASH_LIB}
    rapidgzip::rapidgzip
    dart-damage
)
target_compile_definitions(fqdup PRIVATE HAVE_RAPIDGZIP)

if(ISAL_INCLUDE_DIR AND ISAL_LIB)
    target_include_directories(fqdup PRIVATE ${ISAL_INCLUDE_DIR})
    target_link_libraries(fqdup PRIVATE ${ISAL_LIB})
endif()

if(JEMALLOC_LIB)
    target_link_libraries(fqdup PRIVATE ${JEMALLOC_LIB})
endif()

# Tests
enable_testing()

add_executable(gen_synthetic tests/gen_synthetic.cpp)
target_compile_options(gen_synthetic PRIVATE -O2)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke.sh)
    add_test(NAME smoke
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke.sh $<TARGET_FILE:fqdup>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_damage_auto.sh)
    add_test(NAME damage_auto
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_damage_auto.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_synthetic>)
endif()
