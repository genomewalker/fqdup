cmake_minimum_required(VERSION 3.10)
project(fqdup)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

# Required dependencies
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# zstd — required for .fqcl v2 per-block compression
find_path(ZSTD_INCLUDE_DIR zstd.h)
find_library(ZSTD_LIB NAMES zstd)
if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIB)
    message(FATAL_ERROR "zstd not found. Install libzstd-dev or conda package 'zstd'")
endif()
message(STATUS "Found zstd: include=${ZSTD_INCLUDE_DIR} lib=${ZSTD_LIB}")

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

# htslib (optional - bgzf multi-threaded compressed output)
find_library(HTS_LIB NAMES hts)
find_path(HTS_INCLUDE_DIR htslib/bgzf.h)
if(HTS_LIB AND HTS_INCLUDE_DIR)
    message(STATUS "Found htslib: ${HTS_LIB} (bgzf multi-threaded output enabled)")
else()
    message(STATUS "htslib not found - using single-threaded zlib output")
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
if(ISAL_INCLUDE_DIR AND ISAL_LIB)
    set(LIBRAPIDARCHIVE_WITH_ISAL ON  CACHE BOOL "" FORCE)
else()
    set(LIBRAPIDARCHIVE_WITH_ISAL OFF CACHE BOOL "" FORCE)
endif()
FetchContent_Declare(
    librapidarchive
    GIT_REPOSITORY https://github.com/mxmlnkn/indexed_bzip2.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
    UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(librapidarchive)

# libtaph — taphonomy estimation library (deamination, 8-oxoG, PCR, length/GC,
# mixture model). Provides taph::FrameSelector::{update,finalize}_sample_profile
# and SampleDamageProfile with d_max_combined, all channels, GC bins.
FetchContent_Declare(
    libtaph
    GIT_REPOSITORY https://github.com/genomewalker/libtaph.git
    GIT_TAG        df47bb0fd5271d2bdfd5d8d487bc71536f4fa77c
    GIT_SHALLOW    FALSE
    UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(libtaph)

# BBHash (header-only MPHF — O(1) pass2 lookup for extend, ~3.7 bits/key)
FetchContent_Declare(
    bbhash
    GIT_REPOSITORY https://github.com/rizkg/BBHash.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
    UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(bbhash)

# Build fqdup
add_executable(fqdup
    src/main.cpp
    src/sort.cpp
    src/extend.cpp
    src/derep_pairs.cpp
    src/derep.cpp
    src/damage.cpp
    src/damage_profile.cpp
    src/gen.cpp
    src/fastq_io_backend.cpp
    src/logger.cpp
    src/cluster_format.cpp
    src/view.cpp
)

target_include_directories(fqdup PRIVATE
    ${XXHASH_INCLUDE_DIR}
    ${ZSTD_INCLUDE_DIR}
    ${flat_hash_map_SOURCE_DIR}
    ${bbhash_SOURCE_DIR}
)

target_link_libraries(fqdup PRIVATE
    ZLIB::ZLIB
    Threads::Threads
    ${XXHASH_LIB}
    ${ZSTD_LIB}
    rapidgzip::rapidgzip
    taph
)
target_compile_definitions(fqdup PRIVATE HAVE_RAPIDGZIP HAVE_BBHASH)

if(ISAL_INCLUDE_DIR AND ISAL_LIB)
    target_include_directories(fqdup PRIVATE ${ISAL_INCLUDE_DIR})
    target_link_libraries(fqdup PRIVATE ${ISAL_LIB})
endif()

if(HTS_LIB AND HTS_INCLUDE_DIR)
    target_include_directories(fqdup PRIVATE ${HTS_INCLUDE_DIR})
    target_link_libraries(fqdup PRIVATE ${HTS_LIB})
    target_compile_definitions(fqdup PRIVATE HAVE_BGZF)
endif()

if(JEMALLOC_LIB)
    target_link_libraries(fqdup PRIVATE ${JEMALLOC_LIB})
endif()

# Tests
enable_testing()

add_executable(gen_synthetic tests/gen_synthetic.cpp)
target_compile_options(gen_synthetic PRIVATE -O2)

add_executable(gen_genome_frags tests/gen_genome_frags.cpp)
target_compile_options(gen_genome_frags PRIVATE -O2)

add_executable(test_cluster_format tests/test_cluster_format.cpp src/cluster_format.cpp)
target_include_directories(test_cluster_format PRIVATE include ${ZSTD_INCLUDE_DIR})
target_link_libraries(test_cluster_format PRIVATE ${ZSTD_LIB})
target_compile_options(test_cluster_format PRIVATE -O2)
add_test(NAME cluster_format COMMAND test_cluster_format)

add_executable(test_syncmer_index tests/test_syncmer_index.cpp)
target_include_directories(test_syncmer_index PRIVATE include src ${flat_hash_map_SOURCE_DIR})
target_compile_options(test_syncmer_index PRIVATE -O2)
add_test(NAME syncmer_index COMMAND test_syncmer_index)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke.sh)
    add_test(NAME smoke
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke.sh $<TARGET_FILE:fqdup>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_damage_auto.sh)
    add_test(NAME damage_auto
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_damage_auto.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_synthetic>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_io_failures.sh)
    add_test(NAME io_failures
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_io_failures.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_synthetic>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_math_invariants.sh)
    add_test(NAME math_invariants
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_math_invariants.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_synthetic>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_backends.sh)
    add_test(NAME backends
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_backends.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_synthetic>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_extend_dmgcor.sh)
    add_test(NAME extend_dmgcor
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_extend_dmgcor.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_genome_frags>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_damage.sh)
    add_test(NAME damage
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_damage.sh
                     $<TARGET_FILE:fqdup>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_errcor.sh)
    add_test(NAME errcor
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_errcor.sh
                     $<TARGET_FILE:fqdup>)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_pipeline.sh)
    add_test(NAME pipeline
             COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_pipeline.sh
                     $<TARGET_FILE:fqdup> $<TARGET_FILE:gen_synthetic>)
endif()
