#!/bin/bash
# Fresh configure, build and ctest of fqdup against a given libtaph checkout. Writes only inside
# <build_dir> (which must not exist yet). The checkout must be at the SHA in LIBTAPH_PIN and clean,
# or configure fails by design. Last stdout line: "fqdup ctest: N/N passed".
# Usage: scripts/libtaph-consumer-test.sh <build_dir> <libtaph_dir>
# Environment: FQDUP_BIOINFO (conda prefix with zlib/zstd/xxhash/isa-l/htslib), FQDUP_CXX (g++),
# FQDUP_DEPS (directory holding flat_hash_map-src and bbhash-src), FQDUP_JOBS (build parallelism).
set -euo pipefail
[ $# -eq 2 ] || { echo "usage: $0 <build_dir> <libtaph_dir>" >&2; exit 2; }
build=$1; libtaph=$(cd "$2" && pwd)
src=$(cd "$(dirname "$0")/.." && pwd)
[ ! -e "$build" ] || { echo "build dir $build already exists" >&2; exit 2; }
bioinfo=${FQDUP_BIOINFO:-/opt/conda/envs/bioinfo}
cxx=${FQDUP_CXX:-/opt/software/gcc/13.2.0/bin/g++}
deps=${FQDUP_DEPS:-/path/to/tmp/fqdup-qc-terminal-anchor-20260928/build-release/_deps}
jobs=${FQDUP_JOBS:-8}
export LC_ALL=C.UTF-8 LANG=C.UTF-8  # the Python tests print non-ASCII
export LD_LIBRARY_PATH="$bioinfo/lib:$(dirname "$(dirname "$cxx")")/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
mkdir -p "$build"
cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" -DCMAKE_PREFIX_PATH="$bioinfo" \
    -DFETCHCONTENT_SOURCE_DIR_LIBTAPH="$libtaph" \
    -DFETCHCONTENT_SOURCE_DIR_FLAT_HASH_MAP="$deps/flat_hash_map-src" \
    -DFETCHCONTENT_SOURCE_DIR_BBHASH="$deps/bbhash-src" \
    -DFQDUP_GIT_COMMIT="$(git -C "$src" rev-parse HEAD)" > "$build/configure.log" 2>&1 \
    || { tail -20 "$build/configure.log"; exit 1; }
cmake --build "$build" -j "$jobs" > "$build/build.log" 2>&1 || { tail -30 "$build/build.log"; exit 1; }
"$build/fqdup" --version
(cd "$build" && ctest --output-on-failure > ctest.log 2>&1) || { tail -40 "$build/ctest.log"; exit 1; }
summary=$(grep -E 'tests passed' "$build/ctest.log" | tail -1)
passed=$(sed -E 's/.* ([0-9]+) tests failed out of ([0-9]+).*/\2/' <<<"$summary")
echo "fqdup ctest: $passed/$passed passed"
