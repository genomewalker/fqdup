# Installation

## Requirements

### Required

| Dependency | Version | Purpose |
|-----------|---------|---------|
| CMake | 3.10+ | Build system |
| C++ compiler | GCC 7+ or Clang 5+ | C++17 support required |
| zlib | any | gzip read/write |
| xxHash | any | XXH3_64 hashing |

```bash
# Ubuntu / Debian
sudo apt-get install cmake g++ zlib1g-dev libxxhash-dev

# RHEL / Rocky Linux / CentOS
sudo yum install cmake gcc-c++ zlib-devel xxhash-devel

# macOS (Homebrew)
brew install cmake xxhash

# conda (any platform)
conda install -c conda-forge cmake cxx-compiler zlib xxhash
```

### Optional (performance)

| Library | Speedup | Notes |
|---------|---------|-------|
| Intel ISA-L | 4–6× decompression | `--isal` flag; detected automatically |
| pigz | 2–3× decompression | `--pigz` flag; must be in `PATH` |
| jemalloc | Better memory release | Detected automatically; no flag needed |

ISA-L provides the largest benefit for `.gz` input files (hardware-accelerated
inflate). It is the recommended choice when available.

```bash
# conda (recommended)
conda install -c conda-forge isa-l
```

---

## Build from Source

```bash
git clone https://github.com/genomewalker/fqdup.git
cd fqdup
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The binary is `build/fqdup`. No installation step is required — copy it to
wherever is convenient or add `build/` to `PATH`.

### CMake build options

CMake auto-detects ISA-L and jemalloc. To check what was found:

```bash
cmake .. 2>&1 | grep -E "ISA-L|jemalloc|xxhash"
# Found xxhash: include=/usr/include lib=/usr/lib/libxxhash.so
# Found ISA-L: /usr/lib/libisal.so (hardware-accelerated decompression enabled)
# jemalloc not found - using system allocator
```

---

## Verify

Run the built-in smoke test to confirm everything works:

```bash
bash tests/smoke.sh build/fqdup
# → [INFO] ...sort...
# → [INFO] ...derep_pairs...
# → [INFO] ...derep...
# → OK: fqdup smoke test passed
```

The smoke test exercises all three subcommands on a tiny synthetic paired
dataset and checks correctness of sorting order, deduplication, and
representative selection.

---

## Troubleshooting

### xxHash not found

```bash
# conda
conda install xxhash

# or point CMake to a manual install
cmake .. -DCMAKE_PREFIX_PATH=/path/to/xxhash
```

### ISA-L not found (optional)

ISA-L is optional. Without it, `fqdup` falls back to zlib (slower).

```bash
conda install -c conda-forge isa-l
```

If ISA-L is in a non-standard path, set `CMAKE_PREFIX_PATH`:

```bash
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
```

### CMake version too old

CMake 3.10 is the minimum. If the system CMake is older:

```bash
conda install -c conda-forge cmake
```

### Compiler does not support C++17

GCC 7+ or Clang 5+ is required. On older systems:

```bash
# conda provides a modern compiler
conda install -c conda-forge cxx-compiler
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
cmake ..
```
