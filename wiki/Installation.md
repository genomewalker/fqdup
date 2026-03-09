# Installation

## Requirements

### Required

| Dependency | Version | Purpose |
|-----------|---------|---------|
| CMake | 3.10+ | Build system |
| C++ compiler | GCC 7+ or Clang 5+ | C++17 |
| zlib | any | gzip I/O |
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

| Library | Benefit | Notes |
|---------|---------|-------|
| htslib | Multi-threaded bgzf output compression | auto-detected; no flag needed |
| jemalloc | Better memory release to OS | auto-detected; no flag needed |

Decompression uses rapidgzip (built-in, parallel, no external dependency). Output
compression uses bgzf (htslib) when `--threads > 1` and htslib is detected,
otherwise standard zlib gzip.

```bash
# conda (recommended)
conda install -c bioconda htslib
conda install -c conda-forge jemalloc
```

---

## Build from source

```bash
git clone https://github.com/genomewalker/fqdup.git
cd fqdup
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The binary is `build/fqdup`. Copy it wherever convenient or add `build/` to PATH.

### Check what was detected

```bash
cmake .. 2>&1 | grep -E "htslib|jemalloc|xxhash|rapidgzip"
# Found xxhash: include=/usr/include lib=/usr/lib/libxxhash.so
# Found htslib: /usr/lib/libhts.so (bgzf multi-threaded output enabled)
# jemalloc not found - using system allocator
```

---

## Verify

```bash
bash tests/smoke.sh build/fqdup
# → OK: fqdup smoke test passed
```

The smoke test exercises all three subcommands on a synthetic paired dataset
and checks sorting order, deduplication, and representative selection.

---

## Troubleshooting

### xxHash not found

```bash
# conda
conda install xxhash

# or point CMake to a manual install
cmake .. -DCMAKE_PREFIX_PATH=/path/to/xxhash
```

### htslib not found (optional)

htslib is optional — `fqdup` falls back to single-threaded zlib output without it.

```bash
conda install -c bioconda htslib

# if htslib is in a non-standard path
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
```

### CMake version too old

CMake 3.10+ is required:

```bash
conda install -c conda-forge cmake
```

### Compiler does not support C++17

GCC 7+ or Clang 5+ required:

```bash
conda install -c conda-forge cxx-compiler
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
cmake ..
```
