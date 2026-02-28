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
| Intel ISA-L | 3–5× faster gzip decompression | `--isal` flag; auto-detected |
| pigz | 2–3× faster gzip decompression | `--pigz` flag; must be in PATH |
| jemalloc | Better memory release | auto-detected; no flag needed |

ISA-L gives the largest benefit for `.gz` input on fast storage.

```bash
# conda (recommended)
conda install -c conda-forge isa-l
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
cmake .. 2>&1 | grep -E "ISA-L|jemalloc|xxhash"
# Found xxhash: include=/usr/include lib=/usr/lib/libxxhash.so
# Found ISA-L: /usr/lib/libisal.so (hardware-accelerated decompression enabled)
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

### ISA-L not found (optional)

ISA-L is optional — `fqdup` falls back to zlib without it.

```bash
conda install -c conda-forge isa-l

# if ISA-L is in a non-standard path
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
