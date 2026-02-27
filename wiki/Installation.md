# Installation

## Requirements

**Required:**
- CMake 3.10+
- C++17 compiler (GCC 7+, Clang 5+)
- zlib

```bash
# Ubuntu / Debian
sudo apt-get install cmake g++ zlib1g-dev libxxhash-dev

# RHEL / CentOS / Rocky
sudo yum install cmake gcc-c++ zlib-devel xxhash-devel

# conda
conda install cmake cxx-compiler zlib xxhash
```

**Optional (performance):**

| Library | Benefit | Flag |
|---------|---------|------|
| Intel ISA-L | 4–6× faster `.gz` decompression | `--isal` |
| pigz | 2–3× faster decompression | `--pigz` |
| jemalloc | Better memory return to OS | automatic |

To build with ISA-L support:
```bash
# Install ISA-L (conda recommended)
conda install -c conda-forge isa-l

cmake .. -DHAVE_ISAL=ON
```

## Build

```bash
git clone https://github.com/genomewalker/fqdup.git
cd fqdup
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The binary is `build/fqdup`.

## Verify

```bash
bash tests/smoke.sh build/fqdup
# → OK: fqdup smoke test passed
```
