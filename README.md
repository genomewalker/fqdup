# fqdup

Ultra-fast, memory-efficient paired-end FASTQ deduplication.

## Overview

`fqdup` is a two-subcommand tool for deduplicating paired-end FASTQ files from ancient DNA or other high-duplication sequencing:

- **`fqdup sort`** — External merge sort by read ID (parallel, streaming, low memory)
- **`fqdup derep`** — Two-pass deduplication for sorted paired files (~16 bytes/read)

## Typical Workflow

```bash
# Step 1: Sort both input files by read ID
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G

# Step 2: Deduplicate
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz \
  [-c clusters.tsv.gz]
```

## Memory

| Step | Memory Usage |
|------|-------------|
| `sort` | Configurable via `--max-memory` |
| `derep` | ~16 bytes per input read (e.g. 6.4 GB for 400M reads) |

## Build

```bash
mkdir build && cd build
cmake ..
make
```

### Dependencies

- CMake 3.10+
- C++17 compiler
- zlib (`libz-dev` / `zlib-devel`)
- xxHash (`libxxhash-dev` / `xxhash-devel` / conda `xxhash`)
- Optional: Intel ISA-L for 4-6× faster decompression (`--isal`)
- Optional: pigz for parallel compression/decompression (`--pigz`)
- Optional: jemalloc for better memory release to OS

## Options

### `fqdup sort`

```
fqdup sort -i INPUT -o OUTPUT --max-memory SIZE [options]

  -i FILE            Input FASTQ (.gz supported, or /dev/stdin)
  -o FILE            Output FASTQ (.gz for compression)
  --max-memory SIZE  Memory limit (e.g. 64G, 32G)
  -p THREADS         Sorting threads (default: auto)
  -t DIR             Temp directory (default: .)
  -N                 Natural sort order (numeric suffixes)
  --fast             Uncompressed intermediates (3× faster, more disk)
```

### `fqdup derep`

```
fqdup derep -n NON -e EXT -o-non OUT_NON -o-ext OUT_EXT [options]

  -n FILE     Sorted non-extended FASTQ
  -e FILE     Sorted extended FASTQ
  -o-non FILE Output non-extended FASTQ
  -o-ext FILE Output extended FASTQ
  -c FILE     Cluster statistics (gzipped TSV)
  --no-revcomp  Disable reverse-complement matching
  --pigz      Parallel decompression via pigz
  --isal      Hardware-accelerated decompression (ISA-L)
```

## Algorithm

**Sort**: Single-threaded chunk ingestion → parallel chunk sorting → k-way merge.

**Derep**:
- Pass 1: Stream paired files, store `hash → (record_index, non_len, count)` (~16 bytes/read)
- Pass 2: Stream again, write representative records for each unique hash

The representative is the read pair where the non-extended read is longest. Reverse-complement hashing is applied by default (`min(hash(seq), hash(revcomp(seq)))`).
