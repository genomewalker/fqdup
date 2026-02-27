# fqdup

Ultra-fast, memory-efficient paired-end FASTQ deduplication with ancient DNA damage awareness.

## Overview

`fqdup` deduplicates paired-end FASTQ files produced by ancient DNA (or any high-duplication) sequencing. It is a two-subcommand pipeline:

- **`fqdup sort`** — External merge sort by read ID (parallel, streaming, configurable memory)
- **`fqdup derep`** — Two-pass deduplication for sorted paired files (~16 bytes/read overhead)

Three progressive levels of deduplication are supported:

| Mode | What it collapses |
|------|-------------------|
| Standard | Exact sequence duplicates (+ reverse complements by default) |
| Damage-aware | Standard + C→T / G→A deamination variants at read termini |
| + Error correction | Above + low-count clusters differing by 1 interior substitution (PCR errors) |

## Typical Workflow

```bash
# Step 1: Sort both input files by read ID
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G

# Step 2a: Standard deduplication
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz

# Step 2b: With automatic damage estimation (ancient DNA)
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz \
  --damage-auto

# Step 2c: With damage awareness + PCR error correction
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz \
  --damage-auto \
  --error-correct
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

Required:
- CMake 3.10+
- C++17 compiler (GCC 7+, Clang 5+)
- zlib (`libz-dev` / `zlib-devel`)
- xxHash (`libxxhash-dev` / `xxhash-devel` / conda `xxhash`)

Optional (performance):
- Intel ISA-L — 4–6× faster decompression (`--isal`)
- pigz — parallel decompression (`--pigz`)
- jemalloc — better memory return to OS

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

Required:
  -n FILE     Sorted non-extended FASTQ
  -e FILE     Sorted extended FASTQ
  -o-non FILE Output non-extended FASTQ
  -o-ext FILE Output extended FASTQ

Optional:
  -c FILE          Cluster statistics (gzipped TSV)
  --no-revcomp     Disable reverse-complement matching (default: enabled)
  --pigz           Parallel decompression via pigz
  --isal           Hardware-accelerated decompression (ISA-L)
```

#### Ancient DNA damage-aware deduplication

Reads that differ only by C→T (5') or G→A (3') deamination — the hallmark
of ancient DNA damage — are collapsed into the same cluster. A DART-style
exponential decay model is fitted per-run or supplied manually.

```
  --damage-auto          Estimate damage parameters from input (Pass 0)
  --damage-stride INT    Sampling stride for Pass 0 (default: 1000)
                         Samples every Nth read uniformly across the sorted file.

  -- Manual damage model (alternative to --damage-auto):
  --damage-dmax   FLOAT  d_max for both 5' and 3' ends
  --damage-dmax5  FLOAT  d_max for 5' end only
  --damage-dmax3  FLOAT  d_max for 3' end only
  --damage-lambda FLOAT  Decay constant for both ends
  --damage-lambda5 FLOAT Decay constant for 5' end only
  --damage-lambda3 FLOAT Decay constant for 3' end only
  --damage-bg     FLOAT  Background deamination rate (default: 0.02)
  --mask-threshold FLOAT Mask positions where P(deamination) > T (default: 0.05)

  -- PCR model (informs which terminal positions to protect):
  --pcr-cycles      INT   Number of PCR cycles
  --pcr-efficiency  FLOAT PCR efficiency per cycle, 0–1 (default: 1.0)
  --pcr-error-rate  FLOAT Error rate in sub/base/doubling
                          Q5/HiFi: 5.3e-7 (default), Phusion: 3.9e-6, Taq: 1.5e-4
```

#### PCR error correction (Phase 3)

After standard deduplication, clusters with very low counts that differ from
a high-count neighbour by exactly one substitution **outside** the damage zone
are classified as PCR errors and removed. Uses a count-stratified 3-way
pigeonhole algorithm (O(N) average, AVX2-accelerated Hamming check).

```
  --error-correct        Enable PCR error duplicate removal
  --errcor-ratio  FLOAT  Min count(parent)/count(child) to absorb (default: 50)
  --errcor-max-count INT Only absorb clusters with count ≤ N (default: 5)
  --errcor-bucket-cap INT Max pair-key bucket size (default: 64)
                          Lower = faster on low-complexity regions
```

## Algorithm

### Sort

Single-threaded chunk ingestion → parallel chunk sorting → k-way merge.
Each chunk is compressed (or left plain with `--fast`) before hitting disk.

### Derep

**Pass 0** *(damage-aware only)*: Stride-sample the non-extended file, fit an
exponential decay model `P(deamination, pos) = d_max × exp(-λ × pos) + bg`
independently for the 5' and 3' ends. Positions where excess damage exceeds
`--mask-threshold` define the "damage zone" at each terminus.

**Pass 1**: Stream both sorted files in lockstep. For each read pair:
1. Optionally mask damage-zone bases to a neutral character.
2. Compute canonical hash: `min(XXH3(seq), XXH3(revcomp(seq)))`.
3. Insert into `hash → IndexEntry{record_index, non_len, count}`. On
   collision, keep the entry whose non-extended read is longest (best
   representative); increment count.

**Phase 3** *(error-correct only)*: Build a 3-way pair-key index of all
high-count clusters (count > `--errcor-max-count`). For each low-count
cluster, query three pair-key lookups covering all edit-distance-1 variants
in the interior (non-damage) region. Any child with a parent at count ratio ≥
`--errcor-ratio` and interior Hamming distance ≤ 1 is flagged as a PCR error.

**Pass 2**: Stream both files again. Write the representative record for each
cluster that was not flagged as a PCR error.

### Canonical hashing

```
hash(seq) = min(XXH3_64(seq), XXH3_64(revcomp(seq)))
```

With damage masking, deaminated positions are replaced with a neutral base
before hashing, using `max(P_CT(pos), P_GA(pos))` symmetrically so that
`hash(seq) == hash(revcomp(seq))` is preserved.

## Performance

Tested on a 25.8 M-read ancient DNA library (`a88af16f35`):

| Mode | Unique clusters | vs standard |
|------|----------------|-------------|
| Standard | 5,582,073 | — |
| Damage-aware | 5,547,508 | −34,565 (−0.6%) |
| Damage-aware + error-correct | 5,532,327 | −49,746 (−0.9%) |

Phase 3 error correction added ~2 s on top of a ~30 s total run.

Throughput is typically bottlenecked by I/O (~40,000–50,000 reads/s). Use
`--isal` or `--pigz` to accelerate decompression when reading `.gz` files.
