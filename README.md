# fqdup

Ultra-fast, memory-efficient FASTQ deduplication for paired-end ancient DNA libraries.

## Overview

`fqdup` deduplicates paired-end FASTQ files with three subcommands forming a clean
pipeline:

```
fqdup sort → fqdup derep_pairs → fqdup derep
```

| Step | What it does | Memory |
|------|-------------|--------|
| `sort` | External merge sort by read ID — prerequisite for both dedup steps | Configurable (`--max-memory`) |
| `derep_pairs` | Two-pass paired deduplication: one representative pair per cluster | ~24 bytes/read pair |
| `derep` | Single-file deduplication with damage-aware hashing and PCR error correction | ~16 bytes/read |

The separation is intentional: `derep_pairs` handles structural complexity
(two paired files → one representative pair, selecting the longest non-extended
mate). `derep` handles biological complexity (ancient DNA damage masking, PCR
error correction). Each step can be used independently.

Three levels of deduplication are available:

| Mode | What it collapses |
|------|-------------------|
| Standard | Exact duplicates + reverse complements |
| Damage-aware | Above + reads differing only by terminal C→T / G→A deamination |
| + Error correction | Above + low-count clusters differing by 1 interior substitution (PCR errors) |

## Quick Start

```bash
# 1. Sort both input files by read ID
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G --fast
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G --fast

# 2. Paired deduplication — one representative pair per cluster
fqdup derep_pairs \
  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz

# 3. Single-file dedup — damage estimation + PCR error correction (both on by default)
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz
```

See the [wiki](../../wiki) for detailed tutorials, algorithm description, and
benchmarks.

## Installation

```bash
git clone https://github.com/genomewalker/fqdup.git
cd fqdup
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**Required:** CMake 3.10+, C++17 compiler, zlib, xxHash

```bash
# Ubuntu/Debian
sudo apt-get install cmake g++ zlib1g-dev libxxhash-dev

# RHEL/Rocky
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

Verify the build:

```bash
bash tests/smoke.sh $(pwd)/build/fqdup
# → OK: fqdup smoke test passed
```

> **Note:** Pass the binary as an absolute path — the test script changes directory internally.

## Options

### `fqdup sort`

```
fqdup sort -i INPUT -o OUTPUT --max-memory SIZE [options]

  -i FILE            Input FASTQ (.gz supported, or /dev/stdin)
  -o FILE            Output FASTQ (.gz for compression)
  --max-memory SIZE  Memory budget (e.g. 64G, 16G)
  -p THREADS         Sort threads (default: all cores)
  -t DIR             Temp directory for chunk files (default: .)
  -N                 Natural sort order (numeric read-ID suffixes)
  --fast             Uncompressed chunk intermediates (3× faster I/O, more disk)
```

### `fqdup derep_pairs`

```
fqdup derep_pairs -n NON -e EXT -o-non OUT -o-ext OUT [options]

Required:
  -n FILE      Sorted non-extended FASTQ
  -e FILE      Sorted extended FASTQ
  -o-non FILE  Output non-extended FASTQ
  -o-ext FILE  Output extended FASTQ

Optional:
  -c FILE                 Cluster statistics (gzipped TSV: hash, ext_len, pair_count, non_len)
  --no-revcomp            Disable reverse-complement collapsing (default: enabled)
  --allow-id-mismatch     Warn instead of failing when paired read IDs differ (default: error)
  --pigz                  Parallel decompression via pigz
  --isal                  Hardware-accelerated decompression (ISA-L)
```

Representative selection: the pair whose non-extended read is longest (captures the
most sequence information). No damage modeling or error correction at this step.

### `fqdup derep`

```
fqdup derep -i INPUT -o OUTPUT [options]

Required:
  -i FILE      Sorted input FASTQ (e.g. non-extended output of derep_pairs)
  -o FILE      Output deduplicated FASTQ

Optional:
  -c FILE              Cluster statistics (gzipped TSV: hash, seq_len, count)
  --no-revcomp         Disable reverse-complement collapsing (default: enabled)
  --no-damage          Disable damage estimation and masking (default: auto-enabled)
  --no-error-correct   Disable PCR error correction (default: enabled)
  --pigz               Parallel decompression via pigz
  --isal               Hardware-accelerated decompression (ISA-L)
```

Both `--damage-auto` and `--error-correct` are **on by default** — the primary
use case is ancient DNA where both are always appropriate. Use `--no-damage` or
`--no-error-correct` to disable them for non-aDNA datasets.

#### Damage-aware deduplication (default: on)

Reads differing only by terminal C→T (5') or G→A (3') deamination are collapsed.
A DART-style exponential model is fitted automatically or supplied manually.
If damage is below threshold, standard exact hashing is used automatically.

```
  --no-damage               Disable damage estimation and masking
  --damage-auto             Explicitly enable (already default)
  --damage-dmax   FLOAT     d_max for both 5' and 3' ends (manual model)
  --damage-dmax5  FLOAT     d_max for 5' end only
  --damage-dmax3  FLOAT     d_max for 3' end only
  --damage-lambda FLOAT     Decay constant for both ends
  --damage-lambda5 FLOAT    Decay constant for 5' end only
  --damage-lambda3 FLOAT    Decay constant for 3' end only
  --damage-bg     FLOAT     Background deamination rate (default: 0.02)
  --mask-threshold FLOAT    Mask positions where excess P(deamination) > T (default: 0.05)

  PCR thermocycling model:
  --pcr-cycles      INT     Number of PCR cycles
  --pcr-efficiency  FLOAT   Amplification efficiency per cycle, 0–1 (default: 1.0)
  --pcr-error-rate  FLOAT   Substitution rate in sub/base/doubling
                            Q5/HiFi: 5.3e-7 (default), Phusion: 3.9e-6, Taq: 1.5e-4
```

#### PCR error correction — Phase 3 (default: on)

After deduplication, low-count clusters differing from a high-count neighbour by
exactly one substitution **outside** the damage zone are classified as PCR errors
and removed. Uses a count-stratified 3-way pigeonhole algorithm with packed
2-bit Hamming verification (XOR bit-tricks; no full decode in the hot path).

Sequences containing ambiguous bases (N) are excluded from error correction to
prevent false absorptions. G↔T/C↔A (8-oxoG) substitutions are always protected;
C↔T/G↔A (deamination) additionally protected when damage mode is active.

```
  --no-error-correct      Disable Phase 3
  --error-correct         Explicitly enable (already default)
  --errcor-ratio  FLOAT   count(parent)/count(child) threshold to absorb (default: 50)
  --errcor-max-count INT  Only candidates with count ≤ N are children (default: 5)
  --errcor-bucket-cap INT Max pair-key bucket size (default: 64)
```

## Output Formats

### Cluster statistics (`-c FILE.tsv.gz`)

`derep_pairs` writes a 4-column gzipped TSV:

| Column | Description |
|--------|-------------|
| `hash` | 128-bit XXH3 canonical hash (32 hex chars) |
| `ext_len` | Extended read length of the representative |
| `pair_count` | Number of read pairs in the cluster |
| `non_len` | Non-extended read length of the representative |

`derep` writes a 3-column gzipped TSV:

| Column | Description |
|--------|-------------|
| `hash` | 128-bit XXH3 canonical hash (32 hex chars) |
| `seq_len` | Sequence length of the representative |
| `count` | Number of reads in the cluster |

## Algorithm Summary

```
fqdup sort:        chunk ingestion → parallel sort → k-way merge
fqdup derep_pairs: [Pass 1] stream (ext+non) → build index
                   [Pass 2] stream again → write representative pairs
fqdup derep:       [Pass 0] (opt) stride-sample → fit damage model
                   [Pass 1] stream → build index (with damage masking)
                   [Phase 3] (opt) PCR error correction
                   [Pass 2] stream → write unique representatives
```

**Canonical hash:** `min(XXH3_128(seq), XXH3_128(revcomp(seq)))` — collapses forward
and reverse-complement reads into the same cluster. Uses the full 128-bit hash as the
dedup key (collision probability ~3×10⁻²⁴ at 100 M reads).

**Damage masking:** terminal positions where excess deamination probability exceeds
`--mask-threshold` are replaced with a neutral byte before hashing. Masking is
symmetric (`max(P_CT(i), P_GA(i))`) so `hash(seq) == hash(revcomp(seq))` is
preserved after masking.

**PCR error protection:** substitutions classified as damage-consistent are never
absorbed by Phase 3 error correction. G↔T / C↔A (8-oxoG oxidative damage) is
always protected. C↔T / G↔A (ancient deamination) is additionally protected when
`--damage-auto` or manual damage parameters are active.

See the [wiki](../../wiki) for full algorithmic details.

## Performance

Benchmarked on a 25.8 M-read ancient DNA library (`a88af16f35`, ~65 bp mean length):

| Mode | Unique clusters | Δ vs standard |
|------|----------------|---------------|
| Standard (`derep_pairs` only) | 5,582,073 | — |
| Damage-aware (`derep --no-error-correct`) | 5,547,508 | −34,565 (−0.6%) |
| + Error correction (`derep`, default) | 5,532,327 | −49,746 (−0.9%) |

Wall time (NFS-mounted storage): sort ~20 s, derep_pairs ~25 s, derep ~33 s.
Throughput is I/O-bound (~40,000–50,000 reads/s). Use `--isal` or `--pigz` for
faster decompression.

Memory: ~24 bytes/read pair for `derep_pairs`; ~16 bytes/read for `derep` (plus
~0.25 bytes/bp of unique sequence for error correction — 2-bit packed arena).

## Citation

If you use `fqdup`, please cite:

> [citation placeholder]

## References

- Briggs et al. (2007) Patterns of damage in genomic DNA sequences from a Neandertal. *PNAS*. doi:10.1073/pnas.0704665104
- Jónsson et al. (2013) mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics*. doi:10.1093/bioinformatics/btt193
- Fernandez-Guerra et al. (2025) DART: damage-aware reference-based taxonomic profiling. *bioRxiv*.
- Potapov & Ong (2017) Examining sources of error in PCR by single-molecule sequencing. *PLOS ONE*. doi:10.1371/journal.pone.0169774
