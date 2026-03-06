# fqdup

FASTQ deduplication for paired-end ancient DNA libraries, with damage-aware
hashing and PCR error correction.

---

Ancient DNA deamination turns one original molecule into a cloud of C→T/G→A
sequence variants, so exact-sequence dedup splits true duplicates and inflates
unique read counts. `fqdup` collapses those damage variants into a single cluster
with damage-aware hashing, then optionally merges PCR-error variants with a fast
mismatch-tolerant search — so downstream complexity metrics reflect molecules,
not damage patterns.

## How it works

`fqdup` runs three steps in order, each targeting a distinct layer of the
duplication problem:

**1. `fqdup sort`** — sort both input files by read ID. Required by both
dedup steps.

**2. `fqdup derep_pairs`** — structural deduplication of paired reads.
Uses the Tadpole-extended sequence as the cluster fingerprint rather than the
raw merged read. Short aDNA fragments that share a merged-read core often
diverge in the assembled extension, so false collisions between different
molecules are greatly reduced. The representative kept per cluster is the pair
with the longest merged read.

**3. `fqdup derep`** — biological deduplication of the merged-read output.
Two mechanisms, both on by default:

- **Damage-aware hashing.** Terminal positions where the observed C→T or
  G→A rate exceeds a threshold are replaced with a neutral byte before hashing.
  Two reads from the same molecule that differ only at a deaminated position
  then produce the same hash and collapse. The mask is derived from the data
  itself (Pass 0) using the same exponential decay model as DART and mapDamage2.

- **PCR error correction (Phase 3).** After the index is built, a 3-way
  pigeonhole Hamming search finds clusters with low count that differ from a
  high-count cluster by exactly one interior substitution. These are PCR copying
  errors. Crucially, C↔T and G↔A mismatches are *never* absorbed — they are
  indistinguishable from damage signal and are always kept.

Use `--no-damage` or `--no-error-correct` to disable either mechanism. For
non-aDNA data, use both.

---

## Quick start

### Full ancient DNA pipeline

```bash
# 1. Sort merged and Tadpole-extended reads by read ID
fqdup sort -i merged.fq.gz   -o merged.sorted.fq.gz   --max-memory 64G --fast
fqdup sort -i extended.fq.gz -o extended.sorted.fq.gz --max-memory 64G --fast

# 2. Paired dedup — one representative pair per cluster
fqdup derep_pairs \
  -n merged.sorted.fq.gz \
  -e extended.sorted.fq.gz \
  -o-non merged.deduped.fq.gz \
  -o-ext  extended.deduped.fq.gz

# 3. Damage-aware dedup + PCR error correction (both on by default)
fqdup derep \
  -i merged.deduped.fq.gz \
  -o merged.final.fq.gz
```

### Non-aDNA (modern DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz \
  --no-damage --no-error-correct
```

---

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

**Optional:** jemalloc (better memory return to OS; auto-detected, no flag needed).

Verify:

```bash
bash tests/smoke.sh $(pwd)/build/fqdup
# → OK: fqdup smoke test passed
```

> Pass the binary as an absolute path — the test script changes directory internally.

---

## Options

### `fqdup sort`

```
fqdup sort -i INPUT -o OUTPUT --max-memory SIZE [options]

  -i FILE            Input FASTQ (.gz supported)
  -o FILE            Output sorted FASTQ
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
  -n FILE      Sorted merged (fastp) FASTQ
  -e FILE      Sorted Tadpole-extended FASTQ
  -o-non FILE  Output merged FASTQ (representatives)
  -o-ext FILE  Output extended FASTQ (representatives)

Optional:
  -c FILE        Cluster statistics (gzipped TSV)
  --no-revcomp   Disable reverse-complement collapsing (default: enabled)
```

Representative selection: pair with the longest merged read.

### `fqdup derep`

```
fqdup derep -i INPUT -o OUTPUT [options]

Required:
  -i FILE      Sorted input FASTQ
  -o FILE      Output FASTQ

Optional:
  -c FILE              Cluster statistics (gzipped TSV)
  --no-revcomp         Disable reverse-complement collapsing (default: enabled)
```

Both `--damage-auto` and `--error-correct` are **on by default**. Use the
`--no-*` flags to disable them for non-aDNA datasets.

#### Damage-aware hashing (default: on)

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
  --mask-threshold FLOAT    Mask positions where excess P(deamination) > T
                            (default: 0.05)
```

If both `d_max_5` and `d_max_3` are below 0.02 after estimation, damage is
treated as negligible and standard exact hashing is used — no user intervention
required.

#### PCR error correction — Phase 3 (default: on)

```
  --no-error-correct      Disable Phase 3
  --error-correct         Explicitly enable (already default)
  --errcor-ratio  FLOAT   count(parent) / count(child) threshold (default: 50)
  --errcor-max-count INT  Only absorb clusters with count ≤ N (default: 5)
  --errcor-bucket-cap INT Max pair-key bucket size (default: 64)

  PCR model (for adaptive thresholds):
  --pcr-cycles      INT   Number of PCR cycles
  --pcr-efficiency  FLOAT Amplification efficiency per cycle, 0–1 (default: 1.0)
  --pcr-error-rate  FLOAT Sub/base/doubling — Q5: 5.3e-7 (default),
                          Phusion: 3.9e-6, Taq: 1.5e-4
```

---

## Output formats

### Cluster statistics (`-c FILE.tsv.gz`)

`derep_pairs` writes a 4-column gzipped TSV:

| Column | Description |
|--------|-------------|
| `hash` | 128-bit XXH3 canonical hash (32 hex chars) |
| `ext_len` | Extended read length of the representative |
| `pair_count` | Number of read pairs in the cluster |
| `non_len` | Merged read length of the representative |

`derep` writes a 3-column gzipped TSV:

| Column | Description |
|--------|-------------|
| `hash` | 128-bit XXH3 canonical hash (32 hex chars) |
| `seq_len` | Sequence length of the representative |
| `count` | Number of reads in the cluster |

---

## Performance

Benchmarked on `a88af16f35` — 25.8 M fastp-merged reads, ~91 bp mean length,
`d_max_5 ≈ 0.19`, NFS-mounted storage:

| Step | Unique clusters | Wall time |
|------|----------------|-----------|
| `derep_pairs` (25.8 M in) | 5,582,073 | ~25 s |
| `derep` standard | 3,531,821 | ~22 s |
| `derep --no-error-correct` | 3,511,607 | ~31 s |
| `derep` (default) | 3,510,151 | ~33 s |

Across 31 ancient DNA libraries (2.97 B reads total), `derep` with default
settings absorbed a further 8.6% of reads beyond `derep_pairs`.

**Memory:**

| Component | Formula | Example |
|-----------|---------|---------|
| `derep_pairs` index | ~40 bytes × N_unique pairs | 5 M pairs → ~200 MB |
| `derep` index | ~40 bytes × N_unique clusters | 3.5 M clusters → ~140 MB |
| `derep` SeqArena (2-bit packed) | ~0.25 × L_avg bytes × N_unique | 3.5 M × 91 bp → ~80 MB |

The SeqArena uses 2-bit encoding (A/C/G/T → 2 bits), ~4× more compact than
ASCII. Sequences containing N are stored but skipped during Phase 3.

---

## Algorithm summary

```
fqdup sort:        chunk ingestion → parallel sort → k-way merge
fqdup derep_pairs: [Pass 1] stream (ext + non) → build index
                   [Pass 2] stream again → write representative pairs
fqdup derep:       [Pass 0] stride-sample → fit damage model (--damage-auto)
                   [Pass 1] stream → build index with damage masking
                   [Phase 3] 3-way pigeonhole PCR error correction
                   [Pass 2] stream → write unique representatives
```

**Canonical hash:** `min(XXH3_128(seq), XXH3_128(revcomp(seq)))` — collapses
forward and reverse-complement reads into the same cluster. Collision probability
~3×10⁻²⁴ at 100 M reads.

**Damage masking:** empirical per-position mask derived from observed T/(T+C)
and A/(A+G) frequencies; symmetric masking preserves
`hash(seq) == hash(revcomp(seq))` after masking.

**Phase 3 protection:** G↔T and C↔A (8-oxoG) are always protected. C↔T and
G↔A (deamination) are additionally protected when damage mode is active. Only
A↔T and C↔G transversions are eligible for absorption.

See the [wiki](../../wiki) for detailed algorithm descriptions and benchmarks.

---

## Citation

If you use `fqdup`, please cite:

> [citation placeholder]

## References

- Briggs et al. (2007) Patterns of damage in genomic DNA sequences from a Neandertal. *PNAS*. doi:10.1073/pnas.0704665104
- Jónsson et al. (2013) mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics*. doi:10.1093/bioinformatics/btt193
- Fernandez-Guerra et al. (2025) DART: damage-aware reference-based taxonomic profiling. *bioRxiv*.
- Potapov & Ong (2017) Examining sources of error in PCR by single-molecule sequencing. *PLOS ONE*. doi:10.1371/journal.pone.0169774
