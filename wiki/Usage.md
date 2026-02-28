# Usage

## Prerequisites

Both `derep_pairs` and `derep` require their inputs to be sorted by read ID.
Run `fqdup sort` first. For a paired library:

```bash
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G --fast
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G --fast
```

`--max-memory` controls how much RAM the sort uses for in-memory chunks; set
it to roughly 80% of available RAM. `--fast` keeps the chunk intermediates
uncompressed, which is about 3× faster on disk at the cost of more temporary
space. Omit it if disk is limited.

---

## Full ancient DNA pipeline

The typical workflow runs three commands in sequence.

### Step 1 — sort

Sort both files by read ID, as above.

### Step 2 — paired structural deduplication

```bash
fqdup derep_pairs \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz \
  -o-ext  ext.deduped.fq.gz \
  -c clusters_pairs.tsv.gz
```

`derep_pairs` reads both sorted files in lockstep, hashes the extended read,
and keeps one representative pair per cluster. The representative is the pair
with the longest non-extended read — it preserves the most sequence from reads
trimmed less aggressively during adapter removal. The `-c` flag writes per-cluster
statistics (hash, lengths, and counts for both mates) to a gzipped TSV.

For a 25.8 M-pair library, this step typically finishes in about 25 seconds
and reduces to around 5.6 M unique pairs.

### Step 3 — single-file damage-aware deduplication

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-auto \
  --error-correct \
  -c clusters_derep.tsv.gz
```

`--damage-auto` triggers Pass 0, which scans all reads to estimate the ancient
DNA damage parameters. These are logged before deduplication starts:

```
Pass 0: damage estimation — 5582073 reads processed (5582073 total)
  5'-end: d_max=0.099 lambda=0.246 bg=0.487
  3'-end: d_max=0.020 lambda=0.069 bg=0.509
  Masked positions: 1 (1 bp each end)
  Expected mismatches (L=91): 1.44, 99th-pct tolerance: 5
```

Reads that differ only at masked terminal positions then hash identically and
collapse into the same cluster — recovering pairs split by post-mortem
deamination.

`--error-correct` adds Phase 3 after Pass 1: any cluster with count ≤ 5 that
differs from a high-count cluster (≥ 50× its count) by exactly one substitution
outside the damage zone is absorbed as a PCR copying error.

On the same 5.6 M-read input, the full `--damage-auto --error-correct` run
completes in about 33 seconds and reduces to roughly 3.5 M unique clusters.

---

## Common invocations

### Standard deduplication

For non-ancient DNA, skip the damage and error-correction flags:

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz
```

### Structural deduplication only

If you want to deduplicate paired reads without the biological layer:

```bash
fqdup derep_pairs \
  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz
```

### Manual damage parameters

If damage parameters are already known from a mapDamage2 or DART run, supply
them directly and skip Pass 0:

```bash
fqdup derep \
  -i nonext.deduped.fq.gz -o nonext.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.15 --damage-lambda3 0.30 \
  --mask-threshold 0.05
```

### Strand-specific deduplication

By default, a read and its reverse complement hash to the same cluster
(appropriate for single-stranded libraries). Disable this for protocols where
strand matters:

```bash
fqdup derep_pairs ... --no-revcomp
fqdup derep       ... --no-revcomp
```

---

## All options

### `fqdup sort`

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Input FASTQ (.gz or plain) | required |
| `-o FILE` | Output FASTQ | required |
| `--max-memory SIZE` | Memory budget (e.g. `64G`, `16G`) | required |
| `-p N` | Sort threads | all cores |
| `-t DIR` | Temp directory for chunk files | `.` |
| `-N` | Natural sort (numeric read-ID suffixes) | off |
| `--fast` | Uncompressed chunk intermediates | off |

### `fqdup derep_pairs`

| Flag | Description | Default |
|------|-------------|---------|
| `-n FILE` | Sorted non-extended input FASTQ | required |
| `-e FILE` | Sorted extended input FASTQ | required |
| `-o-non FILE` | Output non-extended FASTQ | required |
| `-o-ext FILE` | Output extended FASTQ | required |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--no-revcomp` | Disable reverse-complement collapsing | off |

### `fqdup derep`

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Sorted input FASTQ | required |
| `-o FILE` | Output FASTQ | required |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--no-revcomp` | Disable reverse-complement collapsing | off |
| `--damage-auto` | Estimate damage parameters from data (Pass 0) | off |
| `--damage-dmax FLOAT` | d_max for both ends | — |
| `--damage-dmax5 FLOAT` | d_max for 5' end | — |
| `--damage-dmax3 FLOAT` | d_max for 3' end | — |
| `--damage-lambda FLOAT` | Decay constant for both ends | — |
| `--damage-lambda5 FLOAT` | Decay constant for 5' end | — |
| `--damage-lambda3 FLOAT` | Decay constant for 3' end | — |
| `--damage-bg FLOAT` | Background substitution rate | 0.02 |
| `--mask-threshold FLOAT` | Mask when excess P(deam) > T | 0.05 |
| `--pcr-cycles INT` | PCR cycles; if omitted, D_eff is estimated from duplication ratio | 0 (auto) |
| `--pcr-efficiency FLOAT` | Efficiency per cycle, 0–1 | 1.0 |
| `--pcr-error-rate FLOAT` | Sub/base/doubling (Q5=5.3e-7, Taq=1.5e-4) | 5.3e-7 |
| `--error-correct` | Enable Phase 3 PCR error correction | off |
| `--errcor-ratio FLOAT` | count(parent)/count(child) threshold | 50.0 |
| `--errcor-max-count INT` | Children must have count ≤ N | 5 |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 64 |
