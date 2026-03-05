# Usage

## Prerequisites

### Upstream steps (outside fqdup)

fqdup operates on reads that have been processed by two upstream tools:

**1. fastp merge** — collapse overlapping R1+R2 pairs into single sequences.
Each merged read represents one complete ancient DNA molecule:

```bash
fastp --merge --merged_out merged.fq.gz \
      --in1 R1.fq.gz --in2 R2.fq.gz \
      --disable_adapter_trimming  # (or with adapter trimming, as appropriate)
```

**2. Tadpole** (BBTools) — extend each merged read from both ends using de Bruijn
graph assembly. The extended reads serve as deduplication fingerprints; the
original merged reads are kept as the actual output:

```bash
tadpole.sh in=merged.fq.gz out=extended.fq.gz mode=extend el=50 er=50
```

### Sort

Both files must be sorted by read ID before deduplication:

```bash
fqdup sort -i merged.fq.gz   -o merged.sorted.fq.gz   --max-memory 64G --fast
fqdup sort -i extended.fq.gz -o extended.sorted.fq.gz --max-memory 64G --fast
```

`--max-memory` controls RAM for in-memory sort chunks; set to ~80% of available
RAM. `--fast` keeps chunk intermediates uncompressed (~3× faster, more disk).
Omit it if disk is limited.

---

## Full ancient DNA pipeline

The typical workflow runs three commands in sequence.

### Step 1 — sort

Sort both files by read ID, as above.

### Step 2 — paired structural deduplication

```bash
fqdup derep_pairs \
  -n merged.sorted.fq.gz \
  -e extended.sorted.fq.gz \
  -o-non merged.deduped.fq.gz \
  -o-ext  extended.deduped.fq.gz \
  -c clusters_pairs.tsv.gz
```

`derep_pairs` reads both sorted files in lockstep and hashes the **extended**
(Tadpole-assembled) read as the cluster fingerprint. Using the extended read
reduces false collisions: two different molecules that happen to share a short
merged sequence will typically diverge in the assembled extension.

The representative kept per cluster is the pair with the longest **merged**
read — preserving the most original ancient DNA sequence. The `-c` flag writes
per-cluster statistics to a gzipped TSV.

For a 25.8 M-pair library, this step typically finishes in about 25 seconds
and reduces to around 5.6 M unique pairs.

### Step 3 — single-file damage-aware deduplication

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  -c clusters_derep.tsv.gz
```

Damage estimation and PCR error correction are **on by default**. Pass 0 scans
all reads to estimate the ancient DNA damage parameters, which are logged before
deduplication starts:

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

Phase 3 (PCR error correction) runs after Pass 1: any cluster with count ≤ 5
that differs from a high-count cluster (≥ 50× its count) by exactly one
substitution outside the damage zone is absorbed as a PCR copying error.

On the same 5.6 M-read input, the full run completes in about 33 seconds and
reduces to roughly 3.5 M unique clusters.

---

## Common invocations

### Standard deduplication (non-aDNA)

For modern DNA, disable damage estimation and error correction explicitly:

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz --no-damage --no-error-correct
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
| `-n FILE` | Sorted merged (fastp) FASTQ | required |
| `-e FILE` | Sorted Tadpole-extended FASTQ | required |
| `-o-non FILE` | Output merged FASTQ (representatives) | required |
| `-o-ext FILE` | Output extended FASTQ (representatives) | required |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--no-revcomp` | Disable reverse-complement collapsing | off |

### `fqdup derep`

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Sorted input FASTQ | required |
| `-o FILE` | Output FASTQ | required |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--no-revcomp` | Disable reverse-complement collapsing | off |
| `--damage-auto` | Estimate damage parameters from data (Pass 0) | **on** |
| `--no-damage` | Disable damage estimation and masking | — |
| `--damage-dmax FLOAT` | d_max for both ends (manual model) | — |
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
| `--error-correct` | Enable Phase 3 PCR error correction | **on** |
| `--no-error-correct` | Disable Phase 3 PCR error correction | — |
| `--errcor-ratio FLOAT` | count(parent)/count(child) threshold | 50.0 |
| `--errcor-max-count INT` | Children must have count ≤ N | 5 |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 64 |
