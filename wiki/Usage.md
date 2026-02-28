# Usage

## Prerequisites

Both `derep_pairs` and `derep` require sorted input. Sort first:

```bash
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G --fast
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G --fast
```

`--fast` keeps chunk intermediates uncompressed during sorting (~3× faster I/O
at the cost of more temporary disk space). Omit it if disk space is limited.

---

## Tutorial: full ancient DNA pipeline

### Step 1 — sort

```bash
fqdup sort \
  -i nonext.fq.gz \
  -o nonext.sorted.fq.gz \
  --max-memory 64G \
  --fast
```

Repeat for the extended file:

```bash
fqdup sort \
  -i ext.fq.gz \
  -o ext.sorted.fq.gz \
  --max-memory 64G \
  --fast
```

Set `--max-memory` to about 80% of available RAM. The sort uses all cores by
default; limit with `-p N` if needed.

### Step 2 — paired deduplication

```bash
fqdup derep_pairs \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz \
  -o-ext  ext.deduped.fq.gz \
  -c clusters_pairs.tsv.gz
```

Both sorted files are read in lockstep. The extended read is hashed; the
representative pair is the one with the longest non-extended read. The `-c`
flag writes cluster statistics:

```
hash            ext_len  ext_count  non_len  non_count
a3f2b1c4...     65       4          58       4
...
```

### Step 3 — single-file damage-aware deduplication

Run `derep` on the non-extended output from Step 2:

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-auto \
  --error-correct \
  -c clusters_derep.tsv.gz
```

`--damage-auto` runs Pass 0: all reads are scanned, per-position C→T and G→A
frequencies are measured, and an exponential decay model is fitted. The
estimated parameters are logged:

```
Pass 0: damage estimation — 5582073 reads processed (5582073 total)
  5'-end: d_max=0.099 lambda=0.246 bg=0.487
  3'-end: d_max=0.020 lambda=0.069 bg=0.509
  --- Damage-Aware Deduplication ---
  5'-end d_max:  0.098842
  Mask threshold:0.050000
  Masked positions: 1 (1 bp each end)
  Expected mismatches (L=91): 1.44, 99th-pct tolerance: 5
```

Reads differing only within the estimated damage zone hash identically.

`--error-correct` adds Phase 3: clusters with count ≤ 5 that differ from a
high-count cluster (count ≥ 50×) by exactly 1 substitution outside the damage
zone are absorbed as PCR errors.

---

## Common invocations

### Standard deduplication (no ancient DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz
```

### Paired deduplication only

```bash
fqdup derep_pairs \
  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz
```

### Damage-aware with manual parameters

If damage parameters are already known from a mapDamage2 or DART run:

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.15 --damage-lambda3 0.30 \
  --damage-bg 0.02 \
  --mask-threshold 0.05
```

### PCR error model

Specify polymerase and cycle count to compute the expected number of PCR-induced
mismatches (informational — does not change the masking, only the reported
tolerance):

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-auto \
  --pcr-cycles 30 \
  --pcr-error-rate 5.3e-7     # Q5 polymerase (default)
```

For Phusion use `--pcr-error-rate 3.9e-6`; for Taq use `1.5e-4`.

### Fast decompression

```bash
# ISA-L (hardware-accelerated, fastest)
fqdup derep_pairs ... --isal
fqdup derep       ... --isal

# pigz (parallel, widely available)
fqdup derep_pairs ... --pigz
fqdup derep       ... --pigz
```

### Reverse-complement collapsing

By default, a sequence and its reverse complement hash to the same cluster.
This is appropriate for single-stranded library protocols. Disable with
`--no-revcomp` for double-stranded protocols where strand matters:

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
| `--pigz` | Parallel decompression via pigz | off |
| `--isal` | Hardware-accelerated decompression (ISA-L) | off |

### `fqdup derep`

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Sorted input FASTQ | required |
| `-o FILE` | Output FASTQ | required |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--no-revcomp` | Disable reverse-complement collapsing | off |
| `--pigz` | Parallel decompression via pigz | off |
| `--isal` | Hardware-accelerated decompression (ISA-L) | off |
| `--damage-auto` | Auto-estimate damage parameters (Pass 0) | off |
| `--damage-dmax FLOAT` | d_max for both ends | — |
| `--damage-dmax5 FLOAT` | d_max for 5' end | — |
| `--damage-dmax3 FLOAT` | d_max for 3' end | — |
| `--damage-lambda FLOAT` | Decay constant for both ends | — |
| `--damage-lambda5 FLOAT` | Decay constant for 5' end | — |
| `--damage-lambda3 FLOAT` | Decay constant for 3' end | — |
| `--damage-bg FLOAT` | Background substitution rate | 0.02 |
| `--mask-threshold FLOAT` | Mask when excess P(deam) > T | 0.05 |
| `--pcr-cycles INT` | PCR cycles (informational) | 0 |
| `--pcr-efficiency FLOAT` | Efficiency per cycle, 0–1 | 1.0 |
| `--pcr-error-rate FLOAT` | Sub/base/doubling (Q5=5.3e-7, Taq=1.5e-4) | 5.3e-7 |
| `--error-correct` | Enable Phase 3 PCR error correction | off |
| `--errcor-ratio FLOAT` | count(parent)/count(child) ≥ ratio | 50.0 |
| `--errcor-max-count INT` | Children must have count ≤ N | 5 |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 64 |
