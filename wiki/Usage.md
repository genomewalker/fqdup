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

## Tutorial: Full ancient DNA pipeline

This tutorial walks through the full three-step workflow on a paired-end library
with expected ancient DNA damage.

### Step 1 — Sort

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

Memory is the main parameter. Set `--max-memory` to about 80% of available RAM
to leave headroom. The sort is parallel (uses all cores by default); limit with
`-p N` if needed.

### Step 2 — Paired deduplication

```bash
fqdup derep_pairs \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz \
  -o-ext  ext.deduped.fq.gz \
  -c clusters_pairs.tsv.gz
```

This reads both sorted files in lockstep, hashes the **extended** read sequence,
and selects one representative pair per cluster. The representative is the pair
with the **longest non-extended read** — this recovers as much non-extended
sequence as possible when duplicate reads differ in trimming length.

The `-c` flag writes cluster statistics (one row per unique cluster):

```
hash            ext_len  ext_count  non_len  non_count
a3f2b1c4...     65       4          58       4
...
```

### Step 3 — Single-file damage-aware deduplication

Run `derep` on the non-extended output from Step 2:

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-auto \
  --error-correct \
  -c clusters_derep.tsv.gz
```

`--damage-auto` runs Pass 0: it stride-samples the input file (every 1000th read
by default), fits an exponential decay model to the observed C→T and G→A
frequencies, and logs the estimated parameters:

```
Pass 0: damage estimation — sampled 100000 reads (every 1000th) ...
  5'-end: d_max=0.071 lambda=0.290 bg=0.301
  3'-end: d_max=0.011 lambda=0.250 bg=0.295
  --- Damage-Aware Deduplication ---
  5'-end d_max:  0.071000
  5'-end lambda: 0.290000
  Mask threshold:0.050000
  Expected mismatches (L=65): 0.51, 99th-pct tolerance: 3
```

Reads whose sequences differ only within the estimated damage zone will now hash
to the same cluster instead of being counted as distinct.

`--error-correct` adds Phase 3: any cluster with count ≤ 5 that differs from a
high-count cluster (count ≥ 50×) by exactly 1 substitution outside the damage
zone is absorbed as a PCR error.

---

## Common Invocations

### Standard deduplication (no ancient DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz
```

### Paired deduplication only (no biological step)

```bash
fqdup derep_pairs \
  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz
```

### Damage-aware with manual parameters

If you know the damage parameters from a previous mapDamage or DART run:

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.15 --damage-lambda3 0.30 \
  --damage-bg 0.02 \
  --mask-threshold 0.05
```

### Damage-aware with PCR model

If you have PCR parameters for your library preparation:

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
# With ISA-L (hardware-accelerated, fastest)
fqdup derep_pairs ... --isal
fqdup derep       ... --isal

# With pigz (parallel, widely available)
fqdup derep_pairs ... --pigz
fqdup derep       ... --pigz
```

### Reverse-complement collapsing disabled

By default, `fqdup` collapses a sequence and its reverse complement into the same
cluster. This is appropriate for single-stranded library protocols. Disable with
`--no-revcomp` for double-stranded protocols where strand matters:

```bash
fqdup derep_pairs ... --no-revcomp
fqdup derep       ... --no-revcomp
```

---

## All Options

### `fqdup sort`

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Input FASTQ (.gz or plain, or `/dev/stdin`) | required |
| `-o FILE` | Output FASTQ (.gz for compression) | required |
| `--max-memory SIZE` | Memory budget (e.g. `64G`, `16G`, `8G`) | required |
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
| `--damage-stride INT` | Sampling stride for Pass 0 | 1000 |
| `--damage-dmax FLOAT` | d_max for both ends | — |
| `--damage-dmax5 FLOAT` | d_max for 5' end | — |
| `--damage-dmax3 FLOAT` | d_max for 3' end | — |
| `--damage-lambda FLOAT` | Decay constant for both ends | — |
| `--damage-lambda5 FLOAT` | Decay constant for 5' end | — |
| `--damage-lambda3 FLOAT` | Decay constant for 3' end | — |
| `--damage-bg FLOAT` | Background deamination rate | 0.02 |
| `--mask-threshold FLOAT` | Mask when excess P(deam) > T | 0.05 |
| `--pcr-cycles INT` | PCR cycles | 0 |
| `--pcr-efficiency FLOAT` | Efficiency per cycle, 0–1 | 1.0 |
| `--pcr-error-rate FLOAT` | Sub/base/doubling (Q5=5.3e-7, Taq=1.5e-4) | 5.3e-7 |
| `--error-correct` | Enable Phase 3 PCR error correction | off |
| `--errcor-ratio FLOAT` | count(parent)/count(child) ≥ ratio | 50.0 |
| `--errcor-max-count INT` | Children must have count ≤ N | 5 |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 64 |
