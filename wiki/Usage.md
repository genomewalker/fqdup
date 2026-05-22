# Usage

## Prerequisites

### Upstream steps

fqdup operates on reads that have been processed by one upstream tool, then
extends them using its built-in assembler:

**1. fastp merge**: collapse overlapping R1+R2 pairs into single sequences.
Each merged read represents one complete ancient DNA molecule:

```bash
fastp --merge --merged_out merged.fq.gz \
      --in1 R1.fq.gz --in2 R2.fq.gz \
      --disable_adapter_trimming  # (or with adapter trimming, as appropriate)
```

**2. fqdup extend**: extend each merged read from both ends using the built-in
de Bruijn graph assembler. The extended reads serve as deduplication fingerprints;
the original merged reads are kept as the actual output:

```bash
fqdup extend -i merged.fq.gz -o extended.fq.gz
```

See [[Extend]] for full options and algorithm details.

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

The typical workflow runs four commands in sequence. Two optional commands can
run beforehand: `fqdup profile` (diagnostic only) and `fqdup trim` (adapter
stub removal — run when profile reports > ~0.5% of reads carry stub remnants).

### Step 0 (optional): trim adapter stubs

```bash
fqdup profile -i merged.fq.gz   # check for adapter stubs first
fqdup trim -i merged.fq.gz -o merged.trimmed.fq.gz
```

`fqdup profile` reports whether enriched terminal hexamers were detected and
what fraction of reads carry them. If stubs are found, run `fqdup trim` to
remove them before extension. `fqdup trim` handles the adapter stubs that end up
at the start of ultra-short collapsed reads that have read into the ligation
junction. fastp trims the 3′ adapter but doesn't catch these on the 5′ end.

**DS libraries only.** Single-stranded (SS) libraries built with the CircLigase
protocol (Gansauge & Meyer) ligate only a 3′ adapter; the 5′ end of each
molecule is native ancient DNA with no P5 ligation junction. No `CTCTTC` stub
can form and `fqdup profile` will not detect 5′ stubs in SS data. Skip this
step for SS libraries.

### Step 1: extend

```bash
fqdup extend -i merged.fq.gz -o extended.fq.gz
```

Extends each merged read outward from both ends using a 3-pass de Bruijn graph
algorithm. Pass 0 estimates damage parameters; Pass 1 builds an oriented k-mer
graph while masking damaged terminals; Pass 2 extends reads via multi-threaded
anchor scan and unitig walk.

Added bases receive quality `#` (Phred 2). Original bases keep their original
quality. For large datasets (≥100 M reads) specify `--threads` to use all
available cores.

### Step 2: sort

Sort both files by read ID, as above.

### Step 3: structural deduplication (extended + merged reads)

```bash
fqdup derep_pairs \
  -n merged.sorted.fq.gz \
  -e extended.sorted.fq.gz \
  -o-non merged.deduped.fq.gz \
  -o-ext  extended.deduped.fq.gz \
  -c clusters_pairs.tsv.gz
```

`derep_pairs` reads both sorted files in lockstep and hashes the **extended**
(`fqdup extend`-assembled) read as the cluster fingerprint. Using the extended
read reduces false collisions: two different molecules that happen to share a
short merged sequence will typically diverge in the assembled extension.

The representative kept per cluster is the entry with the longest **merged**
read, preserving the most original ancient DNA sequence. The `-c` flag writes
per-cluster statistics to a gzipped TSV.

Note: `derep_pairs` is not about R1/R2 paired-end reads. Its input is
already single-end merged reads from `fastp --merge`; `-n` and `-e` are the
same reads in two forms (non-extended and extended), not two sequencing mates.

For a 25.8 M-pair library, this step typically finishes in about 25 seconds
and reduces to around 5.6 M unique pairs.

### Step 4: single-file biological deduplication

```bash
fqdup derep \
  -i merged.deduped.fq.gz \
  -o merged.final.fq.gz \
  -c clusters_derep.tsv.gz
```

PCR error correction is **on by default**. Damage-aware hashing is **off by
default**, if you run DART or mapDamage on the fqdup output, enabling it would
distort per-position damage frequencies (only the most-damaged representative is
retained per cluster). Enable with `--collapse-damage` only when accurate
unique-molecule counting matters more than downstream damage estimation.

When `--collapse-damage` is given, Pass 0 scans all reads and logs the estimated
damage parameters before deduplication starts:

```
Pass 0: damage estimation, 5582073 reads processed (5582073 total)
  5'-end: d_max=0.099 lambda=0.246 bg=0.487
  3'-end: d_max=0.020 lambda=0.069 bg=0.509
  Masked positions: 1 (1 bp each end)
  Expected mismatches (L=91): 1.44, 99th-pct tolerance: 5
```

Reads that differ only at masked terminal positions then hash identically and
collapse into the same cluster, recovering pairs split by post-mortem
deamination.

Phase 3 (PCR error correction) runs after Pass 1: any cluster with count ≤ 5
that differs from a high-count cluster (≥ 50× its count) by exactly one
substitution outside the damage zone is absorbed as a PCR copying error.
C↔T and G↔A mismatches are never absorbed.

On the 5.6 M-read input with `--collapse-damage` enabled, the full run completes
in about 33 seconds and reduces to roughly 3.5 M unique clusters.

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

### `fqdup extend`

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Input merged FASTQ (.gz or plain) | required |
| `-o FILE` | Output extended FASTQ | required |
| `-k N` | k-mer size, odd 5–31 | 17 |
| `--min-count N` | Minimum edge support to include a k-mer | 2 |
| `--max-extend N` | Maximum bases added per side | 100 |
| `--threads N` | Worker threads | all cores |
| `--min-qual N` | Exclude bases below this Phred quality | 20 |
| `--library-type auto\|ds\|ss` | Damage model library type | auto |
| `--no-damage` | Skip damage estimation; no masking | off |
| `--mask-5 N` | Manually mask N bp at 5' end (skips Pass 0) | - |
| `--mask-3 N` | Manually mask N bp at 3' end (skips Pass 0) | - |
| `--mask-threshold FLOAT` | Excess damage threshold for terminal masking | 0.05 |
| `--damage-sample N` | Reads sampled for damage estimation (0=all) | 1000000 |
| `--bbhash` | Build BBHash MPHF for O(1) lookup (saves RAM, slower) | off |
| `--no-damage-qc` | Disable adapter/hexamer QC report | off |
| `--damage-qc-scan-reads N` | Reads sampled for adapter-stub detection (0=all) | 1000000 |
| `--damage-clip-pass MODE` | `off` \| `report` \| `refit` | report |

Added bases receive quality `#` (Phred 2). Reads with no clean interior k-mers
are written unchanged.

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
| `-e FILE` | Sorted fqdup-extended FASTQ | required |
| `-o-non FILE` | Output merged FASTQ (representatives) | required |
| `-o-ext FILE` | Output extended FASTQ (representatives) | required |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--cluster-format FILE` | Write cluster genealogy to `.fqcl` (use with `derep --prior-fqcl`) | off |
| `--no-revcomp` | Disable reverse-complement collapsing | off |
| `--allow-id-mismatch` | Warn instead of failing on read ID mismatches | off |

### `fqdup profile`

Accepts single-end/merged input (`-i`) or raw paired-end reads (`-1`/`-2`).

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Single-end / merged input FASTQ (.gz or plain) | — |
| `-1 FILE -2 FILE` | Paired-end raw reads (R1 → 5' counters, R2 → 3' counters) | — |
| `-p N` | Worker threads | all cores |
| `--library-type auto\|ds\|ss` | Override library-type detection | auto |
| `--mask-threshold FLOAT` | Mask positions where excess P(deam) > T | 0.05 |
| `--tsv FILE` | Write per-position frequency table as TSV | - |
| `--json FILE` | Write full damage profile as JSON | - |
| `--html FILE` | Write interactive damage report as self-contained HTML | - |
| `--length-bins SPEC` | Length-stratified damage: `auto`, `N`, or `e1,e2,...` | off |
| `--adapter-scan-reads N` | Reads sampled for adapter-stub detection (0=all) | 1000000 |

See [[Profile]] for full output description and typical workflow.

### `fqdup split`

Classify reads as damaged or undamaged **without deduplication**. Runs damage
profile estimation (Pass 0), builds the LLR split model, then streams the
input once. The input does **not** need to be sorted. Useful for splitting
already-deduplicated reads or for a quick pre-filter before `fqdup derep`.

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Input FASTQ (raw or .gz); does not need to be sorted | required |
| `--out-damaged FILE` | Write damaged reads here | - |
| `--out-undamaged FILE` | Write undamaged reads here | - |
| `--library-type ds\|ss\|auto` | Override library-type detection | auto |
| `--split-model auto\|bulk\|empirical` | Split model: auto uses empirical when d_max5>0.01, else bulk | auto |
| `--split-threshold FLOAT` | LLR threshold for damaged call | 0.0 |
| `--damage-deam-sample N` | Max reads scanned for Pass 0 damage estimation | 5000000 |
| `-t N` | Threads (I/O capped at 16) | all cores |

At least one of `--out-damaged` / `--out-undamaged` is required.

### `fqdup derep`

`-o` is optional when `--out-damaged`/`--out-undamaged` are given. Use `--help-advanced` for the full errcor/damage/PCR knob list.

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Sorted input FASTQ | required |
| `-o FILE` | Output deduplicated FASTQ | - |
| `-c FILE` | Cluster statistics (gzipped TSV) | off |
| `--out-damaged FILE` | Write LLR-classified damaged reads | - |
| `--out-undamaged FILE` | Write LLR-classified undamaged reads | - |
| `--split-threshold FLOAT` | LLR threshold for damaged/undamaged split | 0.0 |
| `--split-model auto\|bulk\|empirical` | Split model selection | auto |
| `--cluster-format FILE` | Write `.fqcl` cluster genealogy | off |
| `--prior-fqcl FILE` | Load cluster counts from `derep_pairs --cluster-format` output | - |
| `--no-revcomp` | Disable reverse-complement collapsing | off |
| `--damage off\|report\|collapse` | Damage mode (`collapse` = use in clustering; distorts per-pos rates) | report |
| `--library-type auto\|ds\|ss\|unknown` | Library type for damage model | auto |
| `--damage-dmax N[,N]` | Manual d_max override e.g. `0.21,0.13` for 5',3' (skips fit) | - |
| `--damage-dmax5 FLOAT` | Manual d_max 5' end | - |
| `--damage-dmax3 FLOAT` | Manual d_max 3' end | - |
| `--damage-lambda FLOAT` | Manual decay constant both ends | - |
| `--damage-lambda5 FLOAT` | Manual decay constant 5' end | - |
| `--damage-lambda3 FLOAT` | Manual decay constant 3' end | - |
| `--damage-bg FLOAT` | Background substitution rate | 0.02 |
| `--mask-threshold FLOAT` | Mask when excess P(deam) > T | 0.05 |
| `--damage-scan N` | Reads sampled for QC + adapter scan (0=all) | 1000000 |
| `--damage-deam-sample N` | Reads sampled for d_max/lambda fit (0=all) | 5000000 |
| `--damage-clip-pass MODE` | `off` \| `report` \| `refit` | report |
| `--no-error-correct` | Disable Phase 3 PCR error correction | off |
| `--errcor-rescue-indels` | Syncmer-indexed indel rescue (ed≤2) | off |
| `--pcr-cycles INT` | PCR cycles for D_eff estimate | 0 (auto) |
| `--pcr-efficiency FLOAT` | Efficiency per cycle, 0–1 | 1.0 |
| `--pcr-error-rate FLOAT` | Substitution rate per base per doubling | 5.3e-7 |
| `--errcor-min-parent INT` | Min count to index as parent | 3 |
| `--errcor-max-h2-count INT` | Max child count eligible for H=2 absorption | 2 |
| `--errcor-snp-threshold FLOAT` | SNP veto: sig/parent_count threshold | 0.20 |
| `--errcor-snp-min-count INT` | SNP veto: min absolute sig_count | 1 |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 0 (unlimited) |
| `--errcor-empirical` | Empirical posterior-odds model | **on** |
| `--errcor-legacy-veto` | Revert to pre-T5.8 count-ratio veto | off |
| `--errcor-singleton-qual-min INT` | Block absorption when mismatch Phred ≥ N (non-damage) | 25 |

### `fqdup view`

Inspect `.fqcl` cluster genealogy files produced by `derep --cluster-format` or
`derep_pairs --cluster-format`.

| Flag | Description | Default |
|------|-------------|---------|
| `FILE` | Input `.fqcl` file | required |
| *(no flags)* | Print summary header (total clusters, members, edges) | — |
| `--cluster N` | Render ASCII tree for cluster N | — |
| `--staircase N` | Per-node mismatch grid for cluster N | — |
| `--bundle [--end-k K]` | Group clusters by start+end k-mer; K controls anchor length | K=16 |
| `--bundle-staircase HEX` | Render mismatch staircase across one bundle (hex bundle key) | — |
| `--min-bundle-size N` | Skip bundles with fewer than N clusters | 2 |
| `--top-members N` | List top N clusters by member count | — |
| `--top-edges N` | List top N clusters by edge count | — |
| `--dump-members` | Emit TSV: `cluster_id\tmember_id` for all clusters | — |
| `--member-of FASTQ` | Emit TSV: `read_name\tcluster_id` for every read in FASTQ | — |
| `--json` | Emit structured JSON (schema `fqdup.view.v1`) | — |
| `--html PATH` | Write self-contained HTML visualisation (top 50 clusters) | — |
