# fqdup

FASTQ deduplication for ancient DNA libraries, with damage-aware hashing and
PCR error correction.

---

Ancient DNA deamination turns one original molecule into a cloud of C→T/G→A
sequence variants. Exact-sequence dedup splits those true duplicates and inflates
unique read counts. `fqdup` collapses damage variants into a single cluster with
damage-aware hashing, then removes PCR-error variants with a fast
mismatch-tolerant search. Downstream complexity metrics reflect molecules, not
damage patterns.

## How it works

`fqdup` runs four steps in order. An optional diagnostic command, `fqdup profile`,
can be run beforehand to inspect the damage profile and verify library-type
detection.

**0. `fqdup profile` (optional diagnostic)**: standalone multi-threaded damage
profiler. Scans the input reads and reports d_max, lambda, background rate, and
per-position C→T/G→A frequencies. Classifies the library as double-stranded or
single-stranded via a 7-model BIC competition before fitting. Run this before the
full pipeline to confirm that `--collapse-damage` is warranted, check the mask
threshold, or obtain parameters for manual `--damage-dmax`/`--damage-lambda`
overrides. It does not modify reads. When adapter stubs are detected it reports
the fraction of reads affected (e.g. `5'=CTCTTC (1.2% of reads)`).

**0a. `fqdup trim` (optional pre-processing, DS libraries only)**: detect and remove 5′/3′ adapter
stub remnants that upstream tools (fastp, cutadapt) may miss. Typical case: P5
tail hexamers (`CTCTTC`) left at the 5′ end of collapsed reads when fastp trims
only the 3′ adapter from R1 but leaves the P5 tail intact at the 5′ read start.
Uses hexamer frequency analysis on the first `--scan-reads` reads; single-pass —
the scan buffer is replayed into the clip pipeline with no second file open. Run
when `fqdup profile` reports adapter stubs at &gt; ~0.5% of reads.

**DS libraries only.** Single-stranded (SS) libraries prepared with the
CircLigase protocol (Gansauge & Meyer) ligate only a 3′ adapter; the 5′ end of
every molecule is native ancient DNA with no P5 ligation junction. No `CTCTTC`
stub can form. `fqdup profile` will not detect 5′ stubs in SS data; skip
`fqdup trim` for SS libraries.

**1. `fqdup extend`**: extend each merged read outward from both ends using a
built-in de Bruijn graph assembler. The extended reads serve as deduplication
fingerprints; the original merged reads remain the primary output. Extension
also solves QC-trim length variation: PCR copies of the same molecule trimmed
to different lengths by quality/adapter trimming produce different raw sequences
but extend to the same shared interior k-mers, collapsing correctly. Ancient
DNA terminal damage is accounted for automatically, damaged terminal positions
are excluded from k-mer graph construction so they do not create spurious
branches.

**2. `fqdup sort`**: sort both input files by read ID. Required by both
dedup steps.

**3. `fqdup derep_pairs`**: structural deduplication using two files for the same
reads: the original merged reads (`-n`) and their `fqdup extend`-assembled
counterparts (`-e`). The extended sequence is used as the cluster fingerprint.
Short aDNA fragments that share a merged-read core typically diverge in their
extensions, reducing false collisions. The representative kept per cluster is the
entry with the longest merged read. This step is not about R1/R2 paired-end
sequencing; its input is single-end merged reads from `fastp --merge`.

**4. `fqdup derep`**: biological deduplication of the merged-read output.
Two mechanisms, both on by default:

- **PCR error correction (Phase 3, default: on).** A pigeonhole Hamming
  search finds low-count clusters that differ from a high-count cluster by one
  or two interior substitutions (H≤2). An empirical posterior-odds model scores
  each candidate edge using cross-bundle recurrence; reads are absorbed only when
  the evidence favours PCR error over a real variant.
  C↔T and G↔A mismatches are never absorbed; they are indistinguishable from
  damage signal and are always kept.

- **Damage variant collapsing (default: off).** When enabled with `--collapse-damage`,
  terminal positions where the observed C→T or G→A rate exceeds a threshold are
  replaced with a neutral byte before hashing. Reads that differ only by
  deamination collapse into the same cluster. Library type is detected
  automatically via a 7-model BIC competition; override with
  `--library-type ds|ss`.

- **Phase B3: damage-aware H>2 merge (default: on when `--damage collapse/report`).** For
  libraries with d_max ≥ 0.25, PCR copies of the same ancient molecule can
  accumulate 3–5 deamination events at positions just outside the mask zone,
  placing them beyond Phase 3's H≤2 reach. Phase B3 normalizes damage-probable
  interior positions (T→C / A→G), hashes the result, and gates candidate pairs
  by a composite key combining the damage-normalized interior hash with the
  `bundle_key` locus anchor (`start_kmer ⊕ end_kmer`). Only reads sharing
  both the same genomic locus **and** the same normalized interior are compared.
  Within each bucket, all mismatches must be deamination-consistent and pass a
  count-ratio LRT before absorption. See [[Damage-Aware-Deduplication]] for
  details.

  **Leave this off if you plan to run metaDMG or mapDamage on the fqdup output.**
  Here is why. Say one molecule was sequenced 5 times and 3 of those reads
  picked up a C→T at position 1:

  ```
  read_1  TGCATGA...   <- deaminated at pos 1
  read_2  CGCATGA...   <- undamaged
  read_3  TGCATGA...   <- deaminated at pos 1
  read_4  CGCATGA...   <- undamaged
  read_5  TGCATGA...   <- deaminated at pos 1
  ```

  Without `--collapse-damage`, all 5 reads are passed to metaDMG or mapDamage.
  At position 1 they count T=3, C=2, giving a C→T frequency of 0.60 — a
  clear damage signal.

  With `--collapse-damage`, fqdup recognises that these differ only at a
  deaminated terminal position and collapses all 5 into one representative
  (whichever read ID sorted first). Downstream tools see a single read. At
  position 1 they count either T=1 or C=1. The frequency is either 0 or 1 for
  that one molecule, and across the library the per-position frequencies are
  unreliable because coverage is reduced to 1x per molecule.

  Enable `--collapse-damage` only when you need an accurate unique-molecule count
  and are not running downstream damage analysis.

Use `--no-error-correct` to disable PCR error correction. For non-aDNA data,
error correction is the only relevant step, `--no-damage` is already the
default.

---

## Quick start

### Full ancient DNA pipeline

```bash
# 1. Extend merged reads using the built-in de Bruijn graph assembler
fqdup extend -i merged.fq.gz -o extended.fq.gz

# 2. Sort merged and extended reads by read ID
fqdup sort -i merged.fq.gz   -o merged.sorted.fq.gz   --max-memory 64G --fast
fqdup sort -i extended.fq.gz -o extended.sorted.fq.gz --max-memory 64G --fast

# 3. Structural dedup: one representative per cluster, using extended reads as fingerprint
fqdup derep_pairs \
  -n merged.sorted.fq.gz \
  -e extended.sorted.fq.gz \
  -o-non merged.deduped.fq.gz \
  -o-ext  extended.deduped.fq.gz

# 4. PCR error correction (on by default); damage-aware hashing off by default
fqdup derep \
  -i merged.deduped.fq.gz \
  -o merged.final.fq.gz

# Optional: also enable damage-aware hashing (only if NOT running DART/mapDamage downstream)
# fqdup derep -i merged.deduped.fq.gz -o merged.final.fq.gz --collapse-damage
```

### Non-aDNA (modern DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz \
  --no-error-correct
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

> Pass the binary as an absolute path, the test script changes directory internally.

---

## Options

### `fqdup profile`

Standalone damage profiler. Run this before the full pipeline to inspect d_max,
library type, and which positions would be masked.

```
fqdup profile -i INPUT [options]
fqdup profile -1 R1.fq.gz -2 R2.fq.gz [options]   # raw paired-end input

  -i FILE                    Single-end / merged input FASTQ (.gz or plain)
  -1 FILE -2 FILE            Raw paired-end reads (R1 → 5' counters, R2 → 3')
  -p N                       Worker threads (default: all cores)
  --library-type auto|ds|ss  Override library-type auto-detection (default: auto)
  --mask-threshold FLOAT     Flag positions where excess P(deamination) > T
                             (default: 0.05, i.e. 5 pp above background).
                             Output lists which positions exceed the threshold —
                             these are the positions that would be masked in a
                             subsequent fqdup derep --collapse-damage run.
  --tsv FILE                 Write per-position frequency table as TSV
  --json FILE                Write full damage profile as JSON
  --html FILE                Write interactive damage report as self-contained HTML
  --length-bins SPEC         Length-stratified damage curves: auto | N | e1,e2,...
  --adapter-scan-reads N     Reads sampled for adapter-stub detection (default: 1000000; 0=all)
```

When adapter stubs are detected, the report includes a per-stub read fraction
counted over all reads: `adapter stubs: 5'=CTCTTC (1.2% of reads) 3'=TTTCCC (1.6% of reads)`.
The same values appear in `--json` as `adapter_stub5_read_fraction`, `adapter_stub3_read_fraction`,
and `adapter_stub_reads_checked`.

Prints a human-readable report: library type with BIC scores, d_max/lambda per end,
per-position deamination frequencies, and the positions that exceed the mask threshold.
Use the output to decide whether `--collapse-damage` is warranted, verify the library-type
call, or supply parameters directly to `fqdup extend` / `fqdup derep`.

The `--json` output includes the complete machine-readable profile:

- **Deamination** (`deamination`): `d_max_5prime`, `d_max_combined`, `d_metamatch`
  (channel-B anchored blended estimate calibrated against metaDMG, r=0.818 on clean DS
  samples), per-position CT/GA arrays (15 positions each end), BIC scores for
  library-type classification, CpG vs non-CpG amplitude split (`cpg_like`), and
  length-stratified d_max bins (`by_length`, populated when `--length-bins` is given).
- **CpG context** (`deamination.cpg_like`): 5' C→T fitted separately for CpG (`xG`)
  vs non-CpG positions. `cpg_ratio = dmax_cpg / dmax_noncpg`; values below 1.0 indicate
  CpG protection consistent with methylation.
- **Interior C→T clustering** (`interior_ct_cluster`): excess co-occurrence of T at
  non-CpG `{C,T}` sites within the read interior at separations d=1-10 bp.
  `short_asym_log2oe` is CT-track minus AG-control log2 observed/expected; positive values
  indicate clustered interior deamination.
- **8-oxoG asymmetry** (`complement_asymmetry`): `s_oxog` (overall strand asymmetry)
  and `s_oxog_16ctx[16]` (G→T asymmetry across all 16 flanking-dinucleotide contexts).
- **Depurination** (`depurination`): A/G enrichment at 5' and 3' termini; `detected`
  flag with z-score and p-value.
- **Damage-context summary** (`damage_context`): training-free, reference-free summary
  with six damage-process scores (`terminal_deamination_score`, `cpg_context_score`,
  `dipyrimidine_context_score`, `oxidative_context_score`, `fragmentation_context_score`,
  `library_artifact_score`) and a `dominant_process` string indicating the most likely
  damage mechanism. The `evidence` sub-object records the raw inputs used.
- **Damage-type channels** (`damage_types`): array of 8 entries (channels A-H) each
  covering an independent damage mechanism: A=cytosine deamination, B=stop-codon
  deamination validator, C=8-oxoG top strand, D=Chargaff G→T asymmetry, E=purine
  enrichment/depurination, F=8-oxoG bottom strand, G=hydantoin oxidation, H=adenine
  oxidation. Each entry carries `detected`, mechanism-specific rates, and z-scores.
  All eight channels are always present.
- **Preservation score** (`preservation`): composite score (0-1) combining deamination
  evidence, exponential fit quality, CpG age-bias, and oxidation evidence. The
  `authenticity_eff` field is the deamination-only component.
- **Library QC** (`library_qc`): adapter remnant flags, hexamer entropy, composition
  bias z-scores, and a `flags` array listing detected artifacts.

Library-type classification returns `unknown` when no damage signal is detectable
(BIC cannot distinguish DS from SS on a flat profile). This is the correct conservative
call — downstream `fqdup derep` treats `unknown` as DS for masking purposes.

### `fqdup trim`

Detect and remove 5′/3′ adapter stub remnants from collapsed FASTQ. Run after
`fqdup profile` reports adapter stubs at > ~0.5% of reads, before `fqdup extend`.

```
fqdup trim -i INPUT -o OUTPUT [options]

Required:
  -i FILE          Input FASTQ (.gz or plain)
  -o FILE          Output FASTQ (.gz)

Options:
  -p N             Worker threads (default: all cores)
  --scan-reads N   Reads sampled for stub detection (default: 1000000; 0=all)
  --min-length N   Discard trimmed reads shorter than N bp (default: 15)
  --stub5 SEQ      Force 5' stub hexamer (skip auto-detection)
  --stub3 SEQ      Force 3' stub hexamer (skip auto-detection)
```

Single-pass: the scan window is buffered in memory and replayed directly into
the parallel clip pipeline — no second decompression pass. 5′ stubs are clipped
once (first 6 bp); 3′ stubs are clipped iteratively until no match remains.
Reports the fraction of reads clipped at each end.

**DS libraries only.** SS libraries (CircLigase-based prep) have no P5 adapter
at the 5′ end of reads; `fqdup trim` is a no-op on such data.

### `fqdup extend`

```
fqdup extend -i INPUT -o OUTPUT [options]

  -i FILE                Input merged FASTQ (.gz supported)
  -o FILE                Output extended FASTQ

  -k N                   k-mer size (default: 17, max: 31)
  --min-count N          Minimum edge support to include a k-mer (default: 2)
  --max-extend N         Maximum bases added per side (default: 100)
  --threads N            Worker threads (default: all CPU cores)
  --min-qual N           Exclude bases below this Phred quality (default: 20)
  --library-type TYPE    Damage model library type: auto|ds|ss (default: auto)
  --no-damage            Skip damage estimation; no masking
  --mask-5 N             Manually mask N bp at 5' end (skips Pass 0)
  --mask-3 N             Manually mask N bp at 3' end (skips Pass 0)
  --mask-threshold F          Excess damage threshold for terminal masking (default: 0.05)
  --damage-sample N           Estimate damage from first N reads only (default: 1000000; 0=all)
  --bbhash                    Build BBHash MPHF for O(1) k-mer lookup (saves RAM, slower build)
  --no-damage-qc              Disable adapter/hexamer QC report
  --damage-qc-scan-reads N    Reads sampled for adapter-stub detection (default: 1000000; 0=all)
  --damage-clip-pass MODE     off | report (default) | refit
```

Added bases receive quality '#' (Phred 2); original bases keep their original quality scores.

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
  -e FILE      Sorted fqdup-extended FASTQ
  -o-non FILE  Output merged FASTQ (representatives)
  -o-ext FILE  Output extended FASTQ (representatives)

Optional:
  -c FILE               Cluster statistics (gzipped TSV)
  --cluster-format FILE Write a .fqcl cluster genealogy file (wire format v3)
  --no-revcomp          Disable reverse-complement collapsing (default: enabled)
  --allow-id-mismatch   Warn instead of failing on read ID mismatches between -n and -e
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

PCR error correction is **on by default**. Damage variant collapsing is **off by
default**, use `--collapse-damage` only when downstream tools do not
run damage analysis on the fqdup output.

#### Damage variant collapsing (default: off)

```
  --collapse-damage             Enable damage estimation and masking (Pass 0)
  --damage-dmax   FLOAT     d_max for both 5' and 3' ends (manual model)
  --damage-dmax5  FLOAT     d_max for 5' end only
  --damage-dmax3  FLOAT     d_max for 3' end only
  --damage-lambda FLOAT     Decay constant for both ends
  --damage-lambda5 FLOAT    Decay constant for 5' end only
  --damage-lambda3 FLOAT    Decay constant for 3' end only
  --damage-bg     FLOAT     Background deamination rate (default: 0.02)
  --mask-threshold FLOAT    Mask positions where excess P(deamination) > T (default: 0.05)
  --library-type TYPE       auto (default) | ds | ss | unknown
  --damage-scan N           Reads sampled for QC + adapter scan (default: 1000000; 0=all)
  --damage-deam-sample N    Reads sampled for d_max/lambda fit (default: 5000000; 0=all)
  --damage-clip-pass MODE   off | report (default) | refit
```

If both `d_max_5` and `d_max_3` are below 0.02 after estimation, damage is
treated as negligible and standard exact hashing is used, no user intervention
required.

#### Damaged/undamaged split output (optional)

Route each deduplicated read to a separate file based on its per-read
ancient/modern LLR score. Useful for contamination estimation, variant calling
on damaged reads only, or downstream library comparisons.

```bash
fqdup derep -i sorted.fq.gz \
    --out-damaged   ancient.fq.gz \
    --out-undamaged modern.fq.gz \
    --damage collapse           # damage fit required for split scoring
```

`-o` is optional when both `--out-damaged` and `--out-undamaged` are provided.

```
  --out-damaged FILE    Write LLR-classified ancient reads to FILE
  --out-undamaged FILE  Write LLR-classified modern reads to FILE
  --split-threshold F   LLR decision boundary (default: 0.0)
  --split-model MODE    auto (default) | bulk | empirical
                          auto      — per-bin empirical if d_max > 0.01, else bulk
                          bulk      — bulk exponential only, no extra file pass
                          empirical — force per-bin scan regardless of damage level
```

`auto` runs a stripped length-stratified scan after the damage fit, building
empirical per-bin C→T curves via mixture unmixing. Uses precomputed additive
classifier coefficients (no transcendentals in the read loop) and 4 worker
threads. Falls back to the free bulk exponential model on undamaged samples.
See [wiki/Derep](wiki/Derep.md#damagedundamaged-split) for the LLR formula and
per-mode performance numbers.

#### PCR error correction: Phase 3 (default: on)

```
  --cluster-format FILE        Write a .fqcl cluster genealogy file (wire format v3)
  --prior-fqcl FILE            Load cluster_id→n_members from a derep_pairs .fqcl to
                               weight Phase 3 LR count-ratio tests (recommended when
                               using the cascade workflow; see below)

  --no-error-correct           Disable Phase 3
  --error-correct              Explicitly enable (already default)
  --errcor-min-parent    INT   Min count to index as parent (default: 3)
  --errcor-snp-threshold FLOAT SNP veto: sig/parent_count threshold (default: 0.20)
  --errcor-snp-min-count INT   SNP veto: min absolute sig_count (default: 1)
  --errcor-bucket-cap    INT   Max pair-key bucket size (default: 0 = unlimited)
  --errcor-empirical           Empirical posterior-odds model (default: on)
  --errcor-legacy-veto         Revert to pre-T5.8 count-ratio veto
  --errcor-max-h2-count INT    Max child count eligible for H=2 absorption (default: 2)
  --errcor-singleton-qual-min INT  Block absorption when mismatch Phred ≥ N (default: 25)
  --protect-transversions      Protect A↔T / C↔G (Channels H/G) from absorption (default: off)
  --errcor-rescue-indels       Syncmer-indexed indel rescue, ed≤2 (default: off)

  PCR model (for D_eff log estimate only, does not affect absorption):
  --pcr-cycles      INT   Number of PCR cycles
  --pcr-efficiency  FLOAT Amplification efficiency per cycle, 0–1 (default: 1.0)
  --pcr-error-rate  FLOAT Sub/base/doubling, Q5: 5.3e-7 (default),
                          Phusion: 3.9e-6, Taq: 1.5e-4
```

### `fqdup split`

Classify reads as damaged or undamaged **without deduplication**. Runs damage
profile estimation (Pass 0) and an optional length-stratified scan, then streams
the input once. Input does **not** need to be sorted. Useful for splitting
already-deduplicated reads or as a fast pre-filter before `fqdup derep`.

```
fqdup split -i INPUT [options]

Required:
  -i FILE                 Input FASTQ (raw or .gz); does not need to be sorted

Outputs (at least one required):
  --out-damaged  FILE     Write damaged reads here
  --out-undamaged FILE    Write undamaged reads here

Options:
  --library-type STR      ds | ss | auto (default: auto)
  --split-model STR       auto (default) | bulk | empirical
                            auto      — empirical when d_max5 > 0.01, else bulk
                            bulk      — bulk exponential model only
                            empirical — always run length-stratified LSD scan
  --split-threshold F     LLR threshold for damaged call (default: 0.0)
  --damage-deam-sample N  Max reads for Pass 0 damage scan (default: 5000000; 0=all)
  -t N                    Threads (default: all available, capped at 16 for I/O)
```

### `fqdup view`

Inspect `.fqcl` cluster genealogy files produced by `derep --cluster-format` or
`derep_pairs --cluster-format`.

```
fqdup view FILE.fqcl [options]

  (no flags)              Print summary header: total clusters, members, edges
  --cluster N             Render ASCII tree for cluster N
  --staircase N           Per-node mismatch grid for cluster N
  --bundle [--end-k K]    Group clusters by start+end k-mer anchor (default K=16)
  --bundle-staircase HEX  Render mismatch staircase across one bundle (hex key)
  --min-bundle-size N     Skip bundles with fewer than N clusters (default: 2)
  --top-members N         List top N clusters by member count
  --top-edges N           List top N clusters by edge count
  --dump-members          Emit TSV: cluster_id<TAB>member_id for all clusters
  --member-of FASTQ       Emit TSV: read_name<TAB>cluster_id for every read in FASTQ
  --json                  Emit structured JSON (schema fqdup.view.v1)
  --html PATH             Write self-contained HTML visualisation (top 50 clusters)
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

### `.fqcl` cluster genealogy format (wire format v3)

Both `derep_pairs --cluster-format` and `derep --cluster-format` can emit a binary
`.fqcl` file recording the full cluster genealogy. Each cluster stores:

- **nodes** — the root (parent) and all absorbed reads in BFS order
- **edges** — one entry per absorption with `{from, to, pos, ref, alt, n_reads,
  damage_like, score, score_evaluated}`, where `score` is the log-odds LR from the
  empirical model (NaN when not evaluated)

View a `.fqcl` file with:

```bash
fqdup view clusters.fqcl               # terminal summary
fqdup view clusters.fqcl --html out.html  # self-contained HTML visualiser
```

### Cascade workflow (`derep_pairs` → `derep`)

For the best error-correction quality, run `derep_pairs` first and pass its genealogy
to `derep` via `--prior-fqcl`. This seeds Phase 3 count-ratio weights from the
pair-level cluster sizes, preventing over-aggressive merging of independently
captured molecules:

```bash
# Step 1 — pair-aware dedup
fqdup derep_pairs \
  -n sorted.fq.gz -e extended.sorted.fq.gz \
  -o-non derep_pairs.fq.gz -o-ext derep_pairs_ext.fq.gz \
  --cluster-format derep_pairs.fqcl

# Step 2 — damage-aware error correction seeded by pair priors
fqdup sort -i derep_pairs.fq.gz -o derep_pairs.sorted.fq.gz --max-memory 8G
fqdup derep \
  -i derep_pairs.sorted.fq.gz \
  -o final.fq.gz \
  --cluster-format final.fqcl \
  --prior-fqcl derep_pairs.fqcl
```

The `--prior-fqcl` option has no effect when the input was not produced by
`derep_pairs --cluster-format`; it is safe to omit for single-end or pre-merged
libraries.

---

## Performance

### fqdup extend

Benchmarked on DS4, 198 M reads, 16 threads, dandycomp02fl:

| Metric | Value |
|--------|-------|
| Reads extended | 50% |
| Average extension | 3.82 / 3.83 bp per side |
| Wall time | 6:56 (pass1=68 s, finalize=19 s, pass2=5:27) |
| Peak RSS | 41.3 GB |

Comparison: Tadpole achieves 14.75% extension on the same dataset in ~5 min; `fqdup extend` achieves 50% in 7 min.

### fqdup derep_pairs / derep

Benchmarked on `a88af16f35`, 25.8 M fastp-merged reads, ~91 bp mean length,
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
fqdup extend:      [Pass 0] damage estimation (samples 500k reads)
                   [Pass 1] build oriented k-mer graph (masks damaged terminals)
                   [Pass 2] multi-threaded anchor scan + unitig walk → write extended reads
fqdup sort:        chunk ingestion → parallel sort → k-way merge
fqdup derep_pairs: [Pass 1] stream (ext + non) → build index
                   [Pass 2] stream again → write representative pairs
fqdup derep:       [Pass 0] stride-sample → fit damage model (--collapse-damage)
                   [Pass 1] stream → build index with damage masking
                   [Phase 3] 3-way pigeonhole PCR error correction
                   [Pass 2] stream → write unique representatives
```

**Canonical hash:** `min(XXH3_128(seq), XXH3_128(revcomp(seq)))`: collapses
forward and reverse-complement reads into the same cluster. Collision probability
~1.5×10⁻²³ at 100 M reads.

**Damage masking:** empirical per-position mask derived from observed T/(T+C)
and A/(A+G) frequencies; symmetric masking preserves
`hash(seq) == hash(revcomp(seq))` after masking.

**Phase 3 protection:** G↔T and C↔A (8-oxoG, Channel F) are always protected.
C↔T and G↔A (deamination) are additionally protected when damage mode is active.
A↔T (Channel H) and C↔G (Channel G) transversions are eligible for absorption
by default; pass `--protect-transversions` to protect them too — recommended
for high-oxidative-damage libraries where these substitutions are genuine
ancient damage rather than PCR errors.

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
- Rochette NC et al. (2023) On the causes, consequences, and avoidance of PCR duplicates. *Mol Ecol Resour* 23:1299–1318. doi:10.1111/1755-0998.13800
