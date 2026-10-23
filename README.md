# fqdup

> Status as of 2026-10-23: libtaph is pinned to
> `c3b3a5d46e6112e61a461c3a6f10ae9fc30345dd`, the same revision as DART v2.0.0
> (previously `f97613b`). `fqdup profile` now runs libtaph's library profiler
> (`taph::run_library_profile`), the same code as `dart profile`. fqdup's own
> copies of the profiler, HTML report and 8-oxoG code are removed; `--html` and
> `--no-oxog` work through libtaph. Profile JSON changes: each length bin has
> `pi_authentic` (provisional; tau cutoffs come from three labelled
> libraries), `tau_discriminator.point` and `burial_fingerprint.tau_hat` are
> null unless tau is DETECTED, `w_length` is deprecated, and preservation
> f_coh, evidence and score drop by up to 0.074, 0.27 and 0.18 without the
> w_length term. The bulk damage fit uses SQUAREM. The `fqdup merge` source is
> unchanged by this re-pin.

> Status as of 2026-10-22: `fqdup merge` now derives a per-library cutoff
> floor from its own data instead of the fixed 0.05: the smallest grid cutoff
> at which no merged record's certified risk exceeds the bound (never below
> 0.05). A requested cutoff below the floor runs at the floor and is logged;
> `--allow-sub-floor` restores the requested value. At any cutoff at or above
> the floor, including the default 0.316, outputs are byte-identical to the
> previous release on the development benchmark and pilot (every cutoff and
> library checked). `--cutoff auto` is present but refuses to run until its
> error target (epsilon) is declared: it will pick the largest cutoff whose
> truth-free calibrated mate error stays under epsilon without raising
> residual adapter above the 0.316 level, falling back to the floor when
> calibration is bad. `--calibration` writes the held-out family and overlap
> calibration diagnostics. On real libraries the derived floors are higher
> (ds 0.316, ss 0.237) than on simulated ones (0.03–0.06); epsilon will be set
> from real-data measurements and confirmed on a fresh held-out validation.

> Status as of 2026-10-21: `fqdup merge` (native r8qkeep merge; the previous
> pipeline remains reachable only via `--legacy-merge`) now clips confirmed
> adapter from mates it leaves unmerged, accepts 11–19 bp tail/head overlaps
> only when a quality-weighted posterior shows the placement is unique (at
> cutoffs ≥ 0.05 only), rejects 20–30 bp tail/head overlaps whose overlap word
> recurs in a read flank, and places unmerged duplicates at the insert length
> of their merged family siblings. On the simulated development benchmark
> (400,000 pairs with known truth, default cutoff 0.316, against the previous
> release): false merges 10 → 1, missed merges 7,109 → 6,464, exact merges
> 345,198 → 345,589, wrong-length merges 0 → 0 with none added at any cutoff,
> adapter bases left in unmerged mates 25,114 → 1,462, and no true molecule
> base removed by clipping; true positions emitted only wrongly rose by 3
> bases, within an accepted tolerance of 0.001 % of correctly emitted bases.
> Caveats: real libraries with 101 bp reads exercise only the 15–19 bp part of
> the tail/head rescue; terminal-damage and repeat behaviour of rescued reads
> are underpowered on real data.
> Held-out validation (a freshly generated and sealed simulated set of 400,000
> pairs, run once with every tool's settings and the scorers frozen beforehand,
> cutoff 0.316): fqdup emitted the most correct molecule bases (36,222,545) and
> the most exact merges (346,071) of 18 tool configurations, with 0
> wrong-length merges, 3 false merges, 6,506 missed merges and 2,364 adapter
> bases left in unmerged mates. Ranking safety first (false and wrong-length
> merges), then molecule retention, then merge geometry, fqdup is ahead of 12
> of the 14 fair configurations (AdapterRemoval 2 and 3, leeHom, fastp,
> BBMerge with error correction, NGmerge, PEAR, FLASH2, SeqPrep); it is not
> comparable with BBMerge in strict mode (fewer false merges, 1 vs 3, but
> 191,318 missed merges) and the comparison with AdapterRemoval 3 with Mott
> quality trimming is unresolved (one ambiguous pair). PEAR and leeHom miss
> fewer merges and leave no adapter, at the cost of 956 and 52 false merges.

> Historical (undated specifics omitted): across several internal revisions,
> the merge pipeline's default repair-skip, outcomes-file, and
> family-correction behavior changed, and libtaph was repinned without
> changing merge output.

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

Paired-end reads are merged into single sequences first, with `fqdup merge` or
`fastp --merge`. `fqdup` then runs four steps in order. An optional diagnostic
command, `fqdup profile`, can be run beforehand to inspect the damage profile
and verify library-type detection.

**Input. `fqdup merge`**: overlap-merge R1 and RC(R2) into one full-length
ancient DNA molecule, writing into a new output directory (`merged.fq`,
`r1.fq`, `r2.fq`, `orphan.fq`, plus an outcomes file). Tuned for aDNA overlap
detection: high mismatch tolerance in the overlap (`--max-mm-rate` default
0.08) and Bayesian quality consensus at each overlapping base. Adapter
sequences and library type are supplied by the caller — there is no
auto-detection pre-scan — and PCR-duplicate family correction runs by default.
`fastp --merge` remains a supported alternative merger. The previous merge
pipeline (pre-scan auto-detection of library type/adapters/UDG status, a
single merged-file output, and `--json`/`--damage-out`/`--subst-out`
diagnostics) is unchanged but reachable only via `--legacy-merge`; see that
subsection under Options.

**0. `fqdup profile` (optional diagnostic)**: standalone multi-threaded damage
profiler. Scans the input reads and reports d_max, lambda, background rate, and
per-position C→T/G→A frequencies. Classifies the library as double-stranded or
single-stranded via a 7-model BIC competition before fitting. Run this before the
full pipeline to confirm that `--collapse-damage` is warranted, check the mask
threshold, or obtain parameters for manual `--damage-dmax`/`--damage-lambda`
overrides. It does not modify reads. When adapter stubs are detected it reports
the fraction of reads affected (e.g. `5'=CTCTTC (1.2% of reads)`).

**0a. `fqdup trim` (optional pre-processing, DS libraries only)**: remove adapter
stubs that end up at the start of ultra-short collapsed reads that have read into
the ligation junction. fastp trims the 3′ adapter but doesn't catch these on the
5′ end. Uses hexamer frequency analysis on the first `--scan-reads` reads;
single-pass — the scan buffer is replayed into the clip pipeline with no second
file open. Run when `fqdup profile` reports adapter stubs at > ~0.5% of reads.

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
sequencing; its input is single-end merged reads from `fqdup merge` or
`fastp --merge`.

**4. `fqdup derep`**: biological deduplication of the merged-read output.
Three mechanisms:

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
PCR error correction (on by default) is the only relevant step; damage-variant
collapsing is off by default unless `--damage collapse` / `--collapse-damage`
is passed, so no flag is needed to skip it.

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
  -o merged.final.fq.gz \
  --min-length 30            # REQUIRED: caller sets the read-length floor

# Optional: also enable damage-aware hashing (only if NOT running DART/mapDamage downstream)
# fqdup derep -i merged.deduped.fq.gz -o merged.final.fq.gz --collapse-damage
```

### Non-aDNA (modern DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 32G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz \
  --min-length 30
```

---

## Installation

```bash
git clone https://github.com/genomewalker/fqdup.git
cd fqdup
mkdir build && cd build
cmake -DFETCHCONTENT_SOURCE_DIR_LIBTAPH=/path/to/libtaph ..
make -j$(nproc)
```

A bare `cmake ..` fails: CMake refuses to configure unless
`FETCHCONTENT_SOURCE_DIR_LIBTAPH` points at a clean local libtaph checkout at
the pinned commit (see Build provenance below), or `-DFQDUP_ALLOW_REMOTE_LIBTAPH=ON`
is passed to fetch that exact pin from GitHub instead (`CMakeLists.txt`).

**Required:** CMake 3.10+, C++17 compiler, zlib, xxHash, a local or fetchable libtaph checkout

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

### `fqdup merge`

Overlap-merge R1/R2 pairs into single collapsed reads, writing output into a
new directory. `fastp --merge` remains a supported alternative merger.

```
fqdup merge -1 R1.fq[,LANE,...] -2 R2.fq[,LANE,...] -o NEW_OUTPUT_DIR [options]

Required:
  -1 FILE[,...] / -2 FILE[,...]  R1 / R2 input FASTQ, repeatable or
                                  comma-separated; lanes are concatenated into
                                  one logical stream per mate
  -o DIR                          Output directory; must not already exist

Cutoff and threading:
  --cutoff F             Merge-quality cutoff (default: 0.316)
  --cutoffs FILE         Sweep several cutoffs, one per line; each gets its
                         own subdirectory under -o, named by the cutoff text
  -p, --threads N        Worker threads (default: 1, capped at available CPUs)

Merge parameters:
  --adapter1 / --adapter2 SEQ  R1 / R2 adapter sequence (default: fixed
                                adapter sequences; no auto-detection)
  --library-type ds|ss          Library type (default: ds; "auto" not accepted)
  --min-length N                  Minimum merged length (default: 1)
  --min-overlap N                  Minimum overlap length (default: 11)
  --max-mm-rate F                   Max mismatch rate in overlap (default: 0.08)

Family correction (on by default):
  --no-family            Skip PCR-duplicate family correction (changes output)

Prior fit:
  --prior-fit counts|full  counts (default) is the production mode; full is
                           the exact reference mode

Outcomes reporting:
  --outcomes full|summary|none  (default: full). full writes one JSON line
                                 per pair to outcomes.jsonl; summary writes
                                 outcomes.summary.tsv (counts for merged,
                                 unmerged, orphan_r1, orphan_r2, dropped,
                                 total); none writes neither. merged.fq,
                                 r1.fq, r2.fq, and orphan.fq are the same
                                 across all three modes.

Diagnostics:
  --debug-stages   Emit extra timing/memory lines to stderr, plus one
                   debug/ directory directly under -o (not per cutoff) with
                   family, routing, and prior-model intermediate files
  --repair-shadow  Diagnostic repair-routing mode for internal comparison;
                   affects accounting only, not the written reads
  --timings        Emit stage wall-clock lines to stderr

Accepted but currently no-ops: --no-poly-g, --no-internal-panel,
--experimental-short-tail, --experimental-gapped-clipping.
```

Any flag not listed above is rejected, with a hint to use `--legacy-merge`.
Output files are written directly into `-o` (or into each per-cutoff
subdirectory when `--cutoff`/`--cutoffs` names more than one value):
`merged.fq`, `r1.fq`, `r2.fq`, `orphan.fq` (plain, uncompressed FASTQ) and,
unless `--outcomes none`, `outcomes.jsonl` or `outcomes.summary.tsv`.

#### `--legacy-merge`

`fqdup merge --legacy-merge ...` (the flag may appear anywhere in the
arguments) runs the previous merge pipeline instead, with its own CLI and a
single merged-file output rather than a directory:

```
fqdup merge --legacy-merge -1 R1.fq.gz -2 R2.fq.gz -o merged.fq.gz [options]

Required:
  -1 FILE[,...] / -2 FILE[,...]  R1 / R2 input FASTQ; repeatable or comma-separated
  -o FILE                    Output: merged reads (single file)

Library:
  --library-type ds|ss       Override auto-detected library type
  --min-length N              Read-length floor, all emit paths (default: 16; 0=off)
  --max-length N               Drop reads longer than N bp (default: off)

Optional output:
  --r1-out / --r2-out FILE   Unmerged R1 / R2 reads (adapter-trimmed if --adapterN given)
  --orphan-out FILE          Reads where one mate was discarded by QC filters
  --damage-out FILE          Strand-resolved paired damage profile JSON
  --subst-out FILE           Per-position overlap substitution matrix TSV

Overlap:
  --min-overlap N            Minimum overlap length (default: 11)
  --max-mm-rate F            Max mismatch rate in overlap (default: 0.08)
  --clip-r1-5p N             Hard-clip N bases from R1 5' end before overlap
  --min-entropy F            Discard low-complexity merged reads by Shannon entropy (default: off)
  --max-n-rate F             Discard merged reads with N fraction above F (default: 1.0=off)
  --json FILE                Write lossless merge-QC report

Adapter trimming (for unmerged pairs):
  --adapter1 / --adapter2 SEQ  R1 / R2 adapter sequence
  --adapter-fasta FILE         FASTA with adapter pairs (multiple pairs supported)
  --no-internal-panel          Disable built-in aDNA construct read-through trimming (default: on)
```

Unlike the default merge above, library type, adapters, and UDG status are
auto-detected from a 200k-read pre-scan unless overridden with
`--library-type`/`--adapter1`/`--adapter2`, and `--min-length` has a default
rather than being required.

### `fqdup profile`

Standalone damage profiler. Run this before the full pipeline to inspect d_max,
library type, and which positions would be masked.

```
fqdup profile -i INPUT [options]
fqdup profile -1 R1.fq.gz -2 R2.fq.gz [options]   # raw paired-end input
fqdup profile -i MERGED -1 R1.fq.gz -2 R2.fq.gz [options]  # combined: see below

  -i FILE                    Single-end / merged input FASTQ (.gz or plain)
  -1 FILE -2 FILE            Raw paired-end reads (R1 → 5' counters, R2 → 3')
  -p N                       Worker threads (default: all cores)
  --library-type auto|ds|ss  Override library-type auto-detection (default: auto)
  --mask-threshold FLOAT     Flag positions where excess P(deamination) > T
                             (default: 0.05, i.e. 5 pp above background).
                             Output lists which positions exceed the threshold —
                             these are the positions that would be masked in a
                             subsequent fqdup derep --collapse-damage run.
  --min-length N             Drop reads shorter than N bp before profiling
  --max-length N             Drop reads longer than N bp before profiling (default: off)
  --no-oxog                  Skip the 8-oxoG asymmetry channels
  --subst-in FILE[,FILE...]  Read a merge-produced bsubst file; if
                             --library-type is not given, adopt the library
                             type merge explicitly declared (not one it
                             auto-detected)
  --length-bins SPEC         Length-stratified damage curves: auto | off | N | e1,e2,...
                             (default: auto)
  --adapter-scan-reads N     Reads sampled for adapter-stub detection (default: 1000000; 0=all)

Outputs:
  --tsv FILE                 Write per-position frequency table as TSV
  --json FILE                Write the damage profile as JSON
  --json-level summary|standard|full  Detail level of --json (default: standard)
  --count-table-json FILE    Write raw per-position count tables as JSONL
                             (one JSON object per line)
  --damage-summary FILE      Write the taph DamageSummary sidecar (.tdmg) for
                             DART --library-profile (v3 carries pi_damaged_point)
  --damage-json-out FILE     Write a scalar damage model for `fqdup split --damage-json`
  --model-bin-out FILE       Write a full split model for `fqdup split --model-bin`
  --html FILE                Write interactive damage report as self-contained HTML
```

With both `-i` and `-1`/`-2` given (combined mode), the merged input carries
both termini while the unmerged R1/R2 stream adds independent 5'/3' depth; the
two are profiled and finalized separately, not pooled into one set of
estimates.

When adapter stubs are detected, the report includes a per-stub read fraction
counted over all reads: `adapter stubs: 5'=CTCTTC (1.2% of reads) 3'=TTTCCC (1.6% of reads)`.
The same values appear in `--json` as `adapter_stub5_read_fraction`, `adapter_stub3_read_fraction`,
and `adapter_stub_reads_checked`.

Prints a human-readable report: library type with BIC scores, d_max/lambda per end,
per-position deamination frequencies, and the positions that exceed the mask threshold.
Use the output to decide whether `--collapse-damage` is warranted, verify the library-type
call, or supply parameters directly to `fqdup extend` / `fqdup derep`.

`--damage-summary FILE` writes a compact binary DamageSummary sidecar (`.tdmg`) carrying the
finalized per-library damage parameters, for consumption by DART's `--library-profile`. The v3
format adds `pi_damaged_point`, the mixture-based damaged-molecule prevalence used as an
empirical-Bayes prior downstream; it is `-1.0` when the length-stratified joint mixture does not
resolve. This is the same terminal-deamination signal DART would otherwise self-profile, so it is
a prior, not independent evidence.

The `--json` output is libtaph's damage profile, with detail controlled by
`--json-level` (`summary` < `standard` < `full`); `--no-oxog` drops the 8-oxoG
channels from it. The fields themselves (deamination, CpG context, interior
C→T clustering, 8-oxoG asymmetry, depurination, the damage-context summary,
the eight damage-type channels, the ancient-fraction mixture decomposition,
the preservation score, library QC) are libtaph's, not fqdup's — see the
libtaph repository's `docs/fields/` (one page per field) and
`docs/interpretation.md` for what each one means and how to read it. This
README documents only what `fqdup profile` writes and which flags change it.

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
  -n FILE      Sorted merged (fqdup merge / fastp) FASTQ
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
fqdup derep -i INPUT -o OUTPUT --min-length N [options]

Required:
  -i FILE      Sorted input FASTQ
  -o FILE      Output FASTQ
  --min-length N   Drop reads shorter than N bp at ingestion (REQUIRED, >= 1)

Optional:
  --max-length N       Drop reads longer than N bp at ingestion (0 = off; must be >= --min-length)
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
damaged/non-damaged LLR score. Useful for contamination estimation, variant calling
on damaged reads only, or downstream library comparisons.

```bash
fqdup derep -i sorted.fq.gz \
    --min-length 30 \
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
  --min-length N          Drop reads shorter than N bp before classify (default: 16; 0 = off)
  --max-length N          Drop reads longer than N bp before classify (default: off)
  --model-bin FILE        Reuse full split model from `fqdup profile --model-bin-out`
                            (zero estimation: no Pass-0, no LSD scan); digest-validated
  --damage-json FILE      Reuse scalar model from `fqdup profile --damage-json-out`
                            (one LSD pass to build bins); digest-validated
  --allow-model-mismatch  Permit a model whose input digest differs (off by default)
  --damage-deam-sample N  Max reads for Pass 0 damage scan (default: 5000000; 0=all)
  -p N                    Threads (default: all available, capped at 16 for I/O)
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

### `fqdup gen`

Synthetic FASTQ generator for testing: emits IID reads at a configurable GC
content, mixing ancient (damaged) and modern reads at `--f-damaged`. Models
5'/3' deamination, 8-oxoG G→T transversions, and pos-0 purine enrichment from
depurination, plus a constant background substitution rate on all reads.

```
fqdup gen -o OUTPUT [options]

  -o FILE                          Output FASTQ
  -n, --reads N                    Number of reads (default: 1000000)
  --seed N                         PRNG seed (default: 1)
  --read-len N                     Fixed read length (default: 60)
  --gc FLOAT                       GC content fraction (default: 0.45)
  --f-damaged FLOAT                Fraction of ancient reads (default: 0.70)
  --lib-type ds|ss                 Library type (default: ds)
  --dmg5-max / --dmg5-lambda F     5' C→T amplitude / decay rate (default: 0.18 / 0.35)
  --dmg3-ga-max / --dmg3-ga-lambda F  3' G→A amplitude / decay rate (default: same as 5')
  --dmg3-ct-max / --dmg3-ct-lambda F  3' C→T amplitude / decay rate, SS only
  --oxog FLOAT                     Per-G G→T probability (default: 0.0)
  --depur FLOAT                    Pos-0 purine enrichment rate in ancient reads (default: 0.0)
  --bg-ct / --bg-ga FLOAT          Background substitution rate, all reads (default: 0.001)
  --q-score N                      Constant Phred quality (default: 40)
```

### `fqdup consensus` — currently disabled

Designed to emit a per-cluster consensus FASTQ from a `.fqcl` cluster
genealogy (lineage-collapsed Bayesian posterior per parent position, with
damage-zone handling and per-cluster QC flags). Not runnable in this build:
`main.cpp` routes `consensus` to a stub that prints "temporarily disabled
(WIP refactor, task #20)" and exits 2.

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
  --min-length 30 \
  --cluster-format final.fqcl \
  --prior-fqcl derep_pairs.fqcl
```

The `--prior-fqcl` option has no effect when the input was not produced by
`derep_pairs --cluster-format`; it is safe to omit for single-end or pre-merged
libraries.

---

## Performance

No general performance numbers are published here; they depend on library
complexity, depth, and read length. Dataset-specific benchmark figures, where
they exist, are kept in internal task records rather than this README.

### fqdup extend

Extension works outward from both ends of each merged read using the
built-in de Bruijn graph; the fraction of reads extended, the average
per-side extension, wall time, and peak RSS all scale with input size,
library complexity, and thread count.

### fqdup derep_pairs / derep

`derep_pairs` performs the first, structural pass of deduplication; `derep`'s
additional PCR error correction (Phase 3) typically absorbs further reads
beyond what `derep_pairs` alone collapses, at added wall-clock cost. See
Algorithm summary below for the stage breakdown.

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
fqdup extend:      [Pass 0] damage estimation (samples 1M reads by default, --damage-sample)
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
