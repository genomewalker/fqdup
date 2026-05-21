# Changelog

## [2.0.3] - 2026-05-21

### New features

- **`fqdup trim`**: new subcommand to detect and remove 5′/3′ adapter stub
  remnants from collapsed FASTQ. Uses hexamer frequency analysis on the first
  `--scan-reads` (default 1 M) reads to identify enriched terminal 6-mers
  consistent with known adapters, then clips matching reads via a parallel
  work-queue pipeline with ordered BGZF output. Handles multi-stub libraries
  (up to 10 detected candidates per end). Options: `--scan-reads`, `--min-length`,
  `--stub5`/`--stub3` (manual override), `-p` (threads).

- **`fqdup profile` adapter stub fractions**: the human-readable report now
  includes per-stub read fractions over all reads (not just the scan window),
  e.g. `adapter stubs: 5'=CTCTTC (1.2% of reads) 3'=TTTCCC (1.6% of reads)`.
  Same fractions are serialised to JSON as `adapter_stub5_read_fraction`,
  `adapter_stub3_read_fraction`, and `adapter_stub_reads_checked`.

### Performance

- **Single-pass I/O for `fqdup trim`**: the hexamer pre-scan buffers the first
  `--scan-reads` records in memory rather than reading the file twice. After stub
  detection, the buffer is replayed directly into the clip pipeline and the reader
  continues from the current position — no second decompression pass.

- **Single-pass I/O for `fqdup profile`** (SE mode): the adapter-stub pre-scan
  now uses the same multi-threaded reader as the full damage-profiling pass.
  Records are buffered during scanning and fed into the worker queue directly,
  eliminating the previous single-threaded pre-scan pass. On a 94.7 M-read DS
  library: `fqdup profile` 2m07s, `fqdup trim` 1m07s (96 threads, NFS).

## [2.0.2] - 2026-05-14

### New features

- **`--out-damaged` / `--out-undamaged`**: route each deduplicated read to a
  separate FASTQ based on its per-read LLR ancient/modern score. Useful for
  contamination estimation, variant calling on damaged reads only, and library
  comparisons. `-o` is optional when both split outputs are given.

- **`--split-model auto|bulk|empirical`**: controls the LLR classifier used for
  the split.
  - `auto` (default): runs a stripped per-bin empirical scan after the damage
    fit when d_max > 0.01, using mixture-unmixed C→T curves per length bin.
    Falls back to the bulk exponential model (zero extra I/O) on undamaged
    samples.
  - `bulk`: always uses the bulk exponential model already fitted during Pass 0.
    No extra file pass; ~5% less accurate at length-distribution tails.
  - `empirical`: forces the per-bin scan regardless of measured damage level.

- **`DamageSplitModel`** (internal): precomputed per-bin LOD tables built once
  from length-stratified empirical curves. Scoring during Pass 2 is pure array
  lookups with no transcendental calls.

- **`estimate_damage_split_model()`** (internal): stripped accumulation pass
  that collects only T/TC counts per bin via precomputed additive classifier
  coefficients. 7.7× faster than the full `estimate_damage_by_length()` pass
  (150 s vs 1160 s on a 437 M-read SS library), achieved by eliminating CpG
  context tracking, 8-oxoG counters, base-composition profiling, and per-read
  `log()` calls.

- **`DamageEstimate::lsd_hist`**: the bulk damage scan now accumulates a
  log-length histogram, passed as `prebuilt_hist` to
  `estimate_damage_split_model()` to skip its own histogram sub-pass.

- **`--split-threshold F`**: LLR decision boundary (default 0.0). Positive
  values require stronger damage evidence before classifying as ancient.

### No breaking changes

The split flags are purely additive. Existing pipelines that do not use
`--out-damaged` / `--out-undamaged` are unaffected.

## [2.0.1] - 2026-05-11

### New features
- **`--protect-transversions`**: opt-in flag that prevents H=1/H=2 PCR-error
  absorption of A↔T (Channel H) and C↔G (Channel G) singleton reads. Channel F
  (C↔A / G↔T, 8-oxoG) and deamination substitutions were already protected by
  default. Use for high-oxidative-damage libraries where transversion singletons
  are real damage signal rather than PCR errors.

## [2.0.0] - 2026-05-11

### Breaking changes
- **Wire format v3**: `.fqcl` Edge struct gains a `float score` field (NaN = not
  evaluated). `sizeof(Edge)` 16→20 bytes; `kVersion` 2→3. Files written by
  fqdup < 2.0.0 cannot be read by this version (reader hard-fails with a clear
  error message).

### New features
- **LR scores in cluster files**: Each absorbed edge in the `.fqcl` genealogy now
  stores the empirical posterior log-odds score (`score`) from the Phase 3
  error-correction model. `NaN` when the legacy SNP-veto path was used or the
  edge was not evaluated. Visible in `fqdup view` JSON as `score` /
  `score_evaluated`.
- **B3 deamination-aware merge** (`b3_deam_merge.hpp`): For high-damage
  libraries (auto-activated when the profile-derived deamination mass gate
  passes; most effective at `d_max ≥ 0.25`), reads from the same ancient
  molecule that accumulate 3–5 deamination differences outside the mask zone
  now land in the same bucket.
  Interior bases are damage-normalised (T→C at 5′-proximal positions, A→G at
  3′-proximal positions) before hashing, combined with the locus bundle key.
- **HTML cluster visualiser** (`fqdup view --html out.html`): Produces a
  self-contained HTML file (embedded CSS + JS) with interactive cluster
  genealogy trees, LR score histogram, and member sequence staircase. No server
  required.
- **Oxidative channels F/G/H**: Z-scores and prefix-conditioned adapter exclusion
  for complement channels (see `channel_f_z`, `channel_h_z`, `channel_h_z_p2plus`
  in the damage JSON output).

### Improvements
- `json_escape` in `fqdup view` now escapes `<`, `>`, `&` as Unicode escapes to
  prevent `</script>` injection when JSON is embedded in HTML output.
- DART→taph rename in internal phase-timer log messages.

### Bug fixes
- PCR error reads in synthetic test data (`gen_synthetic`, `test_errcor`,
  `test_pipeline`) now carry phred 15 at error positions, correctly modelling
  the singleton quality veto introduced in the empirical model.
