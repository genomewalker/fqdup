# Changelog

## [2.0.4] - 2026-05-29

### New features

- **Damage report HTML redesign**: `fqdup profile --html` now produces a
  tabbed dashboard (Bulk / Fractions / Length Bins / Library / QC) with white
  cards on a light neutral background, a sticky stat strip (d₅, d₃, λ, π,
  read count), and a symmetric smiley plot. The 5′ side is plotted left (0→14)
  and the 3′ side is mirrored right (0→14); a center divider separates the
  two termini. The data blob is emitted as `window.D` for cross-script access.

- **Per-fraction exponential decay fitting**: `ancient_fraction.ancient` and
  `ancient_fraction.modern` in `--json` now include `d_max_5prime_fit`,
  `lambda_5prime`, `d_max_3prime_fit`, `lambda_3prime` — WLS+IRLS exponential
  fits on the per-position hard-call rate arrays. Fraction cards in the HTML
  report display these alongside `n_reads` for each class.

- **Posterior-weighted (soft-EM) ancient-fraction estimation**: π and d₅/d₃
  are now estimated via sigmoid-weighted accumulation rather than hard LLR > 0
  thresholding. A per-read posterior `P(ancient | LLR, π_prior)` weights
  contributions to the T/(T+C) and A/(A+G) soft accumulators over
  `N_SOFT_POS` terminal positions. The π prior is derived from the ratio of
  bulk d_max to the per-fraction LSD d_max. This eliminates PPV collapse at
  low endogenous fractions where hard-call classification is dominated by
  false positives. `ancient_fraction.ancient.fraction` now reflects soft-EM π.

- **DS 3′ pos-0 artifact detection**: the pos-0 artifact flag for the 3′ end
  now also fires on the DS adapter blunting gap pattern (pos-0 G→A depleted
  while pos-1 is elevated), in addition to the existing SS ligation spike
  check. This fixes false-zero `d_max_3prime` in DS libraries where blunting
  suppresses the G→A rate at the terminal base.

- **CircLigase 3′ hexamer-bias confound detection**: `fqdup profile` detects
  when the ligation-site hexamer composition biases the 3′ G→A channel,
  overrides `d_max_combined` to `5prime_only` when the 3′ end is confounded,
  and adjusts `d_metamatch` accordingly.

- **LSD classifier pos-0 exclusion on hexamer artifact**: when a pos-0
  hexamer artifact is detected, the per-read LLR classifier skips that
  terminal position in its damage sum, preventing the artifact spike from
  inflating the ancient-class score for individual reads.

### Performance

- **LSD pass fused into oxoG second pass**: LSD per-read classification runs
  within the same FASTQ iteration as the oxoG Pass 2, eliminating one full
  FASTQ decompression for ancient libraries (~15–30 s per run on NFS).

- **LLR envelope memoized; revcomp allocation eliminated**: per-position
  `addT[p]`/`addC[p]` log-odds increments are precomputed once into a
  300-float table. DS reverse-complement reads are scored in-place via
  `complement_code[decoded[len-1-p]]`, removing one heap alloc + memcpy per
  DS read in the oxoG pass.

### Bug fixes

- **Length-bins 3′ model curve**: the 3′ fit curve in the Length Bins HTML
  panel referenced an undefined `modelCurveRev`. Fixed to `modelCurve3`
  plotted against `pos3x` (mirrored 3′ x-axis), consistent with the smiley
  and fractions panels.

---

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
