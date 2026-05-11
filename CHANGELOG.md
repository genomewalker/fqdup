# Changelog

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
