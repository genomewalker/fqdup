# fqdup split

Classifies reads as damaged/undamaged without deduplicating (`split.cpp:0-2`).
Estimates a damage profile, builds a per-read LLR split model, then streams the
input once routing each read to `--out-damaged` / `--out-undamaged`.

Line references are to `src/split.cpp` on `feat/merge-indel-dimer-qc` at `aa3d27f`.

## Pipeline

**1. Obtain a `DamageProfile`** — three mutually exclusive tiers, selected by flag:

- **Tier 0** `--model-bin FILE` — loads a full split model
  (`fqdup_model_io::read_split_model_bin` :132), produced by
  `fqdup profile --model-bin-out`. Skips both Pass-0 estimation and the LSD
  scan entirely.
- **Tier 1** `--damage-json FILE` — loads a scalar model
  (`load_damage_model_json` :148), produced by `fqdup profile --damage-json-out`.
  Skips Pass-0; one LSD pass still builds the length-stratified bins.
- **Tier 2** (default) — self-estimate via `estimate_damage_with_qc` (:164),
  capped at `--damage-deam-sample` reads (default 5,000,000). Adaptive
  escalation (:171-188): if the profile is inconclusive
  (`pi_state` not `DETECTED`/`LOW_ABUNDANCE`) and the strongest signal
  `d_max > 0.02` at the read cap, rescans at full depth (`max_reads = 0`).

**2. Digest guard** (`guard_digest` :118) — for Tier 0/1, compares the model's
stored input digest against a freshly computed digest of `-i`; refuses a
mismatched pairing unless `--allow-model-mismatch` is set.

**3. Build the empirical LLR scorer** (skipped if Tier 0 already loaded a full
model). Whether to run the length-stratified scan is `--split-model`-gated
(:198-201): `bulk` forces it off, `empirical` forces it on, `auto` (default)
runs it iff `max(d_max_5prime, d_max_3prime) > 0.01`. If run,
`estimate_damage_split_model` (:206) builds a `LengthStratifiedDamageProfile`,
reusing Tier 2's Pass-0 histogram when available. `split_model =
DamageSplitModel::build(lsd_data, profile)` (:214).

**4. Threshold selection** — `split_policy(pi_est)` (:228) gates on the
confidence state of the estimated damage fraction `pi`:

- **not splittable** (`ABSTAIN`/`BELOW_FLOOR`) — threshold set to `+inf`; every
  read routes to undamaged (:240-242).
- **`LOW_ABUNDANCE`** ("yield-locked") — an extra full pass scores every read
  into a 20001-bin LLR histogram over `[-100, 100]`; the cut is placed where
  the upper-tail count matches `pi * n_scan`, i.e. the top-`pi`-fraction
  best-ranked reads are routed as damaged (:243-267).
- **`DETECTED`** (posterior-threshold) — `effective_threshold = split_threshold
  - log_prior_odds` (:268-277).

**5. Streaming pass** (:308-321) — reads the input once, scores each read with
`split_model.score(seq)` (LLR), and writes it to the damaged or undamaged
output by comparing the score to `effective_threshold`.

## Score dump (debug)

`FQDUP_SPLIT_SCORES` env var (not a CLI flag) — if set, dumps
`header<TAB>split_model_llr<TAB>shared_lsd_llr` per read to the given path, for
computing a per-read ROC/AUC against synthetic-data ground truth (:298-313).

## CLI flags

(`usage` :22-48)

| Flag | Effect | Default |
|---|---|---|
| `-i, --input FILE` | Input FASTQ (raw or `.gz`); need not be sorted | required |
| `--out-damaged FILE` | Write damaged reads here | one of the two `--out-*` required |
| `--out-undamaged FILE` | Write undamaged reads here | one of the two `--out-*` required |
| `--library-type STR` | `ds`, `ss`, or `auto` (inferred) | `ds` |
| `--split-model STR` | `auto`, `bulk`, or `empirical` | `auto` |
| `--split-threshold F` | LLR threshold for a damaged call | `0.0` |
| `--damage-deam-sample N` | Max reads for the Pass-0 damage scan | `5000000` |
| `--model-bin FILE` | Reuse a full split model (Tier 0) | — |
| `--damage-json FILE` | Reuse a scalar model (Tier 1) | — |
| `--allow-model-mismatch` | Accept a model with a mismatched input digest | off |
| `-t, --threads N` | Threads | all available, capped at 16 for I/O |

## Output caveat

Every threshold branch logs that the split is a **pi-calibrated enrichment, not a
pure partition** (:266, :276): reference-free per-read damage separation is weak, so
`--out-damaged` is enriched for damage, not a clean partition. To measure the actual
per-read separation for a given run, set `FQDUP_SPLIT_SCORES` on a synthetic dataset
(see Score dump above) — the log no longer quotes a fixed AUC number.

## User-facing

Yes — top-level subcommand (`main.cpp:34,64`), not part of the short
documented workflow list (`main.cpp:42-49`) but listed in the subcommand table.
