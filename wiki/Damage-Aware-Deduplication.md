# Damage-aware deduplication

This page covers the first of two sources of false unique clusters in ancient
DNA deduplication: **post-mortem deamination**. The second — PCR copying
errors — is covered in [[PCR-Error-Correction]]. Both are handled by
`fqdup derep`; see [[Home]] for how they fit together.

---

## Background

Post-mortem DNA degradation produces characteristic patterns of cytosine
deamination. In double-stranded ancient DNA, deamination occurs preferentially
at single-stranded overhangs:

- **5' terminus**: C→T substitutions at positions 0–14, decaying exponentially
  toward the interior
- **3' terminus**: G→A substitutions (the complement of C→T on the opposing
  strand), with the same exponential decay

A read from an ancient molecule and its PCR copy can therefore differ by a C→T
at position 2, making them look non-identical to an exact-hash deduplicator.
Standard deduplication counts them as distinct fragments, inflating unique
cluster counts and creating apparent sequence variation where there is none.

`fqdup` addresses this by masking damage-prone terminal positions before
hashing. Two reads that differ only within the damage zone hash identically
and collapse into the same cluster.

---

## Damage model

`fqdup` uses the same exponential decay model as DART (Fernandez-Guerra et al.
2025) and mapDamage2 (Jónsson et al. 2013):

```
P_CT(pos) = d_max_5 × exp(-lambda_5 × pos) + bg    [C→T at 5' end]
P_GA(pos) = d_max_3 × exp(-lambda_3 × pos) + bg    [G→A at 3' end]
```

where:
- `pos` is the distance (0-indexed) from the relevant terminus
- `d_max_5`, `d_max_3` are the terminal substitution rates
- `lambda_5`, `lambda_3` control the decay rate (typical range 0.2–0.5)
- `bg` is the background substitution rate estimated from the middle of each read

Only the **damage excess** `d_max × exp(-lambda × pos)` is used for masking.
The background rate `bg` reflects the library's base composition, which is the
same for both reads in a duplicate pair and does not cause sequence differences
between them.

---

## Automatic parameter estimation (`--damage-auto`)

Pass 0 scans the entire input to fit the damage model before deduplication
begins.

### Frequency measurement

For each read of length L, the first 15 positions are examined:

```
5' end: freq_T[pos] = count(T at pos) / count(C or T at pos)
3' end: freq_A[pos] = count(A at pos) / count(A or G at 3'-pos)
```

These frequencies approximate `P_CT(pos)` and `P_GA(pos)`.

### Baseline estimation

The background rate comes from the middle third of each read (positions
`L/3` to `2L/3`), matching DART's approach. This avoids contamination from
terminal damage and adapter composition effects.

```
bg_5 = T/(T+C) in the middle third
bg_3 = A/(A+G) in the middle third
```

### Exponential fit

Coverage-weighted ordinary least-squares regression in log-space on positions
1–9 (position 0 is excluded — sorted files place the most common reads first,
which can give atypical terminal composition at position 0). The slope gives
`-lambda`; the intercept gives `log(d_max)`.

```
log(freq[pos] - baseline) ≈ log(d_max) - lambda × pos
```

Lambda is clamped to [0.05, 0.5] to match DART's bounds. The amplitude
`d_max` is initially anchored at position 1, then refined by the OLS intercept.

Damage parameters are logged at the end of Pass 0:

```
Pass 0: damage estimation — 5582073 reads processed (5582073 total)
  5'-end: d_max=0.099 lambda=0.246 bg=0.487
  3'-end: d_max=0.020 lambda=0.069 bg=0.509
```

If both `d_max_5 < 0.02` and `d_max_3 < 0.02`, damage is negligible and
exact hashing is used — Pass 0 adds no masking.

---

## Damage masking

Positions where the fitted damage excess exceeds `--mask-threshold` (default
0.05) are replaced with a neutral byte (`\x01` at 5' C/T positions, `\x02`
at 3' G/A positions) before hashing. The mask is precomputed as a boolean
array over the first 15 positions; no exponential evaluation happens at hash
time.

### Empirical masking

After fitting the model, `fqdup` records which positions actually exceeded the
threshold from the observed frequencies — not just which positions the fitted
curve predicts. This avoids over-masking when the OLS fit overestimates damage
at lightly-affected positions.

```
mask_pos[p] = true  iff  (freq_T[p] - bg_5 > threshold)
                       OR (freq_A[p] - bg_3 > threshold)
```

In practice this typically masks 1–3 terminal positions per end, depending on
library age and preparation protocol.

### Symmetric masking

The mask applies symmetrically to both ends using the same position array:

```
For position i from 5' end:  mask if mask_pos[i]  and base is C or T → \x01
For position i from 3' end:  mask if mask_pos[i]  and base is G or A → \x02
```

This symmetry is required to preserve the canonical hash invariant:

```
canonical_hash(seq) == canonical_hash(revcomp(seq))
```

Without it, `apply_mask(seq)` and `apply_mask(revcomp(seq))` would produce
different masked sequences — a forward-strand read and its reverse complement
would hash to different values, paradoxically *increasing* unique cluster counts.

---

## Manual specification

If damage parameters are known (from a previous mapDamage2 or DART run):

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.20 --damage-lambda3 0.30 \
  --damage-bg 0.02 \
  --mask-threshold 0.05
```

Use `--damage-dmax` to set the same value for both ends when they are similar
(e.g. non-UDG libraries with symmetric damage).

---

## Effect on unique cluster counts

Benchmarked on sample `a88af16f35` — a 25.8 M-read ancient DNA library with
`d_max_5 ≈ 0.099`, `lambda_5 ≈ 0.25`, run through the full pipeline:

| Step | Unique clusters | Change |
|------|----------------|--------|
| `derep_pairs` (structural dedup, 25.8 M pairs in) | 5,582,073 | — |
| `derep` standard (exact hash) | 3,531,821 | — |
| `derep --damage-auto` (empirical masking) | 3,511,607 | −20,214 (−0.6%) |
| `derep --damage-auto --error-correct` | 3,506,272 | −25,549 (−0.7%) |

Damage-aware mode merged 20,214 additional clusters that differed only by
terminal deamination. Error correction absorbed a further 5,335 low-count
PCR-error clusters.

---

## Choosing `--mask-threshold`

The default `0.05` masks positions where the observed damage excess exceeds 5%.
This catches the bulk of deamination events without masking so deep that
unrelated molecules in the interior of a read start colliding.

For heavily damaged libraries (cave sediments, permafrost), a slightly higher
threshold (0.08–0.10) reduces over-masking. For lightly damaged libraries,
the default works well.

The masked positions are reported after Pass 0:

```
--- Damage-Aware Deduplication ---
  5'-end d_max:  0.098842
  3'-end d_max:  0.019759
  5'-end lambda: 0.246287
  3'-end lambda: 0.069124
  Mask threshold:0.050000
  Masked positions: 1 (1 bp each end)
  Expected mismatches (L=91): 1.44, 99th-pct tolerance: 5
```

Cross-check the masked position count against a damage plot from DART or
mapDamage2 to confirm the threshold captures the observed damage zone.
