# Damage-Aware Deduplication

## Background

Post-mortem DNA degradation produces characteristic patterns of cytosine
deamination. In double-stranded ancient DNA, deamination occurs preferentially
at single-stranded overhangs:

- **5' terminus**: C→T substitutions at positions 1–15, decaying exponentially
  toward the interior
- **3' terminus**: G→A substitutions (the complement of C→T on the opposing
  strand), with the same exponential decay

These damage patterns mean that two reads derived from the same original
molecule may differ by a C→T at position 2 on one read and not on another,
making them appear non-identical to an exact-hash deduplicator. Standard
deduplication would count them as distinct fragments, inflating unique cluster
counts and creating apparent sequence variation where there is none.

`fqdup` addresses this by masking damage-prone terminal positions before hashing.
Two reads that differ only within the damage zone hash identically and collapse
into the same cluster.

---

## Damage Model

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
- `bg` is the background deamination rate (default 0.02)

The **damage excess** at position `pos` — the deamination above background —
is `d_max × exp(-lambda × pos)`. Only this excess is used for masking; the
background rate `bg` reflects genomic AT-content, which is the same for both
reads in a duplicate pair and does not cause sequence differences between them.

---

## Automatic Parameter Estimation (`--damage-auto`)

Pass 0 stride-samples the input file to estimate the damage model without
requiring prior knowledge.

### Sampling

Every Nth read is inspected (default N=1000, targeting approximately 100 k
reads). The stride ensures reads are sampled from throughout the sorted file
rather than clustered at the start.

### Frequency measurement

For each sampled read of length L, for the first 15 positions:

```
5' end: freq_T[pos] = count(T at pos) / count(C or T at pos)
3' end: freq_A[pos] = count(A at pos) / count(A or G at 3'-pos)
```

These frequencies approximate `P_CT(pos)` and `P_GA(pos)` respectively.

### Baseline estimation

The background rate is estimated from the middle third of each read (positions
`L/3` to `2L/3`), matching DART's approach. Using the middle third avoids
contamination from terminal damage and adapter composition artifacts.

```
bg_5 = T/(T+C) in the middle third   [5' background]
bg_3 = A/(A+G) in the middle third   [3' background]
```

### Exponential fit

Coverage-weighted ordinary least-squares regression in log-space on positions
1–9 (position 0 is excluded — it can carry adapter or trimming artifacts in
sorted files). The slope gives `-lambda`; the intercept gives `log(d_max)`.

```
log(freq[pos] - baseline) ≈ log(d_max) - lambda × pos
```

Lambda is clamped to [0.05, 0.5] matching DART's bounds. The amplitude
`d_max` is initially estimated from position 1 (a robust anchor), then
refined by the OLS intercept.

Damage parameters are logged at the end of Pass 0:

```
Pass 0: damage estimation — sampled 100000 reads (every 1000th) ...
  5'-end: d_max=0.071 lambda=0.290 bg=0.301
  3'-end: d_max=0.011 lambda=0.250 bg=0.295
```

If both `d_max_5 < 0.02` and `d_max_3 < 0.02`, damage is considered
negligible and exact hashing is used (Pass 0 becomes a no-op).

---

## Damage Masking

Positions where the damage excess exceeds `--mask-threshold` (default 0.05)
are replaced with a neutral byte (`\x01` at 5' C/T positions, `\x02` at 3'
G/A positions) before hashing.

### Symmetric masking

The mask is computed symmetrically:

```
mask(pos) = max(P_CT_excess(pos), P_GA_excess(pos))
           = max(d_max_5 × exp(-lambda_5 × pos),
                 d_max_3 × exp(-lambda_3 × pos))
```

This symmetry is **required** to preserve the canonical hash invariant:

```
canonical_hash(seq) == canonical_hash(revcomp(seq))
```

Without symmetry, `apply_mask(seq)` and `apply_mask(revcomp(seq))` could
produce different masked sequences — causing a forward-strand read and its
reverse complement to hash to different values, and paradoxically *increasing*
unique cluster counts instead of decreasing them.

Because the mask uses the same function at both ends (the max of both damage
contributions at every position), it is inherently symmetric: masking position
`i` in the forward read uses the same threshold as masking position `L-1-i`
in the reverse complement.

---

## Manual Specification

If damage parameters are known (e.g. from a previous mapDamage or DART run):

```bash
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.20 --damage-lambda3 0.30 \
  --damage-bg 0.02 \
  --mask-threshold 0.05
```

Use `--damage-dmax` to set the same value for both 5' and 3' ends when they
are similar (e.g. for non-UDG libraries with symmetric damage).

---

## PCR Thermocycling Damage

PCR itself introduces damage independent of ancient DNA deamination. At
elevated temperatures (~97°C), cytosine deamination occurs at approximately
`1.4 × 10⁻⁶` events per base per cycle, predominantly producing C→T
substitutions. This thermocycling damage accumulates with library amplification
and is distinct from post-mortem deamination.

The total PCR error rate (substitution + thermocycling damage) is modeled as:

```
effective_doublings D = n_cycles × log2(1 + efficiency)
pcr_error_rate = phi × D   [phi = sub/base/doubling for the chosen polymerase]
```

| Polymerase | phi (sub/base/doubling) |
|-----------|------------------------|
| Q5 / HiFi | 5.3 × 10⁻⁷ (default) |
| Phusion / Pfu | 3.9 × 10⁻⁶ |
| KOD | 1.2 × 10⁻⁵ |
| Taq | 1.5 × 10⁻⁴ |

The PCR error rate is used to compute the expected number of mismatches between
two independently-amplified copies of the same molecule, informing the mismatch
tolerance report printed after Pass 0:

```
Expected mismatches (L=65): 0.51, 99th-pct tolerance: 3
```

This is informational — it does not change the masking threshold, only provides
context for interpreting the damage zone boundaries.

---

## Effect on Unique Cluster Counts

Benchmarked on a 25.8 M-read ancient DNA library with estimated
`d_max_5 ≈ 0.071`, `lambda_5 ≈ 0.29`:

| Mode | Unique clusters | Change |
|------|----------------|--------|
| Standard (exact hash) | 5,582,073 | — |
| Damage-aware | 5,547,508 | −34,565 (−0.6%) |

The damage-aware mode correctly merged 34,565 pairs of reads that had
been over-counted due to terminal deamination, reducing apparent unique
molecule count and improving coverage uniformity.

---

## Choosing `--mask-threshold`

The default `0.05` means any position where the damage excess exceeds 5% is
masked. This is conservative enough to catch the bulk of deamination events
while not masking so many positions that false positive collisions arise from
interior sequence similarity.

For highly damaged samples (e.g. cave sediments, permafrost), a slightly
higher threshold (0.08–0.10) can be appropriate to avoid masking too deep
into the read. For lightly damaged samples, the default works well.

The estimated mask boundaries are reported in the log:

```
--- Damage-Aware Deduplication ---
  5'-end d_max:  0.071
  5'-end lambda: 0.290
  Mask threshold: 0.050
  Expected mismatches (L=65): 0.51
```

Check the first few positions of the damage plot (if available from DART or
mapDamage2) to validate that the fitted model captures the observed profile
before proceeding with deduplication.
