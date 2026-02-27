# Damage-Aware Deduplication

## Motivation

Ancient DNA carries characteristic post-mortem damage: cytosine deamination
produces C→T substitutions near the 5' terminus and G→A (the complement on
the reverse strand) near the 3' terminus. Two reads from the same original
molecule can therefore differ only at these terminal positions. Standard
exact-hash deduplication would count them as distinct, inflating unique read
counts and downstream coverage estimates.

## Damage model

`fqdup` uses a DART-style exponential decay model:

```
P_CT(pos) = d_max_5 × exp(−λ_5 × pos) + bg
P_GA(pos) = d_max_3 × exp(−λ_3 × pos) + bg
```

where `pos` is the distance from the respective terminus and `bg` is the
background deamination rate.

Positions where `max(P_CT(pos), P_GA(pos)) > mask_threshold` are considered
part of the **damage zone** and are masked before hashing.

## Symmetric masking

The mask is applied identically at both ends using `max(P_CT, P_GA)` at
each position. This symmetry is required to preserve:

```
canonical_hash(seq) == canonical_hash(revcomp(seq))
```

Without symmetry, a forward-strand read and its reverse complement could
receive different masks and hash to different values, *increasing* unique
counts — the opposite of the intended effect.

## Automatic estimation (`--damage-auto`)

Pass 0 stride-samples the sorted non-extended file (every Nth read,
default N=1000; 100 k reads total). For each of the first 15 terminal
positions it measures the fraction of C→T (5') and G→A (3') events above
the baseline mid-read rate.

Fitting procedure:
- Baseline: mean mismatch rate in the middle third of reads
- OLS fit on positions 1–9 (position 0 excluded — unreliable in sorted files)
- Coverage-weighted (positions with fewer reads contribute less)
- λ clamped to [0.05, 0.5]

Estimated parameters are logged; use `--damage-dmax` / `--damage-lambda` to
override with known values.

## Manual specification

```bash
# Known damage parameters
fqdup derep ... \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.20 --damage-lambda3 0.30 \
  --damage-bg 0.02 \
  --mask-threshold 0.05
```

## PCR thermocycling damage

PCR cycling introduces additional C→T damage at an estimated rate of
`1.4×10⁻⁶ / base / cycle` (97% C→T, polymerase-independent). This is
accounted for in the PCR model alongside the library-preparation deamination.

Relevant flags:
```
--pcr-cycles INT       Number of PCR cycles
--pcr-efficiency FLOAT Efficiency per cycle, 0–1 (default: 1.0)
--pcr-error-rate FLOAT Substitution rate (default: 5.3e-7, Q5/HiFi)
                       Alternatives: Phusion 3.9e-6, KOD 1.2e-5, Taq 1.5e-4
```

## Effect on unique cluster counts

On a 25.8 M-read ancient DNA library (d_max≈0.07, λ≈0.29):

| Mode | Unique clusters |
|------|----------------|
| Standard | 5,582,073 |
| Damage-aware | 5,547,508 (−34,565) |
