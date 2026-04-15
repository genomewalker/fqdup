# PCR error correction (Phase 3)

This page covers the second of two sources of false unique clusters in ancient
DNA deduplication: PCR copying errors. Post-mortem deamination is covered in
[[Damage-Aware-Deduplication]]. The two are complementary, damage masking
handles terminal variation, Phase 3 handles interior variation. See [[Home]]
for the full picture.

`derep` can be run directly on sorted merged reads or after `derep_pairs` in
the full paired pipeline. When following `derep_pairs`, most PCR duplicates
have already been removed; Phase 3 targets the residual errors that survived
because they produced divergent extensions (and were therefore treated as
distinct pairs by `derep_pairs`). When run directly, Phase 3 operates on the
complete amplification depth.

---

## Background

PCR amplification introduces substitution errors at a rate that depends on the
polymerase used. High-fidelity polymerases such as Q5 make roughly
5.3×10⁻⁷ substitutions per base per doubling; Taq makes roughly
1.5×10⁻⁴ (Potapov & Ong 2017). Most errors occur in early cycles, when the
pool of template molecules is small, and therefore propagate into a modest
cluster of copies carrying the wrong base. After deduplication, these clusters
appear as genuine unique sequences, count 2 or 3 rather than count 1, which
is enough to survive any simple singleton filter.

The distinguishing feature of a PCR error cluster is the count ratio: the
correct sequence has been amplified from the original molecule and has a high
read count, while the erroneous copy's cluster is a small fraction of that.
A single interior substitution separating a count=2 cluster from a count=500
cluster is almost never a distinct biological molecule.

---

## Approach

Phase 3 runs entirely in memory on the deduplication index built during Pass 1
— no additional file I/O. It classifies clusters as **parents** (count above
`--errcor-min-parent`, default 3) and **children** (count at or below that
threshold). For each child, Phase 3 finds all neighbours within Hamming
distance ≤ 2 among parents and absorbs them unless the **SNP veto** fires.

Two absorption paths:

- **H=1**: any single non-damage interior substitution. Applies to all children
  regardless of count.
- **H=2**: two non-damage interior substitutions (both must be A↔T or C↔G
  transversions). Only applies to children with count ≤ `--errcor-max-h2-count`
  (default 2). This catches double-error reads from chimera formation, nick
  translation, or clustered oxidative damage that escape H=1 detection.

The SNP veto protects a child when the mismatch appears recurrently across
multiple reads at the same position in the same parent:

```
snp_veto = (sig_count ≥ --errcor-snp-min-count)
        AND (sig_count / parent_count ≥ --errcor-snp-threshold)
```

`sig_count` is the weighted count of all child clusters carrying the same
`(position, alt_base)` mismatch relative to this parent. For H=2 children the
SNP veto is applied independently to both mismatch positions — the child is
protected if either position has population support above the threshold.

---

## Damage-zone handling

Phase 3 includes a special case for C↔T and G↔A mismatches at positions
immediately adjacent to the damage zone (`--mask-threshold` edge). These
positions are not masked (their excess deamination rate is below the threshold)
but they still carry real ancient DNA signal. When `--collapse-damage` is active,
such mismatches bypass the SNP veto and are always absorbed — they are treated
as residual deamination artefacts rather than true SNPs.

This bypass applies only to C↔T and G↔A at the damage-zone edge.

G↔T and C↔A mismatches (consistent with 8-oxoG oxidative damage) are also
unconditionally protected from absorption by `is_damage_sub` — they are never
treated as PCR errors regardless of position, count, or whether
`--collapse-damage` is active. Only A↔T and C↔G transversions can be absorbed
by Phase 3.

---

## PCR error model (logging only)

Phase 3 does not use the PCR error rate to gate absorption decisions; the SNP
veto threshold alone determines whether a child is kept. The PCR model is used
only to compute an expected error rate that is printed in the log, providing
context for interpreting the number of absorptions.

### Error rate formula

Following Potapov & Ong 2017:

```
ε_total = φ × D_eff

where:
  φ     = polymerase error rate (sub/base/doubling)
  D_eff = n × log₂(1 + E)     effective doublings
  n     = number of PCR cycles
  E     = per-cycle efficiency (0–1; 1.0 = ideal doubling)
```

**When `--pcr-cycles` is not given**, D_eff is estimated from the observed
duplication ratio after Pass 1:

```
D_eff_estimated = log₂(total_reads / unique_reads)
```

This estimate is logged:

```
Phase 3: D_eff=2.21 estimated from duplication ratio 25000000/5600000
  (use --pcr-cycles for explicit value)
```

**Caveat when running after `derep_pairs`:** most PCR duplicates are already
removed before `derep` sees the data, so D_eff estimated here is much lower
than the true library D_eff. This affects only the log message; absorption
decisions are unaffected.

### Polymerase error rates (Potapov & Ong 2017)

| Polymerase | φ (sub/base/doubling) | Typical use |
|---|---|---|
| Q5 (NEB) | 5.3×10⁻⁷ | Ancient DNA, high-fidelity |
| Phusion | 3.9×10⁻⁶ | General high-fidelity |
| KOD | 1.2×10⁻⁵ | High-fidelity |
| Taq | 1.5×10⁻⁴ | Older protocols |

All values from Table 3, Potapov & Ong (2017), DOI: 10.1371/journal.pone.0169774.

### Thermocycling deamination

A separate, polymerase-independent source of C→T changes during PCR is
heat-induced cytosine deamination. Potapov & Ong (2017) measured approximately
2.3×10⁻⁵ C→T events per base over 16 cycles (~1.4×10⁻⁶ per base per cycle,
~97% C→T transitions). These are biologically indistinguishable from ancient
deamination; the damage-zone bypass (see above) handles them at the terminal
positions where they concentrate.

---

## The 4-way pigeonhole algorithm

Finding all pairs of sequences that differ by at most two substitutions is
O(n²) with naive pairwise comparison. Phase 3 uses a pigeonhole argument to
reduce this to effectively O(n).

### Principle

Divide the interior region into four equal parts. If two sequences of the same
length differ by **at most two** substitutions, at most two parts contain
mismatches, so at least two parts are identical. Index all C(4,2) = 6 pair
combinations of parts — one of the six pair-keys must fire:

```
interior: [─ part 0 ─ | ─ part 1 ─ | ─ part 2 ─ | ─ part 3 ─]

H=1: mismatch in part k → the other 3 parts all match → 3 pair-keys fire
H=2: mismatches in parts j,k → the other 2 parts match → index (j̄,k̄) fires
H=0: all 6 pair-keys fire
```

The 6 indexes cover all pairs: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).

### Worked example — H=2

Interior (20 bp):

```
Position:  0  1  2  3  4 | 5  6  7  8  9 | 10 11 12 13 14 | 15 16 17 18 19
Parent:    A  C  G  A  A | C  G  A  A  C |  G  A  A  T  G |  C  T  A  G  C
Child:     A  C  G  A  A | C  G  T  A  C |  G  A  A  T  G |  C  T  A  G  C
                                   ^ pos 7 (A→T)     part 2 and 3 identical
```

Wait — one mismatch, so H=1 path fires via any of the 3 pair-keys that
exclude part 1. For a true H=2 example:

```
Child2:    A  C  G  A  A | C  G  T  A  C |  G  A  A  T  G |  C  T  G  G  C
                                   ^ pos 7 (A→T)              ^ pos 17 (A→G)
           part 0 same      part 1 differs   part 2 same       part 3 differs
```

Parts 0 and 2 are identical → pair-key index (0,2) fires → candidate retrieved
→ Hamming check finds 2 mismatches → both are A↔T / A↔G (xr=3 transversions,
non-damage) → child2 count ≤ max_h2_count → absorbed.

### Complexity

Building the 6 indexes is O(n_parents). Each child queries 12 buckets (6 pair
combinations × 2 orientations) — O(1) per bucket. The exact Hamming check runs
only on the small candidate set. In practice Phase 3 completes in 1–3 seconds
for 5 M sequences.

---

## Parameters

`--errcor-min-parent INT` (default 3) — sequences with count above this are
indexed as parents and treated as potential absorption targets. Sequences with
count at or below this are children (candidate PCR errors). Raise this if your
library contains genuine low-count molecules (e.g. single-cell, UMI-tagged).

`--errcor-snp-threshold FLOAT` (default 0.20) — SNP veto ratio: if the
weighted count of reads carrying the same `(position, alt_base)` mismatch
relative to the parent is at least this fraction of the parent count, the child
is protected as a potential true SNP. Lower values absorb more aggressively;
raise to 0.30–0.50 for libraries with genuine heterozygosity or rare variants.

`--errcor-snp-min-count INT` (default 2) — minimum absolute `sig_count`
required for the SNP veto to fire. Raise to 3–5 for libraries where genuine
rare variants are expected at low absolute count.

`--errcor-max-h2-count INT` (default 2) — H=2 absorption (double-error reads)
is only attempted for children with read count ≤ this value. The expected count
of a genuine double-PCR-error is ~φ² × D_eff² × parent_count, which is
astronomically small; a read appearing 3+ times at H=2 from a parent is almost
certainly a real variant, chimera at elevated frequency, or genuine molecule.

`--errcor-bucket-cap INT` (default 64) — prevents low-complexity regions
(poly-A, microsatellites) from producing enormous bucket collisions in the
pair-hash index.

`--pcr-error-rate FLOAT` (default Q5: 5.3e-7), `--pcr-cycles INT`,
`--pcr-efficiency FLOAT` — parameterise the PCR error model used for the
logged D_eff estimate. These values do not affect absorption decisions.

---

## Benchmarks

### Single library

Sample `a88af16f35` (5.58 M reads input to `derep`, 91 bp mean length, Q5 library,
post-`derep_pairs`):

Phase 3 adds approximately 2 seconds to a 31-second `derep` run. It operates
entirely in memory on the index built during Pass 1.

### Batch run: 31 ancient DNA libraries

All libraries processed with `--collapse-damage --error-correct`, Q5 defaults
(φ = 5.3×10⁻⁷), D_eff estimated automatically from duplication ratio.

**Input is post-`derep_pairs`** (`*.derep.non.fq.gz`): fastp-merged reads whose
`fqdup extend`-assembled pairs have already been deduplicated. The duplication
remaining at this stage comes from PCR copies that produced divergent assemblies
and therefore survived `derep_pairs`. D_eff here is the residual D_eff, not the
original PCR amplification depth, hence values are much lower than the true library D_eff.

| Library | Total reads | Unique (P1) | Unique (final) | Dup% | d_max | λ | D_eff | Absorbed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0267130b40 | 15,242,307 | 10,624,058 | 10,624,058 | 30.3 | 0.000 | 0.20 | 0.52 | 0 |
| 0c7e1a3199 | 31,727,991 | 27,011,848 | 27,011,848 | 14.9 | 0.416 | 0.38 | 0.23 | 0 |
| 2a3189abd3 | 28,062,190 | 17,497,629 | 17,496,559 | 37.7 | 0.176 | 0.25 | 0.68 | 1,070 |
| 32cf1bb9a2 | 91,853,574 | 82,727,024 | 82,727,024 | 9.9 | 0.167 | 0.15 | 0.15 | 0 |
| 4940db5c8d | 41,256,425 | 32,439,601 | 32,439,601 | 21.4 | 0.355 | 0.40 | 0.35 | 0 |
| 521588e724 | 94,736,098 | 69,697,677 | 69,697,677 | 26.4 | 0.059 | 0.08 | 0.44 | 0 |
| 55b0793258 | 136,474,361 | 129,837,772 | 129,837,772 | 4.9 | 0.550 | 0.50 | 0.07 | 0 |
| 68825e1df0 | 5,773,793 | 3,782,602 | 3,780,481 | 34.5 | 0.217 | 0.30 | 0.61 | 2,121 |
| 6fb2395c87 | 2,205,562 | 2,076,532 | 2,076,532 | 5.9 | 0.015 | 0.05 | 0.09 | 0 |
| 75c7be7787 | 155,863,098 | 150,386,512 | 150,386,512 | 3.5 | 0.381 | 0.43 | 0.05 | 0 |
| 7646a78f08 | 121,809,247 | 109,821,675 | 109,821,675 | 9.8 | 0.256 | 0.22 | 0.15 | 0 |
| 792a5d1378 | 82,567,031 | 72,228,459 | 72,228,459 | 12.5 | 0.590 | 0.40 | 0.19 | 0 |
| 7abec0fc9d | 23,093,793 | 19,519,181 | 19,519,181 | 15.5 | 0.271 | 0.38 | 0.24 | 0 |
| 83853d2471 | 50,277,392 | 30,919,799 | 30,918,892 | 38.5 | 0.159 | 0.25 | 0.70 | 907 |
| a88af16f35 | 5,582,073 | 3,511,607 | 3,510,151 | 37.1 | 0.193 | 0.25 | 0.67 | 1,456 |
| a97f194584 | 138,638,296 | 132,959,230 | 132,959,230 | 4.1 | 0.366 | 0.41 | 0.06 | 0 |
| abfe73874c | 36,744,399 | 30,825,182 | 30,825,182 | 16.1 | 0.411 | 0.39 | 0.25 | 0 |
| b984d65658 | 41,913,674 | 27,900,335 | 27,899,474 | 33.4 | 0.150 | 0.29 | 0.59 | 861 |
| bd5c178184 | 121,462,889 | 113,302,639 | 113,302,639 | 6.7 | 0.183 | 0.22 | 0.10 | 0 |
| ca85b7a308 | 140,960,298 | 128,936,359 | 128,936,359 | 8.5 | 0.386 | 0.35 | 0.13 | 0 |
| d3d27b4213 | 74,766,739 | 65,848,465 | 65,848,359 | 11.9 | 0.287 | 0.21 | 0.18 | 106 |
| dc0d65f3e9 | 57,976,670 | 49,143,249 | 49,143,249 | 15.2 | 0.431 | 0.36 | 0.24 | 0 |
| dc76c20bfb | 9,718,337 | 5,775,523 | 5,774,048 | 40.6 | 0.123 | 0.20 | 0.75 | 1,475 |
| dfb2272499 | 74,320,027 | 69,044,095 | 69,044,095 | 7.1 | 0.386 | 0.24 | 0.11 | 0 |
| e30685fafc | 97,072,045 | 91,710,984 | 91,710,984 | 5.5 | 0.439 | 0.50 | 0.08 | 0 |
| f024f0d6ad | 72,441,677 | 63,476,566 | 63,476,566 | 12.4 | 0.514 | 0.41 | 0.19 | 0 |
| 3e3ce8e9f7 | 262,729,902 | 251,140,430 | 251,140,430 | 4.4 | 0.249 | 0.303 | 0.065 | 0 |
| 73370e4c6c | 222,366,995 | 214,908,224 | 214,908,205 | 3.4 | 0.271 | 0.311 | 0.049 | 19 |
| 7fe94f44b2 | 266,866,219 | 260,057,150 | 260,057,046 | 2.6 | 0.370 | 0.318 | 0.037 | 104 |
| d12d52cc3a | 203,940,806 | 192,590,764 | 192,590,764 | 5.6 | 0.219 | 0.315 | 0.082 | 0 |
| ed3a67f780 | 263,137,179 | 255,773,761 | 255,773,761 | 2.8 | 0.411 | 0.500 | 0.040 | 0 |

**D_eff** is estimated from duplication ratio (log₂(total/unique_P1)) unless
`--pcr-cycles` is given. **Absorbed** = PCR error sequences removed by Phase 3
(after damage filter). **Unique (P1)** = after damage-aware deduplication but
before Phase 3; **Unique (final)** = output sequence count.

Across all 31 libraries, error correction fires only when D_eff ≥ ~0.6.
The cutoff is not a hard threshold, it reflects the adaptive ceiling formula:
at low D_eff, E[count_child] < 1 for all realistic parent counts, so the
static `--errcor-max-count 5` ceiling dominates and no child qualifies for
absorption. At D_eff ≈ 0.7 and above, the ceiling rises above the static
floor for high-count parents, allowing the algorithm to absorb genuine errors
while ignoring singletons from low-count parents.

---

## References

Potapov V, Ong JL (2017). Examining sources of error in PCR by single-molecule
sequencing. *PLoS ONE* 12(1): e0169774.
DOI: [10.1371/journal.pone.0169774](https://doi.org/10.1371/journal.pone.0169774)

- Table 3: polymerase substitution error rates (φ) for Q5, Phusion, KOD, Taq
- Supplementary: thermocycling deamination experiment (~1.4×10⁻⁶ C→T/base/cycle,
  polymerase-independent)
