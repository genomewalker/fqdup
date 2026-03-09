# PCR error correction (Phase 3)

This page covers the second of two sources of false unique clusters in ancient
DNA deduplication: PCR copying errors. Post-mortem deamination is covered in
[[Damage-Aware-Deduplication]]. The two are complementary — damage masking
handles terminal variation, Phase 3 handles interior variation. See [[Home]]
for the full picture.

In the standard pipeline, `derep` receives the **merged read output of
`derep_pairs`** — reads that were collapsed from R1+R2 by fastp and then
deduplicated against their `fqdup extend`-assembled counterparts. Most PCR
duplicates have already been removed at that stage; `derep`'s Phase 3 targets
the residual PCR errors that survived because they happened to produce divergent
extensions (and were therefore treated as distinct pairs by `derep_pairs`).

---

## Background

PCR amplification introduces substitution errors at a rate that depends on the
polymerase used. High-fidelity polymerases such as Q5 make roughly
5.3×10⁻⁷ substitutions per base per doubling; Taq makes roughly
1.5×10⁻⁴ (Potapov & Ong 2017). Most errors occur in early cycles, when the
pool of template molecules is small, and therefore propagate into a modest
cluster of copies carrying the wrong base. After deduplication, these clusters
appear as genuine unique sequences — count 2 or 3 rather than count 1, which
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
`--errcor-max-count`, default 5) and **children** (count at or below that
threshold). A child is absorbed into a parent when:

1. Their sequences differ by at most one substitution **outside** the damage zone
2. The parent's count is at least `--errcor-ratio` × the child's count (default 50×)
3. The single mismatch (if any) is **not** a damage-consistent substitution (see below)

---

## Damage substitution protection

The damage zone exclusion skips positions already flagged by the masking model,
but sub-threshold positions still carry real C→T and G→A signal. And oxidative
damage (8-oxoG, Channel C in DART) produces G→T transversions uniformly across
the entire read — no position-based masking can protect those.

Phase 3 therefore applies a substitution-type filter: a child is **never**
absorbed if the single mismatch is associated with a known DNA damage mechanism,
regardless of its position:

| Substitution | Damage mechanism |
|---|---|
| C↔T | Ancient deamination (5′ C→T) |
| G↔A | Ancient deamination (3′ G→A, complementary strand) |
| G↔T | Oxidative 8-oxoG (uniform across read) |
| C↔A | Oxidative 8-oxoG (complementary strand) |

Only **A↔T** and **C↔G** mismatches — the two transversion pairs with no known
DNA damage mechanism — are eligible for absorption. In practice this protects
73% of what a naive Phase 3 would absorb on an ancient DNA library with
d_max ≈ 0.10 (sample a88af16f35, see benchmarks below).

---

## PCR error model and adaptive thresholds

### Error rate formula

Following Potapov & Ong 2017, the total substitution error probability per
base position over a full PCR run is:

```
ε_total = φ × D_eff

where:
  φ     = polymerase error rate (sub/base/doubling)
  D_eff = n × log₂(1 + E)     effective doublings
  n     = number of PCR cycles
  E     = per-cycle efficiency (0–1; 1.0 = ideal doubling)
```

The effective doublings formula follows directly from PCR kinetics: if each
cycle amplifies a fraction E of molecules, after n cycles the copy number is
(1+E)ⁿ, giving log₂((1+E)ⁿ) = n × log₂(1+E) effective doublings.

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
heat-induced cytosine deamination. Potapov & Ong (2017) measured this directly
using mock thermocycling without polymerase and found approximately
2.3×10⁻⁵ C→T events per base over 16 cycles, i.e. roughly **1.4×10⁻⁶ per
base per cycle**, of which ~97% were C→T transitions. Because these arise from
heat damage rather than polymerase misincorporation, they are biologically
indistinguishable from ancient deamination and are correctly protected by the
damage substitution filter (C↔T → not absorbed).

### Adaptive child ceiling

The static `--errcor-max-count 5` is calibrated for parents with
~10,000 reads. For a parent with count P, the expected count of a cluster
arising from a specific PCR error (a particular 1-substitution variant) is,
under the equal-rates substitution model:

```
E[count_child] = P × φ × D_eff / 3
```

The /3 term reflects the equal-rates assumption: at any given base, each of
the three alternative bases is equally likely to be substituted in (Potapov &
Ong 2017 note that true substitution spectra show transition/transversion bias,
so this is a first-order approximation). When `--pcr-cycles` is specified,
fqdup computes this per parent and raises the effective child ceiling
accordingly:

```
effective_max(P) = max(--errcor-max-count, ⌈P × φ × D_eff / 3⌉)
```

The global pre-filter (for efficiency) uses the observed maximum parent count.

**When `--pcr-cycles` is not given**, fqdup estimates D_eff automatically from the
observed duplication ratio after Pass 1:

```
D_eff_estimated = log₂(total_reads / unique_reads)
```

This follows directly from PCR kinetics: mean copies per starting molecule after n
cycles with efficiency E is (1+E)ⁿ, so total/unique ≈ (1+E)ⁿ and
D_eff = log₂((1+E)ⁿ) = log₂(total/unique). The estimate is logged:

```
Phase 3: D_eff=2.21 estimated from duplication ratio 25000000/5600000
  (use --pcr-cycles for explicit value)
```

**Caveat when running after `derep_pairs`:** most PCR duplicates are already
removed before `derep` sees the data. The residual duplication that remains
comes from PCR copies whose `fqdup extend` assemblies diverged slightly —
`derep_pairs` treated them as distinct molecules; `derep` collapses them on the original
merged sequence. Because most amplification has already been collapsed, D_eff
estimated here is much lower than the true PCR D_eff, making the adaptive
ceiling conservative. For typical aDNA libraries this does not matter — parent
counts are low enough that the static `--errcor-max-count` dominates regardless.
For mitochondrial-enriched or high-depth modern libraries, supply `--pcr-cycles`
explicitly for an accurate adaptive ceiling.

**Example (Q5, 30 cycles, E=1.0):**

| Parent count | E[count_child] | effective_max |
|---|---|---|
| 1,000 | 0.005 | 5 (static) |
| 10,000 | 0.053 | 5 (static) |
| 100,000 | 0.53 | 5 (static) |
| 1,000,000 | 5.3 | 6 |
| 10,000,000 | 53 | 53 |

For ancient DNA libraries, parent counts rarely exceed 10,000 per cluster, so
the static threshold dominates. For mitochondrial DNA or highly amplified
modern libraries, `--pcr-cycles` becomes important.

---

## The 3-way pigeonhole algorithm

Finding all pairs of sequences that differ by at most one substitution is
O(n²) with naive pairwise comparison. Phase 3 uses a pigeonhole argument to
reduce this to effectively O(n).

### Principle

Divide the interior region into three equal parts. If two sequences of the same
length differ by **exactly one** substitution, that substitution falls in one
of the three parts — so the other two parts must be identical. Three hash
indexes, each covering a different pair of parts, are sufficient to find every
Hamming-distance-1 neighbour:

```
interior: [─── part 0 ─── | ─── part 1 ─── | ─── part 2 ───]

mismatch in part 0 → parts 1 and 2 identical → detected by index (1,2)
mismatch in part 1 → parts 0 and 2 identical → detected by index (0,2)
mismatch in part 2 → parts 0 and 1 identical → detected by index (0,1)
exact match (H=0)  → all three indexes fire
```

### Worked example

Interior (15 bp):

```
Position:  0  1  2  3  4  |  5  6  7  8  9  | 10 11 12 13 14
Parent:    A  C  G  A  A  |  C  G  A  A  C  |  G  A  A  T  G
Child:     A  C  G  A  A  |  C  G  T  A  C  |  G  A  A  T  G
                                    ^ mismatch at pos 7 (A→T, part 1)
```

Hash the three parts:

```
Parent:  h0 = hash(ACGAA)   h1 = hash(CGAAC)   h2 = hash(GAATG)
Child:   h0 = hash(ACGAA)   h1 = hash(CGTAC)   h2 = hash(GAATG)
              ↑ same              ↑ different        ↑ same
```

The pair-key `(h0, h2)` is identical for parent and child — index (0,2) fires.
The candidate is retrieved and the exact Hamming check confirms one mismatch.
A→T is not a damage substitution, so the child is absorbed.

If the mismatch were C→T instead (e.g. position 9: C→T), the same index would
fire, but `is_damage_sub('C','T')` returns true and the child is **kept**.

### Complexity

Building the three indexes is O(n_parents). Each child queries three buckets —
O(1) per bucket under the assumption of few collisions. The full character-level
Hamming check runs only on the small set of candidates returned per child. In
practice Phase 3 completes in 1–3 seconds for 5M sequences.

---

## Parameters

`--errcor-ratio` (default 50) controls how large the count gap must be. A ratio
of 50 means a count=1 child needs a count≥50 parent. Raise it (100–200) if
your library contains genuine rare variants at low coverage, such as in
single-cell or low-input sequencing.

`--errcor-max-count` (default 5) sets the static child ceiling. This is always
a hard lower bound; with `--pcr-cycles` the effective ceiling may be higher for
large parents (see adaptive thresholds above).

`--errcor-bucket-cap` (default 64) prevents low-complexity regions — poly-A
stretches, microsatellites — from producing enormous bucket collisions in the
pair-hash index.

`--pcr-cycles`, `--pcr-efficiency`, `--pcr-error-rate` (default Q5: 5.3e-7)
enable adaptive thresholds. Without `--pcr-cycles`, the PCR model is used only
for the expected-mismatch log line; the static `--errcor-max-count` applies.

---

## Benchmarks

### Damage filter: single library

Sample `a88af16f35` (5.58 M reads input to `derep`, 91 bp mean length, Q5 library):

| Mode | Absorbed | Output sequences |
|------|----------|-----------------|
| `--error-correct` (naive, no damage protection) | 5,335 | 3,506,272 |
| `--error-correct` (with damage substitution filter) | 1,456 | 3,510,151 |

73% of what the naive version absorbed were C↔T or G↔A — real damage signal,
not PCR errors. The 1,456 absorptions are A↔T and C↔G transversions only.

Phase 3 adds approximately 2 seconds to a 31-second `derep` run. It operates
entirely in memory on the index built during Pass 1.

### Batch run: 31 ancient DNA libraries

All libraries processed with `--damage-auto --error-correct`, Q5 defaults
(φ = 5.3×10⁻⁷), D_eff estimated automatically from duplication ratio.

**Input is post-`derep_pairs`** (`*.derep.non.fq.gz`): fastp-merged reads whose
`fqdup extend`-assembled pairs have already been deduplicated. The duplication
remaining at this stage comes from PCR copies that produced divergent assemblies
and therefore survived `derep_pairs`. D_eff here is the residual D_eff, not the
original PCR amplification depth — hence values are much lower than the true library D_eff.

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
The cutoff is not a hard threshold — it reflects the adaptive ceiling formula:
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
