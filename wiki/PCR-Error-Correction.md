# PCR error correction (Phase 3)

This page covers the second of two sources of false unique clusters in ancient
DNA deduplication: PCR copying errors. Post-mortem deamination is covered in
[[Damage-Aware-Deduplication]]. The two are complementary — damage masking
handles terminal variation, Phase 3 handles interior variation. See [[Home]]
for the full picture.

---

## The problem

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

**Caveat when running after `derep_pairs`:** the input to `derep` has already had
most PCR duplicates removed. The residual duplication reflects structural reads
(same non-extended sequence, different extended) rather than pure PCR amplification.
The estimated D_eff will be lower than the true PCR D_eff, making the adaptive
ceiling more conservative. For aDNA libraries this does not matter — parent counts
are low enough that the static `--errcor-max-count` dominates regardless. For
mitochondrial-enriched or high-depth modern libraries, specify `--pcr-cycles`
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

Sample `a88af16f35` (5.58 M reads input to `derep`, 91 bp mean length, Q5 library):

| Mode | Absorbed | Output sequences |
|------|----------|-----------------|
| `--error-correct` (naive, no damage protection) | 5,335 | 3,506,272 |
| `--error-correct` (with damage substitution filter) | 1,456 | 3,510,151 |

73% of what the naive version absorbed were C↔T or G↔A differences —
real damage-carrying reads, not PCR errors. The 1,456 retained absorptions
are A↔T and C↔G transversions only.

Phase 3 adds approximately 2 seconds to a 31-second `derep` run. It operates
entirely in memory on the index built during Pass 1.

---

## References

Potapov V, Ong JL (2017). Examining sources of error in PCR by single-molecule
sequencing. *PLoS ONE* 12(1): e0169774.
DOI: [10.1371/journal.pone.0169774](https://doi.org/10.1371/journal.pone.0169774)

- Table 3: polymerase substitution error rates (φ) for Q5, Phusion, KOD, Taq
- Supplementary: thermocycling deamination experiment (~1.4×10⁻⁶ C→T/base/cycle,
  polymerase-independent)
