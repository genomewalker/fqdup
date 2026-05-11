# Damage-aware deduplication

This page covers the first of two sources of false unique clusters in ancient
DNA deduplication: post-mortem deamination. PCR copying errors are covered in
[[PCR-Error-Correction]]. Both are handled by `fqdup derep`; see [[Home]] for
how they fit together.

`derep` operates on fastp-merged reads: the original collapsed R1+R2 molecules
on which the damage signal is present and should be modelled. In the full
paired pipeline these reads arrive from `derep_pairs`; `derep` can also be
run directly on sorted merged reads without a preceding `derep_pairs` step.

---

## Background

Ancient DNA carries characteristic substitution patterns at the ends of
sequenced fragments. Cytosine residues in single-stranded overhangs deaminate
to uracil over time, producing C→T substitutions at the 5' terminus. The same
process on the complementary strand appears as G→A at the 3' terminus. Both
signals decay exponentially from the ends toward the read interior.

When two reads come from different original ancient molecules at the same genomic
locus — one deaminated at the 5' terminus, one not — exact-match deduplication
treats them as different fragments. Deamination is a post-mortem chemical
modification fixed in the ancient template before library preparation; every PCR
copy of a deaminated molecule carries the C→T change, so reads with and without
the substitution at the same position must come from distinct original molecules.
The unique molecule count is inflated by however many locus-copies were sequenced
both with and without terminal damage, which in a heavily-amplified library
can be substantial.

---

## The model

`fqdup` uses the same exponential decay model as DART (Fernandez-Guerra et al.
2025) and mapDamage2 (Jónsson et al. 2013). The expected C→T frequency at
position `pos` from the 5' end is:

```
P_CT(pos) = d_max_5 × exp(-lambda_5 × pos) + bg
```

and analogously for G→A from the 3' end with `d_max_3` and `lambda_3`.
The parameters `d_max` (terminal rate), `lambda` (decay rate, typically 0.2–0.5),
and `bg` (background from the library's base composition) are either estimated
automatically or supplied manually.

Only the damage excess, `d_max × exp(-lambda × pos)`, drives masking. The
background `bg` is the same for both copies of a duplicate pair and does not
cause sequence differences between them.

---

## Automatic estimation (`--collapse-damage`)

Pass 0 scans all reads before deduplication begins. The per-thread
accumulation, background estimation, exponential fit, and library-type
classification all run inside [libtaph](https://github.com/genomewalker/libtaph)
via `taph::FrameSelector`; fqdup only owns the scan loop and the subsequent
masking decision.

For each read, libtaph tallies the T/(T+C) frequency at each of the first 15
positions from the 5' end, and A/(A+G) from the 3' end. The background rate is
taken from the middle third of the read, avoiding terminal damage and adapter
composition bias, consistent with DART's approach.

After scanning, a coverage-weighted ordinary least-squares fit in log-space
(positions 1-9; position 0 is excluded because it can carry first-cycle or
ligation artifacts) gives the exponential parameters. Lambda is clamped to
[0.05, 0.5]. The excess T/(T+C) at each position is divided by `(1 - bg)`, the
C fraction of the T+C pool, to convert from raw frequency excess into estimated
P(C->T) deamination rate, matching DART's convention and making `d_max`
directly comparable to metaDMG output. The estimated parameters are
logged:

```
Pass 0: damage estimation, 5582073 reads processed (5582073 total)
  5'-end: d_max=0.193 lambda=0.246 bg=0.487
  3'-end: d_max=0.040 lambda=0.069 bg=0.509
```

If both `d_max_5` and `d_max_3` are below 0.02, damage is treated as
negligible and no masking is applied, standard exact hashing is used.

---

## Masking

After libtaph returns the finalized `SampleDamageProfile`, fqdup records which
positions actually exceeded the `--mask-threshold` (default 0.05) in the
observed frequencies, not just which positions the fitted curve predicts. This empirical approach avoids over-masking
when OLS overestimates damage at lightly-affected positions. In practice, 1–3
terminal positions per end are masked.

At hash time, each masked position is replaced with a neutral byte (`\x01` at
C or T in the 5' zone, `\x02` at G or A in the 3' zone). Two reads from the
same molecule that differ only within the masked zone then produce identical
hashed sequences and collapse into the same cluster.

The mask uses the same position indices from both ends, position `i` from the
5' end and position `i` from the 3' end share the same mask flag. This symmetry
is required to preserve the canonical hash invariant:

```
canonical_hash(seq) == canonical_hash(revcomp(seq))
```

Without it, a forward-strand read and its reverse complement would hash
differently after masking, paradoxically increasing unique cluster counts
instead of reducing them.

The masked positions and expected mismatch count are reported after Pass 0:

```
--- Damage-Aware Deduplication ---
  5'-end d_max:  0.192754
  3'-end d_max:  0.040259
  5'-end lambda: 0.246287
  Mask threshold:0.050000
  Masked positions: 1 (1 bp each end)
  Expected mismatches (L=91): 1.44, 99th-pct tolerance: 5
```

---

## Manual parameters

If damage parameters are already known from a previous mapDamage2 or DART run,
skip Pass 0 and supply them directly:

```bash
fqdup derep \
  -i merged.deduped.fq.gz -o merged.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.20 --damage-lambda3 0.30 \
  --mask-threshold 0.05
```

Use `--damage-dmax` to set the same value for both ends (e.g. non-UDG libraries
with symmetric damage).

---

## Effect on cluster counts

Benchmarked on sample `a88af16f35`, 5.58 M merged reads input to `derep`
(from a 25.8 M read-pair library after `derep_pairs`), 91 bp mean length,
`d_max_5 ≈ 0.193`, `lambda_5 ≈ 0.25`:

| Step | Unique clusters | Change |
|------|----------------|--------|
| `derep_pairs` output (merged reads in) | 5,582,073 | - |
| `derep` standard (exact hash) | 3,531,821 | - |
| `derep --collapse-damage` | 3,511,607 | −20,214 (−0.6%) |
| `derep --collapse-damage --error-correct` | 3,510,151 | −21,670 (−0.6%) |

Damage-aware mode merged 20,214 clusters split by terminal deamination. Error
correction (with damage substitution protection) absorbed a further 1,456
PCR-error clusters. A↔T (Channel H) and C↔G (Channel G) transversions are eligible for absorption by default; pass `--protect-transversions` to protect them too.
C↔T and G↔A at **terminal positions** (within 5 bp of the read ends) are handled
by the damage-aware bypass: they are preferentially absorbed as residual
deamination rather than protected as SNPs (see below). G↔T and C↔A (8-oxoG
oxidative damage) are always protected.

---

## Library-type detection (DS vs SS)

Ancient DNA libraries prepared with different protocols show distinct damage
patterns:

- **Double-stranded (DS) libraries** carry symmetric C→T signal at the 5' end
  and G→A at the 3' end; both signals decay at similar rates.
- **Single-stranded (SS) libraries** carry only C→T at the 5' end (original
  strand) or only G→A at the 3' end (complementary strand), depending on
  ligation orientation. Some SS protocols produce a sharp G→A spike at the last
  3' position from ligation artifacts, without the gradual 3' decay seen in DS.

`fqdup extend` and `fqdup derep` both auto-detect library type through a
7-model BIC competition (powered by
[libtaph](https://github.com/genomewalker/libtaph)) before
fitting the decay parameters. The models jointly evaluate four biological
channels, 5' C→T, 3' G→A decay, 3' G→A spike at position 0, and 3' C→T
(SS-original reads only). The winning model determines whether the damage
profile is treated as DS, SS complementary, SS original, or negligible.

The detected library type is printed in the log:

```
Library type: DS (d_max_5=0.193, d_max_3=0.040)
Library type: SS (complement-only; d_max_3=0.221, spike at pos 0)
```

Override with `--library-type ds|ss` if the auto-detection is incorrect (e.g.
mixed libraries or very low-coverage samples where the BIC test is
underpowered).

### Damage-aware bypass

The detected library type also controls which terminal mismatches are treated as
residual deamination rather than candidate SNPs during Phase 3 error correction:

- **DS libraries**: C→T within 5 bp of the 5′ end and G→A within 5 bp of the
  3′ end are flagged as damage-channel mismatches. These are absorbed more
  aggressively (exempt from the singleton high-quality veto).
- **SS libraries**: C→T within 5 bp of **either** end is flagged. G→A is not a
  primary damage pattern in SS protocols and is treated like any other
  substitution.

Interior C↔T and G↔A mismatches (beyond 5 bp from the ends) are scored by the
empirical posterior-odds model like any other transition. They are absorbed if
the cross-bundle recurrence evidence is low (consistent with PCR error) and
protected if the same signature recurs across multiple independent bundles
(consistent with a real variant).

### Paired-end inputs

When `-1/-2` is given, the classifier runs in native PE mode rather than
SE-on-each-mate. `R1[i]` measures top-strand 5'-end position `i`; the
complement of `R2[i]` measures the same fragment's 3'-end position `i`. Pairs
whose insert is shorter than the read length read into the adapter — those
adapter bases would otherwise pollute ct5 / ga3 and (in v8) misclassify
clearly-DS libraries as `SS_orig`. v9 detects short-insert pairs by R1/R2
overlap (15 bp window, ≤3 mismatches) and skips them; the count is reported
as `pe_short_insert_skipped`.

### Chemistry-aware Briggs fit

The exponential model `d(p) = bg + d_max·exp(-p/λ)` uses a tail-anchored
background estimated from positions 20..49 (denom ≥ 100), where the damage
component has decayed by 1–2 orders of magnitude. A global-mean baseline
leaks deamination signal into `bg` and biases `d_max` downward in damaged
libraries. libtaph also reports a non-parametric **area-excess** statistic
(sum of `rate − bg` over the first 10 positions) and a likelihood-ratio
score; the BIC classifier consumes both, so libraries that fit the
exponential form poorly but show clear terminal excess are still classified
correctly.

---

## Damage estimation in `fqdup extend`

`fqdup extend` runs the same Pass 0 damage estimation before building its k-mer
graph. The estimated mask lengths determine which terminal positions are excluded
from k-mer indexing, damaged terminal k-mers would create spurious branches in
the de Bruijn graph because some reads carry C→T and others do not. Excluding
them from graph construction lets the walk traverse damaged positions via k-mers
contributed by overlapping longer reads, which carry the same terminal sequence
without damage.

The `--damage-sample` parameter (default 500 000 reads) controls the estimation
sample size for `fqdup extend`. For most libraries, 500k reads gives stable
estimates; for very low-complexity inputs the full file can be used with
`--damage-sample 0`.

---

## Choosing `--mask-threshold`

The default of 0.05 catches the bulk of deamination events without masking
so far into the read interior that unrelated molecules start colliding. For
heavily damaged libraries (cave sediments, permafrost), raising the threshold
to 0.08–0.10 reduces over-masking. Cross-check the logged masked position count
against a damage plot from DART or mapDamage2 to confirm the threshold
captures the actual damage zone.

---

## Phase B3: damage-aware H>2 merge

The terminal masking (Phase 1) and H≤2 pigeonhole correction (Phase 3) together handle the majority of duplicate fragmentation caused by deamination. However, for libraries with very high terminal damage (d_max ≥ 0.25), PCR copies of the same ancient molecule can accumulate 3–5 deamination events at positions just outside the mask zone (where P(damage) is 1–5%). These reads are too distant for Phase 3's H≤2 pigeonhole to ever pair them.

Phase B3 handles this case with a dedicated damage-normalised bucket hash.

### How it works

For each eligible read, Phase B3 computes two values:

**1. `bundle_key` (locus anchor)**
The same `start_kmer ⊕ end_kmer` used by Phase 3. Two reads sharing a `bundle_key` originated from the same genomic locus (same fragment endpoints). This is the reference-free equivalent of mapping coordinates.

**2. `kdamage_hash` (damage-normalized interior)**
The interior sequence (between the damage mask zones) is extracted into a packed 2-bit buffer. Positions where `P(damage | profile, L) > b3_deam_threshold` (default 0.01) are normalized in-place: T→C at 5'-proximal positions, A→G at 3'-proximal positions. The result is hashed with XXH3. PCR copies of the same molecule that differ only in deamination state at these positions yield identical `kdamage_hash` values.

**3. Composite bucket key `b3k`**
The two values are combined into a single 64-bit key:
```
b3k = XXH3_64bits_withSeed(&kdamage_hash, 8, bundle_key)
```
Two reads land in the same B3 bucket if and only if they share both their locus anchor **and** their damage-normalized interior. This prevents cross-locus collisions — a critical safeguard for high-depth capture data where unrelated molecules at different loci can coincidentally share similar interior sequences.

### Within-bucket filtering

For each pair within a bucket (sorted by descending count, parent-child monotone ascent):

1. **Count ratio**: `parent_count / child_count ≥ b3_count_ratio` (default 5×)
2. **Hamming distance**: H ∈ [3, `b3_max_hamming`] (default max 5)
3. **Deamination consistency**: every mismatch must be C↔T at 5'-proximal or G↔A at 3'-proximal positions with `P(damage) > b3_deam_threshold`
4. **LRT**: the count-ratio likelihood-ratio test must favour PCR error over independent origin (`LRT > log(10)`)

All four conditions must pass for absorption.

### Parameters

| Flag | Default | Effect |
|------|---------|--------|
| `--b3-disable` | — | Disable Phase B3 entirely |
| `--b3-min-dmax` | 0.25 | Skip B3 if d_max_avg below this |
| `--b3-deam-threshold` | 0.01 | P(damage) threshold for normalization and consistency check |
| `--b3-max-hamming` | 5 | Maximum H for B3 absorption |
| `--b3-count-ratio` | 5.0 | Minimum parent/child count ratio |

### Shotgun vs capture

For shotgun data, Phase 3 and terminal masking already handle the bulk of damage-driven fragmentation. B3 fires only when `d_max_avg ≥ 0.25` (2+ Myr libraries), and the bundle_key locus gate has negligible precision cost because cross-locus interior collisions are rare over a large genome.

For capture data, the locus gate is essential. Many independent molecules land on the same narrow target region, so `kdamage_hash` alone would group unrelated molecules into the same bucket. The `bundle_key` gate restricts comparisons to reads sharing identical fragment endpoints — the same constraint that keeps Phase 3 precise.

### Log output

```
Phase B3: running (d_max_avg=0.411)
Phase B3: candidates=12483 absorbed=183 protected=41
```

`candidates` = pairs passing the count ratio and Hamming filter before the deamination-consistency and LRT checks. `absorbed` = pairs fully merged. `protected` = pairs that passed Hamming + count ratio but failed deamination-consistency or LRT (kept separate).
