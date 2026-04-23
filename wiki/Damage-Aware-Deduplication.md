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
PCR-error clusters. Only A↔T and C↔G transversions are eligible for absorption;
C↔T and G↔A (deamination) and G↔T and C↔A (8-oxoG oxidative damage) are all
protected as potential ancient DNA damage signal.

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
underpowered). The library type only affects which channels are used for the
damage model fit and masking; the downstream deduplication steps are identical.

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
