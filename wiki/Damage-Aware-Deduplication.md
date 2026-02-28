# Damage-aware deduplication

This page covers the first of two sources of false unique clusters in ancient
DNA deduplication: post-mortem deamination. PCR copying errors are covered in
[[PCR-Error-Correction]]. Both are handled by `fqdup derep`; see [[Home]] for
how they fit together.

---

## The problem

Ancient DNA carries characteristic substitution patterns at the ends of
sequenced fragments. Cytosine residues in single-stranded overhangs deaminate
to uracil over time, producing C→T substitutions at the 5' terminus. The same
process on the complementary strand appears as G→A at the 3' terminus. Both
signals decay exponentially from the ends toward the read interior.

When two reads come from the same original molecule but one carries a deaminated
C→T at position 1, exact-match deduplication treats them as different fragments.
The unique molecule count is inflated by however many molecules were sequenced
both with and without terminal damage — which in a heavily-amplified library
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

Only the damage excess — `d_max × exp(-lambda × pos)` — drives masking. The
background `bg` is the same for both copies of a duplicate pair and does not
cause sequence differences between them.

---

## Automatic estimation (`--damage-auto`)

Pass 0 scans all reads before deduplication begins. For each read, `fqdup`
tallies the T/(T+C) frequency at each of the first 15 positions from the 5' end,
and A/(A+G) from the 3' end. The background rate is taken from the middle third
of the read — avoiding terminal damage and adapter composition bias, consistent
with DART's approach.

After scanning, a coverage-weighted ordinary least-squares fit in log-space
(positions 1–9; position 0 is excluded because sorted files tend to put the
most abundant reads first, which can skew terminal composition) gives the
exponential parameters. Lambda is clamped to [0.05, 0.5]. The estimated
parameters are logged:

```
Pass 0: damage estimation — 5582073 reads processed (5582073 total)
  5'-end: d_max=0.099 lambda=0.246 bg=0.487
  3'-end: d_max=0.020 lambda=0.069 bg=0.509
```

If both `d_max_5` and `d_max_3` are below 0.02, damage is treated as
negligible and no masking is applied — standard exact hashing is used.

---

## Masking

After fitting, `fqdup` records which positions actually exceeded the
`--mask-threshold` (default 0.05) in the observed frequencies, not just which
positions the fitted curve predicts. This empirical approach avoids over-masking
when OLS overestimates damage at lightly-affected positions. In practice, 1–3
terminal positions per end are masked.

At hash time, each masked position is replaced with a neutral byte (`\x01` at
C or T in the 5' zone, `\x02` at G or A in the 3' zone). Two reads from the
same molecule that differ only within the masked zone then produce identical
hashed sequences and collapse into the same cluster.

The mask uses the same position indices from both ends — position `i` from the
5' end and position `i` from the 3' end share the same mask flag. This symmetry
is required to preserve the canonical hash invariant:

```
canonical_hash(seq) == canonical_hash(revcomp(seq))
```

Without it, a forward-strand read and its reverse complement would hash
differently after masking — paradoxically increasing unique cluster counts
instead of reducing them.

The masked positions and expected mismatch count are reported after Pass 0:

```
--- Damage-Aware Deduplication ---
  5'-end d_max:  0.098842
  3'-end d_max:  0.019759
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
  -i nonext.deduped.fq.gz -o nonext.final.fq.gz \
  --damage-dmax5 0.35 --damage-lambda5 0.40 \
  --damage-dmax3 0.20 --damage-lambda3 0.30 \
  --mask-threshold 0.05
```

Use `--damage-dmax` to set the same value for both ends (e.g. non-UDG libraries
with symmetric damage).

---

## Effect on cluster counts

Benchmarked on sample `a88af16f35` — 25.8 M read pairs, 91 bp mean length,
`d_max_5 ≈ 0.099`, `lambda_5 ≈ 0.25`:

| Step | Unique clusters | Change |
|------|----------------|--------|
| `derep_pairs` (25.8 M pairs in) | 5,582,073 | — |
| `derep` standard (exact hash) | 3,531,821 | — |
| `derep --damage-auto` | 3,511,607 | −20,214 (−0.6%) |
| `derep --damage-auto --error-correct` | 3,506,272 | −25,549 (−0.7%) |

Damage-aware mode merged 20,214 clusters split by terminal deamination. Error
correction absorbed a further 5,335 low-count PCR-error clusters.

---

## Choosing `--mask-threshold`

The default of 0.05 catches the bulk of deamination events without masking
so far into the read interior that unrelated molecules start colliding. For
heavily damaged libraries (cave sediments, permafrost), raising the threshold
to 0.08–0.10 reduces over-masking. Cross-check the logged masked position count
against a damage plot from DART or mapDamage2 to confirm the threshold
captures the actual damage zone.
