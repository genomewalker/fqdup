# fqdup

FASTQ deduplication for paired-end ancient DNA libraries.

---

## The problem

Standard deduplication identifies duplicate reads by exact sequence match. Two
reads are duplicates if and only if they are identical. This works for modern
DNA, but ancient DNA has two properties that break this assumption.

**Post-mortem deamination.** Cytosine residues in ancient DNA spontaneously
deaminate to uracil over time, producing C→T substitutions at the 5' end and
G→A substitutions at the 3' end (the complement on the reverse strand). The
rate is highest at the terminal positions and decays exponentially toward the
interior. Two reads from the same original molecule may therefore differ at
position 1 — one carrying the original C, the other carrying a deaminated T —
and an exact-match deduplicator counts them as distinct. The true unique molecule
count is inflated.

**PCR errors.** Library amplification introduces substitution errors at a low
but nonzero rate (~10⁻⁷ per base per doubling for high-fidelity polymerases).
A molecule amplified 30 cycles has thousands of copies; a single early copying
error propagates into a cluster of slightly-wrong copies. These clusters tend
to have low read counts relative to the correct sequence, but they survive
exact-match deduplication as apparent unique molecules.

Both effects push in the same direction: more apparent unique sequences than
there are true unique molecules. For ancient DNA studies that count unique
molecules to estimate coverage or population size, this matters.

---

## How fqdup addresses this

`fqdup` adds two steps on top of standard deduplication, each targeting one
source of false positives.

**Damage-aware hashing.** Before hashing, positions where the empirically
observed deamination rate exceeds a threshold are replaced with a neutral byte.
Two reads that differ only at a masked terminal position then produce the same
hash and collapse into the same cluster. The mask is derived from the data
itself (Pass 0) using the same exponential decay model as DART and mapDamage2.
See [[Damage-Aware-Deduplication]].

**PCR error correction.** After the deduplication index is built (Pass 1),
Phase 3 identifies clusters with low read count that differ from a high-count
cluster by exactly one substitution in the interior (outside the damage zone).
These are almost certainly PCR copying errors. The count ratio threshold
(default 50×) ensures only unambiguous cases are absorbed. See
[[PCR-Error-Correction]].

The two steps are complementary: damage masking handles terminal variation from
deamination; error correction handles interior variation from PCR. Neither
alone is sufficient.

---

## Pipeline overview

```
fqdup sort        fqdup derep_pairs        fqdup derep
──────────────    ──────────────────────   ──────────────────────────────────────
Sort reads by  →  One representative    →  Damage-aware hashing + PCR error
read ID           pair per cluster         correction (single-file)
(prerequisite)    (structural dedup)        (biological dedup)
```

`derep_pairs` handles paired-file complexity: two files read in lockstep,
one representative pair per cluster, chosen by longest non-extended read.
No biology here — just hashing and representative selection.

`derep` handles the biological layer: optional damage masking, optional PCR
error correction, single-file input and output.

The split matters because the paired structure (extended vs non-extended reads)
is a different concern from the damage biology. `derep_pairs` can be used
without `derep` (if damage correction is not needed), and `derep` can be used
on single-end data without `derep_pairs`.

---

## Wiki pages

| Page | Contents |
|------|----------|
| [[Installation]] | Build from source, dependencies |
| [[Usage]] | Tutorials and options for all three subcommands |
| [[Algorithm]] | Internal mechanics of sort, derep_pairs, and derep |
| [[Damage-Aware-Deduplication]] | The deamination model, empirical masking, and symmetry |
| [[PCR-Error-Correction]] | Phase 3: 3-way pigeonhole Hamming search |
| [[Performance]] | Benchmarks and memory usage |

---

## Quick reference

### Full ancient DNA pipeline

```bash
# Sort
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G --fast
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G --fast

# Structural dedup (paired)
fqdup derep_pairs \
  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz

# Biological dedup (damage + error correction)
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-auto --error-correct
```

### Standard (no ancient DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 64G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz
```
