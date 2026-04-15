# fqdup

FASTQ deduplication for ancient DNA libraries.

---

## Background

Standard deduplication identifies duplicate reads by exact sequence match. Two
reads are duplicates if and only if they are identical. This works for modern
DNA, but ancient DNA has two properties that break this assumption.

**Post-mortem deamination.** Cytosine residues in ancient DNA spontaneously
deaminate to uracil over time, producing C→T substitutions at the 5' end and
G→A substitutions at the 3' end (the complement on the reverse strand). The
rate is highest at the terminal positions and decays exponentially toward the
interior. Two reads from different original ancient molecules at the same genomic locus
may therefore differ at position 1 — one molecule was deaminated at that site,
the other was not — and an exact-match deduplicator counts them as distinct.
The true unique molecule count is inflated.

**PCR errors.** Library amplification introduces substitution errors at a low
but nonzero rate (~10⁻⁷ per base per doubling for high-fidelity polymerases).
A molecule amplified 30 cycles has ~10⁹ copies (2³⁰); a single early copying
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

The upstream steps that produce fqdup's inputs are:

1. **fastp merge**, collapse overlapping R1+R2 into a single sequence representing the full ancient DNA molecule
2. **`fqdup extend`**, extend each merged read outward from both ends using the built-in de Bruijn graph assembler, recovering sequence that may be present in other reads in the library

This produces two files for the same set of molecules:
- **merged** (`-n`): the original fastp-merged read, actual ancient DNA sequence
- **extended** (`-e`): the same read after `fqdup extend` assembly, longer, assembly-assisted fingerprint

fqdup then runs in three steps:

```
fqdup extend      fqdup sort        fqdup derep_pairs          fqdup derep
──────────────    ──────────────    ────────────────────────   ─────────────────────
Extend merged  →  Sort merged   →   One representative      →  Damage-aware hashing
reads via         and extended      pair per cluster           + PCR error correction
de Bruijn         reads             (uses extended as           on merged reps
graph                               fingerprint; keeps          (biological dedup)
                                    merged as output)
```

`derep_pairs` deduplicates using the `fqdup extend`-assembled sequence as the
cluster fingerprint, it is longer and more unique than the merged read alone,
reducing false collisions between reads from different molecules. The
representative kept is the pair with the longest merged read (most original
ancient DNA sequence).

`derep` then handles the biological layer on the merged output: damage masking
collapses reads that differ only due to deamination; Phase 3 absorbs low-count
clusters that differ from a high-count cluster by a single non-damage substitution.

Any residual duplicates `derep` sees are merged reads that were identical but whose
extensions diverged slightly, `derep_pairs` kept both; `derep` collapses
them on the original merged sequence.

---

## Wiki pages

| Page | Contents |
|------|----------|
| [[Installation]] | Build from source, dependencies |
| [[Usage]] | Tutorials and options for all subcommands |
| [[Extend]] | fqdup extend: de Bruijn graph extension algorithm, CLI, benchmarks |
| [[Damage]] | fqdup damage: standalone damage profiler, inspect d_max, library type, mask positions |
| [[Algorithm]] | Internal mechanics of sort, derep_pairs, and derep |
| [[Damage-Aware-Deduplication]] | The deamination model, empirical masking, and symmetry |
| [[PCR-Error-Correction]] | Phase 3: 3-way pigeonhole Hamming search |
| [[Performance]] | Benchmarks and memory usage |

---

## Quick reference

### Full ancient DNA pipeline

```bash
# Extend merged reads using the built-in de Bruijn graph assembler
fqdup extend -i merged.fq.gz -o extended.fq.gz

# Sort merged and extended reads
fqdup sort -i merged.fq.gz   -o merged.sorted.fq.gz   --max-memory 64G --fast
fqdup sort -i extended.fq.gz -o extended.sorted.fq.gz --max-memory 64G --fast

# Structural dedup: one representative pair per cluster (extended = fingerprint)
fqdup derep_pairs \
  -n merged.sorted.fq.gz \
  -e extended.sorted.fq.gz \
  -o-non merged.deduped.fq.gz \
  -o-ext  extended.deduped.fq.gz

# Biological dedup: PCR error correction (on by default)
# Add --collapse-damage only if NOT running DART/mapDamage on this output
fqdup derep \
  -i merged.deduped.fq.gz \
  -o merged.final.fq.gz
```

### Standard (no ancient DNA)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 64G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz --no-error-correct
```
