# fqdup

FASTQ deduplication for paired-end ancient DNA libraries.

## Pipeline overview

```
fqdup sort        fqdup derep_pairs        fqdup derep
──────────────    ──────────────────────   ──────────────────────────────────────
Sort reads by  →  One representative    →  Damage-aware hashing + PCR error
read ID           pair per cluster         correction (single-file)
(prerequisite)    (structural dedup)        (biological dedup)
```

The design separates two concerns:

- **Structural deduplication** (`derep_pairs`): two paired files, read in
  lockstep, one representative pair per cluster (longest non-extended mate
  wins). No biology — just hashing and representative selection.

- **Biological deduplication** (`derep`): single file, optional ancient DNA
  damage masking, optional PCR error correction.

Either step can be run independently. For single-end data, skip `derep_pairs`
and run `derep` directly after `sort`.

## Wiki pages

| Page | Contents |
|------|----------|
| [[Installation]] | Build from source, dependencies |
| [[Usage]] | Tutorials and options for all three subcommands |
| [[Algorithm]] | How sort, derep_pairs, and derep work internally |
| [[Damage-Aware-Deduplication]] | The C→T / G→A damage model and masking |
| [[PCR-Error-Correction]] | Phase 3 count-stratified 3-way pigeonhole algorithm |
| [[Performance]] | Benchmarks and memory usage |

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

# Biological dedup (single-file, damage + error correction)
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
