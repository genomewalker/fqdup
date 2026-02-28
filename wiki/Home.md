# fqdup Wiki

Ultra-fast, memory-efficient FASTQ deduplication for paired-end ancient DNA libraries.

## Pipeline Overview

```
fqdup sort        fqdup derep_pairs        fqdup derep
──────────────    ──────────────────────   ──────────────────────────────────────
Sort reads by  →  One representative    →  Damage-aware hashing + PCR error
read ID           pair per cluster         correction (single-file)
(prerequisite)    (structural dedup)        (biological dedup)
```

The design separates two orthogonal concerns:

- **Structural deduplication** (`derep_pairs`): two paired files, read in
  lockstep, one representative pair chosen per cluster (longest non-extended
  mate wins). No biology — just hashing and representative selection.

- **Biological deduplication** (`derep`): single file, optional ancient DNA
  damage masking, optional PCR error correction. Does not touch the paired
  structure at all.

Either step can be run independently. If you have single-end data, skip
`derep_pairs` and run `derep` directly after `sort`.

## Wiki Pages

| Page | Contents |
|------|----------|
| [[Installation]] | Build from source, dependencies, Conda |
| [[Usage]] | Tutorials and common invocations for all three steps |
| [[Algorithm]] | How sort, derep_pairs, and derep work internally |
| [[Damage-Aware-Deduplication]] | The C→T / G→A damage model and masking strategy |
| [[PCR-Error-Correction]] | Phase 3 count-stratified 3-way pigeonhole algorithm |
| [[Performance]] | Benchmarks, memory usage, I/O tips |

## Quick Reference

### Full ancient DNA pipeline

```bash
# Sort
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G --fast
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G --fast

# Paired deduplication (structural)
fqdup derep_pairs \
  -n nonext.sorted.fq.gz -e ext.sorted.fq.gz \
  -o-non nonext.deduped.fq.gz -o-ext ext.deduped.fq.gz

# Single-file dedup (biological — damage + error correction)
fqdup derep \
  -i nonext.deduped.fq.gz \
  -o nonext.final.fq.gz \
  --damage-auto --error-correct
```

### Standard (no damage, no error correction)

```bash
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 64G
fqdup derep -i reads.sorted.fq.gz -o reads.deduped.fq.gz
```
