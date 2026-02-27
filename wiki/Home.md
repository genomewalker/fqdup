# fqdup Wiki

Ultra-fast, memory-efficient paired-end FASTQ deduplication with ancient DNA damage awareness.

## Pages

- [[Installation]] — build from source, dependencies
- [[Usage]] — quick-start workflow and common invocations
- [[Algorithm]] — detailed description of the three-pass algorithm
- [[Damage-Aware-Deduplication]] — ancient DNA C→T / G→A damage model
- [[PCR-Error-Correction]] — Phase 3 count-stratified error correction
- [[Performance]] — benchmarks, memory, throughput

## Quick Start

```bash
# Sort
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 64G

# Deduplicate (standard)
fqdup derep -n non.sorted.fq.gz -e ext.sorted.fq.gz \
            -o-non out.non.fq.gz -o-ext out.ext.fq.gz

# Deduplicate (ancient DNA — automatic damage estimation + PCR error correction)
fqdup derep -n non.sorted.fq.gz -e ext.sorted.fq.gz \
            -o-non out.non.fq.gz -o-ext out.ext.fq.gz \
            --damage-auto --error-correct
```
