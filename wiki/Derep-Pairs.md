# fqdup derep_pairs

## Purpose

`fqdup derep_pairs` is the structural deduplication step for paired reads.
It deduplicates using the `fqdup extend`-assembled sequence as the cluster
fingerprint rather than the raw merged read.

Using the extended sequence as fingerprint reduces false collisions: two
different molecules that happen to share a short merged-read sequence will
typically diverge in their assembly-extended versions. Conversely, PCR
duplicates that differ only because one was captured in a slightly shifted
window extend to the same flanking sequence and collapse correctly.

The representative kept per cluster is the pair with the **longest merged
read**, preserving the most original ancient DNA sequence for downstream
analysis.

`fqdup derep_pairs` does not perform damage-aware hashing or PCR error
correction. Those are handled by the subsequent `fqdup derep` step.

---

## Algorithm

`fqdup derep_pairs` is a two-pass streaming algorithm.

**Pass 1: index construction**
Both sorted files are streamed in lockstep. For each pair, the canonical hash
of the **extended** read is computed: `min(XXH3_128(ext), XXH3_128(revcomp(ext)))`.
The index stores `{record_offset, count, merged_len}`, ~40 bytes per unique
pair. When the same extended hash appears again, the count is incremented; if
the new merged read is longer than the stored representative, the representative
is updated.

**Pass 2: output**
Both files are re-streamed. For each pair, if this is the representative
(longest merged read of its cluster), it is written to both output files.

---

## CLI reference

```
fqdup derep_pairs -n NON -e EXT -o-non OUT -o-ext OUT [options]

Required:
  -n FILE       Sorted merged (fastp) FASTQ
  -e FILE       Sorted fqdup-extended FASTQ
  -o-non FILE   Output merged FASTQ (representatives)
  -o-ext FILE   Output extended FASTQ (representatives)

Optional:
  -c FILE        Cluster statistics (gzipped TSV)
  --no-revcomp   Disable reverse-complement collapsing (default: enabled)
```

Both input files must be sorted by read ID (`fqdup sort`). The two files must
contain reads in the same order, the same read ID at the same position in
both files.

---

## Output

### Reads

Two FASTQ files:
- `-o-non`: merged reads (representatives), use this for downstream analysis
- `-o-ext`: extended reads (representatives), only needed if further
  processing requires the extended fingerprints

### Cluster statistics (`-c FILE.tsv.gz`)

Gzipped TSV, one row per cluster:

| Column | Description |
|--------|-------------|
| `hash` | 128-bit XXH3 canonical hash of extended read (32 hex chars) |
| `ext_len` | Extended read length of the representative |
| `pair_count` | Number of read pairs in the cluster |
| `non_len` | Merged read length of the representative |

---

## Memory

The index stores ~40 bytes per unique pair (record offset + count + merged
length). For a 25 M-pair library with ~5.6 M unique pairs: ~220 MB. Memory
scales linearly with unique pair count, not total reads.

---

## Benchmarks

On sample `a88af16f35`, 25.8 M read pairs in, ~91 bp merged-read mean
length:

| Metric | Value |
|--------|-------|
| Input pairs | 25,800,000 |
| Unique clusters | 5,582,073 |
| Wall time | ~25 s |
| Index memory | ~220 MB |

See [[Performance]] for full benchmarks.
