# Performance

## Benchmarks

All benchmarks were run on sample `a88af16f35`, 25.8 M fastp-merged reads
(from paired-end sequencing), mean read length ~91 bp, ancient DNA with
`d_max_5 ≈ 0.193`, `lambda_5 ≈ 0.25`. Storage: NFS-mounted.

### Full pipeline: deduplication quality

The full pipeline runs three steps. Numbers shown are for the merged read file
at each stage (the `-n` / output file of `derep_pairs`).

| Step | Unique clusters | Wall time |
|------|----------------|-----------|
| `derep_pairs` (25.8 M merged reads in) | 5,582,073 | ~25 s |
| `derep` standard (exact hash) | 3,531,821 | ~22 s |
| `derep --damage-auto` | 3,511,607 | ~31 s |
| `derep --damage-auto --error-correct` | 3,510,151 | ~33 s |

- `derep_pairs` reduces 25.8 M merged reads to 5.58 M unique (78.4% structural dedup)
- `derep --damage-auto` merges a further 20,214 clusters split by terminal
  deamination (−0.6% vs standard)
- `--error-correct` absorbs 1,456 PCR-error clusters (−0.04% vs damage-auto only;
  C↔T and G↔A mismatches are protected as potential damage signal)

Pass 0 (full-scan damage estimation) adds ~9 seconds. Phase 3 (error
correction) adds ~2 seconds, it runs in memory on the Pass 1 index with no
additional I/O.

---

## Memory

### `fqdup derep_pairs`

The index stores one entry per unique cluster. With hash map overhead, expect
approximately 40 bytes per unique pair.

```
Memory ≈ 40 bytes × N_unique_pairs
```

| Unique pairs | Memory |
|-------------|--------|
| 5 M | ~200 MB |
| 25 M | ~1 GB |
| 100 M | ~4 GB |

### `fqdup derep` (without `--error-correct`)

Same structure, one entry per unique cluster.

```
Memory ≈ 40 bytes × N_unique_clusters
```

| Unique clusters | Memory |
|----------------|--------|
| 3.5 M | ~140 MB |
| 25 M | ~1 GB |
| 100 M | ~4 GB |
| 200 M | ~8 GB |

The index size from our benchmark: 134 MB for 3.53 M unique clusters from
5.58 M reads (38 bytes/unique cluster, consistent with the formula).

### `fqdup derep` (with error correction, default on)

Error correction adds the `SeqArena`, a 2-bit packed array storing each base
in 2 bits rather than 8, approximately 4× more compact than ASCII:

```
Memory ≈ 40 bytes × N_unique              (Pass 1 index)
       + 0.25 × L_avg bytes × N_unique    (SeqArena, 2-bit packed)
```

| Unique clusters | Mean length | SeqArena |
|----------------|------------|---------|
| 3.5 M | 91 bp | ~80 MB |
| 25 M | 65 bp | ~410 MB |
| 100 M | 65 bp | ~1.6 GB |
| 260 M | 65 bp | ~4.2 GB |

For a 25 M-pair library with 50% duplication (12.5 M unique) and default error
correction at 65 bp mean length:

```
Index:    40 × 12.5 M        = 500 MB
SeqArena: 0.25 × 65 × 12.5 M = ~203 MB
Total:   ~700 MB
```

Use `--no-error-correct` to skip Phase 3 entirely and avoid allocating the
SeqArena. For most aDNA libraries the SeqArena is a small fraction of total
memory, the index dominates.

---

## fqdup extend

Benchmarked on DS4, 198 M reads, 16 threads, dandycomp02fl:

| Metric | Value |
|--------|-------|
| Reads extended | 50% |
| Average extension | 3.82 / 3.83 bp per side |
| Pass 1 (build) | 68 s |
| Finalize | 19 s |
| Pass 2 (extend + write) | 5:27 |
| Total wall time | 6:56 |
| Peak RSS | 41.3 GB |

Memory breakdown: 6.7 B k-mer observations × 8 bytes during accumulation; after
finalize 1.644 B distinct k-mers × 16 bytes = 26 GB; peak RSS 41 GB (k-mer store
+ read buffers + thread overhead).

Comparison with Tadpole on the same DS4 dataset: Tadpole achieves 14.75%
extension in ~5 min; `fqdup extend` achieves 50% in 7 min.

---

## Throughput and I/O

Deduplication is I/O-bound on network storage. On the NFS-mounted server:

| Step | Observed throughput |
|------|-------------------|
| Pass 0 (damage estimation) | ~600 k reads/s |
| Pass 1 (index build) | ~400–500 k reads/s |
| Pass 2 (output) | ~200–300 k reads/s |

On local NVMe, throughput is higher and CPU becomes the bottleneck.

### Decompression and compression

Decompression uses rapidgzip, a built-in parallel gzip decoder. It is always
active and requires no flags. On `.gz` input it uses multiple threads
automatically.

Output compression uses bgzf (htslib) when `--threads > 1`, otherwise standard
zlib gzip. No flags are required.

### Sort acceleration

```bash
# Uncompressed temporaries: ~3× faster I/O, more temporary disk space
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 64G --fast

# Limit parallelism if needed
fqdup sort ... -p 16
```

The sort scales well with available cores for the in-memory chunk sort phase.

### Two-pass overhead

Both `derep_pairs` and `derep` read their input twice (Pass 1 + Pass 2).
Plan for 2× the read I/O of the input file. Output is written once.

---

## Scaling to larger datasets

For libraries above 400 M read pairs or 100 M unique clusters:

1. **Memory**: run `derep_pairs` first to reduce unique count, then run
   `derep` on the smaller merged read output. The merged file after
   `derep_pairs` is typically 60–80% smaller.

2. **I/O**: use `--fast` for the sort step. Place temp files on fast local
   storage with `-t /local/scratch`. Decompression uses rapidgzip automatically.

3. **Error correction**: on by default. For very large datasets with 260 M
   unique clusters the 2-bit SeqArena adds ~4 GB; the quality gain (typically
   < 0.2% additional merges) may not justify this for extremely large
   libraries. Use `--no-error-correct` to skip it.
