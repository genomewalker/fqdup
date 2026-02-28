# Performance

## Benchmarks

All benchmarks were run on sample `a88af16f35` — 25.8 M paired-end reads,
mean read length ~91 bp, ancient DNA with `d_max_5 ≈ 0.099`, `lambda_5 ≈ 0.25`.
Storage: NFS-mounted.

### Full pipeline: deduplication quality

The full pipeline runs three steps. Numbers shown are from the non-extended file
at each stage.

| Step | Unique clusters | Wall time |
|------|----------------|-----------|
| `derep_pairs` (25.8 M pairs in) | 5,582,073 | ~25 s |
| `derep` standard (exact hash) | 3,531,821 | ~22 s |
| `derep --damage-auto` | 3,511,607 | ~31 s |
| `derep --damage-auto --error-correct` | 3,506,272 | ~33 s |

- `derep_pairs` reduces 25.8 M pairs to 5.58 M unique (78.4% structural dedup)
- `derep --damage-auto` merges a further 20,214 clusters split by terminal
  deamination (−0.6% vs standard)
- `--error-correct` absorbs 5,335 PCR-error clusters (−0.15% vs damage-auto only)

Pass 0 (full-scan damage estimation) adds ~9 seconds. Phase 3 (error
correction) adds ~2 seconds — it runs in memory on the Pass 1 index with no
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

Same structure — one entry per unique cluster.

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

### `fqdup derep` (with `--error-correct`)

Error correction adds the `SeqArena` — a contiguous byte array storing the
full sequence of every unique cluster:

```
Memory ≈ 40 bytes × N_unique    (Pass 1 index)
       + L_avg bytes × N_unique  (SeqArena)
```

| Unique clusters | Mean length | SeqArena |
|----------------|------------|---------|
| 3.5 M | 91 bp | ~320 MB |
| 25 M | 65 bp | ~1.6 GB |
| 100 M | 65 bp | ~6.5 GB |

For a 25 M-pair library with 50% duplication (12.5 M unique) and
`--error-correct` at 65 bp mean length:

```
Index:    40 × 12.5 M = 500 MB
SeqArena: 65 × 12.5 M = 813 MB
Total:   ~1.3 GB
```

---

## Throughput and I/O

Deduplication is I/O-bound on network storage. On the NFS-mounted server:

| Step | Observed throughput |
|------|-------------------|
| Pass 0 (damage estimation) | ~600 k reads/s |
| Pass 1 (index build) | ~400–500 k reads/s |
| Pass 2 (output) | ~200–300 k reads/s |

On local NVMe, throughput is higher and CPU becomes the bottleneck.

### Decompression acceleration

| Flag | Speedup | Requirement |
|------|---------|-------------|
| `--isal` | 3–5× decompression | Intel ISA-L library |
| `--pigz` | 2–3× decompression | pigz in PATH |

Use `--isal` when reading `.gz` files from fast storage. Both flags affect
decompression only; gzip output uses pigz automatically when available.

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
   `derep` on the smaller non-extended output. The non-extended file after
   `derep_pairs` is typically 60–80% smaller.

2. **I/O**: use `--isal` or `--pigz` and `--fast` for the sort step. Place
   temp files on fast local storage with `-t /local/scratch`.

3. **Error correction**: `--error-correct` is optional. For very large
   datasets the SeqArena can be substantial; the quality gain (typically
   < 0.2% additional merges) may not justify the memory cost.
