# Performance

## Benchmarks

All benchmarks were run on a 96-core server with NFS-mounted storage.

### Deduplication quality

Tested on sample `a88af16f35` — a 25.8 M-read ancient DNA library with
mean read length ~65 bp and estimated `d_max_5 ≈ 0.071`, `lambda_5 ≈ 0.29`:

| Pipeline | Unique clusters | Change vs standard |
|----------|-----------------|--------------------|
| `derep_pairs` only (standard) | 5,582,073 | — |
| `derep --damage-auto` | 5,547,508 | −34,565 (−0.6%) |
| `derep --damage-auto --error-correct` | 5,532,327 | −49,746 (−0.9%) |

The damage-aware step correctly merged 34,565 read pairs over-split by
terminal deamination. Phase 3 error correction absorbed a further 15,181
PCR-error clusters.

### Wall time

| Step | Time |
|------|------|
| `sort` (25.8 M reads, gzip input) | ~20 s |
| `derep_pairs` (25.8 M pairs, gzip input) | ~25 s |
| `derep` (standard, gzip input) | ~30 s |
| `derep --damage-auto` | ~31 s |
| `derep --damage-auto --error-correct` | ~33 s |

Phase 3 error correction adds only ~2–3 seconds because it operates entirely
in memory on the index built during Pass 1 — no additional I/O.

---

## Memory

### `fqdup derep_pairs`

The Pass 1 index stores one `IndexEntry` per **input read pair** — not per
unique cluster — because we need to track which record index is the
representative for each cluster. `IndexEntry` is approximately 24 bytes.

```
Memory ≈ 24 bytes × N_input_pairs
```

| Input size | Memory |
|------------|--------|
| 25 M pairs | ~600 MB |
| 100 M pairs | ~2.4 GB |
| 400 M pairs | ~9.6 GB |

### `fqdup derep` (without `--error-correct`)

The Pass 1 index stores one `IndexEntry` (~16 bytes) per **input read**:

```
Memory ≈ 16 bytes × N_input_reads
```

| Input size | Memory |
|------------|--------|
| 25 M reads | ~400 MB |
| 100 M reads | ~1.6 GB |
| 400 M reads | ~6.4 GB |

### `fqdup derep` (with `--error-correct`)

Error correction adds the `SeqArena` — a contiguous byte array storing the
full sequence of every unique cluster (one sequence per cluster, not per read):

```
Memory ≈ 16 bytes × N_reads           (Pass 1 index)
       + L_avg bytes × N_unique        (SeqArena)
```

| Unique clusters | Mean length | SeqArena |
|----------------|------------|---------|
| 5 M | 65 bp | ~325 MB |
| 50 M | 65 bp | ~3.25 GB |
| 200 M | 65 bp | ~13 GB |

For a 400 M-read library with 50% duplication rate and `--error-correct`:

```
Pass 1 index:  16 × 400 M = 6.4 GB
SeqArena:      65 × 200 M = 13 GB
Total:        ~19–20 GB
```

Without `--error-correct` the SeqArena is never allocated: 400 M reads at
16 bytes/read = 6.4 GB.

---

## Throughput and I/O

Deduplication is primarily I/O-bound, not CPU-bound. Reads are processed at
approximately **40,000–50,000 reads/s** on NFS. On local NVMe this rises
significantly and CPU becomes the bottleneck.

### Decompression acceleration

| Flag | Speedup | Requirement |
|------|---------|-------------|
| `--isal` | 4–6× decompression | Intel ISA-L library |
| `--pigz` | 2–3× decompression | pigz in PATH |

Use `--isal` when reading `.gz` files from fast storage. Both flags affect
decompression only; gzip **output** uses pigz automatically when available.

### Sort acceleration

```bash
# Faster sort: uncompressed temporaries (3× faster I/O, more disk)
fqdup sort -i reads.fq.gz -o reads.sorted.fq.gz --max-memory 64G --fast

# Control parallelism
fqdup sort ... -p 16   # use 16 threads for sorting
```

The sort is the most parallelism-friendly step — it scales linearly with
available cores for the in-memory chunk sort phase.

### Two-pass overhead

Both `derep_pairs` and `derep` read their input **twice** (Pass 1 + Pass 2).
Plan for 2× the read I/O of the input file size. Output is always written
once (Pass 2 only).

---

## Scaling to Larger Datasets

For very large libraries (> 400 M read pairs, > 100 M unique clusters):

1. **Memory**: use `derep_pairs` first to reduce unique clusters, then run
   `derep` on the smaller non-extended output. The non-extended file after
   `derep_pairs` will typically be substantially smaller than the input.

2. **I/O**: use `--isal` or `--pigz` and `--fast` (for sort). If possible,
   place temp files on fast local storage (`-t /local/scratch`).

3. **Error correction**: if memory is tight, `--error-correct` is optional.
   The quality gain (typically < 1% additional merges) may not justify the
   SeqArena overhead for extremely large datasets.
