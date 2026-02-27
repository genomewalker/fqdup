# Algorithm

## Sort

External merge sort with three phases:

1. **Chunk ingestion** (single-threaded): Stream input, accumulate reads into
   memory-limited chunks.
2. **Parallel sort** (multi-threaded): Each chunk is sorted in-place by read
   ID using `std::sort` across all available cores.
3. **k-way merge** (single-threaded): Sorted chunks are merged with a
   priority queue; output is written as a single sorted stream.

`--fast` keeps chunk intermediates as uncompressed FASTQ (3× faster I/O,
uses more disk). Without `--fast`, chunks are gzip-compressed on disk.

## Derep

### Overview

```
[Pass 0]  (damage-aware only)  Estimate deamination model
    ↓
[Pass 1]  Build hash → IndexEntry map
    ↓
[Phase 3] (error-correct only) Absorb PCR 1-error variants
    ↓
[Pass 2]  Stream again, write representative records
```

### Pass 0 — Damage estimation

Stride-samples every Nth non-extended read (default N=1000, 100 k reads
scanned). For each of the first 15 positions at the 5' and 3' ends, the
fraction of C→T (5') and G→A (3') substitutions relative to the baseline
mid-read rate is measured. An exponential decay model is fitted:

```
P(deamination, pos) = d_max × exp(−λ × pos) + bg
```

Fit: coverage-weighted OLS on positions 1–9 (not 0, which is unreliable in
sorted files). λ is clamped to [0.05, 0.5].

### Pass 1 — Index construction

For each read pair `(non, ext)`:
1. Compute canonical hash of the non-extended sequence:
   ```
   h = min(XXH3_64(seq), XXH3_64(revcomp(seq)))
   ```
   With damage masking, positions in the damage zone are replaced with a
   neutral character before hashing. Masking is **symmetric**:
   `mask(pos) = max(P_CT(pos), P_GA(pos))`, so
   `canonical_hash(seq) == canonical_hash(revcomp(seq))`.
2. Insert `(h, seq_len) → IndexEntry{record_idx, non_len, count}` into the
   hash map. On collision:
   - Increment `count`.
   - Keep the entry with the longest non-extended read as representative.

Memory: ~16 bytes per **input** read (not per unique cluster).

### Phase 3 — PCR error correction

See [[PCR-Error-Correction]] for full details.

### Pass 2 — Output

The index maps each unique `(hash, length)` fingerprint to a record index
(byte offset in the sorted input). The sorted files are streamed a second
time; each record is written if and only if its record index matches the
representative for its fingerprint and it was not flagged as a PCR error in
Phase 3.

## Hashing

[XXH3_64](https://github.com/Cyan4973/xxHash) is used throughout — it is
non-cryptographic, extremely fast (≥ 20 GB/s on modern CPUs), and has
excellent avalanche properties suitable for deduplication.

`ska::flat_hash_map` (open-addressing with robin-hood probing) provides
cache-friendly O(1) average-case lookups.
