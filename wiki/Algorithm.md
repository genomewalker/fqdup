# Algorithm

## fqdup sort

External merge sort in three phases:

```
Phase 1: Chunk ingestion (single-threaded)
  Stream input → accumulate reads into memory-limited chunks
  Each full chunk: sort in-place (parallel, std::sort) → write to disk

Phase 2: Parallel sort (already done per-chunk in Phase 1)

Phase 3: k-way merge (single-threaded)
  Priority queue over all sorted chunk files → single sorted output stream
```

Chunks are gzip-compressed by default. With `--fast`, chunks stay
uncompressed (~3× faster I/O at the cost of more temporary disk space).

Reads are sorted lexicographically by read ID (everything before the first
space in the `@` header line). Natural-numeric sort (`-N`) is also available
for IDs with numeric suffixes.

---

## fqdup derep_pairs

Two-pass paired deduplication on sorted files.

### Overview

```
[Pass 1]  Stream ext+non in lockstep → build hash → IndexEntry map
             ↓
[Pass 2]  Stream ext+non again → write representative pairs
```

### Pass 1 — index construction

Both sorted files are read simultaneously, one pair at a time. Reads at each
position are expected to share the same read ID (a warning is logged if they
differ). For each pair `(ext, non)`:

1. **Hash the extended read:**

   ```
   h = canonical_hash(ext.seq) = min(XXH3_64(ext.seq), XXH3_64(revcomp(ext.seq)))
   ```

   Taking the min over the sequence and its reverse complement ensures that a
   molecule sequenced from either strand lands in the same cluster.

2. **Form a fingerprint:** `(h, len(ext.seq))` — the sequence length is
   included to reduce false collisions between reads of different lengths.

3. **Insert into the index:**
   - New fingerprint: `IndexEntry{ record_idx=current, non_len=len(non.seq), count=1 }`
   - Existing fingerprint:
     - Increment `count`
     - If `len(non.seq) > stored non_len`: update `record_idx` and `non_len`

   **Representative selection:** the pair with the longest non-extended read
   is kept. Longer non-extended reads retain more sequence — they were trimmed
   less aggressively during adapter removal.

Memory: the index has one entry per unique cluster. With hash map overhead,
expect approximately 40 bytes per unique cluster.

### Pass 2 — output

The index maps each fingerprint to the record index of its representative.
A reverse map `record_idx → fingerprint` is built, then both files are
streamed again:

```
for each pair (ext, non) at position record_idx:
    if record_idx in records_to_write:
        write ext to out_ext
        write non to out_non
        (optionally) write cluster stats to -c file
```

No damage masking or error correction happens here.

---

## fqdup derep

Single-file deduplication with optional damage masking and PCR error
correction. Input is a sorted FASTQ (typically the non-extended output
of `derep_pairs`).

### Overview

```
[Pass 0]  (damage-aware only)   Full scan → fit damage model
              ↓
[Pass 1]  Stream input → build hash → IndexEntry map
              ↓
[Phase 3] (error-correct only)  Absorb PCR 1-error variants
              ↓
[Pass 2]  Stream again → write representative records
```

### Pass 0 — damage estimation

See [[Damage-Aware-Deduplication]] for the full model description.

All reads are scanned to measure per-position C→T (5' end) and G→A (3' end)
frequencies. The background rate comes from the middle third of each read. An
exponential decay model is fitted to each end:

```
P_CT(pos) = d_max_5 × exp(-lambda_5 × pos) + bg
P_GA(pos) = d_max_3 × exp(-lambda_3 × pos) + bg
```

The fit uses coverage-weighted OLS on positions 1–9. Lambda is clamped to
[0.05, 0.5]. After fitting, positions where the observed excess exceeds
`--mask-threshold` are stored in a boolean mask array — no model evaluation
happens at hash time.

Damage is treated as significant when `d_max_5 > 0.02` or `d_max_3 > 0.02`.
Below that, Pass 0 still runs but no masking is applied.

### Pass 1 — index construction

For each read `rec`:

1. **Compute hash** (with or without damage masking):

   Without damage:
   ```
   h = min(XXH3_64(rec.seq), XXH3_64(revcomp(rec.seq)))
   ```

   With damage masking:
   ```
   masked = apply_damage_mask(rec.seq)
   h      = min(XXH3_64(masked), XXH3_64(apply_damage_mask(revcomp(rec.seq))))
   ```

   The masking substitutes a neutral byte at each position where the damage
   excess exceeds `--mask-threshold`. Independently-deaminated copies of the
   same molecule then hash identically. See [[Damage-Aware-Deduplication]] for
   why the mask must be symmetric.

2. **Form fingerprint:** `(h, len(rec.seq))`

3. **Insert into index:**
   - New fingerprint: record index and (if `--error-correct`) append sequence
     to SeqArena
   - Existing fingerprint: increment `count`

### Phase 3 — PCR error correction

See [[PCR-Error-Correction]] for full details.

### Pass 2 — output

```
for each read at position record_idx:
    if record_idx in records_to_write (and not flagged as PCR error):
        write rec to output
        (optionally) write cluster stats to -c file
```

---

## Canonical hashing

```
canonical_hash(seq) = min(XXH3_64(seq), XXH3_64(revcomp(seq)))
```

[XXH3_64](https://github.com/Cyan4973/xxHash) is a non-cryptographic hash
with throughput above 20 GB/s and excellent distribution. `ska::flat_hash_map`
(open-addressing, robin-hood probing) provides cache-friendly O(1) lookup and
insert. The composite key `(hash, seq_len)` reduces false-positive collisions
between reads of different lengths.

---

## Memory model

| Component | Formula | Example (5.6 M reads, 63% dup) |
|-----------|---------|-------------------------------|
| `derep_pairs` Pass 1 index | ~40 bytes × N_unique_pairs | ~90 MB (2.2 M unique) |
| `derep` Pass 1 index | ~40 bytes × N_unique_clusters | ~140 MB (3.5 M unique) |
| SeqArena (error-correct) | ~1 byte/bp × N_unique × L_avg | ~225 MB (3.5 M × 65 bp) |

The index grows with unique clusters, not total reads. For large libraries with
low duplication rates the distinction matters less; for highly-amplified
libraries (50–80% duplication) the index is substantially smaller than the
input.

The SeqArena stores one full sequence per unique cluster and is only allocated
when `--error-correct` is active.
