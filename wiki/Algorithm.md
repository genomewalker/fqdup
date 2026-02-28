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

Chunks are gzip-compressed on disk by default. With `--fast`, chunks stay
uncompressed (~3× faster I/O, ~3× more temporary disk usage).

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

### Pass 1 — Index construction

Both sorted files are read simultaneously, one pair at a time. The reads at
each position are expected to share the same read ID (a warning is logged if
they differ). For each pair `(ext, non)`:

1. **Hash the extended read:**

   ```
   h = canonical_hash(ext.seq) = min(XXH3_64(ext.seq), XXH3_64(revcomp(ext.seq)))
   ```

   The `min()` over the sequence and its reverse complement ensures that a
   molecule sequenced from either strand lands in the same cluster.

2. **Form a fingerprint:** `(h, len(ext.seq))` — the sequence length is
   included in the key to reduce false collisions between reads of different
   lengths that happen to share a hash prefix.

3. **Insert into the index:**
   - New fingerprint: `IndexEntry{ record_idx=current, non_len=len(non.seq), count=1 }`
   - Existing fingerprint:
     - Increment `count`
     - If `len(non.seq) > stored non_len`: update `record_idx` and `non_len`
       (the current record becomes the new representative)

   **Representative selection rule:** the pair with the longest non-extended
   read is chosen as representative. Longer non-extended reads contain more
   sequence information — they were trimmed less aggressively during adapter
   removal.

Memory usage: approximately `24 × N` bytes where N is the total number of
read pairs (regardless of duplication level, because we store one entry per
input record in the index map).

### Pass 2 — Output

The index tells us which record index holds the representative for each
fingerprint. We build a reverse map `record_idx → fingerprint`, then stream
both files again:

```
for each pair (ext, non) at position record_idx:
    if record_idx in records_to_write:
        write ext to out_ext
        write non to out_non
        (optionally) write cluster stats to -c file
```

No damage masking or error correction occurs in this step.

---

## fqdup derep

Single-file deduplication with optional damage masking and PCR error
correction. Takes sorted FASTQ as input (typically the non-extended output
of `derep_pairs`).

### Overview

```
[Pass 0]  (damage-aware only)   Stride-sample → fit damage model
              ↓
[Pass 1]  Stream input → build hash → IndexEntry map
              ↓
[Phase 3] (error-correct only)  Absorb PCR 1-error variants
              ↓
[Pass 2]  Stream again → write representative records
```

### Pass 0 — Damage estimation

See [[Damage-Aware-Deduplication]] for the full model description.

In brief: every Nth read (default N=1000, targeting 100 k sampled reads)
is inspected. For each of the first 15 terminal positions, T/(T+C) (5' end)
and A/(A+G) (3' end) frequencies are measured. Baseline rates come from the
middle third of reads. An exponential decay model is fitted independently to
each end:

```
P_CT(pos) = d_max_5 × exp(-lambda_5 × pos) + bg
P_GA(pos) = d_max_3 × exp(-lambda_3 × pos) + bg
```

The fit uses coverage-weighted OLS on positions 1–9 (position 0 is excluded
— sorted files place the most common reads first, which can give atypical
terminal composition at position 0). Lambda is clamped to [0.05, 0.5].

Damage is considered significant when `d_max_5 > 0.02` or `d_max_3 > 0.02`.
If damage is below threshold, Pass 0 is a no-op and standard exact hashing
is used.

### Pass 1 — Index construction

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
   excess `max(P_CT(i), P_GA(i))` exceeds `--mask-threshold`. This makes
   independently-deaminated copies of the same molecule hash identically.
   See [[Damage-Aware-Deduplication]] for why symmetry matters.

2. **Form fingerprint:** `(h, len(rec.seq))`

3. **Insert into index:**
   - New fingerprint: `IndexEntry{ record_idx=current, seq_id=arena.append(seq), count=1 }`
     (`seq_id` and arena append only when `--error-correct` is active)
   - Existing fingerprint: increment `count` (representative stays as the
     first occurrence — position selection is irrelevant here since all reads
     of a single-file cluster are biologically equivalent)

### Phase 3 — PCR error correction

See [[PCR-Error-Correction]] for full details.

### Pass 2 — Output

Same structure as `derep_pairs` Pass 2, but for a single file:

```
for each read at position record_idx:
    if record_idx in records_to_write (and not flagged as PCR error):
        write rec to output
        (optionally) write cluster stats to -c file
```

---

## Canonical Hashing

```
canonical_hash(seq) = min(XXH3_64(seq), XXH3_64(revcomp(seq)))
```

[XXH3_64](https://github.com/Cyan4973/xxHash) is a non-cryptographic hash
with excellent speed (≥ 20 GB/s) and distribution. It is used exclusively for
deduplication — no cryptographic guarantees are assumed.

`ska::flat_hash_map` (open-addressing with robin-hood probing) provides
cache-friendly O(1) average-case lookup and insert. The composite key
`(hash, seq_len)` is used to reduce false-positive hash collisions between
reads of different lengths.

---

## Memory Model

| Component | Formula | Example (400 M reads, 50% dup) |
|-----------|---------|-------------------------------|
| `derep_pairs` Pass 1 index | ~24 bytes × N_pairs | ~9.6 GB |
| `derep` Pass 1 index | ~16 bytes × N_reads | ~3.2 GB (200 M unique) |
| SeqArena (error-correct) | ~1 byte/bp × N_unique × L_avg | ~13 GB (200 M × 65 bp) |
| Phase 3 pair-key index | ~O(N_parents) | small |

The `derep_pairs` index is proportional to **total** input read pairs (not unique
clusters), because we must track every record's position to know which one is the
representative. The `derep` index grows with unique clusters only if `--error-correct`
is used (because the SeqArena only stores one sequence per unique cluster).

Without `--error-correct`, `derep` uses ~16 bytes × total reads, regardless of
duplication level.
