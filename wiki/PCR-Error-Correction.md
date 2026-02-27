# PCR Error Correction (Phase 3)

## Motivation

PCR amplification introduces substitution errors at a low rate. A cluster
with count=1 that differs from a cluster with count=500 by a single
substitution in a non-damaged position is almost certainly a PCR error, not
a distinct molecule. Phase 3 identifies and removes these error clusters.

## Algorithm

### Definitions

- **Parent**: a cluster with `count > --errcor-max-count`
- **Child**: a cluster with `count ≤ --errcor-max-count`
- **Absorbed**: a child whose parent satisfies `count(parent) ≥ ratio × count(child)`
  and whose interior sequence differs from the parent by at most 1 substitution

### Interior region

Only positions **outside** the damage zone are used for error correction.
The damage zone boundaries `(k5, k3)` are computed from the fitted damage
model. The **interior** is `seq[k5 : L-k3]`.

This prevents deaminated terminal positions from masking real PCR errors or
creating false positives.

### 3-way pigeonhole principle

If two sequences of the same length differ by exactly 1 substitution in the
interior of length `ilen`, at least 2 of the 3 equal-width parts of the
interior must match exactly.

`fqdup` builds three pair-key indices over the parent set, one per pair of
parts `(p0,p1)`, `(p0,p2)`, `(p1,p2)`. For each child, three hash lookups
cover all possible edit-distance-1 parents. Candidates are verified with an
AVX2-accelerated interior Hamming check.

```
interior = [k5 ─────────────────── L-k3)
              p0      p1      p2
pair-key index 0: hash(p0) ⊕ hash(p1)
pair-key index 1: hash(p0) ⊕ hash(p2)
pair-key index 2: hash(p1) ⊕ hash(p2)
```

Sequences are sharded by interior length — prevents cross-length false
candidates (edit distance 1 cannot change length).

### Epoch deduplication

Each child makes 3 pair-key queries that may return overlapping candidate
sets. A `seen_epoch` vector (uint64_t, one per unique sequence) ensures each
candidate is checked at most once per child, without clearing the vector
between children.

### Complexity

- Space: O(N) for epoch vector + O(parents) for pair-key index
- Time: O(N) pair-key lookups + O(candidates × ilen / 32) AVX2 checks
- In practice: adds ~2 s to a 30 s run on 25 M reads

## Parameters

```
--error-correct         Enable Phase 3
--errcor-ratio  FLOAT   count(parent)/count(child) ≥ ratio to absorb
                        Default: 50.0 — a 1-in-50 chance the child is real
--errcor-max-count INT  Only children with count ≤ N are candidates
                        Default: 5
--errcor-bucket-cap INT Max bucket size in pair-key index
                        Default: 64 — caps low-complexity / repeat regions
```

## Effect on unique cluster counts

On a 25.8 M-read ancient DNA library:

| Mode | Unique clusters |
|------|----------------|
| Damage-aware | 5,547,508 |
| + Error correction (ratio=50, max-count=5) | 5,532,327 (−15,181) |

Phase 3 absorbed 15,181 PCR error clusters in ~2 seconds.

## Implementation notes

- `pair_key(ha, hb, tag, ilen)` is **truly commutative**: `ha` and `hb` are
  sorted before mixing so `pair_key(h0, h1, ...) == pair_key(h1, h0, ...)`.
- `epoch` and `seen_epoch` are `uint64_t` — no wrap-around risk.
- `IndexEntry::count` is `uint64_t` — handles datasets with > 4 G PCR
  copies of the same molecule without overflow.
- All arena sizes and pool sizes are bounds-checked before uint32_t casts
  with descriptive `std::runtime_error` throws.
