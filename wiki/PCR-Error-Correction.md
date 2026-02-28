# PCR Error Correction (Phase 3)

## Motivation

PCR amplification introduces substitution errors at a low but nonzero rate.
A read cluster with count=1 that differs from a cluster with count=500 by a
single substitution in a non-damaged interior position almost certainly
represents a PCR copying error, not a distinct original molecule.

Phase 3 identifies these clusters and removes them, reducing the count of
spurious unique sequences in the output. It operates on the deduplication
index (after Pass 1), so it adds no extra I/O.

---

## Definitions

| Term | Meaning |
|------|---------|
| **Parent** | A cluster with `count > --errcor-max-count` |
| **Child** | A cluster with `count ≤ --errcor-max-count` |
| **Interior** | Positions `[k5, L - k3)` — outside the damage zone |
| **Absorbed** | A child matched to a parent at interior Hamming distance ≤ 1 with count ratio ≥ `--errcor-ratio` |

---

## Interior Region

Only positions **outside** the damage zone are examined. The damage zone
boundaries `k5` (from the 5' end) and `k3` (from the 3' end) are derived from
the fitted damage model: `k5` is the largest position where the 5' damage
excess still exceeds `--mask-threshold`, and similarly for `k3`.

```
seq: [ damaged zone 5' | ←── interior ──→ | damaged zone 3' ]
          k5 positions                         k3 positions
```

The interior length is `ilen = L - k5 - k3`.

This separation is essential for two reasons:

1. **Prevents false negatives**: a deaminated terminal C→T would look like a
   single substitution from its parent, but it is not a PCR error. Excluding
   the damage zone avoids absorbing genuinely distinct reads that merely carry
   terminal damage.

2. **Prevents false positives**: two reads from different molecules that happen
   to share terminal sequence but differ in the interior would not be falsely
   merged.

When `--damage-auto` or manual damage parameters are not specified, `k5 = k3 = 0`
and Phase 3 uses the full sequence interior.

---

## Algorithm: 3-way Pigeonhole

### Principle

If sequences A and B of the same length differ by exactly one substitution in
the interior of length `ilen`, then at least 2 of the 3 equal-width parts of
the interior must be identical:

```
interior: [─── part 0 ─── | ─── part 1 ─── | ─── part 2 ───]

If mismatch falls in part 0: parts 1 and 2 match → detected by index (p1, p2)
If mismatch falls in part 1: parts 0 and 2 match → detected by index (p0, p2)
If mismatch falls in part 2: parts 0 and 1 match → detected by index (p0, p1)
```

Therefore, three pair-key indices covering `(p0,p1)`, `(p0,p2)`, `(p1,p2)`
are sufficient to detect every edit-distance-1 parent for every child. No
approximation — this is exact for Hamming distance 1 in the interior.

### Index construction (Phase 3a)

For each parent sequence (count > `--errcor-max-count`):

1. Compute the three part hashes: `h0 = XXH3(interior[0:s0])`, etc.
2. Insert into three pair-key maps, sharded by `ilen` (sequences of different
   lengths cannot be Hamming-distance-1 neighbours):

```
map_01[pair_key(h0, h1, tag=0, ilen)] → parent_id
map_02[pair_key(h0, h2, tag=1, ilen)] → parent_id
map_12[pair_key(h1, h2, tag=2, ilen)] → parent_id
```

`pair_key` is commutative (sorts `ha`, `hb` before mixing), so
`pair_key(h0, h1, ...) == pair_key(h1, h0, ...)` — the index is symmetric.

Bucket sizes are capped at `--errcor-bucket-cap` (default 64) to prevent
low-complexity or repeat regions from creating extremely large buckets that
would slow down Phase 3b.

### Parent lookup (Phase 3b)

For each child sequence (count ≤ `--errcor-max-count`):

1. Compute the child's three part hashes.
2. Query all three pair-key maps:

```
candidates ← map_01[pair_key(h0, h1)] ∪ map_02[pair_key(h0, h2)] ∪ map_12[pair_key(h1, h2)]
```

3. For each candidate parent, verify:
   - Same interior length (should always hold due to length sharding, but
     checked for safety)
   - Not already flagged as a PCR error
   - `count(parent) ≥ --errcor-ratio × count(child)`
   - `interior_hamming(parent, child) ≤ 1`

   The Hamming check uses AVX2 SIMD (32 bases at a time) when available,
   falling back to a scalar loop.

4. If any parent passes all checks, the child is marked as absorbed.

### Epoch deduplication

The three pair-key queries for a single child may return overlapping candidate
sets (a parent could appear in multiple maps). Each candidate is checked at
most once per child using a `seen_epoch` vector: each child increments a
global epoch counter, and candidates are skipped if their recorded epoch
matches the current one.

This avoids the cost of clearing a visited-set between children.

---

## Parameters

| Flag | Description | Default |
|------|-------------|---------|
| `--errcor-ratio FLOAT` | Minimum `count(parent)/count(child)` to absorb | 50.0 |
| `--errcor-max-count INT` | Children must have `count ≤ N` | 5 |
| `--errcor-bucket-cap INT` | Max candidates per pair-key bucket | 64 |

### Choosing `--errcor-ratio`

The default ratio of 50 means a child with count=1 requires a parent with
count ≥ 50, a child with count=2 requires count ≥ 100, etc. This is
conservative — a 1-in-50 abundance suggests a PCR error is far more likely
than a rare biological variant.

For highly sensitive applications (e.g. single-cell or ultra-deep sequencing)
where low-count unique sequences are biologically meaningful, increase the
ratio (100–200) or disable Phase 3 entirely.

### Choosing `--errcor-max-count`

Only clusters with count ≤ `--errcor-max-count` are candidates for absorption.
The default of 5 is appropriate for typical libraries. Raising it risks
absorbing true biological variants that happen to have low coverage; lowering
it is safe but reduces sensitivity.

### Choosing `--errcor-bucket-cap`

This prevents low-complexity regions (e.g. poly-A stretches, microsatellites)
from dominating the pair-key index. When a bucket exceeds the cap, additional
insertions are silently discarded. Children in low-complexity regions may
therefore not find all of their parents, but this is acceptable — those regions
are also likely to have high false positive rates regardless.

---

## Complexity

| Operation | Cost |
|-----------|------|
| Phase 3a index | O(N_parents × ilen / 3) hashes |
| Phase 3b lookup | O(N_children × 3) pair-key queries |
| Phase 3b verify | O(N_candidates × ilen / 32) AVX2 Hamming checks |

In practice Phase 3 adds ~2–3 seconds to a 30-second deduplication run on
25 M reads. The pair-key index is small compared to the main index (only
parents are stored).

---

## Benchmarks

On a 25.8 M-read ancient DNA library:

| Mode | Unique clusters | Change |
|------|----------------|--------|
| Damage-aware only | 5,547,508 | — |
| + Error correction (ratio=50, max-count=5) | 5,532,327 | −15,181 |

15,181 PCR error clusters absorbed in ~2 seconds.

---

## Implementation Notes

- `pair_key(ha, hb, tag, ilen)` normalises `(ha, hb)` by sorting before
  mixing, ensuring commutativity. The `tag` (0, 1, or 2) and `ilen` are
  folded in to prevent collisions between the three maps.

- `IndexEntry::count` is `uint64_t` — supports datasets with > 4 G PCR
  copies of the same molecule without silent overflow.

- `SeqArena` stores sequences as byte arrays with `uint32_t` offsets. It
  is bounds-checked before every append and raises `std::runtime_error` on
  overflow. Sequences longer than 65,535 bp are rejected (stored length is
  `uint16_t`).

- The `is_error_` vector is indexed by `seq_id` (arena index), not by hash
  map slot, so its iteration order during Phase 3b is cache-friendly.
