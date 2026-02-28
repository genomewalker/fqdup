# PCR error correction (Phase 3)

## Motivation

PCR amplification introduces substitution errors at a low but nonzero rate.
A cluster with count=1 that differs from a cluster with count=500 by a single
interior substitution almost certainly represents a copying error, not a
distinct original molecule.

Phase 3 identifies these clusters and flags them for removal. It runs
entirely in memory on the index built during Pass 1, so it adds no
extra file I/O.

---

## Definitions

| Term | Meaning |
|------|---------|
| **Parent** | A cluster with `count > --errcor-max-count` |
| **Child** | A cluster with `count ≤ --errcor-max-count` |
| **Interior** | Positions `[k5, L - k3)` — outside the damage zone |
| **Absorbed** | A child matched to a parent at interior Hamming distance ≤ 1 with count ratio ≥ `--errcor-ratio` |

---

## Interior region

Only positions **outside** the damage zone are examined. The damage zone
boundaries `k5` (5' end) and `k3` (3' end) come from the fitted damage model:
`k5` is the number of masked positions from the 5' end, `k3` from the 3' end.

```
seq: [ damaged zone 5' | ←── interior ──→ | damaged zone 3' ]
          k5 positions                         k3 positions
```

Interior length: `ilen = L - k5 - k3`.

This separation matters because a deaminated terminal C→T would look exactly
like a single-substitution PCR error. Excluding the damage zone prevents
legitimate damage variants from being absorbed as errors.

When damage parameters are not specified, `k5 = k3 = 0` and Phase 3 uses
the full read as the interior.

---

## Algorithm: 3-way pigeonhole

### Principle

If sequences A and B of the same length differ by exactly one substitution in
the interior of length `ilen`, then at least 2 of the 3 equal-width parts of
the interior are identical:

```
interior: [─── part 0 ─── | ─── part 1 ─── | ─── part 2 ───]

mismatch in part 0 → parts 1 and 2 match → detected by index (p1, p2)
mismatch in part 1 → parts 0 and 2 match → detected by index (p0, p2)
mismatch in part 2 → parts 0 and 1 match → detected by index (p0, p1)
```

Three pair-key indices over `(p0,p1)`, `(p0,p2)`, `(p1,p2)` are sufficient
to find every Hamming-distance-1 parent for every child. This is exact — no
approximation.

### Index construction (Phase 3a)

For each parent (count > `--errcor-max-count`):

1. Compute three part hashes: `h0 = XXH3(interior[0:s0])`, etc.
2. Insert into three pair-key maps, sharded by `ilen` (sequences of different
   lengths cannot be Hamming-distance-1 neighbours):

```
map_01[pair_key(h0, h1, tag=0, ilen)] → parent_id
map_02[pair_key(h0, h2, tag=1, ilen)] → parent_id
map_12[pair_key(h1, h2, tag=2, ilen)] → parent_id
```

`pair_key` normalises `(ha, hb)` by sorting before mixing, so the index is
symmetric. Bucket sizes are capped at `--errcor-bucket-cap` (default 64)
to prevent low-complexity regions from producing pathologically large buckets.

### Parent lookup (Phase 3b)

For each child (count ≤ `--errcor-max-count`):

1. Compute the child's three part hashes.
2. Query all three maps:

```
candidates ← map_01[pair_key(h0, h1)]
           ∪ map_02[pair_key(h0, h2)]
           ∪ map_12[pair_key(h1, h2)]
```

3. For each candidate parent, verify:
   - Not already flagged as a PCR error
   - `count(parent) ≥ --errcor-ratio × count(child)`
   - `interior_hamming(parent, child) ≤ 1`

   The Hamming check uses AVX2 SIMD (32 bases at a time) when available,
   falling back to a scalar loop.

4. If any parent passes, the child is marked absorbed.

### Epoch deduplication

The three queries for a single child may return overlapping candidate sets.
Each candidate is checked at most once per child via a `seen_epoch` vector:
each child increments a global epoch counter, and candidates with a matching
epoch are skipped. This avoids clearing a visited-set between children.

---

## Parameters

| Flag | Description | Default |
|------|-------------|---------|
| `--errcor-ratio FLOAT` | Minimum `count(parent)/count(child)` to absorb | 50.0 |
| `--errcor-max-count INT` | Children must have `count ≤ N` | 5 |
| `--errcor-bucket-cap INT` | Max candidates per pair-key bucket | 64 |

### Choosing `--errcor-ratio`

A ratio of 50 means a child with count=1 needs a parent with count ≥ 50. This
is conservative for typical libraries — a 1:50 ratio strongly suggests a PCR
copying error rather than a rare biological variant. For ultra-deep sequencing
where low-count variants are biologically meaningful, raise the ratio to 100–200
or disable Phase 3 entirely.

### Choosing `--errcor-max-count`

Only clusters with count ≤ `--errcor-max-count` are candidates for absorption.
The default of 5 works well for most libraries. Raising it risks absorbing
true low-coverage variants; lowering it reduces sensitivity.

### Choosing `--errcor-bucket-cap`

Caps bucket size in the pair-key index to prevent low-complexity regions
(poly-A, microsatellites) from dominating lookup time. Insertions beyond the
cap are silently discarded; affected children may miss some parents, but
low-complexity regions already have elevated false-positive rates.

---

## Complexity

| Operation | Cost |
|-----------|------|
| Phase 3a index | O(N_parents × ilen / 3) hashes |
| Phase 3b lookup | O(N_children × 3) pair-key queries |
| Phase 3b verify | O(N_candidates × ilen / 32) AVX2 Hamming checks |

Phase 3 typically adds 2–4 seconds to a run that otherwise takes 30 seconds.

---

## Benchmarks

On sample `a88af16f35` (5.58 M reads input to `derep`, 91 bp mean length):

| Mode | Unique clusters | Change |
|------|----------------|--------|
| `derep --damage-auto` | 3,511,607 | — |
| `derep --damage-auto --error-correct` | 3,506,272 | −5,335 (−0.15%) |

5,335 PCR-error clusters absorbed in ~2 seconds.

---

## Implementation notes

- `pair_key(ha, hb, tag, ilen)` sorts `(ha, hb)` before mixing to ensure
  commutativity. The `tag` (0, 1, or 2) and `ilen` are folded in to prevent
  cross-map collisions.

- `IndexEntry::count` is `uint64_t` — supports datasets with more than 4 G
  PCR copies of the same molecule without overflow.

- `SeqArena` stores sequences as a contiguous byte array with `uint32_t`
  offsets, bounds-checked on every append. Sequences longer than 65,535 bp
  are rejected (stored length is `uint16_t`).

- The `is_error_` vector is indexed by `seq_id` (arena index), not hash map
  slot, keeping Phase 3b iteration cache-friendly.
