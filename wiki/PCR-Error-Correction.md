# PCR error correction (Phase 3)

This page covers the second of two sources of false unique clusters in ancient
DNA deduplication: PCR copying errors. Post-mortem deamination is covered in
[[Damage-Aware-Deduplication]]. The two are complementary — damage masking
handles terminal variation, Phase 3 handles interior variation. See [[Home]]
for the full picture.

---

## The problem

PCR amplification introduces substitution errors at roughly 10⁻⁷ per base per
doubling for high-fidelity polymerases (Potapov & Ong 2017). Most errors occur
in early cycles, when the pool of template molecules is small, and therefore
propagate into a modest cluster of copies carrying the wrong base. After
deduplication, these clusters appear as genuine unique sequences — count 2 or 3
rather than count 1, which is enough to survive any simple singleton filter.

The distinguishing feature of a PCR error cluster is the count ratio: the
correct sequence has been amplified from the original molecule and has a high
read count, while the erroneous copy's cluster is a small fraction of that.
A single interior substitution separating a count=2 cluster from a count=500
cluster is almost never a distinct biological molecule.

---

## Approach

Phase 3 runs entirely in memory on the deduplication index built during Pass 1
— no additional file I/O. It classifies clusters as **parents** (count above
`--errcor-max-count`, default 5) and **children** (count at or below that
threshold). A child is absorbed into a parent when:

1. Their sequences differ by at most one substitution **outside** the damage zone
2. The parent's count is at least `--errcor-ratio` × the child's count (default 50×)

The damage zone exclusion is critical: a deaminated terminal C→T looks exactly
like a one-substitution PCR error. By restricting to the interior, Phase 3
avoids absorbing reads that differ legitimately due to ancient damage.

---

## The 3-way pigeonhole algorithm

Finding all pairs of sequences that differ by at most one substitution is
expensive with naive pairwise comparison. Phase 3 uses a pigeonhole argument
to make this efficient.

Divide the interior region into three equal parts. If two sequences of the same
length differ by exactly one substitution, that substitution falls in one of the
three parts — so the other two parts must be identical. Three hash indexes, each
covering a different pair of parts, are therefore sufficient to find every
Hamming-distance-1 neighbour without comparing all pairs:

```
interior: [─── part 0 ─── | ─── part 1 ─── | ─── part 2 ───]

mismatch in part 0 → parts 1 and 2 are identical → detected by index (1,2)
mismatch in part 1 → parts 0 and 2 are identical → detected by index (0,2)
mismatch in part 2 → parts 0 and 1 are identical → detected by index (0,1)
```

Phase 3a builds the three pair-hash indexes from all parent sequences. Phase 3b
queries them for each child. Candidates returned by the index are verified with
an exact Hamming check (AVX2 SIMD when available, scalar fallback), and verified
parents are checked against the count ratio. Sequences of different lengths are
automatically separated because the indexes are sharded by interior length.

---

## Parameters

`--errcor-ratio` (default 50) controls how large the count gap must be. A ratio
of 50 means a count=1 child needs a count≥50 parent. Raise it (100–200) if
your library contains genuine rare variants at low coverage, such as in
single-cell or low-input sequencing, or disable Phase 3 entirely with
`--error-correct` omitted.

`--errcor-max-count` (default 5) sets the child ceiling. Clusters with count
above 5 are never absorbed as PCR errors, regardless of their abundance relative
to any neighbour. Raising it risks absorbing true biological variants at low
coverage; the default is appropriate for typical libraries.

`--errcor-bucket-cap` (default 64) prevents low-complexity regions — poly-A
stretches, microsatellites — from producing enormous bucket collisions in the
pair-hash index. Insertions beyond the cap are silently dropped; affected
children may miss some parents, but these regions are also high false-positive
territory.

---

## Benchmarks

On sample `a88af16f35` (5.58 M reads input to `derep`, 91 bp mean length):

| Mode | Unique clusters | Change |
|------|----------------|--------|
| `derep --damage-auto` | 3,511,607 | — |
| `derep --damage-auto --error-correct` | 3,506,272 | −5,335 (−0.15%) |

5,335 PCR-error clusters absorbed in approximately 2 seconds. Phase 3 is cheap
because it works entirely on the in-memory index with no additional file reads.
