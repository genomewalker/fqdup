# Algorithm

## fqdup sort

The sort is a standard external merge sort. The input is streamed and
accumulated into memory-limited chunks; each chunk is sorted in-place using
parallel `std::sort` and written to disk. Once all chunks are flushed, a
priority-queue k-way merge produces the single sorted output stream.

Chunks are gzip-compressed on disk by default. With `--fast`, chunks stay
uncompressed — roughly 3× faster I/O at the cost of 3× more temporary disk
space. For a 25 M-read library on NFS, the sort completes in about 20 seconds
with `--fast`; without it, closer to 60 seconds.

Reads are sorted lexicographically by read ID (the portion of the `@` header
line before the first space). Natural-numeric sort (`-N`) is also supported for
IDs with numeric suffixes.

---

## fqdup derep_pairs

`derep_pairs` is a two-pass algorithm over two sorted FASTQ files read in
lockstep. The extended file drives the hash; the non-extended file provides
the mate.

**Pass 1** streams both files simultaneously, one pair at a time. For each
pair, the extended read's sequence is hashed to a canonical fingerprint
`(min(XXH3(seq), XXH3(revcomp(seq))), seq_len)`. If the fingerprint is new,
a new index entry is created recording the current record position and the
non-extended read's length. If it already exists, the count is incremented
and, if the current non-extended read is longer than the stored one, the
representative is updated to this pair. The representative selection rule —
longest non-extended mate — preserves the most sequence from adapter trimming.

**Pass 2** streams both files a second time. For each pair, if its record
position matches the representative stored in the index, it is written to the
output files. Cluster statistics (if `-c` is given) are written here as well.

No damage masking or error correction occurs in `derep_pairs`. It is a
structural step.

Memory: the index has one entry per unique cluster, approximately 40 bytes each
(20-byte IndexEntry plus hash-map overhead). For a 25 M-pair library with 78%
duplication, that is about 5.6 M entries — roughly 220 MB.

---

## fqdup derep

`derep` takes a single sorted FASTQ and runs up to four phases depending on
which options are active.

**Pass 0** (with `--damage-auto`) scans all reads to estimate the ancient DNA
damage model. It measures the T/(T+C) frequency at each of the first 15 5'-end
positions and the A/(A+G) frequency at the corresponding 3'-end positions.
Background rates come from the middle third of each read. A coverage-weighted
OLS fit in log-space (positions 1–9) gives the exponential decay parameters
`d_max` and `lambda` for each end. Position 0 is excluded because in sorted
files the most abundant reads appear first, which can skew terminal composition.
Lambda is clamped to [0.05, 0.5].

After fitting, positions where the observed damage excess exceeds
`--mask-threshold` are stored in a boolean array. All subsequent hashing uses
this precomputed array — no exponential evaluation at hash time.

**Pass 1** streams the input and builds the index. For each read, the sequence
is optionally masked (neutral bytes at C/T in masked 5' positions, G/A in
masked 3' positions), then hashed to a canonical fingerprint. New fingerprints
get a new index entry; existing fingerprints increment the count. If
`--error-correct` is active, the full sequence of each unique cluster is also
stored in a contiguous sequence arena for Phase 3.

**Phase 3** (with `--error-correct`) runs entirely in memory on the index. It
finds child clusters (count ≤ 5) that differ from a parent cluster (count ≥ 50×
the child's) by exactly one interior substitution and marks them as PCR errors.
See [[PCR-Error-Correction]] for how the 3-way pigeonhole search works.

**Pass 2** streams the input again and writes each read whose record position
is the representative of a non-absorbed cluster.

---

## Canonical hashing

All three subcommands use the same canonical hash:

```
canonical_hash(seq) = min(XXH3_64(seq), XXH3_64(revcomp(seq)))
```

Taking the minimum of the two hashes means that a molecule sequenced from
either strand maps to the same fingerprint — forward and reverse-complement
reads collapse into the same cluster. The fingerprint key also includes the
sequence length to reduce false collisions between reads of different lengths
that happen to share a hash value.

[XXH3_64](https://github.com/Cyan4973/xxHash) exceeds 20 GB/s on modern
hardware with excellent distribution. The hash map is `ska::flat_hash_map`
(open-addressing, robin-hood probing), giving cache-friendly O(1) average-case
lookup and insert.

---

## Memory model

The index grows with unique clusters, not total reads. For libraries with high
duplication rates this makes a substantial difference.

| Component | Size | Example (5.6 M reads, 37% dup) |
|-----------|------|-------------------------------|
| `derep` index | ~40 bytes × N_unique | ~140 MB (3.5 M unique) |
| SeqArena (`--error-correct`) | ~L_avg bytes × N_unique | ~320 MB (3.5 M × 91 bp) |

The SeqArena stores one full sequence per unique cluster and is only allocated
when `--error-correct` is active. Without it, `derep` on 5.6 M reads uses
roughly 140 MB regardless of the duplication rate.

For `derep_pairs`, the index similarly grows with unique pairs, not total pairs.
