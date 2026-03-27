# Algorithm

## fqdup extend

`fqdup extend` builds a self-referential de Bruijn k-mer graph from the merged
reads and extends each read outward from both ends. It runs three passes:

**Pass 0: damage estimation.** Samples the first 500k reads (configurable via
`--damage-sample`) to fit an exponential decay model to C→T and G→A frequencies.
Identifies terminal positions where excess damage exceeds `--mask-threshold`
(default 0.05). These positions are masked in Pass 1 to exclude damaged k-mers
from the graph. Pass 0 is skipped when `--no-damage` or `--mask-5`/`--mask-3` are
given.

**Pass 1: k-mer graph construction.** Streams all reads and inserts every
interior k-mer (excluding masked terminal positions) into an oriented k-mer
store. Oriented storage: for each canonical k-mer K, both K and RC(K) are
inserted so that left-extension is equivalent to a right-walk on the reverse
complement. This avoids separate forward/reverse traversal logic.

Key properties of the graph:
- **Exact 2-bit k-mer keys** (no hashing, no collision risk). k≤31 fits in a
  single `uint64_t`. Two separate stores are NOT kept, one canonical store
  covers all orientations.
- **Low-complexity filter**: k-mers with fewer than 3 distinct 2-mers (e.g.
  homopolymers) are discarded. This prevents spurious extensions in repetitive
  regions.
- **Minimum count filter**: k-mers with edge support below `--min-count`
  (default 2) are discarded after counting.
- **KMC3-in-RAM approach**: k-mer observations are accumulated in sharded
  append buffers, then sorted and reduced in place. On DS4 (198 M reads):
  6.7 B observations × 8 bytes demand-paged during accumulation; after
  finalize: 1.644 B distinct k-mers × 16 bytes = 26 GB. Peak RSS 41 GB.
- **Per-shard prefix index**: 2^12 bucket directory reduces binary search
  depth from 19 to 7 during Pass 2 lookups.

**Pass 2: extension and output.** Extends each read via a 3-stage pipeline:
reader → work queue → N worker threads → done queue → writer with reorder
buffer (to preserve input order).

Each worker applies the following algorithm to one read:
1. Find the rightmost clean (unmasked) k-mer as the right anchor; find the
   leftmost clean k-mer as the left anchor.
2. Walk right from the right anchor along the k-mer graph. A step is taken only
   when the transition is **reciprocally unique**: K→K' is unique AND the
   back-edge K'→K is also unique (unitig walk). Walk stops at `--max-extend`
   bases added.
3. Repeat in the reverse direction for left extension (implemented as a
   right-walk on the reverse complement).
4. n_skip logic: anchoring at position p re-derives L−(p+k) in-read bases
   before producing new sequence, anchoring never extends past the read end.
5. Reads with no clean interior k-mer are written unchanged.

Added bases receive quality `#` (Phred 2). Original base qualities are
preserved.

---

## fqdup sort

The sort is a standard external merge sort. The input is streamed and
accumulated into memory-limited chunks; each chunk is sorted in-place using
parallel `std::sort` and written to disk. Once all chunks are flushed, a
priority-queue k-way merge produces the single sorted output stream.

Chunks are gzip-compressed on disk by default. With `--fast`, chunks stay
uncompressed, roughly 3× faster I/O at the cost of 3× more temporary disk
space. For a 25 M-read library on NFS, the sort completes in about 20 seconds
with `--fast`; without it, closer to 60 seconds.

Reads are sorted lexicographically by read ID (the portion of the `@` header
line before the first space). Natural-numeric sort (`-N`) is also supported for
IDs with numeric suffixes.

---

## fqdup derep_pairs

`derep_pairs` is a two-pass algorithm over two sorted FASTQ files read in
lockstep:

- **`-n` (merged)**: fastp-merged reads, the original R1+R2 collapsed
  sequence, representing the actual ancient DNA molecule
- **`-e` (extended)**: the same reads after `fqdup extend` de Bruijn graph
  extension from both ends, longer, assembly-assisted sequences used as
  deduplication fingerprints

The extended file drives the hash. Using the `fqdup extend`-assembled sequence
as the fingerprint reduces false collisions: two different molecules that happen
to share a short merged core will typically diverge in the assembled extension.

**Pass 1** streams both files simultaneously, one pair at a time. For each
pair, the extended read's sequence is hashed to a canonical fingerprint
`(min(XXH3(seq), XXH3(revcomp(seq))), seq_len)`. If the fingerprint is new,
a new index entry is created recording the current record position and the
merged read's length. If it already exists, the count is incremented and, if
the current merged read is longer than the stored one, the representative is
updated to this pair. The representative selection rule, longest merged mate
- preserves the most original ancient DNA sequence.

**Pass 2** streams both files a second time. For each pair, if its record
position matches the representative stored in the index, it is written to the
output files. Cluster statistics (if `-c` is given) are written here as well.

No damage masking or error correction occurs in `derep_pairs`. It is a
structural step.

Memory: the index has one entry per unique cluster, approximately 40 bytes each
(20-byte IndexEntry plus hash-map overhead). For a 25 M-pair library with 78%
duplication, that is about 5.6 M entries, roughly 220 MB.

---

## fqdup derep

`derep` takes a single sorted FASTQ (either directly from `fqdup sort` or
the merged-read output of `derep_pairs`) and runs up to four phases depending
on which options are active.

**Pass 0** (with `--collapse-damage`) scans all reads to estimate the ancient DNA
damage model. It measures the T/(T+C) frequency at each of the first 15 5'-end
positions and the A/(A+G) frequency at the corresponding 3'-end positions.
Background rates come from the middle third of each read. A coverage-weighted
OLS fit in log-space (positions 1–9) gives the exponential decay parameters
`d_max` and `lambda` for each end. Position 0 is excluded because in sorted
files the most abundant reads appear first, which can skew terminal composition.
Lambda is clamped to [0.05, 0.5].

After fitting, positions where the observed damage excess exceeds
`--mask-threshold` are stored in a boolean array. All subsequent hashing uses
this precomputed array, no exponential evaluation at hash time.

**Pass 1** streams the input and builds the index. For each read, the sequence
is optionally masked (neutral bytes at C/T in masked 5' positions, G/A in
masked 3' positions), then hashed to a canonical fingerprint. New fingerprints
get a new index entry; existing fingerprints increment the count. Each unique
cluster's sequence is also stored in a 2-bit packed arena for Phase 3 (since
error correction is on by default).

**Phase 3** (PCR error correction, on by default) runs entirely in memory on the index. It
finds child clusters (count ≤ 5) that differ from a parent cluster (count ≥ 50×
the child's) by exactly one interior substitution and marks them as PCR errors.
See [[PCR-Error-Correction]] for how the 3-way pigeonhole search works.

**Pass 2** streams the input again and writes each read whose record position
is the representative of a non-absorbed cluster.

---

## Canonical hashing

All three subcommands use the same canonical hash:

```
canonical_hash(seq) = min(XXH3_128(seq), XXH3_128(revcomp(seq)))
```

Taking the minimum of the two 128-bit hashes means that a molecule sequenced
from either strand maps to the same fingerprint, forward and reverse-complement
reads collapse into the same cluster. The fingerprint key also includes the
sequence length to reduce false collisions between reads of different lengths
that happen to share a hash value. Using the full 128-bit hash reduces the
collision probability to ~3×10⁻²⁴ at 100 M reads.

[XXH3](https://github.com/Cyan4973/xxHash) exceeds 20 GB/s on modern
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
| SeqArena (default on) | ~0.25 × L_avg bytes × N_unique | ~80 MB (3.5 M × 91 bp) |

The SeqArena uses 2-bit packed encoding (A=00, C=01, G=10, T=11), storing
each base in 2 bits rather than 8, approximately 4× more compact than ASCII.
Sequences containing ambiguous bases (N) are stored but flagged ineligible;
they participate in deduplication but are skipped during Phase 3 error correction.
Use `--no-error-correct` to skip Phase 3 and avoid allocating the arena entirely.

For `derep_pairs`, the index similarly grows with unique pairs, not total pairs.
