# fqdup extend

## Purpose

`fqdup extend` is the built-in de Bruijn graph read extender that replaces the
external Tadpole (BBTools) step in the fqdup pipeline. It takes fastp-merged
reads and extends each read outward from both ends, producing a longer
assembly-assisted fingerprint for the `derep_pairs` deduplication step.

### Why extension helps

**Shifted-window duplicates.** Two PCR duplicates from the same original
molecule may be captured in slightly different "windows", one read starts 2 bp
earlier or ends 3 bp later. After fastp merge, these reads have different
lengths and different sequences, so an exact-match deduplicator counts them as
distinct molecules. Extension recovers the shared flanking sequence present in
other reads of the library and grows both reads toward the same
assembly-assisted representation, so they hash identically and collapse
correctly. This "shifted window" problem is described in Rochette et al. (2023).

**QC-trim recovery.** Adapter and quality trimming produces PCR copies of the
same molecule at a range of lengths determined by where adapter contamination
and base-quality dropoff happen to fall in each copy. Because the copies are
trimmed at different positions, their sequences differ and `fqdup derep` (like
any exact-match deduplicator) treats them as distinct molecules. Extension
re-grows each trimmed copy outward along interior k-mers shared with longer
copies of the same molecule until they all reach the same fingerprint and
collapse.

The severity of this problem grows with length variability. In benchmarks
using synthetic data with Gaussian-distributed trim lengths
(read-length 65 bp, σ=15 bp, 30× coverage, no damage):

| Tool | Unique reported | Error vs truth |
|------|----------------|----------------|
| `fqdup derep` (no extend) | 4 035 | +1 245 % |
| `fqdup extend` → `fqdup derep` | 410 | +37 % |

With σ=5 bp (mild trimming), `fqdup extend` reaches +4.3 %, essentially
perfect. The residual error at σ=15 bp comes from reads trimmed so short
(<~20 bp) that no valid k-mer anchor exists; those are written unchanged by
`fqdup extend` and remain as separate clusters.

Note: QC-trim recovery works best without ancient DNA damage. If un-trimmed
copies carry intact terminal damage (C→T, G→A), those damaged k-mers appear in
the graph alongside the clean k-mers from trimmed copies, creating branches that
prevent consistent extension. For pre-trimmed input in a modern-DNA or
fully-trimmed aDNA workflow, use `--no-damage`.

---

## Algorithm

`fqdup extend` is a 3-pass algorithm.

### Pass 0: damage estimation

Samples the first N reads (default: 500k, set with `--damage-sample`) to fit an
exponential decay model of ancient DNA damage:

- Measures T/(T+C) at 5' positions 0–14 and A/(A+G) at 3' positions 0–14
- Background rates come from the middle third of each read
- A coverage-weighted OLS fit in log-space (positions 1–9) yields `d_max` and
  `lambda` for each end
- Positions where excess damage exceeds `--mask-threshold` (default 0.05) are
  marked for masking

These masked terminal positions are excluded from k-mer graph construction in
Pass 1. The damaged zone remains traversable in Pass 2 via k-mers contributed
by overlapping longer reads.

Pass 0 is skipped when `--no-damage` is given, or when `--mask-5`/`--mask-3`
are provided (explicit mask lengths skip damage estimation entirely).

### Pass 1: k-mer graph construction

Streams all reads and builds an oriented k-mer store from interior k-mers
(masked terminals excluded).

**Oriented storage.** For each canonical k-mer K, both K and RC(K) are
inserted. This means left-extension is equivalent to a right-walk on the
reverse complement, no separate forward/reverse traversal is needed.

**Exact 2-bit k-mer encoding.** Each base is packed into 2 bits (A=00, C=01,
G=10, T=11). For k≤31 the entire k-mer fits in a `uint64_t`. There is no
hashing, keys are exact with zero collision risk.

**Low-complexity filter.** K-mers with fewer than 3 distinct 2-mers are
discarded before insertion. This eliminates homopolymers and other repetitive
k-mers that would create spurious graph paths.

**Minimum count filter.** After accumulation, k-mers with edge support below
`--min-count` (default 2) are removed. This discards sequencing errors that
appear only once.

**KMC3-in-RAM approach.** K-mer observations are accumulated into sharded
append buffers (2^12 shards). After all reads are processed, each shard is
sorted and reduced in-place (counting unique k-mers). A per-shard prefix index
(2^12 bucket directory) reduces binary search depth from 19 to 7 during Pass 2.

Memory model (DS4, 198 M reads, k=17):
- Accumulation: 6.7 B observations × 8 bytes, demand-paged
- After finalize: 1.644 B distinct k-mers × 16 bytes = 26 GB
- Peak RSS: 41 GB (k-mer store + read buffers + thread overhead)

### Pass 2: extension and output

A 3-stage pipeline processes reads in parallel while preserving input order:

```
reader → work_queue → N workers → done_queue → writer (reorder buffer)
```

Each worker extends one read:

1. Find the rightmost clean (unmasked, in-range) k-mer as the right anchor;
   find the leftmost clean k-mer as the left anchor.
2. **Right extension**: walk right from the right anchor. A step K→K' is taken
   only when the edge is **reciprocally unique**: K→K' is the only forward edge
   AND K'→K is the only back-edge (unitig walk condition). Stop at `--max-extend`
   bases added.
3. **Left extension**: equivalent right-walk on the reverse complement of the
   read.
4. **n_skip logic**: when the anchor is at position p inside the read,
   L−(p+k) in-read bases are re-derived before new bases are produced. This
   ensures the extension never extends past the original read end.
5. Reads with no clean interior k-mer are written unchanged.

Added bases receive quality `#` (Phred 2). Original base qualities are
unchanged.

---

## Output format

`fqdup extend` produces a standard FASTQ file. The FASTQ header lines and
original sequence/quality are preserved; added bases are appended with quality
`#`. The output can be `.gz` (auto-detected from file extension) or plain text.

---

## CLI reference

```
fqdup extend -i INPUT -o OUTPUT [options]

Required:
  -i FILE                Input merged FASTQ (.gz or plain)
  -o FILE                Output extended FASTQ

K-mer graph:
  -k N                   K-mer size (default: 17, max: 31)
  --min-count N          Minimum edge support (default: 2)
  --max-extend N         Maximum bases added per side (default: 100)
  --threads N            Worker threads (default: all CPU cores)
  --min-qual N           Exclude bases below this Phred quality (default: 20)

Damage handling:
  --library-type TYPE    Library type for damage model: auto|ds|ss (default: auto)
  --no-damage            Skip damage estimation; no terminal masking
  --mask-5 N             Manually mask N bp at 5' end (skips Pass 0)
  --mask-3 N             Manually mask N bp at 3' end (skips Pass 0)
  --mask-threshold F     Excess damage threshold for masking (default: 0.05)
  --damage-sample N      Use first N reads for damage estimation (default: 500000; 0=all)
```

---

## Performance

Benchmarked on DS4, 198 M reads, 16 threads, dandycomp02fl:

| Metric | Value |
|--------|-------|
| Reads extended | 50% |
| Average extension | 3.82 / 3.83 bp per side |
| Pass 1 (build) | 68 s |
| Finalize | 19 s |
| Pass 2 (extend + write) | 5:27 |
| Total wall time | 6:56 |
| Peak RSS | 41.3 GB |

Comparison with Tadpole on the same DS4 dataset:

| Tool | Extension rate | Wall time |
|------|---------------|-----------|
| Tadpole (k=17, el=100, trimends=9) | 14.75% | ~5 min |
| `fqdup extend` (k=17, --max-extend 100, 16 threads) | 50% | 7 min |

`fqdup extend` achieves significantly higher extension rates because it uses
all reads in the library to build the graph, reads that carry flanking context
for other reads contribute their k-mers regardless of their own extension
status.

---

## References

- Rochette NC et al. (2023) On the causes, consequences, and avoidance of PCR
  duplicates. *Mol Ecol Resour* 23:1299–1318. doi:10.1111/1755-0998.13800
- Briggs et al. (2007) Patterns of damage in genomic DNA sequences from a
  Neandertal. *PNAS*. doi:10.1073/pnas.0704665104
