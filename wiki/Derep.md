# fqdup derep

## Purpose

`fqdup derep` is the biological deduplication step. It collapses reads that
represent the same original molecule into a single representative, handling two
sources of false unique clusters specific to ancient DNA:

- **Post-mortem deamination**, C→T / G→A substitutions at fragment termini
  split true duplicates into apparent distinct sequences. Damage-aware hashing
  replaces affected terminal positions with a neutral byte before hashing,
  collapsing them. See [[Damage-Aware-Deduplication]].

- **PCR copying errors**, low-count clusters that differ from a high-count
  cluster by a single interior substitution. Phase 3 absorbs these as PCR
  artefacts. See [[PCR-Error-Correction]].

`fqdup derep` can be run directly on sorted merged reads or after
`fqdup derep_pairs` in the full paired pipeline. When following `derep_pairs`,
it handles residual duplicates: reads whose extended fingerprints diverged
slightly (keeping them through `derep_pairs`) but whose original merged
sequences are identical.

---

## Algorithm

`fqdup derep` runs in up to four stages:

**Pass 0: damage estimation** (only with `--collapse-damage`)
Scans all reads to fit an exponential decay model of terminal deamination
(same model as [[Damage-Aware-Deduplication]]). Identifies which positions
exceed the mask threshold. Skipped if `--collapse-damage` is not given (default), or manual
parameters are supplied.

**Pass 1: index construction**
Streams all reads and builds a hash index. With damage masking active, masked
terminal positions are replaced by neutral bytes before hashing. The canonical
hash `min(XXH3_128(seq), XXH3_128(revcomp(seq)))` ensures forward and
reverse-complement reads land in the same cluster. The representative kept per
cluster is the first occurrence.

**Phase 3: PCR error correction**
After Pass 1, a 3-way pigeonhole Hamming search identifies clusters with low
count that differ from a high-count cluster by exactly one interior
substitution. These are absorbed as PCR copying errors. C↔T and G↔A
mismatches are never absorbed (indistinguishable from damage). Enabled by
default; disable with `--no-error-correct`.

**Pass 2: output**
Streams reads again, writing one representative per cluster.

---

## CLI reference

```
fqdup derep -i INPUT -o OUTPUT [options]

Required:
  -i FILE              Sorted input FASTQ (.gz or plain)
  -o FILE              Output FASTQ

Optional:
  -c FILE              Cluster statistics (gzipped TSV)
  --no-revcomp         Disable reverse-complement collapsing (default: enabled)
```

### Damage variant collapsing (default: off)

| Flag | Description | Default |
|------|-------------|---------|
| `--collapse-damage` | Estimate damage from data (Pass 0) | off |
| `--damage-dmax FLOAT` | d_max for both ends (manual) | - |
| `--damage-dmax5 FLOAT` | d_max for 5' end | - |
| `--damage-dmax3 FLOAT` | d_max for 3' end | - |
| `--damage-lambda FLOAT` | Decay constant for both ends | - |
| `--damage-lambda5 FLOAT` | Decay constant for 5' end | - |
| `--damage-lambda3 FLOAT` | Decay constant for 3' end | - |
| `--damage-bg FLOAT` | Background substitution rate | 0.02 |
| `--mask-threshold FLOAT` | Mask when excess P(deam) > T | 0.05 |
| `--library-type auto\|ds\|ss` | Override library-type detection | auto |

Run `fqdup damage -i FILE` first to inspect d_max and which positions would
be masked before committing to `--collapse-damage`. See [[Damage]].

### PCR error correction (default: on)

| Flag | Description | Default |
|------|-------------|---------|
| `--error-correct` | Enable Phase 3 | **default** |
| `--no-error-correct` | Disable Phase 3 | - |
| `--errcor-mode capture\|shotgun` | Coverage regime | shotgun |
| `--errcor-min-parent INT` | Min count to index as parent | 3 |
| `--errcor-snp-threshold FLOAT` | SNP veto: sig/parent_count threshold | 0.20 |
| `--errcor-snp-min-count INT` | SNP veto: min absolute sig_count | 2 (capture), 1 (shotgun) |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 64 |
| `--pcr-cycles INT` | PCR cycles for D_eff log estimate | 0 (auto) |
| `--pcr-efficiency FLOAT` | Amplification efficiency per cycle | 1.0 |
| `--pcr-error-rate FLOAT` | Substitutions per base per doubling (log only) | 5.3e-7 |

**Mode:**
- `shotgun` (default): conservative — `snp_min_count=1`, Phase B3 disabled;
  appropriate when true molecules may appear only once.
- `capture`: aggressive — `snp_min_count=2`, Phase B3 (iterative small-cluster
  pass) enabled; appropriate for deep-coverage capture libraries.

---

## Output

### Reads

Standard FASTQ. One read per unique cluster (the representative, first
occurrence in sorted input).

### Cluster statistics (`-c FILE.tsv.gz`)

Gzipped TSV, one row per cluster:

| Column | Description |
|--------|-------------|
| `hash` | 128-bit XXH3 canonical hash (32 hex chars) |
| `seq_len` | Sequence length of the representative |
| `count` | Number of reads in the cluster |

---

## Benchmarks

On sample `a88af16f35`, 5.58 M reads in, 91 bp mean length,
`d_max_5 ≈ 0.193`:

| Mode | Unique clusters | Wall time |
|------|----------------|-----------|
| Standard (exact hash, EC on) | 3,531,821 | ~22 s |
| `--collapse-damage` only | 3,511,607 | ~31 s |
| `--collapse-damage --error-correct` | 3,510,151 | ~33 s |

Across 31 ancient DNA libraries (2.97 B reads total) run in the full paired
pipeline, `fqdup derep` absorbed a further 8.6% of reads beyond what
`fqdup derep_pairs` kept.

See [[Performance]] for full benchmarks and memory usage.
