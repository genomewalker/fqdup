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
After Pass 1, a 4-way pigeonhole Hamming search identifies clusters with low
count that differ from a high-count cluster by one or two interior
substitutions (H≤2). These are absorbed as PCR copying errors. C↔T, G↔A,
G↔T, and C↔A mismatches are never absorbed (damage-consistent). A↔T and C↔G
transversions are eligible by default; use `--protect-transversions` to protect
them too (recommended for high-oxidative-damage libraries). H=2 absorption
requires both mismatches to be eligible (non-protected) transversions and the
child count to be ≤ `--errcor-max-h2-count` (default 2). Enabled by default;
disable with `--no-error-correct`.

**Pass 2: output**
Streams reads again, writing one representative per cluster. When
`--out-damaged` / `--out-undamaged` are given, each read is also scored by
the damage split model (see [Damaged/undamaged split](#damagedundamaged-split))
and routed to the appropriate output file.

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

Run `fqdup profile -i FILE` first to inspect d_max and which positions would
be masked before committing to `--collapse-damage`. See [[Damage]].

### Damaged/undamaged split

After deduplication, `derep` can route each representative read to one of two output
files based on its per-read LLR ancient/modern score:

```bash
fqdup derep -i sorted.fq.gz -o dedup.fq.gz \
    --out-damaged   ancient.fq.gz \
    --out-undamaged modern.fq.gz \
    --damage collapse
```

`-o` is optional when both split outputs are given (defaults to `/dev/null`).

| Flag | Description | Default |
|------|-------------|---------|
| `--out-damaged FILE` | Write reads with LLR ≥ threshold here | - |
| `--out-undamaged FILE` | Write reads with LLR < threshold here | - |
| `--split-threshold F` | LLR decision boundary | 0.0 |
| `--split-model MODE` | `auto` \| `bulk` \| `empirical` | `auto` |

#### `--split-model` modes

| Mode | Extra file pass | Accuracy | When to use |
|------|----------------|----------|-------------|
| `auto` | only if d_max > 0.01 | empirical per-bin | **default** — optimal for most aDNA libraries |
| `bulk` | never | bulk exponential | quick splits, exploratory work, very large files |
| `empirical` | always | empirical per-bin | force even when damage is borderline |

**`auto` (default):** runs a stripped length-stratified scan when damage is
detected (d_max > 0.01), building per-bin empirical C→T rate curves via mixture
unmixing. Falls back to the bulk exponential model (zero extra I/O) when no
damage is detected. The empirical scan uses precomputed additive classifier
coefficients (no log() calls in the read loop) and is capped at 4 worker threads
to match the single FASTQ reader throughput.

**`bulk`:** uses the exponential decay model already fitted during Pass 0
(d_max, λ, bg). No extra file read. Slightly less accurate at very short
(<50 bp) and very long (>150 bp) reads where the exponential approximation
deviates from the empirical curves.

**`empirical`:** same as `auto` but forces the per-bin scan even if damage
appears low (e.g. heavily contaminated libraries where contamination suppresses
the bulk d_max).

#### LLR scoring

For each read of length L, the model finds the matching length bin and computes:

```
LLR = Σ_{i=1}^{n_pos}  T_i·log(p_dam[i]/p_und[i])  +  C_i·log((1-p_dam[i])/(1-p_und[i]))
```

where `p_dam[i]` and `p_und[i]` are the empirical C→T rates at position i in
the damaged and undamaged read classes respectively (recovered by mixture
unmixing from bulk + per-class rates). For SS libraries the 3′ end is scored
symmetrically. LLR > 0 → classified as damaged.

Performance on a 437 M-read SS library (d_max ≈ 0.16, 3 length bins):

| Phase | Time |
|-------|------|
| Bulk damage fit (Pass 0) | ~50 s |
| Split model scan (`auto`) | ~150 s |
| Derep pass 1 (index) | ~1650 s |

### PCR error correction (default: on)

| Flag | Description | Default |
|------|-------------|---------|
| `--error-correct` | Enable Phase 3 | **default** |
| `--no-error-correct` | Disable Phase 3 | - |
| `--errcor-min-parent INT` | Min count to index as parent | 3 |
| `--errcor-max-h2-count INT` | Max child count eligible for H=2 absorption | 2 |
| `--errcor-snp-threshold FLOAT` | SNP veto: sig/parent_count threshold | 0.20 |
| `--errcor-snp-min-count INT` | SNP veto: min absolute sig_count | 2 |
| `--errcor-bucket-cap INT` | Pair-key bucket size cap | 0 (unlimited) |
| `--protect-transversions` | Protect A↔T / C↔G (Channels H/G) from absorption | off |
| `--pcr-cycles INT` | PCR cycles for D_eff log estimate | 0 (auto) |
| `--pcr-efficiency FLOAT` | Amplification efficiency per cycle | 1.0 |
| `--pcr-error-rate FLOAT` | Substitutions per base per doubling (log only) | 5.3e-7 |

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
