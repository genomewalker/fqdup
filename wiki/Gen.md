# fqdup gen

Line references are to `src/gen.cpp` at `fb62967`.

## Status

Internal/dev tool: a synthetic FASTQ generator for testing and benchmarking
against a known ground truth. It is listed as a subcommand
(`fqdup --help`, `src/main.cpp:36`: "Generate synthetic FASTQ with
configurable damage patterns") but does not appear in the typical pipeline
workflow (`src/main.cpp:42-49`) — it produces test input, not pipeline
output.

## Purpose

Generates IID reads at a fixed length and GC content, split into an "ancient"
and a "modern" class (`--f-damaged` sets the ancient fraction), then applies
three independent damage channels to the ancient reads plus background
substitution noise to all reads (file header comment, `gen.cpp:1-9`):

- **Channel A — deamination**: 5' C→T, exponential decay from position 0
  (amplitude `--dmg5-max`, rate `--dmg5-lambda`). For `ds` libraries also 3'
  G→A with its own amplitude/rate; for `ss` libraries additionally 3' C→T.
- **Channel D — 8-oxoG**: uniform per-position G→T transversion probability
  (`--oxog`), independent of read position.
- **Channel E — depurination**: with probability `--depur`, a pyrimidine
  (C/T) at position 0 of an ancient read is replaced by a purine (A or G,
  drawn in proportion to base composition) (`gen_main` :221-225).

Every read's header records its ground-truth class:
`@SYN:<index>:damaged` or `@SYN:<index>:undamaged` (:248-250).

## CLI flags

`print_usage` :21-50. All arguments take a value; `-h`/`--help` prints this
and exits 0.

| Flag | Effect | Default |
|---|---|---|
| `-o, --output FILE` | Output FASTQ path (`.gz` or plain, by extension) | required |
| `-n, --reads N` | Number of reads | 1000000 |
| `--seed N` | PRNG seed (`mt19937_64`) | 1 |
| `--read-len N` | Fixed read length | 60 |
| `--gc FLOAT` | GC content fraction | 0.45 |
| `--f-damaged FLOAT` | Fraction of reads drawn as "ancient" | 0.70 |
| `--lib-type ds\|ss` | Library type; gates 3' C→T and default `--dmg3-ct-max` | ds |
| `--dmg5-max FLOAT` | 5' C→T amplitude at position 0 | 0.18 |
| `--dmg5-lambda FLOAT` | 5' C→T decay rate | 0.35 |
| `--dmg3-ga-max FLOAT` | 3' G→A amplitude | `--dmg5-max` |
| `--dmg3-ga-lambda FLOAT` | 3' G→A decay rate | `--dmg5-lambda` |
| `--dmg3-ct-max FLOAT` | 3' C→T amplitude (adds to Channel A on `ss`) | 0 on `ds`, `--dmg5-max` on `ss` |
| `--dmg3-ct-lambda FLOAT` | 3' C→T decay rate | `--dmg5-lambda` |
| `--oxog FLOAT` | Per-G G→T probability (Channel D) | 0.0 |
| `--depur FLOAT` | Position-0 purine-enrichment rate on ancient reads (Channel E) | 0.0 |
| `--bg-ct FLOAT` | Background C→T rate, all reads, all positions | 0.001 |
| `--bg-ga FLOAT` | Background G→A rate, all reads, all positions | 0.001 |
| `--q-score N` | Constant Phred quality written for every base | 40 |

Validation: `--read-len` must be ≥1, `--q-score` in `[0,93]`, all
probability flags (`--gc`, `--f-damaged`, `--bg-ct`, `--bg-ga`, `--oxog`,
`--depur`, `--dmg5-max`, and any explicitly-set `--dmg3-*-max`) in `[0,1]`,
and all resolved lambda rates `>0` (:119-165).

## Input / output

No input file — sequence and damage are drawn from the seeded PRNG. Output
is a single FASTQ file (gzip level 6 if `-o` ends in `.gz`, plain text
otherwise); a one-line run summary (fraction damaged, library type, damage
parameters) is printed to stdout on completion (:266-273).
