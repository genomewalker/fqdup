# Usage

## Standard workflow

```bash
# Step 1: Sort both files by read ID
fqdup sort -i nonext.fq.gz -o nonext.sorted.fq.gz --max-memory 64G --fast
fqdup sort -i ext.fq.gz    -o ext.sorted.fq.gz    --max-memory 64G --fast

# Step 2: Deduplicate
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz
```

`--fast` keeps intermediates uncompressed during sorting (3× faster, more disk).

## Ancient DNA workflow

```bash
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz \
  --damage-auto           # auto-estimate C→T / G→A damage
```

## With PCR error correction

```bash
fqdup derep \
  -n nonext.sorted.fq.gz \
  -e ext.sorted.fq.gz \
  -o-non out.non.fq.gz \
  -o-ext out.ext.fq.gz \
  --damage-auto \
  --error-correct         # remove 1-error PCR variants
```

## Cluster statistics

Add `-c clusters.tsv.gz` to write a gzipped TSV with per-cluster counts:

```bash
fqdup derep ... -c clusters.tsv.gz
```

## fqdup sort options

```
-i FILE            Input FASTQ (.gz or plain, or /dev/stdin)
-o FILE            Output FASTQ (.gz for compression)
--max-memory SIZE  Memory limit (e.g. 64G, 32G, 8G)
-p THREADS         Sorting threads (default: all cores)
-t DIR             Temp directory for chunks (default: .)
-N                 Natural sort order (numeric read-ID suffixes)
--fast             Uncompressed chunk intermediates (faster, more disk)
```

## fqdup derep options

```
-n FILE            Sorted non-extended FASTQ
-e FILE            Sorted extended FASTQ
-o-non FILE        Output non-extended FASTQ
-o-ext FILE        Output extended FASTQ
-c FILE            Cluster statistics TSV (gzipped)
--no-revcomp       Disable reverse-complement collapsing
--pigz             Parallel decompression (pigz)
--isal             Hardware-accelerated decompression (ISA-L)

Damage-aware:
--damage-auto          Auto-estimate damage parameters (Pass 0)
--damage-stride INT    Sampling stride for Pass 0 (default: 1000)
--damage-dmax   FLOAT  d_max for both ends
--damage-dmax5  FLOAT  d_max for 5' end
--damage-dmax3  FLOAT  d_max for 3' end
--damage-lambda FLOAT  Decay constant for both ends
--mask-threshold FLOAT Mask positions with P(deam) > T (default: 0.05)
--pcr-cycles INT       PCR cycles (informs damage model)
--pcr-error-rate FLOAT Sub/base/doubling: Q5=5.3e-7 (default), Taq=1.5e-4

Error correction:
--error-correct        Enable Phase 3 PCR error removal
--errcor-ratio  FLOAT  count(parent)/count(child) threshold (default: 50)
--errcor-max-count INT Absorb clusters with count ≤ N (default: 5)
--errcor-bucket-cap INT Pair-key bucket cap (default: 64)
```
