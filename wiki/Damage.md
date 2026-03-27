# fqdup damage

## Purpose

`fqdup damage` profiles ancient DNA deamination damage and computes the
per-position mask that would be applied by `fqdup extend` and `fqdup derep
--collapse-damage`. It is a standalone diagnostic command, run it before the
full pipeline to inspect the damage profile, verify that library-type
auto-detection is correct, and decide whether `--collapse-damage` is warranted.

---

## Usage

```bash
fqdup damage -i FILE [options]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i FILE` | Input FASTQ (.gz or plain) | required |
| `-p N` | Worker threads | all cores |
| `--library-type auto\|ds\|ss` | Override library-type auto-detection | auto |
| `--mask-threshold FLOAT` | Mask positions where excess P(deam) > T | 0.05 |
| `--tsv FILE` | Write per-position frequency table as TSV | - |

---

## Algorithm

`fqdup damage` is a fully multi-threaded single-pass profiler. The producer
streams the FASTQ into batches of 8,192 sequences; N worker threads each
own their own `SampleDamageProfile` (no locking). After all reads are
processed, the per-thread profiles are merged and finalized in a single call.

The profiling and finalization logic is provided by DART's
[libdart-damage](https://github.com/genomewalker/libdart-damage). It measures
four biological channels per position:

- **5' C→T**, T/(T+C) at positions 0–14 from the 5' end
- **3' G→A decay**, A/(A+G) at positions 0–14 from the 3' end
- **3' G→A spike**, amplitude at 3' position 0 (ligation-site signal in SS libraries)
- **3' C→T**, T/(T+C) at 3' positions 0–14 (SS-original reads only)

Background rates are estimated from the middle third of reads (consistent with
DART and mapDamage2).

### Library-type detection

A 7-model BIC competition determines whether the library is double-stranded
(DS) or single-stranded (SS, complement-only or original-strand). The models
jointly evaluate which combination of channels explains the data better than
the null (no damage). The winning model is reported along with individual BIC
scores and channel amplitudes.

See [[Damage-Aware-Deduplication]] for a full description of the classifier.

### Mask computation

After finalization, positions where the observed excess deamination rate
exceeds `--mask-threshold` are flagged for masking. The mask is reported in
the human-readable output and, if `--tsv` is specified, in the per-position
table.

---

## Output

### Human-readable report

```
=== fqdup damage ===
Input:   merged.fq.gz
Threads: 16
Reads:   5,582,073 scanned
Length:  min=30  mean=91.2  max=150

Library: DS (auto-detected)
  BIC  bias=0.0  DS=125432.1  SS=42.0  SS_full=89.3
  fit  CT5_amp=0.1928  ΔBIC=125432.1  GA3_amp=0.0403  ΔBIC=1248.7  GA0_amp=0.0201  ΔBIC=22.3  CT3_amp=0.0012  ΔBIC=0.1
  5' terminal shift: +0.0193  (z=18.4)
  3' terminal shift: +0.0040  (z=3.9)

5'-end   d_max=0.1928  lambda=0.246  bg=0.4872
3'-end   d_max=0.0403  lambda=0.069  bg=0.5091
combined d_max=0.1928  (source=5prime) [validated]

Mask threshold: 0.05 → 1 position masked (pos 0)

pos  5'_CT   5'_GA   3'_GA
---  ------  ------  ------
  0  0.2194  0.4895  0.0521  *
  1  0.1571  0.4888  0.0458
  2  0.1048  0.4891  0.0432
  ...
```

Fields:

- **Library**, DS or SS, auto-detected or forced. Composition-bias warnings are printed if the read-end base composition makes library-type detection unreliable.
- **BIC**, BIC scores for the null model (bias), DS, SS, and SS-full (diagnostic). Higher means better fit.
- **fit**, Amplitude and ΔBIC for each of the four biological channels.
- **terminal shift / z-score**, Observed minus background T/(T+C) or A/(A+G) at the terminal positions; the z-score indicates significance.
- **d_max / lambda / bg**, Exponential decay parameters for each end, and the combined d_max (the larger of the two, used as the aggregate damage indicator).
- **\[validated\] / \[mixture\] / \[ARTIFACT\]**, Additional flags from the libdart-damage classifier. `[ARTIFACT]` signals the damage pattern is inconsistent with genuine ancient DNA.
- **Mask threshold**, Which positions exceed the threshold and would be masked in `fqdup extend` / `fqdup derep --collapse-damage`.
- **Per-position table**, Observed T/(T+C) and A/(A+G) frequencies at each of the first 15 positions; `*` marks masked positions.

### TSV output (`--tsv`)

```
pos	end5	freq5	end3	freq3	mask	cov5	cov3
0	5prime	0.219400	3prime	0.052100	1	5582073	5582073
1	5prime	0.157100	3prime	0.045800	0	5581842	5581842
...
```

`freq` values are −1 when coverage at that position is below 100 reads.

---

## Typical workflow

Run `fqdup damage` first to inspect the profile before committing to
`--collapse-damage`:

```bash
fqdup damage -i merged.fq.gz --tsv damage_profile.tsv
```

Then use the output to decide:

1. **d_max_5 < 0.02 and d_max_3 < 0.02**, no meaningful damage; skip
   `--collapse-damage` in `fqdup derep` entirely.

2. **Library type looks wrong**, override with `--library-type ds|ss` in
   `fqdup damage` first to confirm, then pass the same flag to `fqdup extend`
   and `fqdup derep`.

3. **mask-threshold produces too many / too few masked positions**, adjust
   `--mask-threshold` and re-run `fqdup damage` until the masked zone matches
   the damage plot. Then pass the same threshold to `fqdup extend` and
   `fqdup derep`.

4. **Supply manual parameters**, use the printed d_max and lambda values
   directly:

```bash
fqdup derep \
  -i nonext.deduped.fq.gz -o nonext.final.fq.gz \
  --damage-dmax5 0.193 --damage-lambda5 0.246 \
  --damage-dmax3 0.040 --damage-lambda3 0.069 \
  --mask-threshold 0.05
```

This skips Pass 0 in `fqdup derep` and uses the pre-computed parameters.
