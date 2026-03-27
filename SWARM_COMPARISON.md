# fqdup vs SWARM: detailed comparison

fqdup is designed to replace SWARM+vsearch in ancient DNA deduplication pipelines.
This document covers the algorithmic differences, failure modes of each tool, and
benchmark results on synthetic libraries with controlled ground truth.

---

## What SWARM does

SWARM clusters sequences by edit distance (d=1 by default) using breadth-first search.
vsearch first collapses exact duplicates, then SWARM links remaining unique sequences
whose edit distance is ≤ 1.

```
Input reads
    │
    ▼
vsearch --derep_fulllength        (exact dedup, different lengths never merge)
    │
    ▼
SWARM -d1 --fastidious            (BFS: link nodes within 1 edit, then merge
    │                              low-abundance clusters into large neighbours)
    ▼
OTU centroids
```

The critical property: SWARM uses **transitive closure**. If A and B are within d=1,
and B and C are within d=1, then A, B, and C all end up in the same cluster -
regardless of counts, regardless of whether A and C are related.

---

## What fqdup does

fqdup evaluates every potential merge **independently** with a count-ratio guard.
An edge A→B is only accepted if the minority cluster is small enough relative to
the majority that it is more likely a PCR error than a real variant.

```
Input reads (sorted)
    │
    ▼
Phase 1: exact canonical hash              (damage-masked XXH3_128 → clusters)
    │
    ▼
Phase 2: pigeonhole Hamming d≤2            (link candidate pairs)
    │
    ▼
Phase 3 (B3): per-edge ratio guard         (accept edge only if count(B)/count(A) < threshold)
    │                                       independent evaluation, no transitive closure
    ▼
Output clusters
```

---

## Failure mode 1: SWARM collapses real SNP alleles at any depth

The simplest case: two genuine molecules A and B differ by 1 substitution.
Both are at high coverage (200x each). They are real genetic variants.

```
         200x        100x
          A    ─d1─   B

SWARM BFS:
  A (large) ──absorbs──► B (also large, but d=1 away)
  Result: 1 cluster     ← WRONG, should be 2

fqdup ratio guard:
  edge A→B:  count(B)/count(A) = 100/200 = 0.50
  threshold: 0.20
  0.50 ≥ 0.20 → SNP veto: reject edge
  Result: 2 clusters    ← CORRECT
```

SWARM has no count-ratio test. A sequence within d=1 of any larger cluster is always
absorbed, regardless of its own abundance. High-depth minority alleles are silently lost.

**Benchmark (F1–F3, 2 SNP molecules):**

| Scenario | True | fqdup | SWARM |
|---|---|---|---|
| A=200x, B=100x (ratio 50%) | 2 | **2 (+0.0%)** | 1 (−50.0%) |
| A=200x, B=50x  (ratio 25%) | 2 | **2 (+0.0%)** | 1 (−50.0%) |
| A=100x, B=99x  (ratio 99%) | 2 | **2 (+0.0%)** | 1 (−50.0%) |

---

## Failure mode 2: SWARM transitive closure collapses multi-allele chains

Three genuine alleles A, B, C, each differs from its neighbour by 1 SNP.

```
     A ─d1─ B ─d1─ C

     All at 100x each.
```

SWARM BFS starts from the largest node, follows all d=1 edges, and expands:

```
SWARM BFS from A:
  visit B (d=1 from A) → add to cluster
  visit C (d=1 from B) → add to cluster   ← C never compared to A!
  Result: {A, B, C} = 1 cluster            ← WRONG

A vs C directly:  they may differ by 2 edits, but SWARM never checks.
The chain A─B─C is enough.
```

fqdup evaluates each directed edge independently:

```
fqdup edge evaluation:
  A→B:  count(B)/count(A) = 100/100 = 1.0 ≥ 0.20 → veto  (SNP)
  B→C:  count(C)/count(B) = 100/100 = 1.0 ≥ 0.20 → veto  (SNP)
  Result: {A}, {B}, {C} = 3 clusters                ← CORRECT
```

**Benchmark (F4–F6, 3–4 allele structures):**

```
F4: linear chain        F5: star              F6: 4-allele same position
   A─B─C                   B                     A─B
  200 100 100              /                      │╲
                        A─╌                       │ ╲
                          \                       C─D
                           C                   all pairs d=1
                        200 100 100             all 100x
```

| Scenario | True | fqdup (capture) | SWARM |
|---|---|---|---|
| F4: 3-allele chain A=200x,B=C=100x | 3 | **3 (+0.0%)** | 1 (−66.7%) |
| F5: 3-allele star  A=200x,B=C=100x | 3 | **3 (+0.0%)** | 1 (−66.7%) |
| F6: 4-allele, all pairs d=1, all 100x | 4 | **4 (+0.0%)** | 1 (−75.0%) |

---

## Failure mode 3: SWARM absorbs singletons and low-coverage alleles

In shotgun sequencing, genuine molecules appear only 1–2 times. SWARM's fastidious mode
is specifically designed to rescue low-abundance sequences, but it cannot distinguish
a PCR error singleton from a real singleton allele.

```
Singleton SNP scenario (shotgun):

    A (1x) ─d1─ B (1x)      ← two real molecules, each seen once

SWARM fastidious:
  B is low-abundance, A is low-abundance
  fastidious links them (both below --boundary 3)
  Result: 1 cluster          ← WRONG: collapses 2 real molecules

fqdup shotgun mode (snp_min_count = 1):
  edge A→B: count(B)/count(A) = 1/1 = 1.0 ≥ 0.20 → veto
  Result: 2 clusters         ← CORRECT
```

Multi-hop singleton chain (SSH5–SSH6):

```
SSH5: A─B─C, each 1x              SSH6: all 4 alleles differ pairwise, each 1x

SWARM BFS:                              A
  A absorbs B (d=1)                    /│\
  B absorbs C (d=1)                   B─╌─C
  {A,B,C} = 1 cluster ← WRONG         \│/
                                        D
fqdup shotgun:                   SWARM → 1 cluster  ← WRONG
  every edge: ratio = 1.0        fqdup → 4 clusters ← CORRECT
  every edge: vetoed
  {A},{B},{C} = 3 clusters ← CORRECT
```

**Benchmark (SSH scenarios):**

| Scenario | True | fqdup-shotgun | fqdup-capture | SWARM |
|---|---|---|---|---|
| SSH1: 2 SNPs × 1x | 2 | **2 (+0.0%)** | 1 (−50%) | 1 (−50%) |
| SSH2: 2 SNPs × 2x | 2 | **2 (+0.0%)** | 2 (+0.0%) | 1 (−50%) |
| SSH3: A=5x, B=1x | 2 | **2 (+0.0%)** | 1 (−50%) | 1 (−50%) |
| SSH4: A=10x, B=1x* | 2 | 1 (−50%) | 1 (−50%) | 1 (−50%) |
| SSH5: 3-chain × 1x | 3 | **3 (+0.0%)** | 1 (−66.7%) | 1 (−66.7%) |
| SSH6: 4-allele × 1x | 4 | **4 (+0.0%)** | 1 (−75%) | 1 (−75%) |

\* SSH4 is the resolution limit: a 10:1 depth ratio with a single minority read is
genuinely ambiguous, cannot distinguish real SNP from PCR error at any ratio threshold.
Both fqdup and SWARM fail here; this is a fundamental limit, not a bug.

---

## Typical capture library accuracy (C1–C7)

For standard PCR libraries without SNP alleles, both tools perform similarly.
fqdup capture mode is designed for high-PCR-depth libraries where singletons are errors.

```
Error %: (reported_unique - true_unique) / true_unique
Negative = over-merged (false positives)
Positive = under-merged (false negatives)
```

| Scenario | True | fqdup | SWARM |
|---|---|---|---|
| C1: 300 mol ×30x, 40bp, PCR | 300 | +0.7% | +0.7% |
| C2: 300 mol ×30x, 40bp, damage+PCR | 300 | +4.7% | +4.7% |
| C3: 500 mol ×10x, 65bp ±15bp, damage+PCR | 500 | +11.0% | **+6.8%** |
| C4: 300 mol ×30x, 65bp ±15bp, high damage | 300 | **+1.3%** | +1.7% |
| C5: 1000 mol ×5x, low PCR | 1000 | −0.9% | −0.9% |
| C6: 2 SNP mols ×200x | 2 | **+0.0%** | +0.0% |
| C7: 2 SNP mols ×5x | 2 | **+0.0%** | +0.0% |

**C3** is the only case where SWARM wins (+6.8% vs +11.0%). The difference comes from
variable-length reads (±15bp from adapter trimming): SWARM's d=1 edit distance chains
reads of different lengths, while fqdup requires exact sequence match. However, for
real ancient DNA, variable-length reads can also represent genuinely different molecules
at overlapping positions, both tools are making an assumption here. See the C3 note
below.

> **C3 caveat:** In the synthetic benchmark, all length variation is defined as adapter
> trimming (PCR copies of the same molecule). In real aDNA data, variable-length reads
> can also be genuinely different molecules that share a locus. fqdup's conservative
> behaviour (no length-variant merging by default) is the safer choice when molecule
> identity cannot be confirmed from sequence alone.

---

## Throughput

```
fqdup  = sort (4 threads) + derep --error-correct --damage-auto (capture mode)
SWARM  = seqkit | vsearch --derep_fulllength | swarm -d1 --fastidious --boundary 3

Dataset (65bp, damage+PCR)        reads    fqdup r/s   SWARM r/s   speedup
─────────────────────────────   ───────   ──────────   ─────────   ───────
1k mol × 100x  (high redundancy) 100,000    168,000     150,000      1.1×
5k mol × 30x   (moderate)        150,000    182,000      98,000      1.9×
10k mol × 20x  (high diversity)  200,000    177,000      82,000      2.2×
```

SWARM throughput degrades with diversity because the d=1 BFS scales with the number of
unique sequences. fqdup throughput is flat, the hash lookup is O(1) per read.

---

## Summary

| Property | fqdup | SWARM |
|---|---|---|
| SNP allele protection | Per-edge ratio guard, protects any allele above count threshold | None, any d=1 neighbour is absorbed by the larger cluster |
| Transitive closure | No, each edge evaluated independently | Yes, BFS chains through unlimited hops |
| Low-coverage alleles (shotgun) | Protected (shotgun mode: singleton ratio test) | Absorbed by fastidious mode |
| Damage awareness | Yes, C→T / G→A masked before hashing; library type auto-detected | No |
| Variable-length reads | Conservative, no merging across lengths | Chains via d=1 edit distance |
| Throughput | Flat with diversity, O(1) hash lookup | Degrades, BFS over unique sequences |
| Speed | 1.1–2.2× faster | Slower at high diversity |

fqdup is strictly better for ancient DNA metagenomics and population genomics where
preserving true genetic diversity matters. SWARM was designed for amplicon OTU clustering -
a different problem where the goal is to collapse sequencing noise rather than to
distinguish genuine variants.
