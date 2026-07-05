# Fable design: library-type-aware clip / trim / merge engine (fqdup merge)

Author: Fable 5 (system architect). Scope: **read-preprocessing only** — adapter/
linker DETECTION, TRIMMING, CLIPPING, pair MERGING, EMIT/orphan routing in
`fqdup merge` (`src/merge.cpp`). **No estimator content.** Design only, no code
edits. Handles **both ss and ds** library preps; every section states its ss vs
ds behavior. Unverifiable claims tagged `[ASSUMPTION]`.

Grounding (line refs = `src/merge.cpp`, worktree `feat/unmerged-clip`):
`detect_merge_params` :541 (geometry/damage/linker :633-713, k-mer accum
:624-630), `find_adapter_in` :717, `consensus_merge` :342, `trim_adapter_5p`
:388, `trim_adapter` :402, `trim_polybase` :426 (3'-only), `find_adapter_tail`
:459, `passes_qc` :768, `clip_unmerged` :785, `emit_unmerged` :822, merged path
:1051-1070, Phase 0/0b/0c/1/2 :875-1043, `MergeOpts` :501, `DetectedParams` :524.
Memory realm=ellesmere: prep-table authority over filename/geometry, quality-blind
tally, ss/ds co-occurrence prep rules.

---

## 0. Thesis and the ss/ds axis

No aligner runs downstream — a reference-free deamination estimator reads the
**raw base at molecule position 0..~9**. So the governing invariant is:

> **Damage-preserving terminus fidelity:** strip every adapter/linker/artifact
> base and *not one biological base*, and never emit a fabricated terminus. One
> over-trimmed, frame-shifted, or cross-strand-"corrected" terminal base
> silently destroys the C→T estimand. Generic mergers (cutadapt/fastp/
> AdapterRemoval/leeHom) tolerate a fuzzy terminus because an aligner re-frames
> later; we cannot.

**The ss/ds split is not cosmetic — it changes which termini are biological:**

| | **ds prep** | **ss prep** |
|---|---|---|
| Library structure | insert flanked by TruSeq P5/P7, both strands | single-strand; extra 5' linker `CGCAATGCTCATGGACTCAA` + `CTCTTCCGATCT` complement bleedthrough |
| 5' terminus | biological, C→T | biological, C→T (behind the linker) |
| 3' terminus | **biological, G→A** — protect as carefully as 5' | **prep artifact (dA-tail)** — not biological signal |
| Termini to preserve | **both** | **5' only** |
| 3' handling | trim adapter only, keep biological bases | trim adapter **and** treat terminal A-run as prep, exclude from 3' signal |
| Geometry signal | R1[0:4] vs RC(R2) prefix **agrees** (ds) | prefix **disagrees** (`prefix_agree_rate<0.3` ⇒ ss, :638) |

The current code already computes `is_ss` from prefix agreement (:637-638) — this
design **elevates that flag to a first-class branch** through detection, clip,
merge, emit, and the JSON contract, and hardens the detection itself
(§1.1) because memory shows filename/geometry can both be wrong on a real
library (a named-`ss` lib was proven ds by 3' geometry + BIC tournament).

---

## 1. Architecture + damage-preserving invariants

### 1.1 Library-type detection (auto, with authoritative override)
```
PRESCAN (once/library, CPU):
  1. prefix-agreement geometry     : R1[0:4] vs RC(R2) insert prefix (:592-600)
                                      agree-rate >=0.7 → ds ; <0.3 → ss ; else ambiguous
  2. 3' terminus geometry cross-check: ds shows A-excess/G-depletion (G→A) at 3';
                                      ss shows a pure dA-tail (A-run independent of insert)
  3. 5' linker presence            : a dominant non-adapter 5' k-mer (>=5%, :709)
                                      is an ss signature (ds has no extra linker)
  4. RESOLVE: geometry + linker + 3' agree → set library_type with confidence.
             disagree → library_type=ambiguous, needs param/prep-table.
  OVERRIDE: --library-type {ss,ds,auto} ALWAYS wins over detection.
            record library_type_source = {declared, prep_table, detected}.
```
[ASSUMPTION grounded in memory: geometry alone misclassified a real ds lib named
"ss"; therefore an authoritative `--library-type` / prep-table value must
override, and the JSON must carry `library_type_source` so downstream never
trusts a bad auto-call. The three-signal vote (geometry + linker + 3') is more
robust than the single prefix-agreement test at :638.]

### 1.2 Stages (both types share the skeleton; branches marked ss/ds)
```
 PRESCAN → library_type + adapters + linker(ss) + poly-G + dark-cycle rate
   │
 per pair (stream, SIMD/GPU):
   5' GUARDS ─→ OVERLAP (Phase 0/0b/0c/1/2) ─merged→ CLIP(merged, type-aware) ─┐
      │                    │                                                   │
      │                    └─ no overlap → emit_unmerged (type-aware)          │
      ▼                                                                        │
  dark-cycle / linker-boundary flags                                          │
                                                                              QC + ROUTE
                                            ┌────────┬──────┬──────────┬────────┐
                                          merged  unm-pair orphan(R1|R2) dropped
                                                                              ▼
                              PROFILE TALLY (strand-native termini; §1.3 I3)
                                          ss: 5' only biological ; ds: 5'+3'
                                                                              ▼
                                       clipped FASTQ + profile.json (§6)
```

### 1.3 Invariants (upheld by every stage)
- **I1 — no biological trim.** A terminal base is removed only if it matches a
  confirmed construct at an anchored position. `trim_adapter_5p` (:388) anchors
  at pos 0, ≤1 mm — keep that; extend to the full learned linker (ss, §2.2).
  *ds:* applies to **both** ends. *ss:* applies to 5'; the 3' dA-tail is prep, so
  its removal is allowed (not biological).
- **I2 — never fabricate a terminus.** If post-clip pos 0 (or the ds 3' terminus)
  cannot be guaranteed the true molecule end (dark-cycle run, ambiguous linker
  boundary, indel-uncertain cut), flag `terminus_untrusted` and **exclude that
  end from its composition channel** — keep the read for length/derep, never
  clip-and-keep a fabricated terminus.
- **I3 — terminal bases strand-native, not consensus.** In a merged short insert
  both strands cover both termini. `consensus_merge` (:368-375) picks the
  higher-quality base on disagreement, so it *can* overwrite a damage T at the 5'
  with the undamaged partner's C. **Magnitude, honestly bounded:** the partner
  base at the merged **5'** comes from **R2's 3' end** — Illumina's lowest-quality
  region — so `q1 >= q2` usually holds and `consensus_merge` keeps R1's damaged
  base. So the pos-0 corruption is *small* (a minority of reads where q2 happens
  to exceed q1), and it **grows with position** deeper into the read as R2's
  quality recovers. The fix is still worth applying — it is correct and free — but
  the claim is "removes a position-dependent, mostly-interior bias", **not**
  "erases the pos-0 signal". Fix: within `skip_terminal` (:505/:656), the
  **composition tally** uses the base from the strand that natively reads that
  terminus:
  - 5' (both types): R1 top strand (reads C→T).
  - 3' **ds**: RC(R2) reads the bottom-strand 3' (G→A frame) — tally strand-
    native, **protect it**. Symmetric bound: the partner at the merged 3' comes
    from R1's 3' end (R1's own lowest-Q region), so `q2 >= q1` usually holds and
    RC(R2)'s native 3' base is kept anyway; corruption again grows toward the
    interior. Applying the strand-native tally removes the residual bias.
  - 3' **ss**: it is a dA-tail prep artifact — **not** tallied as biological.
  The emitted merged SEQ may remain consensus; only the tally is strand-native.
  *A correct, free cleanup that generic mappability-tuned mergers skip; invisible to a
  mappability-tuned tool.* [ASSUMPTION: validate strand-native vs consensus
  comp5/comp3 against clay-test mapped truth.]
- **I4 — length preserved.** Clipping removes only adapter / ss-linker /
  ss-dA-tail / low-complexity artifact; the surviving span is the biological
  insert. No hidden hard-clip (`clip_5p` :506 stays 0; if used, record it).

---

## 2. Data-driven detection & learning

### 2.1 Robust adapter-FAMILY detection (Defect 2) — both types
**Flaw:** the prescan finds overlaps with a *hardcoded* `TRUSEQ_R1` seed (:579),
so `a1_freq/a2_freq` (:625-629) only accumulate k-mers adjacent to
TruSeq-found overlaps, then `adapter1 = revcomp(adapter2)` blindly (:694).
A spurious Nextera k-mer (`GGAAGAGCGTCG`) beat the true TruSeq
(`AGATCGGAAGAGCACAC`) in the vote; band-aided with the hardcoded universal
fallback (:792).

**Fix — panel-anchored family vote** (adapter seqs are *protocol* constants,
explicitly PI-allowed):
1. Panel of read-through families as constants: TruSeq `AGATCGGAAGAGC`, Nextera
   `CTGTCTCTTATACACATCT`, MGI/BGI `AAGTCGGAGGCCAAGCGGTC`, small-RNA
   `TGGAATTCTCGGGTGCCAAGG`. [ASSUMPTION: covers libraries in play; extensible via
   `--adapter-fasta`, exists :513/:1133.]
2. For each family, count supported read-through hits at a valid insert boundary
   (indel-tolerant, §3.2). Score = hits.
3. Winner = max hits above a **statistical** support floor (Poisson over buffer
   size, not a damage constant). Tie/none → `adapter_family=ambiguous`, keep the
   universal fallback.
4. Refine the exact adapter from the aligned read-through **consensus** (captures
   library variants); derive `adapter1/adapter2` from the winning family's known
   P5/P7 geometry, not blind revcomp.

*ss vs ds:* same panel and vote. ds typically resolves to TruSeq on both mates;
ss additionally carries the 5' linker (§2.2), which the family vote must **not**
mistake for a read-through adapter (it is not in the panel; the is_known check
:703-707 already excludes adapters — extend it to exclude the learned linker).

### 2.2 Variable-length 5' linker learning (Defect 1) — **ss only**
**Flaw:** linker captured as a fixed 12-mer (`r1_5p_freq` :630, chosen :701-711).
Real construct is 20 bp `CGCAATGCTCATGGACTCAA`; the 8 bp tail `GGACTCAA` survives.
If a linker precedes a *real* insert, that tail sits on the biological 5' = I1
violation.

**Fix — greedy self-calibrated extension to the entropy boundary:**
1. Seed = dominant 5' 12-mer `S12` (≥5%, not a known adapter, :709).
2. Reads whose 5' matches `S12` ≤1 mm: walk `j=12,13,…`, compute majority-base
   fraction `m(j)`.
3. Extend while `m(j)` stays fixed (construct ≈0.95–1.0); **stop when `m(j)`
   drops to the library's own interior positional purity** `m_int` (positions
   15–25 of the same reads, the biological null ≈0.25–0.45). Crossing = linker→
   insert boundary. Self-calibrated; no magic constant. [ASSUMPTION: ligated
   linker is per-library near-invariant; must recover the 20 bp construct
   exactly on the known case.]
4. Learn full var-length linker + the **length distribution of what follows**:
   dimer-only (followed by adapter/zero insert) → clip whole, body drops below
   min-length; insert-leading → clip **exactly** to learned length, preserving
   insert pos 0.
5. **Clip the learned linker on the MERGED 5' too** — *current bug:* the merged
   path (:1060-1061) trims `adapter2` (5') and `adapter1` (3') but **never the
   linker**, so a merged linker-led insert keeps its tail. Add
   `trim_adapter_5p(merged, linker, …)` (full linker) before the adapter2 trim.

*ds:* **no linker stage** — ds has no extra 5' construct; skip entirely. Running
linker learning on ds would risk clipping a real 5' terminus, so it is gated on
`library_type==ss`.

---

## 3. Matching, merge consensus, SIMD/GPU

### 3.1 Overlap/merge cascade (keep, harden) — both types
Phase 0 (adapter-anchored short insert) → 0b (8/6 bp anchors) → 0c (1-5 bp tail)
→ 1 (quality-weighted Hamming, long insert) → 2 (tail-head) (:875-1043) is sound
and short-insert-first (the aDNA case). Hardenings: route terminal tally through
I3 (type-aware: ss 5'-only, ds 5'+3'); feed panel + indel-tolerant matcher into
every `find_adapter_in` call site.

### 3.2 Indel-tolerant anchor (Defect 3) — both types
**Flaw:** `find_adapter_in` (:717) is fixed-offset Hamming; residual
`GATCGGAAGAGC` = `AGATCGGAAGAGC` minus one 5' base (1 bp indel), missed (~0.37%).
**Fix — seed + ±1-shift confirm:** at each candidate anchor, also test the
adapter shifted −1 and +1 and take the min-mismatch of the three (one indel in a
12 bp anchor = 3 Hamming passes, O(3·12), invoked **only** on exact-anchor
misses → amortized ~0). Full banded (±2) DP only behind a flag for pathological
libraries. [ASSUMPTION: single-indel covers the 0.37%; verify residual <0.05% on
Med11.] Type-agnostic: same adapter chemistry both preps.

### 3.3 Merge consensus (leeHom-equivalent + I3) — type-aware termini
Keep the Bayesian consensus (:335-378): agree → `min(q1+q2,60)`; disagree →
higher-qual base, `q=|q1−q2|`. Only the **tally** within `skip_terminal` bypasses
consensus for the terminus-native strand (I3). ss: protect 5' only. ds: protect
5' **and** 3' (RC(R2) native, G→A frame). `skip_terminal` is already set
type-sensitively by the damage/UDG logic (:654-666); keep that, and additionally
force the ds 3' into the protected/tallied set.

### 3.4 SIMD plan (CPU hot path)
`pack_2bit` (:199) exists; overlap Hamming (`mismatch_prefix` :212) = XOR of
packed 2-bit vectors + popcount (8–16 bases/instr AVX2); `find_best_ov` (:261)
inner loop vectorized. Adapter anchor = broadcast 12 bp pattern,
`_mm_cmpeq_epi8`+`movemask`+popcount per position; ±1 indel = two shifted
compares. Quality-weighted mismatch (`min(q1,q2)/30.f`, :916/:944) already f32 →
SIMD min + FMA.

### 3.5 GPU plan (f32, 2-GPU) — what maps, what stays CPU
Per-pair independent → data-parallel across 2 GPUs (shard pairs, host reduce).
Maps: panel matching (pair × family), packed overlap score vector `mm[ov]`,
`consensus_merge` (per-position, coalesced). Run all phase scores on GPU,
**select the phase on the host** to avoid warp divergence from the branchy
cascade. Stays CPU: prescan (hash maps, linker greedy-extend, family vote —
once/library, trivial), gzip decode (libdeflate, N threads), variable-length
string surgery, orphan routing, phase selection. f32 is correct precision — all
scores are rates/qualities in [0,1]·small ints. **Throughput** [ASSUMPTION,
hardware-dependent]: I/O + derep-sort bound, not merge-compute bound; ~1–3 GB/s
aggregate at decompression line rate, GPU near-idle on the merge scan — its real
payoff is fusing the derep-key emit (§5), not the scan. Type-agnostic (same
kernels; the type branch is a cheap host-side flag).

---

## 4. The five known defects — fix table (with ss/ds behavior)

| # | Defect (loc) | Mechanism | ss | ds |
|---|---|---|---|---|
| 1 | Fixed-12 linker → 8 bp tail (:630,:701); merged path never clips linker (:1060) | Greedy-extend to interior-purity boundary; learn full var-length linker; clip on merged 5' too; dimer-only vs insert-leading by following-length (§2.2) | **Active** | **Skip** (no linker) |
| 2 | Wrong family (Nextera vs TruSeq); circular seed (:579,:694) | Panel family vote by supported hits; refine exact seq from aligned consensus; a1/a2 from geometry (§2.1) | Active (+exclude linker from vote) | Active |
| 3 | 1 bp indel missed (`GATCGGAAGAGC`, :717) | Seed + ±1-shift confirm, only on exact-anchor misses (§3.2) | Active | Active |
| 4 | `trim_polybase` 3'-only (:426); R2 5' ~28% dark-cycle poly-G | 5' dark-cycle guard: low-Q leading G-run ⇒ terminus lost ⇒ `terminus_untrusted`, exclude from comp, don't clip-and-keep (§4.1) | Active | Active |
| 5 | Orphans 43k→6.25M; mate-of-origin lost (:834) | Preserve which mate survived; route by mate; balance identity; emit all buckets (§4.2) | Active (orphan-R2 = ss dA-tail 3', low value) | Active (orphan-R2 = biological 3') |

### 4.1 Defect 4 — 5' dark-cycle guard (both types)
2-colour dark cycles emit `G` at near-floor quality; a biological G-rich 5' has
normal quality. **Discriminator:** leading `G`-run ≥ `min_run` **AND** mean run
quality < (interior mean quality − margin), margin from the read's own quality
distribution (self-calibrated). If both: true 5' overwritten, position unknowable
⇒ `terminus_untrusted`, exclude from the 5' composition (I2), count
`dropped_dark_cycle`; the body after the run may still feed length/derep if it
passes QC. **Exclude, never clip-to-fabricate.** Applied to R1 5' (rare) and R2
5' (common, ~28%). Note R2 5' = RC(molecule 3'): a lost R2 5' loses a molecule-3'
terminus — **material for ds** (3' is biological), **immaterial for ss** (3' is
the dA-tail artifact anyway).

### 4.2 Defect 5 — orphan accounting + balance (both types)
`emit_unmerged` (:822-838): both QC-pass → unmerged pair; one → orphan carrying
only `mr.merged` (mate identity lost); none → drop. The 6.25M orphan surge is
**expected/correct** — after clipping, dimer/adapter mates collapse below
min-length while the real-insert mate survives; that survivor is a genuine
terminus and must reach the estimator. Fixes:
- **Preserve mate-of-origin:** tag `ORPHAN_R1` (survivor = molecule 5', C→T
  frame) vs `ORPHAN_R2` (survivor = RC(molecule 3'), complement like the R2
  handling). Conflation corrupts the channel. *ds:* orphan-R2 is a **biological
  3'** terminus — valuable. *ss:* orphan-R2 is the **dA-tail** — retained for
  count/length balance but flagged non-biological for the 3' channel.
- **Balance identity (emitted + asserted, no silent loss):**
  `2·pairs_in = 2·merged + 2·unmerged_pairs + orphan_r1 + orphan_r2
  + dropped_qc + dropped_dark_cycle`. Emit `balanced: true/false`; false = hard
  bug.
- **Where written:** define explicit output streams — `*.merged.fq.gz`,
  `*.unmR1.fq.gz`, `*.unmR2.fq.gz`, `*.orphanR1.fq.gz`, `*.orphanR2.fq.gz` (the
  current run already emits merged/unmR1/unmR2; add the two orphan streams so
  orphans are never silently folded into unmerged or dropped).

---

## 5. Single-pass dataflow + derep ordering (both types)

### 5.1 Clip-before-derep (reasoned)
PCR/optical duplicates of one molecule have **different raw sequences** —
different adapter read-through lengths, read lengths, adapter errors — and become
**byte-identical only after clipping to the insert**. Therefore:
- Derep-before-clip is wrong: misses true duplicates (raw reads differ), can
  falsely collapse distinct molecules sharing a prefix.
- Clip-before-derep is correct: clip to insert, then **exact** full-length derep
  collapses true duplicates. Deamination is pre-amplification, so all copies
  carry the *same* C→T (and, for ds, the same 3' G→A) and collapse correctly.
  Derep must be **exact**, never fuzzy — fuzzy would merge molecules differing by
  a real C→T and erase signal.
- **"Derep first, then profile"** = first relative to *profiling*, not clipping.
  Total order: **detect → clip/merge → exact-derep → profile.**

Type note: the derep key must be **type-aware at the 3'** — ds uses the full
insert (both termini biological); ss should key on 5'+length and treat the
dA-tail cautiously (a variable-length dA-tail must be normalized before it enters
the key, or identical molecules with different tail lengths won't collapse).
[ASSUMPTION: ss dA-tail length varies; normalize/trim it pre-key.]

### 5.2 Single-pass fusion
Per read: detect-frozen params → clip/merge → emit clipped record **and its
exact derep key**. Derep = sort/hash on the key (GPU radix sort across 2 GPUs =
the real hot loop on billions of reads). Profile = tally each unique **once**,
strand-native and type-aware at termini. Minimizers only *bucket*; final equality
is exact bytes. Merge arithmetic + tally are cheap; the pass is I/O/sort-bound.

---

## 6. Output / JSON contract (library-type-aware) + current gaps

```json
{
  "schema_version": "1.0",
  "provenance": {
    "tool_version": "...", "input_sha256_r1": "...", "input_sha256_r2": "...",
    "min_length": 15, "min_ov": 11, "max_mm_rate": 0.08, "clip_5p": 0
  },
  "library": {
    "name": "Med11_ss",
    "library_type": "ss|ds|ambiguous",
    "library_type_source": "declared|prep_table|detected",
    "library_type_confidence": 0.98,
    "is_udg": false, "is_half_udg": false, "skip_terminal": 4,
    "biological_termini": ["5p"]          // ds → ["5p","3p"] ; ss → ["5p"]
  },
  "detection": {
    "adapter_family": "TruSeq|Nextera|MGI|ambiguous",
    "adapter1": "...", "adapter2": "...",
    "adapter_5p_linker": "CGCAATGCTCATGGACTCAA", "linker_len": 20,   // ss only; ds → null
    "linker_confidence": 0.98, "linker_mode": "dimer_only|insert_leading",
    "poly_g_detected": true, "dark_cycle_5p_rate_r2": 0.28,
    "damage_5p": 0.11, "damage_3p": 0.02, "prefix_agree_rate": 0.12,
    "da_tail_3p_rate": 0.31               // ss diagnostic; ds → low
  },
  "ledges": [0,35,50,75,100,150],
  "composition": {
    "note": "termini within skip_terminal are STRAND-NATIVE (I3); only biological_termini are valid",
    "merged_comp5":  "int[6][30][4]",     // both types
    "merged_comp3":  "int[6][30][4]",     // ds: biological G->A ; ss: dA-tail (flagged non-bio)
    "unmerged_r1_5p":"int[30][4]",        // molecule 5', C->T frame
    "unmerged_r2_5p":"int[30][4]",        // RC(molecule 3'); ds biological, ss dA-tail
    "orphan_r1_5p":  "int[30][4]",        // mate-preserved 5'
    "orphan_r2_5p":  "int[30][4]"         // mate-preserved RC 3'
  },
  "length_hist": { "merged_N": "int[6]", "unmerged_r1_N": "int[...]", "fine_hist": "int[...]" },
  "counts": {
    "pairs_in": 30900000, "reads_in": 61800000,
    "merged": 0, "unmerged_pairs": 0, "orphan_r1": 0, "orphan_r2": 0,
    "dropped_qc": 0, "dropped_dark_cycle": 0, "terminus_untrusted": 0,
    "balanced": true
  }
}
```
The `biological_termini` field is the key type-aware handshake: the estimator
reads only the termini listed there, so it applies ds-3' (G→A) logic or skips the
ss dA-tail automatically. **Validation:** ship a JSON Schema; refuse a profile
with `schema_version` mismatch, `balanced:false`, missing composition arrays, or
a `library_type` that contradicts `biological_termini`.

**Gaps in the current profile** (today: a pooled, quality-blind
`SampleDamageProfile` :812/:1053 + subst TSV/`.bsubst`). Add:
1. Length-binned per-channel `comp5/comp3` (`[6][30][4]`).
2. `unmerged_r1_5p`/`unmerged_r2_5p` composition (only clipped FASTQ emitted now).
3. Mate-preserved `orphan_r1_5p`/`orphan_r2_5p` (orphans lose mate identity :834).
4. Full var-length linker + `linker_mode` (today only 12-mer).
5. `adapter_family` + refined exact adapters (today one blind revcomp).
6. **`library_type` + `library_type_source` + `biological_termini`** — the
   type-aware contract; not emitted today.
7. `da_tail_3p_rate`, dark-cycle / `terminus_untrusted` counts.
8. Read-count balance identity (`balanced`).
9. Strand-native + quality-floor tally flags (memory: tally is quality-blind, so
   low-quality dark-cycle G pollutes it).

---

## 7. Build priority

- **P0 — library-type detection hardening + `library_type`/`biological_termini`
  in JSON.** Everything branches on this; also cheapest safety. (§1.1)
- **P1 — orphan mate-preservation + balance identity + orphan output streams.**
  (Defect 5)
- **P2 — ss linker learning + merged-5' linker clip** (gated `library_type==ss`).
  Highest terminus-corruption risk. (Defect 1)
- **P3 — 5' dark-cycle guard + `terminus_untrusted` exclusion** (both types;
  material for ds 3'). (Defect 4)
- **P4 — panel family detection + exact-seq refine.** (Defect 2)
- **P5 — indel-tolerant anchor (±1 shift).** (Defect 3)
- **P6 — per-channel comp5/comp3 + unmerged/orphan channels + I3 strand-native,
  type-aware tally + quality floor in JSON.**
- **P7 — SIMD overlap/anchor, then GPU fuse** (panel + scores + derep radix
  sort), validated bit-identical to the CPU reference first.

**Do NOT build:** fuzzy derep (erases C→T / ds G→A); default full-DP adapter
alignment (±1 shift suffices); GPU prescan or GPU phase-selection (divergence/
triviality); consensus-corrected terminal tallies (violates I3); ss linker
clipping on ds libraries (would corrupt a real 5'); treating the ds 3' as a
dA-tail (it is biological G→A).

---

### Assumptions ledger
- A1: I3 strand-native tally improves fidelity — validate vs clay-test mapped
  truth (5' both types; 3' ds).
- A2: ss linker is per-library near-invariant; greedy-extend recovers 20 bp.
- A3: adapter panel covers libraries in play (`--adapter-fasta` extensible).
- A4: single-indel model covers the 0.37% residual; ±2 DP only if not.
- A5: throughput is I/O/sort-bound, not merge-compute-bound; numbers
  hardware-dependent, unbenchmarked.
- A6: auto library-type detection can be wrong (memory: named-`ss` lib proven ds
  by geometry + BIC); therefore `--library-type`/prep-table overrides and the
  JSON carries `library_type_source`.
- A7: ss dA-tail length varies and must be normalized before the derep key or
  identical molecules won't collapse.
