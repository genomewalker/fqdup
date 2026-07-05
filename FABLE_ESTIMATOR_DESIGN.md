# Fable design: fast + accurate reference-free ss aDNA deamination estimator

Author: Fable 5 (system architect). Scope: **the downstream damage estimator**
that consumes the fqdup-merge profile (see companion `FABLE_MERGE_ENGINE_DESIGN.md`
for the preprocessing engine). Design only, no source edits.

Grounding: `estimator_final.py` (locked rho estimator), `estimator.py`
(length-decomposition ceiling), `estimator2.py`, `compos.py`, `analyze_reffree.py`
(all in `clay_estimator/reffree/`); chitta realm=ellesmere (rho lock,
length-decomp ceiling result, first-moment non-identifiability, ss-3' death).
Unverifiable claims tagged `[ASSUMPTION]`.

---

## 0. The one hard truth this design is built around

From a **single terminal composition profile** the pair (pi, d_max) is
**not separately identifiable** — only the product enters. Code-proven and
triply-converged in memory (`[insight][FINAL]`): the ss terminal observable is
`E[Tfrac(i)] = baseline + pi·d_max·geom(beta,i)`; `(pi,d_max)` couple **only**
through the amplitude `A = pi·d_max`. No curve-fitting on one profile separates
them.

So the estimator's job is not to "fit harder". It is to **add a second
observable that carries pi and d_max with different weights**, and to state
precisely when that second observable identifies d_max — and when it does not,
to say NO-GO. This design supplies two axes (fragment **length**, primary;
within-read **co-occurrence**, research arm) and uses the now-clean
**unmerged-long modern channel** as the empirical undamaged template that both
anchors the length deconvolution and drives the false-positive guard.

**Deliverable in one line:** a per-library, calibration-free 5'-only estimator
that (a) reports `pi_peak·d_max` as an honest *lower bound* always, (b) lifts it
to point-identified `d_max` when the ancient and modern fragment-length
distributions are separable and the ancient component is above a self-measured
noise floor, and (c) issues a self-calibrated NO-GO otherwise — the correct
call at ultra-dilution.

---

## 1. System architecture + dataflow

```
 raw R1/R2.fq.gz  ──▶ fqdup merge (companion design) ──▶ profile.json + channels
                                                              │
   MERGED (short, both termini)   UNM_R1 (mol 5' long)   UNM_R2 (RC mol 3')
                    │                     │                     │
                    ▼                     ▼                     ▼
 ┌──────────────────────────────────────────────────────────────────┐
 │ DEREP-FIRST (exact) → GPU PROFILE PASS (f32, 2-GPU): one scan      │
 │   comp5[Lbin,pos,base], comp3[Lbin,pos,base], N[Lbin],            │
 │   modern templates umod5/umod3 from UNM channels                  │
 └───────┬──────────────────────────────────────────────────────────┘
         ▼  host reduce (720-cell histograms, trivial)
 ┌──────────────────────────────────────────────────────────────────┐
 │ CPU ESTIMATOR (tiny): per-Lbin rho, q_int, Delta_rho, exp fit     │
 │   → M(b)=pi(b)·d_max ; modern-anchored length deconvolution       │
 │   → d_max (or lower bound) ; 4-gate self-calibrated GO/NO-GO      │
 └──────────────────────────────────────────────────────────────────┘
```

Channel roles: **MERGED** = short-insert, ancient-enriched, core signal (5' for
ss; 5'+3' for ds — see companion `biological_termini`). **UNM_R1 5'** =
molecule-5' of long inserts (modern-enriched, pi≈0) = empirical undamaged
template + deconvolution anchor. **UNM_R2 5'** = RC(molecule 3'), complemented
into the C→T frame. For **ss** the 3' merged channel is a proven-dead dA-tail
prep artifact (chitta; task#13) → QC/abstain only. For **ds** the 3' G→A channel
is a valid second damage estimator and doubles the signal.

---

## 2. Accuracy method — the estimator math

### 2.1 Observable (locked)
Per length bin `b` (edges `[0,35,50,75,100,150,∞)`), per 5' position `i∈0..29`,
from merged unique reads:
```
fC(b,i)=comp5[b,i,C]/N,  fT(b,i)=comp5[b,i,T]/N,  rho(b,i)=fT/(fC+fT)
```
`rho` is immune to the universal 5' purine-fragmentation artifact (depurination
inflates fC,fT equally → cancels in the ratio; C→T raises rho). **⇒ pos0 usable
directly.** Per-library, per-bin, self-measured:
```
q_int(b)=Σ_{8..14}fC / Σ_{8..14}(fC+fT)     rho_int(b)=mean_{8..14}rho
sig_int(b)=std_{8..14}rho                    Delta_rho(b,i)=rho(b,i)-rho_int(b)
```
Locked decomposition: `Delta_rho(b,i)=pi(b)·d_max·q_int(b)·exp(-i/beta(b))`.
Fit `A(b)·exp(-i/beta(b))` to `Delta_rho(b,0..9)`; **per-length mixed amplitude**
```
M(b) := A(b)/q_int(b) = pi(b)·d_max
```

### 2.2 Tier 1 — peak estimator (always reported, honest lower bound)
`d_hat_lb = max_b M(b) = pi_peak·d_max ≤ d_max`. Equality iff a length bin is
ancient-pure. Measured recoveries ARE `pi_peak`: Med11 0.356/0.423=**84%**, Med25
0.131/0.178=**74%**. 198A recovers only ~20% because no 198A length bin exceeds
~20% ancient (overlapping ancient/modern lengths) — a library property, not
estimator weakness.

### 2.3 Tier 2 — modern-anchored length deconvolution (primary new method)
Identifying idea: **d_max is one scalar shared across all length bins; pi(b) is
not.** Two-component mixture:
```
N(b)       = Na·ga(b;θ) + Nm·gm(b)              (observed length hist)
M(b)·N(b)  = d_max · Na · ga(b;θ)               (observed damage-weighted counts)
```
- `ga(b;θ)` = ancient length shape, 1–2 param (truncated geometric). [ASSUMPTION:
  aDNA length ≈ geometric over 25–150nt; verify by fit residual.]
- `gm(b)` = modern short-fragment shape; **long tail pinned by UNM_R1** length
  distribution (modern extends long). [ASSUMPTION: modern short = smooth
  continuation of the observed modern-long shape.]

Both observed vectors share `Na·ga(b;θ)`. Joint NLS over `{θ,Na,Nm,P}` with
`P:=d_max·Na`:
```
min Σ_b [N(b)−Na·ga−Nm·gm]²/N(b) + Σ_b [M(b)N(b)−P·ga]²/var
then d_max = P / Na
```
`d_max` recovered **without any bin reaching pi=1**: the ratio of the shape-
matched ancient-damage curve to the ancient-count curve. This lifts 198A-class
libraries toward the true ceiling. Fallback to Tier-1 lower bound + `pi-limited`
flag when shapes don't separate.

### 2.4 Tier 3 — within-5' co-occurrence (research arm, gated OFF, not shipped)
A **second moment** separates pi from d_max: for two 5' C sites,
`Cov(i,j) ≈ pi(1−pi)·d_max²·e^{-(i+j)/beta}`, so
`Cov / [Delta_rho(i)Delta_rho(j)/q²] ≈ (1−pi)/pi` → pi (hence d_max) from one
length bin, no pure-ancient regime needed. **Status: UNVERIFIED for ss** (prior
test invalid — raw counts, task#14; ss-3' cross-terminal variant dead, task#13).
Ship only behind a **pre-registered ds-control gate**: reproduce known ds
co-occurrence separation first. Diagnostics-only until then; never overrides
Tier 1/2.

---

## 3. Identifiability analysis (adversarial, honest bounds)

**Theorem (first moment, proven).** A single 5' `rho` profile identifies only
`A=pi·d_max·q_int` and shape `beta`. `(pi,d_max)` non-identifiable from it alone.

**Length lifts it — conditionally.** `M(b)=pi(b)·d_max` is K equations in K+1
unknowns (K pi(b) + d_max): under-determined by one. Tier 2's shape prior
`pi(b)=Na·ga(b;θ)/N(b)` supplies it. Therefore:

> **d_max is point-identified iff ancient and modern fragment-length distributions
> are distinguishable AND the ancient component is large enough to estimate its
> shape θ from the mixture.**

**Quantitative floor (order-of-magnitude, from the truth table).**
`A_short ≈ pi_short·d_max·q_int` (q_int≈0.5); interior noise
`sig_int ≈ sqrt(0.25/N)`. GO requires `A_short > K·sig_int`, K=3:

| lib | pi_short | d_max | A_short | N | sig | A/sig | call |
|---|---|---|---|---|---|---|---|
| Med11 | ~0.9 | 0.42 | ~0.19 | 1e7 | 1.6e-4 | ~1200 | GO |
| 198A | ~0.2 | 0.315 | ~0.031 | 1e7 | 1.6e-4 | ~190 | GO; Tier-2 lifts 0.063→~0.3 |
| Med25 | ~0.75 | 0.178 | ~0.067 | 1e7 | 1.6e-4 | ~420 | GO |
| IS80 | small | 0.070 | ~0.005 | 1e7 | 1.6e-4 | ~30 | weak / floor |
| FLB45m | ~0.002 | 0.41(anc) | ~4e-4 | 1e7 | 1.6e-4 | ~2.5 | **NO-GO (<3σ)** |

FLB45m (ancient = 0.2% of mapped) sits at ~2.5σ, **below 3σ** → correct answer is
**NO-GO**: d_max is not identifiable from a 0.2% component statistically
indistinguishable from Poisson noise. The design proves the limit, not a guess.
(Only Tier-3's second moment could break this, and only past the ds-control
gate.)

### 3.1 Self-calibrated GO/NO-GO (4 gates; only free constant K=3, a σ-guard)
1. **Amplitude vs self-noise:** `A_short > K·sig_int`. Kills blanks/floor.
2. **Shape:** `Delta_rho(i)` monotone decreasing i=0→4, physical `0.5<beta<8`.
3. **Modern-anchor consistency:** UNM_R1 template `rho_mod(i)` must be **flat**
   (`A_mod < K·sig_mod`) and free of residual period-2. If the *modern* channel
   shows terminal decay, the "damage" is a shared artifact → abstain. (Single
   most valuable use of the newly-clean channel.)
4. **Length concentration:** damage-weighted counts `M(b)·N(b)` concentrate in
   short bins beyond the modern template's expectation. Length-flat → flag.

**Output contract:** GO+Tier-2 converges → `d_max` + bootstrap CI. GO+Tier-2
fails → Tier-1 lower bound + `pi-limited`. Any gate fails → NO-GO + which gate.

**ds note:** ds gives a second independent estimate from the 3' G→A channel
(same math, comp3). Agreement between 5' and 3' d_max is an extra ds-only
consistency gate; ss cannot use it (3' dead).

---

## 4. Speed design (derep-first, GPU f32, 2-GPU)

**Derep-first, exact.** Derep must be **exact full-length** (or 5'-anchored key),
never fuzzy: PCR copies of one ancient molecule are identical (deamination is
pre-amplification) → collapsing is correct and de-biases the tally; fuzzy derep
would merge molecules differing by a damage C→T and **erase signal**. Minimizers
only bucket; final equality is exact bytes.

**GPU kernels (f32 for ratios; counts stay int).** Kernel A: 2-bit pack + Lbin +
atomic histogram into shared-mem tiles `comp5/comp3[Lbin,30,4]` (1440 cells,
negligible traffic) + `N[Lbin]` + modern templates (R2 complemented in-kernel).
Kernel B: derep key + GPU radix sort (the real hot loop on billions of reads).
2-GPU: shard reads, host reduce. The 12-param exp fits stay **CPU** (log-linear
seed + one LM polish) — do NOT build a GPU fit kernel (pure overhead).

**What stays CPU:** gzip decode (libdeflate, N threads), the fits, the Tier-2
joint NLS, the 4-gate logic. Bottleneck is decompression + derep sort, not the
tiny histogram — fuse derep+tally into one scan. **Throughput** [ASSUMPTION,
hardware-dependent]: ~1–3 GB/s effective, billions of reads in tens of minutes.

---

## 5. Validation plan (validate, never fit)

Truth table (metaDMG mapped d_max): 198A 0.315, Med11 0.423, Med25 0.178, IS80
0.070; FLB45m floor ~0, bamdam ancient-subset ceiling 0.413 (0.2% of mapped).
1. **Ceiling recovery:** Tier-2 ≥ Tier-1 monotone on Med11/Med25/198A; CI covers
   metaDMG on GO libraries; 198A lifts from ~0.063 toward 0.315 without touching
   truth values.
2. **Pooling-independence:** per-library estimate must not change with batch
   composition (only K is global, a σ-guard) — any dependence is a bug.
3. **Blank safety (decisive):** FLB45m → NO-GO via Gate 1 (~2.5σ) and Gate 4;
   IS80 → NO-GO/weak lower bound. A false GO on FLB45m fails the design.
4. **Modern-anchor falsification:** UNM_R1 `rho_mod(i)` flat on all libs (Gate 3);
   inject pre-fix contaminated unmerged data → Gate 3 must flip to abstain
   (regression guard for the merge-engine clip fix).
5. **Tier-3 gate (if pursued):** reproduce ds co-occurrence pi-separation before
   trusting any ss second-moment number.
6. **ds cross-terminal:** on ds libs, 5' and 3' d_max estimates must agree.

---

## 6. Build plan — additions to fqdup/libtaph (priority order)

- **P0** — profile-in-derep-pass plumbing: length-binned `comp5/comp3`, `N[Lbin]`,
  modern `umod5/umod3` from UNM channels (port `compos.py`, incl. R2 complement).
  CPU first, correctness reference. (Depends on the merge-engine JSON contract.)
- **P1** — the CPU estimator: `rho/q_int/Delta_rho`, per-bin exp fit `M(b)`,
  Tier-1 peak + lower bound, 4-gate GO/NO-GO. Reproduces current best + honest
  FLB45m NO-GO — ship before Tier 2.
- **P2** — Tier-2 modern-anchored length deconvolution + bootstrap CI +
  `pi-limited` fallback. The 198A accuracy win.
- **P3** — GPU fuse (Kernel A/B, f32, 2-GPU); validate GPU histograms
  bit-identical to the P0 CPU reference.
- **P4** — Tier-3 co-occurrence (research, gated behind ds-control).
- **ds path** — enable the 3' G→A second estimate + 5'/3' agreement gate on ds
  libraries (reads `biological_termini` from the merge profile).

**Do NOT build:** a GPU curve-fit kernel, a fuzzy-derep path, or an ss 3' damage
estimator (proven dead).

---

### Assumptions ledger
- A1: aDNA length ≈ truncated geometric over 25–150nt (Tier-2 `ga`).
- A2: modern short shape = smooth continuation of UNM_R1 modern-long shape.
- A3: `pi_short` in §3 is order-of-magnitude, back-derived from recoveries +
  truth table, not independently measured per-lib.
- A4: throughput numbers hardware-dependent, unbenchmarked.
- A5: Tier-3 second-moment separation sound but **unverified for ss**; gated
  behind ds-control, never the primary number.
