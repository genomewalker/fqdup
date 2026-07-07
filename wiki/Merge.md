# fqdup merge

Reconstructs full aDNA molecules from paired-end reads by overlap-merging R1/RC(R2),
cleaning adapter / library-construct / poly-G artifacts, QC-gating, and emitting a
merged stream (+ unmerged pairs, orphans). This is the input stage for reference-free
damage: the merged insert is the molecule the estimator profiles.

Line references are to `src/merge.cpp` on `feat/merge-complexity-merged` at `aa3d27f`.

## Pipeline (per pair, `merge_worker` :1382)

1. **Pre-scan auto-detection** (`detect_merge_params` :862) — first *n_scan* pairs.
2. **Overlap detection** — phases 0 / 0-extra / 0b / 0c / 1 / 2, in order.
3. **Consensus merge** (`consensus_merge` :350) — Bayesian per-base over the overlap.
4. **Trimming** — construct 5', adapter 5'/3', poly-G.
5. **QC gate** — different tests on merged vs unmerged (see below).
6. **Emit** — merged, or unmerged pair, or single orphan; zero-insert dimers dropped.

## Auto-detection (pre-scan, no user constants)

- **Adapters** — `adapter1` (P7_RC suffix past overlap on R1), `adapter2` (P5_RC on RC(R2)).
- **Library geometry** — 5' prefix agreement R1[0:8] vs RC(R2)[0:8]; `is_ss = prefix_agree_rate < 0.3` (:1063).
- **ss/ds construct-panel stem vote** (Ellesmere supp §4.3.2) — the SCR splint core
  `GGAAGAGCGTCG` (`SS_SPLINT_CORE` :898) reads through R1 3'. In ds it sits behind the
  TruSeq stem (`AGATC` upstream); in ss it does not. Decisive votes only (undeterminable
  / 2-mm-ambiguous reads abstain; symmetric tolerance avoids a ds→ss bias). Overrides
  geometry when positive (`type_from_panel` :1072-1075) — fixes SCR-splint ss libs
  miscalled as ds.
- **UDG / half-UDG status** and a **5'/3' damage pre-estimate** (`damage_5p` :1087).
- **5' linker** — a fixed library construct (adapter-dimer prefix) spikes in the R1 5'
  12-mer frequency; learned (`adapter_5p_linker` :1264) and clipped from unmerged reads.
- **Technical constructs** — every distinct construct with ≥ `ADAPT_MINN` = 100 supporting
  pairs is consensus-learned, up to `MAX_TECH` = 6 (:878-879, learn loop :1153-1186):
  P7/P5/splint/index in both orientations. The dominant one becomes `adapter1`; the full
  set (`tech_seqs`) drives merged-read technical QC.
- **Merged-mate complexity reference** (:901-905) — worst-window entropy + base-dominance
  of every overlap-verified insert, per mate. These trusted real fragments derive the
  unmerged low-complexity gate (entropy < 1st pct OR dominance > 99th pct; `derive` :1056);
  no hardcoded entropy/dominance constant. Below `COMPLEXITY_REF_FLOOR` = 500 inserts
  (:538) the gate stays disabled — a tiny/blank library cannot derive a degenerate gate.

## Overlap detection (phases, in order)

- **Phase 0** — adapter anchor in R1 (12bp seed, ≤2 mm), cross-validated against adapter2 in RC(R2) (`find_adapter_in` :1472).
- **Phase 0-extra** — extra adapters from `--adapter-fasta` or the internal panel (:1495), tried only when phase 0 finds no overlap.
- **Phase 0b** — progressively shorter anchors for near-read-length inserts (:1523).
- **Phase 0c** — adapter tail scan (`find_adapter_tail` :1552) for inserts 96-100bp (1-5bp of adapter visible).
- **Phase 1** — quality-weighted d=0 Hamming (`find_best_ov_qwt` :1582) for long inserts (insert ≥ read length).
- **Phase 2** — tail-head scan (`find_tail_head_ov` :1593) for inserts *longer* than read length.

## Trimming

- `trim_construct_5p` (:414) — learned variable-length 5' construct.
- `trim_adapter_5p` (:396) / `trim_adapter` (:431) — 5' / 3' adapter.
- `trim_polybase` (:455) + `derive_polyg_k` (:475) — 3' poly-G run-through (NovaSeq 2-color dark cycles).
- `clip_unmerged` (:1352) — full adapter + linker + poly-G clip on the unmerged mates.

## QC gate

`passes_qc` (:1312) checks `min_length`, `max_n_rate`, `min_entropy`, and the worst-window
low-complexity gate (entropy floor + base-dominance ceiling). Merged and unmerged reads
take different paths:

- **Unmerged mates** (:1403-1407) — full `passes_qc` (complexity gate derived from the
  merged-mate reference) plus two 5' adapter-dimer rejects: `is_adapter_dimer_5p` (:593,
  rigid Hamming frame) and `is_adapter_dimer_5p_indel` (:686, Myers ±2 anchor). The indel
  test catches zero-insert dimers carrying a single interior indel in the read-through
  adapter (the `AGATCGGAAGAG` G-run site) that the rigid frame misses.
- **Merged reads** (:1634-1635, :1672) — a merged read is `[insert][3'adapter]`, so its 5'
  is the insert. It is technical only if the *whole* read is adapter/construct: a zero-insert
  dimer that overlap-merged into pure adapter. `is_technical_read` (:735) tests this across
  all learned constructs via `is_adapter_fragment` (:651, Myers, indel-tolerant, exact-k-mer
  fast-reject) and drops it. The low-complexity gate is **not** applied to merged reads:
  it is derived from the merged inserts it would gate (circular — rejects ~1-2% of genuine
  inserts by construction), and overlap verification is itself the complexity check. Merged
  reads keep only the absolute floors, `passes_qc(mr.merged, opts, 0.f, 1.0f)`
  (min_length / max_n_rate / min_entropy). The 5'-anchor dimer checks are intentionally not
  run on merged reads — they are calibrated for unmerged mates and false-drop real inserts
  whose 5' resembles the adapter start (measured: 98% of their drops were real DNA).

## Adapter matcher

- **Hamming** `find_adapter_in` (:848, :1272) — anchored, ≤2 mismatch, seed 12, no indels.
  Used for all overlap/adapter placement.
- **Myers** `BitMatcher` (:614) — bit-parallel semi-global edit distance (Myers 1999),
  ~12% edit budget, indel-tolerant. Two uses: `is_adapter_fragment` (:651) drops a merged
  read whose whole sequence aligns as an adapter substring within budget;
  `is_adapter_dimer_5p_indel` (:686) rejects an unmerged mate whose 5' is a single-indel
  adapter dimer. This closes the 1-bp-indel dimer blind spot the rigid Hamming frame left.
  The k-budget gates the drop, so widening the seed cannot false-drop a real insert (which
  will not align full-length within ~12%).

The derep tree has a separate banded matcher (`derep_detail/banded_ed.hpp`) — a different
implementation, not the merge one.

## Multi-construct

Two mechanisms, both about libraries whose read-through goes beyond the standard Illumina
adapters:

- **Learned constructs (`tech_seqs`, default, always on)** — auto-detection learns every
  distinct construct with ≥100-pair support (see Auto-detection). Merged-read technical QC
  scans all of them in both orientations. No flag; nothing to configure.
- **Internal panel (`use_internal_panel` :811, default on; `--no-internal-panel` to disable)**
  — a built-in ss/ds-type-aware table of known aDNA construct read-throughs, added to the
  phase-0-extra adapter set (:1872) only when `--adapter-fasta` supplied none. It is a
  fixed top-up for constructs the per-library scan underlearns (e.g. an ss library where
  auto-detection only recovers the 12bp splint core, not the full 34bp universal). Validated
  3'-terminus-neutral before being defaulted on: on the 16-lib panel it changes 3' C→T by
  ≤0.02% and does not shave damage.

`--adapter-fasta` (odd lines = R1, even = R2; `load_adapter_fasta` :1745) supplies explicit
construct pairs and takes precedence over the internal panel.

## Output streams & optional profiles

- **merged** (primary), **unmerged pairs**, **orphans** (`--orphan-out`; one mate survived QC),
  **dropped** (zero-insert dimers, sub-min-length, QC-failed).
- `--damage-out` — paired damage profile JSON. `--subst-out` — overlap substitution matrix TSV.

## CLI flags

`--adapter1 --adapter2 --adapter-fasta --no-internal-panel --clip-r1-5p --min-overlap
--max-mm-rate --min-length --min-entropy --max-n-rate --r1-out --r2-out --orphan-out
--damage-out --subst-out`

`--internal-panel` is kept as a no-op for back-compat (the panel is now on by default).

## Commit state

- **389bea7** (= main): unmerged-mate QC gate only (low-complexity + adapter-dimer). Merged
  reads carry no technical drop. No Myers.
- **feat/merge-complexity-merged** (`aa3d27f`, three commits ahead of 389bea7): x86-64-v3
  build default (`b7e0c4c`), internal panel default on (`c412a48`), and the merged clean-data
  pipeline (`aa3d27f`) — Myers `is_adapter_fragment` merged-read drop, `is_adapter_dimer_5p_indel`
  on unmerged mates, the embedded-adapter path removed, and the non-circular merged gate.
  Validated fixed-vs-prefix on the 16-lib panel: 0 state-verdict flips, inert on the 15 real
  libraries (|Δestimate| ≤ 0.0125), and on the 11S-NC negative-control blank it only moves the
  spurious estimate toward truth. Committed, not yet pushed or merged to main.

Build with `-DFQDUP_MARCH=x86-64-v3` for cluster-portable binaries (`native` bakes in AVX-512
opcodes that SIGILL on non-AVX-512 nodes).
