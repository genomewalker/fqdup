# Reconciliation: integration ↔ feat/derep-phase3-split

fable+sol fusion (room fqdup_merge_reconcile), grounded in real diffs, with my own
verification of the two load-bearing "need to verify" items. WIRE-AND-HOLD: local only.

## Verdict
Trunk = **feat/derep-phase3-split** (carries the streaming stack + C1 refactor that must
survive). Do NOT naive `git merge integration` — integration edits the pre-refactor
phase3 body that phase3-split moved into the `.inc`, so a 3-way merge yields modify/delete
conflicts across the whole Phase-3 region. Instead **replay integration's intent as
commits on a reconciliation branch, then record ancestry with `git merge --no-ff -s ours
integration`** (sol) — the reconciled tree is already authoritative, so `-s ours` avoids
meaningless conflict-resolution and can't silently clobber C1 (fable's concern).

## --min-length / --max-length union
- `--min-length` REQUIRED, `>= 1` (integration): CLI errors if unset. **Adjustment
  (commit 2, executed):** dropped the `-1` engine sentinel — kept `min_length_ = 0`
  default and enforce required-ness at the CLI (`if (min_length < 1) error`). The
  sentinel bought nothing once the CLI gate is the single source of truth, and 0 is
  the honest "no lower cap" identity for the gate `Lrec < min_length_`. Keep
  phase3-split's `set_length_filter(mn,mx)` superset; drop integration's `set_min_length`.
- `--max-length` optional, 0 = disabled (phase3-split only).
- **Accounting bug both caught:** integration's `n_input_reads = total_reads_ +
  n_below_min_length_` UNDERCOUNTS once `--max-length` fires. Fix: separate
  `n_above_max_length_` counter; `n_input_reads = indexed + below + above`. Keep both logs.
- Add CLI validation `max_length == 0 || max_length >= min_length` (no diff has it — new).
- Gate advances `record_idx` in BOTH ingest paths (serial `derep.cpp:353`, parallel
  post-join merge `derep.cpp:497`) + emit pass. VERIFIED: the parallel gate is in the
  single-threaded post-`join()` merge loop → counter increment needs NO atomic.

## Short-interior EC policy (the true conflict)
Adopt integration's **unconditional short-read EC-skip**, folded INTO the `.inc` (not by
reverting C1): remove `kBruteforceMinLen` + `short_parents` indexing/logging + the
brute-force child branch; add the single `ec_eligible(id)` predicate (interior >=
kMinInteriorLen) and apply it at the two probes integration guards (adjacent-length +
rescue). Keep `Phase3Runner` + the alias block.

## Dead members — VERIFIED schema-facing, so DO NOT delete
`short_interior_skipped`, `short_brute_evaluated`, `short_brute_found` are serialized in
the `.fqcl` JSON (`cluster_format.cpp:327-329`). Keep them **zeroed** in the merge
(integration already assigns the stats to 0). Removing them is a SEPARATE, schema-versioned
decision — not part of this merge. `loss_short_interior_skipped_` now counts RETAINED reads
(kept, not dropped) — its name lies; rename only under a schema bump, else keep zeroed.

## Validation (old md5 guards die BY DESIGN — output legitimately changes)
Replace byte-identical guard with:
1. **Determinism invariant** (the thing the guard actually protected): identical output at
   -p 1/2/8 for each new corpus. This still holds and is the acceptance gate.
2. **Structural invariant** (sol, stronger than fable's byte-diff): NO Phase-3 correction
   edge involves an interior-<20 read; eligible-only (interior>=20) correction behavior
   unchanged. Expect equal-OR-MORE representatives (not strictly more); do NOT require
   long-read output bytes untouched (short-cluster changes can alter emitted metadata).
3. Boundary tests: lengths min-1 / min / max / max+1; pass1/pass2 record_idx alignment;
   metadata conservation (indexed+below+above == physical input).
4. metaDMG/bamdam = frozen downstream diagnostic ONLY. Never tune thresholds/expected
   outputs to improve agreement.

## The real blocker (fable, terra-verified): required-flag turns CI red
~9 test scripts + README omit `--min-length`. The CLI guard reddens CI the instant it
lands. Test-script + doc updates MUST be in the SAME commit as the guard (atomic) —
otherwise steps between are un-guardable. **Executed (commit 2):** 9 test scripts taken
wholesale from integration (`git checkout integration -- <9>`, verified pure migration:
every changed line = old derep call → `--min-length 1`, behavior-preserving since no
test read < 1 bp). README taken from integration (strict superset — documents required
`--min-length` across all subcommands; neither side documented `--max-length`), then one
`--max-length` synopsis line added for the phase3-split-only flag. `test_merge_complexity.sh`
has NO derep call (its `--min-length 25` is `fqdup merge`'s own flag) — nothing to migrate.

## Commit sequence (each independently guardable)
1. Record phase3-split guard baseline (no code change) — already have af3ac4/f71967/9a0bed.
2. **Atomic** (EXECUTED, branch feat/derep-reconcile, off C1 tip 4f5b824): required-min +
   optional-max semantics + both counters (n_below_min_length_/n_above_max_length_) +
   physical-input accounting (`n_input_reads = total_reads_ + below + above`) + `max>=min`
   validation + 9 test-script migrations + README + usage text. Build clean; required-flag
   and `max<min` CLI errors verified firing. Guard: `--min-length 1` guard variant across
   default/bucketcap/kappa × -p 1/2/8 — must stay byte-identical (corpus floor = 16 bp, so
   `--min-length 1` drops nothing). [guard_c2.sh: **C2_GUARD_ALL_OK**, 9/9 byte-identical.]
3. `.inc` EC-policy replacement (EXECUTED). Folded integration's unconditional short-read
   EC-skip into `phase3_error_correct.inc` (7 surgical edits): deleted kBruteforceMinLen,
   short_parents index + its log, the ~55-line child brute-force branch; short-interior
   children now just `short_interior_skipped++; continue;`. Kept Phase3Runner + alias block.
   short_brute_evaluated/found/too_small_skipped kept declared/merged/persisted as **0**
   (schema-facing, cluster_format.cpp:327-329); loss_short_interior_skipped_ now = kept-count.
   Structural invariant holds BY CONSTRUCTION (edge-generating short path deleted; grep CLEAN).
   Guard C3_ALL_DETERMINISTIC — identical md5 across -p 1/2/8, all corpora CHANGED vs old
   (EC-skip fired), +243 representatives each (strict superset, 0 lost). NEW baselines:
   default e6635f57, bucketcap 97924cc1, kappa 5d257c56 (guard_c3.sh).
4. Cherry-pick wiki docs commit `8078031` (no code collision).
5. `git merge --no-ff -s ours integration` — record ancestry after replay ledger audit.
6. (Deferred, separate) schema-versioned removal/rename of dead short_brute_* fields.
