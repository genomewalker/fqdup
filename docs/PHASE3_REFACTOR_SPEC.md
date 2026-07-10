# Phase 3 (`phase3_error_correct`) modular refactor — execution spec

Source: fusion room `derep_phase3_refactor` (fable + sol), grounded in verified
source facts. Goal: fully modular/professional decomposition of the 2088-line
`DerepEngine::phase3_error_correct` **while producing byte-identical output**.

## Hard guard (non-negotiable)
- derep output on the 3M-read slice (`/maps/projects/caeg/scratch/kbd606/tmp/derep_refactor/slice.fq.gz`,
  md5 `164ce4b3…`) must stay md5 **`af3ac4e14f12240b1ce1af8a400169b1`** at `-p 1/2/8`.
- Two distinct assertions every commit: (a) md5 == baseline at each thread count;
  (b) cross-thread `cmp` byte-identical.
- Never touch the deferred-debt comparator `derep.cpp:770-771` (key-only, non-total;
  an id tie-breaker changes the drop set → changes md5 → behavioral).

## Codegen-safety rules
- Stay in **derep.cpp's translation unit**: all new files are `.hpp`, `#include`d at the
  current splice site. No `.cpp`, no cross-TU move, no LTO, no `always_inline`/`flatten`
  (each is its own codegen change — sol, accepted).
- `Phase3Runner` re-parents the body: alias every touched private at top of `run()`
  (`auto& index_ = e_.index_;` …) so the body stays **textually byte-identical**.
- B2 stays **monolithic + explicitly `inline`**; do NOT factor κ arithmetic or LR/SNP-veto
  into helpers. Do not replace `[&]` worker captures. Preserve every declaration, loop,
  sort, comparator, container type, reserve point, capture, thread partition, per-thread
  merge order, and float expression order.

## Architecture (2 files, cohesive not shattered)
- `src/derep_detail/phase3_workspace.hpp` — `struct Phase3Workspace`, cross-stage locals only:
  - identity/count: `N, id_count, id_fwd_count, seq_entry, acc_count`
  - occupancy: normal/rescue bundle keys + occupancy maps
  - routing/index: `layout_cache, build_map, short_parents, shards`, interior slab/offsets,
    ordered IDs, length segments, bucket histogram
  - B1 output: `per_thread_mm`, per-thread stats, merged `mismatches`
  - ascent/model: `parent_chain`, empirical model, κ counts/gate, joint-adj counters
  - shared stats/genealogy accumulators
  - (stage-only atomics, worker lambdas, clocks, scratch stay method-local)
- `src/derep_detail/phase3_runner.hpp` — `class Phase3Runner { DerepEngine& e_; Phase3Workspace w_; }`,
  all 8 stage methods inline in-class in execution order:
  `explicit Phase3Runner(DerepEngine&); void run();` + private
  `setup_occupancy(), build_phase_a(), collect_b1(), apply_b2_directed_ascent(),
   merge_b3(), probe_adjacent_lengths(), rescue_indels(), commit()`.
  No-param (read/write `e_`,`w_`) EXCEPT mutated accumulators (`errcor_absorbed`, `loss_*`)
  passed as explicit `uint64_t&` out-params so the write is visible to review.
- `DerepEngine`: forward-declare `Phase3Runner` before the class, add `friend class Phase3Runner;`,
  include `phase3_runner.hpp` AFTER the complete engine definition, thin wrapper
  `void phase3_error_correct(){ Phase3Runner(*this).run(); }`.

## Verified coarse-stage map (original derep.cpp line ranges)
setup/occupancy 535-681 · PhaseA index+routing+short_parents+build_shard 682-797 ·
B1 collect (∥) 798-1260 · B2 ascent+κ+empirical+H2-veto 1261-1813 ·
B3 H>2 merge 1814-2024 · T5.6 adjacent-len 2025-2175 ·
T8.7 indel-rescue (DEFAULT OFF) 2176-2578 · commit acc_count+genealogy 2580-2623.

## Commit sequence (each a guarded md5 checkpoint)
- **C0 (MANDATORY FIRST): characterization corpora.** Build inputs that activate the two
  dormant channels — one with `bucket_overflow_drops>0`, one with interior `κ≥2.0` — and
  capture their own baseline md5s. The 3M guard slice is blind to PhaseA-sort/bucket and
  B2-κ regressions (drops=0, κ=0.091). Without C0, PhaseA/B2 extraction is refactoring blind.
- **C1:** forward-decl + `friend` + `Phase3Runner` scaffold + move body into `run()` via the
  alias block, body textually unchanged (collapses room's C1+C2). Guard must pass trivially;
  if red, wrapper leaked a scope/capture → revert, migrate aliases individually.
- **C2:** hoist function-locals into `Phase3Workspace` one field-group at a time.
- **C3…:** carve one coarse stage per commit. With C0 corpora landed, PhaseA/B2 order stops
  mattering; else dormant-last (commit, indel-rescue, adjacent-len first; B2 and PhaseA last).
- **Final:** delete `phase3_error_correct.inc`.
- Fail handling: stage extraction reddens → revert that extraction, keep block inline in
  `run()`, localize via `.o` disasm diff. NEVER rebaseline or "repair" the algorithm.

## Status
- Step 0 (guard characterization): DONE — baseline current, thread-stable, codegen-robust for guard corpus.
- Step 1 (.inc quarantine): DONE + verified, commit 9559638, derep.cpp 3634→1547.
- C0 corpora: TODO (next).
