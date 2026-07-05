# P0 — separate merged vs unmerged stream finalize (fqdup profile)

## Goal
In combined mode (`profile -i merged -1 unmR1 -2 unmR2`), the merged stream (short fragments,
damage-enriched) and the unmerged R1/R2 stream (LONG fragments, length-enriched for modern DNA)
are currently pooled into ONE `states` pool and finalized ONCE — destroying the short-vs-long
contrast that is the entire point. P0 finalizes the two streams SEPARATELY and emits their contrast.

## Verified hook points (src/profile.cpp, this worktree — re-confirm line numbers before editing)
- `states` pool declared ~:365
- PE (paired) workers bind to the same `states` ~:446  ← must rebind to a new `states_pe`
- single merge of per-thread states ~:496
- single finalize ~:568
- JSON input assembled ~:1116 (`ProfileJsonInput`), serialized in `profile_to_json` ~:1263
- ss 3′ marked non-biological ~:1249

## P0 tasks
1. Add `states_pe(n_threads)` beside the `states` declaration (~:365).
2. Rebind the paired/PE worker path (~:446) to write into `states_pe`, NOT `states`.
3. Build `dp_merged` from `states` and `dp_unmerged` from `states_pe`; finalize BOTH through the
   existing finalize path (~:568) — do not pool them.
4. Emit into the JSON a `streams` object: `streams.merged`, `streams.unmerged`, and
   `streams.length_contrast` with:
   - per-stream terminal amplitudes (d5/d3 raw) and per-stream π,
   - `delta = excess_merged - b_long` (b_long = unmerged terminal excess). **delta MUST be allowed
     to be zero or NEGATIVE** — `b_long` is the long-fragment NULL stratum, not "modern by
     definition": in an all-ancient core b_long>0 is real residual ancient damage, so also emit the
     RAW per-stream amplitudes alongside delta, never only the subtraction.

## Constraints (from room review — non-negotiable)
- Calibration-free: NO new hardcoded damage/threshold constants. No min-N constant for the unmerged
  stream — instead propagate binomial variance on b_long and, when its CI overlaps the merged
  interior baseline, emit `null_unresolved` (an abstention, itself constant-free). Wiring the CI is
  P1; for P0 just carry the per-position denominators through so P1 can compute it.
- Per-library only; never pooled.
- While in here, CONFIRM d5 and π are independently estimated (not d5 ≡ 0.39·π by a fit anchor).
  Empirically d_anc = d5/π ≈ 0.390 across 11/13 validation libs — if that constant is baked into the
  fit rather than emergent, flag it (do not silently "fix").

## Out of scope for P0 (P1/P2 — do NOT attempt now)
- d_max/π confidence intervals, ss-π rescue via R2-5′, joint length-stratified fit.
- NOTE for whoever does P1: "append unmerged as top length stratum before `estimate_damage_by_length`"
  does NOT compose — that call takes the merged `in_path`; unmerged lives in `states_pe`/`dp_unmerged`,
  unreachable via the reader. Inject `dp_unmerged` per-position counts as a synthetic top bin into
  `lsd_prebuilt`, or run the joint fit outside `estimate_damage_by_length`.

## Build & smoke test
- Build: `cmake --build /maps/projects/caeg/scratch/kbd606/tmp/fqdup_unmerged_clip/build --parallel`
- Smoke: run `profile` combined mode on one clay ds lib's merged+unmR1+unmR2 and confirm the JSON now
  carries `streams.{merged,unmerged,length_contrast}` with plausible per-stream numbers and that the
  legacy top-level fields are unchanged (regression guard).
- Worktree only. NOTHING to main.
