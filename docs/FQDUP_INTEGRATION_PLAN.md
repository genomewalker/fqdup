# fqdup integration + worktree cleanup plan

Target: land everything real into shared `fernandezguerra/apps/repos/fqdup`
(`integration` branch, currently `d86278a`). libtaph UNTOUCHED (frozen; refactor
is fqdup-only). All merges/deletes are WIRE-AND-HOLD — proposed here, executed
only on explicit approval. `main` is NEVER modified.

## Key finding
`feat/derep-phase3-split` (tip 9559638, 21 ahead of integration) is a *stack*
that already absorbed almost all outstanding work. Verified `git merge-base
--is-ancestor` — these branches are ALREADY CONTAINED in it (land for free when
it merges): merge-max-length, integrated, ss-shortfrag-dev, damage-summary-out,
profile-earlystop, profile-single-pass, short-fragment-base, merge-auto-entropy,
merge-complexity-merged, derep-phase3-streaming, derep-sortscan-discovery,
derep-pass1-fast, derep-drop-all-edges. (14 branches.)

13 more branches are ALREADY IN integration (ancestors of d86278a) → pure delete:
bulk-damage-design-b, clean-data-qc-main, lsd-accumulator-to-libtaph,
phase3-csr-payload, phase3-gcfilter, phase3-shardrelayout, purine-joint-classifier,
unmerged-clip, perf-finalize-threads, perf-parallel-ingest, phase3-csr,
wip/fqdup-join-preserve, integration-candidate.

## The only real decisions: 5 branches with unique unmerged commits
| branch | unique commit | verdict |
|---|---|---|
| feat/derep-entropy-cap-restore (=main tip 36d285b) | `derep: restore entropy-guarded bucket cap on streaming engine` | **LAND** — real feature; directly relevant to the bucket-cap path (C0 corpus exercises it). Cherry-pick from the feature branch, NOT main. |
| feat/phase3-hugepage (7317ff8) | `derep: HugeSlab (mmap+MADV_HUGEPAGE) for Phase-3 interior slab` | **LAND (evaluate)** — perf-only, should be byte-identical; verify against the 3 guard corpora. |
| feat/short-frag-floor-wip (4d57bbd) | `wip: --short-fragment-floor routing + CLI flag` | **EVALUATE** — WIP; may be superseded by the shipped `--min-length` short-floor. Confirm before landing. |
| feat/derep-win3-radix-candidates (4b1a044) | `WIN3 S1: drop FlatPairIndex directory, lower_bound on CSR keys` | **DEFER** — experimental Phase-3 index variant; conflicts with streaming path. Bench first. |
| t5-empirical-priors (1cbc9ba, 8 ahead/150 behind) | per-occ-bin absorption counters + walk/complexity counters | **ARCHIVE** — very stale; salvage individual ideas later, do not merge. |

WIP-preserve snapshots (discard after confirming nothing unique): derep-gc-prefilter
(0ed34de), multichannel-perread-scorer (a63a322), merge-indel-dimer-qc (9360da3,
its log-fix already landed as e9fa8d0 in the stack).

## Execution order
1. **Finish phase3 refactor** on `feat/derep-phase3-split` (C0 done; C1→C3 in progress),
   byte-identical vs all 3 guard corpora (`af3ac4…`, `f71967…`, `9a0bed…`) at -p 1/2/8.
2. **Merge `feat/derep-phase3-split` → `integration`** in shared fqdup. One merge lands
   the stack + the refactor. Rebuild + run the 3 guard corpora post-merge.
3. **Cherry-pick the LAND set** onto integration one at a time, guard-checking each:
   entropy-cap-restore (36d285b), then phase3-hugepage (7317ff8). Resolve short-frag-floor
   after confirming vs shipped `--min-length`.
4. **Delete** the 13 already-integrated branches + the 14 now-contained branches
   (after step 2) + their worktrees. Archive t5 + WIP snapshots as tags, then delete.
5. **Prune worktrees**: keep `fqdup_sf` (active), shared `fqdup` [integration],
   `fqdup_base`. `git worktree remove` the ~17 duplicates (4×a4a8eae, 2×8a87177,
   2×50f0f27, 2×36d285b, dead experiments).

## Guard (unchanged)
Every merge/cherry-pick into integration: rebuild, then md5 the 3 corpora at -p 1/2/8:
default `af3ac4e14f12240b1ce1af8a400169b1`, bucketcap `f71967f38008088b6c6f72befbd45701`,
kappa `9a0bed46a456e6ee0e0d6ca902569bd3`. A cherry-pick that changes derep output must be
justified (feature) or rejected (regression).
