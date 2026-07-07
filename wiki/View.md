# fqdup view

Line references are to `src/view.cpp` at `fb62967`.

## Purpose

`fqdup view` inspects `.fqcl` cluster genealogy files — the binary format
written by the derep clustering stage (read via `cf::Reader` /
`fqdup/cluster_format.hpp`). It is a read-only debug/reporting tool: header
summary, single-cluster trees, cross-cluster bundling, and a FASTQ-to-cluster
lookup, in plain text, JSON, or a static HTML viz.

User-facing (listed as a subcommand in `fqdup --help`, `src/main.cpp:37`).

## Modes

`usage()` at `view.cpp:899`. Exactly one mode applies per invocation; default
is the header summary.

- **(no flags)** — header summary: `n_clusters` and the file's `metadata`
  JSON blob (`print_header` :117).
- **`--cluster N`** — ASCII tree of cluster `N`: walks the edge list from the
  root (parent) outward, printing each node's position, from→to base, and
  `n_reads` (subtree size) (`print_tree` :122).
- **`--staircase N`** — per-node mismatch grid for cluster `N`: one row per
  node, `.` where a base matches the parent, the differing base otherwise
  (`print_staircase` :154, edits applied via `apply_edits` :96).
- **`--bundle [--end-k K]`** — groups all clusters by a hash of their first
  and last K bases (`fqdup::bundlekey::from_cluster`, `bundle_key.hpp`),
  default `K=16`. Prints one row per bundle: key, cluster count, member
  cluster IDs (`print_bundles` :187).
- **`--bundle-staircase HEX`** — renders `--staircase` for every cluster
  sharing bundle key `HEX` (`print_bundle_staircase` :210).
- **`--min-bundle-size N`** (default `2`) — with `--bundle`, hides bundles
  with fewer than `N` member clusters.
- **`--top-members N`** / **`--top-edges N`** — top-`N` clusters by member
  count or edge count, via a bounded heap over all clusters (`print_top`
  :224).
- **`--dump-members`** — TSV `cluster_id\tmember_id`, one row per member,
  for every cluster.
- **`--member-of FASTQ`** — for every record in `FASTQ`, emits TSV
  `read_name\tcluster_id` (`cluster_id=-1` if unmapped). Matching is by exact
  canonical (strand-agnostic) sequence fingerprint (xxhash3-128 over the
  sequence and its reverse complement, `mo_canonical_fp` :783); it indexes
  every reconstructed node sequence of every cluster (`emit_member_of` :803).
  Damage-mask-merged twins whose exact base sequence differs from any indexed
  node report `-1` — this is an exact-match lookup, not damage-aware.
- **`--json`** — emit structured JSON instead of plain text. Envelope:
  `{"schema":"fqdup.view.v1","fqcl":{path,n_clusters,metadata},"query":{mode,args},"data":...}`
  (`emit_envelope_open` :76). Exception: `--dump-members --json` emits NDJSON
  (one cluster object per line, no envelope).
- **`--html PATH`** — writes a self-contained HTML visualization of the top
  50 clusters by member count to `PATH` (`emit_html_static` :538). Combine
  with `--top-members N` / `--top-edges N` to change the count and ranking.

## Input / output

Input: one `.fqcl` file (positional argument). Output: stdout (plain text,
JSON, or NDJSON depending on mode/flags) or the file at `--html PATH`.
