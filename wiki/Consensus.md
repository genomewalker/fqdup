# fqdup consensus

Line references are to `src/consensus.cpp` at `fb62967` unless noted.

## Status: disabled, not built

`fqdup consensus` is currently a stub. `src/main.cpp:21-24` defines
`consensus_main` as:

```cpp
// consensus_main: TEMP DISABLED — WIP refactor (task #20). Stubbed to keep linker happy.
static int consensus_main(int, char**) {
    std::cerr << "fqdup consensus: temporarily disabled (WIP refactor, task #20)\n";
    return 2;
}
```

`src/consensus.cpp` — the real implementation — is excluded from the build
(`CMakeLists.txt:139`, commented out: `# src/consensus.cpp  # TEMP DISABLED:
WIP refactor — task #20`). Its own commit is tagged WIP (`155af34`: "...
consensus/syncmer WIP").

As checked out, `consensus.cpp` would not compile if re-added to the build:
`consensus_main` calls `terminal_load(r.damage.term_5, lambda_5)` (:540,
:596) and `damage_qc(r.damage.term_5, r.damage.term_3, bg, damage_enabled)`
(:604), but no `terminal_load` or `damage_qc` with those signatures is
defined anywhere in the file — only `terminal_load_from_parent` (:199,
`std::string` parent + lambda + end + base) and `damage_qc_bin` (:347,
takes `InteriorComp[8]` arrays, not the raw `DamageTrack` `uint8_t[8]`
arrays passed at the call sites) exist. Neither of those is called from
`consensus_main`.

## Intended design (from the file's header comment, not a working build)

Per the design note at `consensus.cpp:1-20`:

- Alignment-free per-cluster consensus from a `.fqcl` v2 genealogy.
- Lineage-collapsed Bayesian posterior per parent position: every unique
  node in the cluster's edge tree is one lineage/vote; symmetric Dir(1)
  prior + Multinomial likelihood gives `p_top = (count_top+1)/(N_lineages+4)`,
  `Q = clamp(-10·log10(1-p_top), 0, 40)`.
- Damage-zone positions (route 2) emit the parent base with a clamped Q and
  set a cluster-level `damage_compatible` flag, rather than folding damage
  into the posterior.
- Background base composition for damage QC is aggregated per log-length
  (LSD) bin from parent-only interior positions, excluding jackpot
  (>P99-size) clusters.
- Per-cluster `damage_excess_z` = `(terminal_load - bin_mean) / bin_sigma_robust`,
  a signed, bias-aware terminal-damage-load z-score, not a probability.
- QC state `damage_model ∈ {pass, warn, fail}` from monotonicity + T-excess/
  C-depletion consistency + sample support; a separate `ligation_bias_state
  ∈ {none, suspect, strong}`.
- Intended outputs: `<prefix>.fq[.gz]`, `<prefix>.clusters.jsonl`,
  `<prefix>.summary.json`.

## CLI flags (as written in the disabled source)

`usage()` :443-462:

| Flag | Effect | Default |
|---|---|---|
| `-i, --in PATH` | Input `.fqcl` file | required |
| `-o, --out PREFIX` | Output prefix | required |
| `--min-lineages N` | Skip clusters with fewer lineages than `N` | 1 |
| `--no-gzip` | Write plain FASTQ instead of gzip | gzip on |
| `--threads N` | bgzf threads for gzip output | 1 |
| `-h, --help` | Print usage and exit | — |

This flag table describes source that is not part of the current build; it
does not reflect a runnable command.
