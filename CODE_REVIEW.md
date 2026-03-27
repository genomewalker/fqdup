# Code Review

Date: 2026-03-27

## Scope

Recheck of the current worktree against the previously reported correctness, math/validation, and performance findings.

## Validation

- `ctest --test-dir build --output-on-failure`
- Result: `6/6 tests passed`

## Findings

No confirmed open findings from the previous review remain in the current worktree.

The previously reported issues now appear fixed, including:

- truncated `.gz` FASTQ acceptance on the `rapidgzip` path
- unsigned CLI parsing that accepted negative literals via wraparound
- invalid `derep` float tuning values
- damage-profile "typical length" logging bias
- output paths that ignored final write or close failures

## Spot Checks

The earlier output-integrity repros now fail correctly instead of reporting success after data loss:

- `fqdup derep -i /tmp/trunc.out.fq -o /dev/full`
- `fqdup derep_pairs -n /tmp/trunc.out.fq -e /tmp/trunc.out.fq -o-non /dev/full -o-ext /tmp/ok.ext`
- `fqdup derep -i /tmp/trunc.out.fq -o /tmp/ok.out -c /dev/full`

These now exit non-zero and report write or close failures as expected.

## Recommended Next Checks

These are not confirmed bugs in the current worktree, but they are the next areas worth checking.

### 1. I/O Failure Injection

Exercise every command against:

- short writes
- disk-full targets such as `/dev/full`
- broken pipes
- unreadable inputs
- truncated gzip members
- temp-file rename failures

The recent fixes covered several final write and close paths, but this class of bug tends to recur in less-traveled output code.

### 2. Math and Model Invariants

Add property-style tests around `gen`, `damage`, and Phase 3 error correction to verify:

- probabilities remain in `[0, 1]`
- stricter thresholds behave monotonically
- masking never increases mismatch counts
- expected mismatch and tolerance summaries move sensibly with read length and damage rates

This is the best next step for math review beyond the current example-based tests.

### 3. Performance on Realistic Large Inputs

Profile these paths on large and adversarial datasets:

- `derep` Phase 3 on low-complexity and high-depth buckets
- `sort` on large compressed inputs with many chunks
- `rapidgzip` under multi-threaded callers

Measure wall time, peak RSS, and temp-disk usage, not just CPU time.

### 4. Cross-Backend Consistency

Run the same datasets through plain, `gzip`, `bgzf`, and `rapidgzip` paths and compare outputs.

The main remaining risk is not obvious crashes, but subtle behavioral drift between readers and writers.
