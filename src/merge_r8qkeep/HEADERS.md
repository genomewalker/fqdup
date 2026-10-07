# Frozen header provenance

`native_headers/` contains the 25-header quoted-include closure from
`native-quality-v1/production-pair-quality-v5/include/`, falling back to the
terminal-anchor tree at commit `9d2997a` exactly as the frozen native build did.
`census_headers/` contains the seven-header closure from the frozen
`realistic-dev-staged/support/` census sources.

Header algorithms and original notices are preserved. Translation-unit-only
namespace macros isolate native and census models from each other and from
legacy merge. Include paths are source-specific, not target-wide. Scratch
artifacts were read only; these local snapshots remove their build-time include
dependency. `repair_headers/` is the separately gated repair snapshot.
