# Design brief (SCOPE CORRECTED): smart clip / trim / merge engine for ss aDNA

You are Fable 5, system architect. SCOPE = the **read-preprocessing engine only**: adapter/linker
detection, trimming, clipping, and pair MERGING in `fqdup merge` (src/merge.cpp). NOT the downstream
damage estimator. Design a **smart, fast, accurate, damage-preserving, decisively superior**
clip/trim/merge system. Deliverable = concrete buildable design written to
`FABLE_DESIGN_OUTPUT.md` in this dir. Do not edit code. Mark unverifiable claims as "assumption".

## Why this matters (the domain constraint that makes it non-generic)
This preprocessing feeds a reference-free ancient-DNA **deamination** estimator that reads the 5' (and
3') terminal base composition. The C->T damage signal lives in the FIRST FEW BASES of the molecule
terminus. Therefore the engine must be **damage-preserving**: it must remove adapter/linker/artifact
bases and NOTHING biological — a single over-trimmed or frame-shifted terminal base corrupts the
estimand. Generic trimmers (cutadapt/fastp/AdapterRemoval/leeHom) are tuned for mappability, not for
preserving a pristine, correctly-framed terminus on ultra-short (25-75nt) ancient fragments. "Superior"
= better than those FOR THIS: zero terminus corruption + maximal true-artifact removal + fast.

## Current implementation (this worktree, branch feat/unmerged-clip) — what works, what doesn't
fqdup merge collapses overlapping R1/R2 into merged reads (short inserts = ancient-enriched); non-
overlapping pairs go out as unmerged R1/R2 (long inserts = modern-enriched). Read src/merge.cpp:
- `detect_merge_params` (~:539): prescans a buffer, auto-detects adapter1 (R1 3' readthrough),
  adapter2 (R1 5' complement bleedthrough), poly-G, and NEW: `adapter_5p_linker` (top 5' 12-mer that
  is not a known adapter, freq>=5%).
- `clip_unmerged` (~:785): whole-read scan for earliest of {detected adapter, universal AGATCGGAAGAG}
  -> resize; then 5' strips of linker, adapter_5p, adapter_3p, universal; then 3' poly-G.
- `emit_unmerged` (~:810): clips both mates, QC, routes drop/orphan/pair.
- Merged path (~:1059): consensus_merge + trim_adapter_5p(adapter2) + trim_adapter(adapter1) + QC.

VERIFIED on clay-test Med11 (30.9M pairs): merged count identical pre/post (merged path untouched);
unmerged R1 5' adapter contamination 13% -> 0.27%; linker CGCAATGCTCAT auto-detected; 92% of unmerged
R1 now clean full-length 101bp.

KNOWN DEFECTS you must design fixes for:
1. **Linker length**: real library linker is 20bp `CGCAATGCTCATGGACTCAA` but detection captures only
   the first 12bp (fixed-width substr). Universal-adapter scan masks it for dimers, but a linker
   preceding a real insert would leave an 8bp `GGACTCAA` tail on the biological 5' = TERMINUS CORRUPTION.
   Need full-length linker learning (variable length, data-driven).
2. **Wrong adapter auto-detected**: detection locked onto Nextera `GGAAGAGCGTCG` while the true
   readthrough is TruSeq `AGATCGGAAGAGCACAC`. We patched with a hardcoded universal fallback; the
   DETECTION itself should robustly identify the correct adapter family per library.
3. **Indel-tolerant matching**: residual `GATCGGAAGAGC` (~0.37%) = adapter with a 1bp 5' indel, missed
   by the fixed-offset 2-mismatch scan. Need an alignment tolerant to a single indel without blowing up cost.
4. **R2 5' poly-G**: unmerged R2 5' is ~28% poly-G (2-colour dark-cycle no-signal at read START).
   trim_polybase only trims 3'. Need a 5' low-complexity/dark-cycle guard (and decide: clip vs drop).
5. **Orphan accounting**: after the fix, unmerged 8.35M->1.83M, orphans 43k->6.25M (one mate was pure
   adapter, correctly clipped to empty, good mate retained). Validate this is correct, define the
   contract for where orphans are written, and ensure no silent read loss.

## Speed dimension (PI directive, high priority)
"we need speed for derep, then derep and run profile on it"; "use the f32 with 2 GPU".
- The engine runs before/with dereplication on billions of reads. Design for throughput: single-pass,
  SIMD/GPU-amenable adapter matching, minimal per-read allocation.
- Where can detection + clipping + merge overlap-scan be vectorised (SIMD) or offloaded (GPU, f32)?
  What is inherently branchy/CPU? Give a concrete kernel/dataflow plan and expected throughput shape.
- Derep-first vs clip-first ordering: argue the correct order (clipping changes sequence identity, so
  it likely must precede exact-duplicate collapse — reason it out).

## Output-contract dimension (PI: "validate json structure, provide all we need in fqdup profile")
The estimator consumes fqdup's profile output. Part of this design: specify the JSON/emit contract the
merge+profile stage must produce so the estimator has EVERYTHING it needs and nothing ambiguous:
- Per-channel composition (merged comp5/comp3 [Lbin,pos,base]; unmerged R1 5'; unmerged R2 5').
- Provenance: detected adapters/linker, counts (merged/unmerged/orphan/dropped), min-length, per-read-length histograms.
- A schema you'd validate against. Note any field the estimator needs that the current profile likely omits.

## Deliverable sections (write to FABLE_DESIGN_OUTPUT.md)
1. Architecture of the smart clip/trim/merge engine: detection -> clip -> merge -> emit, damage-preserving invariants.
2. Data-driven adapter+linker LEARNING (variable-length, correct-family, per-library, calibration-free).
3. Matching algorithm: indel-tolerant, fast, SIMD/GPU plan (f32), the merge overlap-consensus.
4. The 5 known-defect fixes above, each with a concrete mechanism.
5. Speed design: single-pass dataflow, derep ordering argument, GPU/CPU split, throughput.
6. Output/JSON contract + schema the estimator needs; gaps in current fqdup profile.

Ground it in src/merge.cpp (read clip_unmerged, emit_unmerged, detect_merge_params, trim_adapter,
trim_adapter_5p, find_adapter_in, consensus_merge, passes_qc). You may call mcp__chitta__recall
(realm=ellesmere) for library-prep specifics. Mark assumptions explicitly.
