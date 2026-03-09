#!/usr/bin/env bash
# test_extend_dmgcor.sh — Verify damage-correction improves extension through
# damaged terminals.
#
# Design:
#   Reads are random fragments from a shared reference genome (2000 bp).
#   Different fragment lengths overlap at each locus → longer reads provide
#   k-mers that extend beyond shorter reads' termini (the realistic aDNA model).
#
#   With dmax5=0.40: ~40% of reads have C→T at position 0.
#   At coverage ~25x and min_count=2, BOTH C and T edges appear in the graph
#   (both meet min_count threshold) → the walk sees n_supported==2 and stops.
#
#   Damage correction distinguishes damaged branches from real SNPs by checking
#   that the minor-allele fraction matches the expected deamination model.
#
# Expected:
#   --no-damage (no masking):         highest extension (upper bound)
#   --mask-5 4 --mask-3 4:            masking active, correction DISABLED
#                                     (manual mask → profile.enabled=false)
#   auto damage + correction:         correction ENABLED, should exceed masked-only
#
set -euo pipefail

FQDUP="${1:-/maps/projects/fernandezguerra/apps/repos/fqdup/build/fqdup}"
GEN="${2:-/maps/projects/fernandezguerra/apps/repos/fqdup/build/gen_genome_frags}"
TMPDIR_LOCAL="${TMPDIR:-/scratch/tmp}"
GENOME_LEN=10000
N_READS=2000     # ~11x coverage on 10kb genome for reads of mean length 55bp

SYNTH_NODMG="$TMPDIR_LOCAL/dmgcor_nodmg.fq"
SYNTH_DMG="$TMPDIR_LOCAL/dmgcor_dmg.fq"
EXT_NODMG="$TMPDIR_LOCAL/dmgcor_ext_nodmg.fq"
EXT_MASKED="$TMPDIR_LOCAL/dmgcor_ext_masked.fq"
EXT_CORR="$TMPDIR_LOCAL/dmgcor_ext_corr.fq"

# ── Generate genomic fragments ─────────────────────────────────────────────
echo "Generating fragments from ${GENOME_LEN}bp reference genome..."
"$GEN" \
    --genome-len "$GENOME_LEN" \
    --n-reads    "$N_READS"   \
    --min-len    30           \
    --max-len    80           \
    --dmax5      0.0          \
    --dmax3      0.0          \
    --seed       42           \
    > "$SYNTH_NODMG"

"$GEN" \
    --genome-len "$GENOME_LEN" \
    --n-reads    "$N_READS"   \
    --min-len    30           \
    --max-len    80           \
    --dmax5      0.40         \
    --dmax3      0.35         \
    --lambda5    0.35         \
    --lambda3    0.35         \
    --seed       42           \
    > "$SYNTH_DMG"

N=$(grep -c '^@' "$SYNTH_DMG")
echo "  Generated $N reads"

ext_rate() {
    local fqdup="$1" input="$2" out="$3" opts="$4"
    "$fqdup" extend -i "$input" -o "$out" $opts --min-count 2 2>&1 \
        | grep "extend pass2:.*extended="
}

# ── Run 1: no-damage input, no masking ────────────────────────────────────
echo ""
echo "=== Run 1: undamaged reads, --no-damage (upper bound) ==="
ext_rate "$FQDUP" "$SYNTH_NODMG" "$EXT_NODMG" "--no-damage"

# ── Run 2: damaged reads, --no-damage (masking disabled, graph has terminal kmers) ──
echo ""
echo "=== Run 2: damaged reads, --no-damage (no masking — shows damage impact on graph) ==="
ext_rate "$FQDUP" "$SYNTH_DMG" "$EXT_NODMG" "--no-damage"

# ── Run 3: damaged reads, manual mask, NO correction ──────────────────────
# Manual --mask-5/--mask-3 sets mask lengths but leaves profile_.enabled=false
# → correction never fires inside walk_right.
# mask-6 matches what auto-estimation produces at dmax=0.40, threshold=0.05
# (excess damage at pos 5: 0.40*exp(-1.75)≈0.069 > 0.05 → masked)
echo ""
echo "=== Run 3: damaged reads, --mask-5 6 --mask-3 6 (masking, NO correction) ==="
ext_rate "$FQDUP" "$SYNTH_DMG" "$EXT_MASKED" "--mask-5 6 --mask-3 6"

# ── Run 4: damaged reads, auto damage estimation + correction ─────────────
echo ""
echo "=== Run 4: damaged reads, auto damage estimation (masking + correction) ==="
"$FQDUP" extend -i "$SYNTH_DMG" -o "$EXT_CORR" \
    --damage-sample 0 --mask-threshold 0.05 --min-count 2 2>&1 \
    | grep -E "(d_max|Masked positions|mask length|extend pass2:.*extended=)"

echo ""
echo "=== Key comparison ==="
echo "  (higher % = more reads fingerprinted, better dedup)"
grep "extend pass2:.*extended=" "$TMPDIR_LOCAL"/dmgcor_ext_*.fq 2>/dev/null || true

# ── Cleanup ───────────────────────────────────────────────────────────────
rm -f "$SYNTH_NODMG" "$SYNTH_DMG" "$EXT_NODMG" "$EXT_MASKED" "$EXT_CORR"
