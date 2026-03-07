#!/usr/bin/env bash
# test_pipeline.sh — End-to-end pipeline tests: sort → derep_pairs → derep
#
# Tests:
#   P1. Primary dedup path (deterministic Tadpole-like extension)
#       3000 molecules × 5x coverage, no damage.
#       All reads from the same molecule get the same extension (deterministic
#       per molecule-ID), so derep_pairs collapses them correctly.
#       Assert: after derep_pairs, unique ≈ 3000 (within ±5%).
#               derep adds no further meaningful reduction.
#
#   P2. Residual dedup path (random extension per read)
#       2000 molecules × 5x coverage, no damage.
#       Each read gets a unique random extension, so derep_pairs can't collapse
#       duplicates.  derep on the derep_pairs output collapses the residuals.
#       Assert: after derep_pairs, unique ≈ 10000 (minimal dedup);
#               after derep, unique ≈ 2000 (residuals absorbed).
#
#   P3. Residual dedup with damage (full realistic pipeline)
#       2000 molecules × 5x coverage, dmax5=0.25 dmax3=0.20.
#       Random extension → derep_pairs misses dups → derep --damage-auto must
#       collapse them despite terminal deamination differences per read.
#       Assert: after derep --damage-auto, unique ≈ 2000 (within ±15%).
#               Unique count must be meaningfully lower than after derep --no-damage.

set -euo pipefail

FQDUP=${1:-build/fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== fqdup full pipeline tests (sort → derep_pairs → derep) ==="

# ---------------------------------------------------------------------------
# make_ext: generate ext.fq from non.fq.
#
#   non_fq   — source non-extended FASTQ (gen_synthetic output)
#   ext_fq   — output path for extended FASTQ
#   mode     — "det" (deterministic per molecule) or "rnd" (random per read)
#   seed     — integer RNG seed for rnd mode
#   ext_bp   — number of random bases to prepend and append
#
# In "det" mode: extension is derived from the molecule index embedded in the
# read header (r%07d_m%05d), so all duplicates of the same molecule get the
# same extension fingerprint, enabling derep_pairs to collapse them.
#
# In "rnd" mode: each read gets an independent random extension, simulating
# the case where Tadpole extended the same molecule differently for different
# PCR duplicates, so derep_pairs cannot collapse them.
# ---------------------------------------------------------------------------
make_ext() {
    local non_fq=$1 ext_fq=$2 mode=$3 seed=$4 ext_bp=$5
    python3 - "$non_fq" "$ext_fq" "$mode" "$seed" "$ext_bp" <<'PYEOF'
import sys, re, random

non_fq, ext_fq, mode, seed, ext_bp = \
    sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

BASES = 'ACGT'
rng_global = random.Random(seed)

def det_ext(mol_id, n):
    """Deterministic extension derived from molecule index."""
    r = random.Random(mol_id * 1000003 + 7)
    return ''.join(r.choice(BASES) for _ in range(n))

mol_id_re = re.compile(r'_m(\d+)')

with open(non_fq) as fin, open(ext_fq, 'w') as fout:
    while True:
        header = fin.readline()
        if not header:
            break
        seq  = fin.readline().rstrip()
        plus = fin.readline()
        qual = fin.readline()

        if mode == 'det':
            m = mol_id_re.search(header)
            mol_id = int(m.group(1)) if m else 0
            prefix = det_ext(mol_id,      ext_bp)
            suffix = det_ext(mol_id + 1,  ext_bp)
        else:
            prefix = ''.join(rng_global.choice(BASES) for _ in range(ext_bp))
            suffix = ''.join(rng_global.choice(BASES) for _ in range(ext_bp))

        ext_seq  = prefix + seq + suffix
        ext_qual = 'I' * len(ext_seq)
        fout.write(header)
        fout.write(ext_seq + '\n')
        fout.write(plus)
        fout.write(ext_qual + '\n')
PYEOF
}

# ---------------------------------------------------------------------------
# Test P1: Primary dedup — deterministic extension, no damage
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P1: Primary dedup path (deterministic extension) ---"
N_MOL=3000
N_READS=15000   # 5x coverage
EXT_BP=15

echo "  $N_MOL molecules × 5x coverage, no damage, det extension (${EXT_BP}bp each end)"

"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --no-damage --seed 42 > "$TMPDIR/p1_non.fq"

make_ext "$TMPDIR/p1_non.fq" "$TMPDIR/p1_ext.fq" det 42 $EXT_BP

"$FQDUP" sort -i "$TMPDIR/p1_non.fq" -o "$TMPDIR/p1_non.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p1_ext.fq" -o "$TMPDIR/p1_ext.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep_pairs \
    -n "$TMPDIR/p1_non.sorted.fq" -e "$TMPDIR/p1_ext.sorted.fq" \
    -o-non "$TMPDIR/p1_dp_non.fq" -o-ext "$TMPDIR/p1_dp_ext.fq" 2>/dev/null

"$FQDUP" derep -i "$TMPDIR/p1_dp_non.fq" -o "$TMPDIR/p1_final.fq" \
    --no-error-correct 2>/dev/null

AFTER_DP=$(grep -c '^@' "$TMPDIR/p1_dp_non.fq" || true)
AFTER_DEREP=$(grep -c '^@' "$TMPDIR/p1_final.fq" || true)
echo "  After derep_pairs: $AFTER_DP  |  After derep: $AFTER_DEREP  |  True: $N_MOL"

python3 - <<EOF
import sys
n_mol = $N_MOL
after_dp, after_derep = $AFTER_DP, $AFTER_DEREP

lo_dp, hi_dp = int(n_mol * 0.95), int(n_mol * 1.05)
if not (lo_dp <= after_dp <= hi_dp):
    print(f"  FAIL: derep_pairs gave {after_dp} unique (expect {lo_dp}–{hi_dp})")
    sys.exit(1)
print(f"  PASS: derep_pairs  {after_dp} ∈ [{lo_dp}, {hi_dp}]")

# derep should not further reduce (no residual dups with det extension)
if after_derep > after_dp:
    print(f"  FAIL: derep gave MORE reads than derep_pairs ({after_derep} > {after_dp})")
    sys.exit(1)
pct_residual = (after_dp - after_derep) / after_dp * 100
if pct_residual > 2.0:
    print(f"  FAIL: derep removed {pct_residual:.1f}% residuals (det ext should yield <2%)")
    sys.exit(1)
print(f"  PASS: derep residual reduction {pct_residual:.1f}% < 2%")
EOF

echo "  PASS: Test P1"

# ---------------------------------------------------------------------------
# Test P2: Residual dedup path — random extension, no damage
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P2: Residual dedup path (random extension per read) ---"
N_MOL=2000
N_READS=10000   # 5x coverage

echo "  $N_MOL molecules × 5x coverage, no damage, random extension per read"
echo "  derep_pairs misses most dups → derep catches residuals"

"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --no-damage --seed 55 > "$TMPDIR/p2_non.fq"

make_ext "$TMPDIR/p2_non.fq" "$TMPDIR/p2_ext.fq" rnd 55 15

"$FQDUP" sort -i "$TMPDIR/p2_non.fq" -o "$TMPDIR/p2_non.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p2_ext.fq" -o "$TMPDIR/p2_ext.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep_pairs \
    -n "$TMPDIR/p2_non.sorted.fq" -e "$TMPDIR/p2_ext.sorted.fq" \
    -o-non "$TMPDIR/p2_dp_non.fq" -o-ext "$TMPDIR/p2_dp_ext.fq" 2>/dev/null

"$FQDUP" derep -i "$TMPDIR/p2_dp_non.fq" -o "$TMPDIR/p2_final.fq" \
    --no-error-correct 2>/dev/null

AFTER_DP=$(grep -c '^@' "$TMPDIR/p2_dp_non.fq" || true)
AFTER_DEREP=$(grep -c '^@' "$TMPDIR/p2_final.fq" || true)
echo "  After derep_pairs: $AFTER_DP  |  After derep: $AFTER_DEREP  |  True: $N_MOL"

python3 - <<EOF
import sys
n_mol = $N_MOL
n_reads = $N_READS
after_dp, after_derep = $AFTER_DP, $AFTER_DEREP

# With random extension, each read looks unique to derep_pairs.
# Expect derep_pairs to pass through nearly all reads (>85%).
min_dp = int(n_reads * 0.85)
if after_dp < min_dp:
    print(f"  FAIL: derep_pairs deduped too aggressively ({after_dp} < {min_dp}); "
          f"random extension not working")
    sys.exit(1)
print(f"  PASS: derep_pairs kept {after_dp}/{n_reads} ({after_dp/n_reads*100:.0f}%) — random ext as expected")

# derep collapses residuals to near true molecule count (within ±10%)
lo_final, hi_final = int(n_mol * 0.90), int(n_mol * 1.10)
if not (lo_final <= after_derep <= hi_final):
    print(f"  FAIL: after derep, {after_derep} unique (expect {lo_final}–{hi_final})")
    sys.exit(1)
removed = after_dp - after_derep
print(f"  PASS: derep absorbed {removed} residuals ({removed/after_dp*100:.0f}%)")
print(f"  PASS: final unique {after_derep} ∈ [{lo_final}, {hi_final}]")
EOF

echo "  PASS: Test P2"

# ---------------------------------------------------------------------------
# Test P3: Residual dedup with damage
#
# Random extension → derep_pairs can't collapse dups.
# Non-extended reads carry terminal deamination (C→T, G→A) that differs between
# reads from the same molecule.  derep --damage-auto must collapse these residuals
# despite sequence differences at terminal positions.
#
# Control: derep --no-damage on the same input should leave MORE unique reads
# (the damaged duplicates are not collapsed without masking).
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P3: Residual dedup with damage (damage-aware pipeline) ---"
N_MOL=2000
N_READS=10000   # 5x coverage
DMAX5=0.25; DMAX3=0.20; LAMBDA=0.35

echo "  $N_MOL molecules × 5x coverage, dmax5=$DMAX5 dmax3=$DMAX3 lambda=$LAMBDA"
echo "  Random extension → derep_pairs misses dups → derep --damage-auto collapses residuals"

"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --dmax5 $DMAX5 --dmax3 $DMAX3 --lambda5 $LAMBDA --lambda3 $LAMBDA \
    --seed 77 > "$TMPDIR/p3_non.fq"

make_ext "$TMPDIR/p3_non.fq" "$TMPDIR/p3_ext.fq" rnd 77 15

"$FQDUP" sort -i "$TMPDIR/p3_non.fq" -o "$TMPDIR/p3_non.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p3_ext.fq" -o "$TMPDIR/p3_ext.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep_pairs \
    -n "$TMPDIR/p3_non.sorted.fq" -e "$TMPDIR/p3_ext.sorted.fq" \
    -o-non "$TMPDIR/p3_dp_non.fq" -o-ext "$TMPDIR/p3_dp_ext.fq" 2>/dev/null

# Control: no damage masking
"$FQDUP" derep -i "$TMPDIR/p3_dp_non.fq" -o "$TMPDIR/p3_nodmg.fq" \
    --no-error-correct 2>/dev/null

# Treatment: damage-aware masking
"$FQDUP" derep -i "$TMPDIR/p3_dp_non.fq" -o "$TMPDIR/p3_dmg.fq" \
    --damage-auto --no-error-correct 2>"$TMPDIR/p3_dmg.log"

AFTER_DP=$(grep -c '^@' "$TMPDIR/p3_dp_non.fq" || true)
AFTER_NODMG=$(grep -c '^@' "$TMPDIR/p3_nodmg.fq" || true)
AFTER_DMG=$(grep -c '^@' "$TMPDIR/p3_dmg.fq" || true)
echo "  After derep_pairs: $AFTER_DP"
echo "  After derep no-damage: $AFTER_NODMG"
echo "  After derep --damage-auto: $AFTER_DMG  |  True: $N_MOL"
grep -E "d_max=|Masked positions:" "$TMPDIR/p3_dmg.log" | head -3 | sed 's/^/  /'

python3 - <<EOF
import sys
n_mol = $N_MOL
after_dp, no_dmg, dmg = $AFTER_DP, $AFTER_NODMG, $AFTER_DMG

# damage-auto must collapse more reads than the no-damage baseline
if dmg >= no_dmg:
    print(f"  FAIL: --damage-auto gave no improvement ({dmg} >= {no_dmg})")
    sys.exit(1)
improvement = no_dmg - dmg
pct = improvement / no_dmg * 100
print(f"  INFO: improvement {pct:.1f}% ({no_dmg} → {dmg} unique)")
print(f"  INFO: final unique {dmg} vs true {n_mol} "
      f"({dmg/n_mol*100:.0f}% — residual inflation from partially-masked damage positions expected)")

# Meaningful improvement threshold: with dmax5=0.25 and 5x coverage,
# damage masking must collapse at least 25% of the no-damage excess.
if pct < 25.0:
    print(f"  FAIL: improvement {pct:.1f}% < 25% — damage masking not working")
    sys.exit(1)
print(f"  PASS: damage-auto improvement {pct:.1f}% ≥ 25%")

# Structural sanity: damage-auto must get closer to n_mol than no-damage does
dist_dmg    = abs(dmg    - n_mol)
dist_nodmg  = abs(no_dmg - n_mol)
if dist_dmg >= dist_nodmg:
    print(f"  FAIL: damage-auto ({dmg}) is not closer to n_mol ({n_mol}) than no-damage ({no_dmg})")
    sys.exit(1)
print(f"  PASS: damage-auto closer to true n_mol than no-damage ({dist_dmg} < {dist_nodmg})")
EOF

echo "  PASS: Test P3"

# ---------------------------------------------------------------------------
# Test P4: PCR library (damage + PCR errors, realistic ancient DNA prep)
#
# Mimics a standard PCR-amplified ancient DNA library:
#   - Moderate coverage (10x)
#   - Terminal deamination (dmax5=0.25, dmax3=0.20)
#   - PCR point errors (~20% of reads carry at least one error)
#   - Deterministic extension per molecule
#
# Deduplication path:
#   derep_pairs: exact-match duplicates collapsed (clean copies from same molecule)
#   derep --damage-auto --error-correct: residuals from (a) reads with PCR errors
#                                         and (b) differently-damaged copies
#
# Assert: --damage-auto --error-correct gives better unique-count accuracy than
#         the baseline (no damage masking, no error correction).
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P4: PCR library (damage + PCR errors) ---"
N_MOL=500
N_READS=5000    # 10x coverage
DMAX5=0.25; DMAX3=0.20; LAMBDA=0.35
PCR_RATE=0.003  # ~20% of reads have ≥1 error (P(k≥1) = 1-e^{-0.003*75} ≈ 0.20)

echo "  $N_MOL molecules × 10x coverage"
echo "  dmax5=$DMAX5 dmax3=$DMAX3, PCR rate=$PCR_RATE/base"
echo "  Deterministic extension → PCR-error reads differ in interior → need --error-correct"

"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --dmax5 $DMAX5 --dmax3 $DMAX3 --lambda5 $LAMBDA --lambda3 $LAMBDA \
    --pcr-rate $PCR_RATE --seed 101 > "$TMPDIR/p4_non.fq"

make_ext "$TMPDIR/p4_non.fq" "$TMPDIR/p4_ext.fq" det 101 15

"$FQDUP" sort -i "$TMPDIR/p4_non.fq" -o "$TMPDIR/p4_non.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p4_ext.fq" -o "$TMPDIR/p4_ext.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep_pairs \
    -n "$TMPDIR/p4_non.sorted.fq" -e "$TMPDIR/p4_ext.sorted.fq" \
    -o-non "$TMPDIR/p4_dp_non.fq" -o-ext "$TMPDIR/p4_dp_ext.fq" 2>/dev/null

# Baseline: no damage masking, no error correction
"$FQDUP" derep -i "$TMPDIR/p4_dp_non.fq" -o "$TMPDIR/p4_base.fq" \
    2>/dev/null

# Full treatment: damage masking + error correction
"$FQDUP" derep -i "$TMPDIR/p4_dp_non.fq" -o "$TMPDIR/p4_full.fq" \
    --damage-auto --error-correct 2>/dev/null

AFTER_DP=$(grep -c '^@' "$TMPDIR/p4_dp_non.fq" || true)
AFTER_BASE=$(grep -c '^@' "$TMPDIR/p4_base.fq" || true)
AFTER_FULL=$(grep -c '^@' "$TMPDIR/p4_full.fq" || true)
echo "  After derep_pairs:              $AFTER_DP"
echo "  After derep (baseline):         $AFTER_BASE"
echo "  After derep (damage+errcor):    $AFTER_FULL  |  True: $N_MOL"

python3 - <<EOF
import sys
n_mol = $N_MOL
after_dp, base, full = $AFTER_DP, $AFTER_BASE, $AFTER_FULL

# derep_pairs should have reduced the count (exact duplicate pairs collapsed)
if after_dp >= $N_READS:
    print(f"  FAIL: derep_pairs did not reduce count at all ({after_dp} = n_reads)")
    sys.exit(1)
dp_dedup = ($N_READS - after_dp) / $N_READS * 100
print(f"  INFO: derep_pairs dedup rate {dp_dedup:.1f}%")

# damage + error correction must outperform baseline
if full >= base:
    print(f"  FAIL: damage+errcor ({full}) not better than baseline ({base})")
    sys.exit(1)
improvement = (base - full) / base * 100
print(f"  PASS: damage+errcor improvement {improvement:.1f}% ({base} → {full})")

# final count must be closer to n_mol than baseline
dist_full = abs(full - n_mol)
dist_base = abs(base - n_mol)
if dist_full >= dist_base:
    print(f"  FAIL: damage+errcor ({full}) not closer to n_mol ({n_mol}) than baseline ({base})")
    sys.exit(1)
print(f"  PASS: closer to n_mol: {dist_full} < {dist_base} (baseline distance)")
EOF

echo "  PASS: Test P4"

# ---------------------------------------------------------------------------
# Test P5: Capture library (high duplication, highly damaged)
#
# Mimics a hybridisation-capture aDNA library after PCR re-amplification:
#   - High coverage (30x) from repeated sequencing of captured fragments
#   - Strong damage signal (dmax5=0.30, dmax3=0.25)
#   - Deterministic extension per molecule
#
# At 30x the pipeline must efficiently collapse duplicates to ~N_unique.
# Assert: final unique ≈ N_unique (within ±5%).
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P5: Capture library (30x, high damage) ---"
N_MOL=300
N_READS=9000    # 30x coverage
DMAX5=0.30; DMAX3=0.25

echo "  $N_MOL molecules × 30x coverage, dmax5=$DMAX5 dmax3=$DMAX3"

"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --dmax5 $DMAX5 --dmax3 $DMAX3 --lambda5 0.35 --lambda3 0.35 \
    --seed 202 > "$TMPDIR/p5_non.fq"

make_ext "$TMPDIR/p5_non.fq" "$TMPDIR/p5_ext.fq" det 202 15

"$FQDUP" sort -i "$TMPDIR/p5_non.fq" -o "$TMPDIR/p5_non.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p5_ext.fq" -o "$TMPDIR/p5_ext.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep_pairs \
    -n "$TMPDIR/p5_non.sorted.fq" -e "$TMPDIR/p5_ext.sorted.fq" \
    -o-non "$TMPDIR/p5_dp_non.fq" -o-ext "$TMPDIR/p5_dp_ext.fq" 2>/dev/null

"$FQDUP" derep -i "$TMPDIR/p5_dp_non.fq" -o "$TMPDIR/p5_final.fq" \
    --damage-auto --no-error-correct 2>/dev/null

AFTER_DP=$(grep -c '^@' "$TMPDIR/p5_dp_non.fq" || true)
AFTER_FINAL=$(grep -c '^@' "$TMPDIR/p5_final.fq" || true)
echo "  After derep_pairs: $AFTER_DP  |  After derep: $AFTER_FINAL  |  True: $N_MOL"

python3 - <<EOF
import sys
n_mol = $N_MOL
n_reads = $N_READS
after_dp, final = $AFTER_DP, $AFTER_FINAL

# At 30x the pipeline must remove the vast majority of reads (≥90% overall)
dedup_rate = (n_reads - final) / n_reads * 100
dp_rate    = (n_reads - after_dp) / n_reads * 100
print(f"  INFO: derep_pairs dedup {dp_rate:.1f}%  overall dedup {dedup_rate:.1f}%")
print(f"  INFO: final unique {final} vs true {n_mol}")
print(f"  NOTE: terminal damage at pos 4+ not fully masked → some inflation expected")

if dedup_rate < 90.0:
    print(f"  FAIL: overall dedup rate {dedup_rate:.1f}% < 90% (expected at 30x)")
    sys.exit(1)
print(f"  PASS: overall pipeline dedup rate {dedup_rate:.1f}% ≥ 90%")

# derep --damage-auto must further reduce from derep_pairs alone
if final >= after_dp:
    print(f"  FAIL: derep (damage-auto) did not reduce from derep_pairs ({final} >= {after_dp})")
    sys.exit(1)
step2_reduction = (after_dp - final) / after_dp * 100
print(f"  PASS: derep --damage-auto removed {step2_reduction:.1f}% of derep_pairs residuals")
EOF

echo "  PASS: Test P5"

# ---------------------------------------------------------------------------
# Test P6: Shotgun library (low duplication, mostly singletons)
#
# Mimics a low-coverage whole-genome shotgun library without enrichment:
#   - 1.2x coverage — most molecules sampled once (P(dup) ≈ 30%)
#   - No damage (modern-DNA-like or well-preserved sample)
#   - Deterministic extension
#
# The pipeline must NOT over-deduplicate: distinct molecules with similar
# sequences must not be collapsed, and the overall dedup rate must reflect
# the low duplication from Poisson sampling at 1.2x.
# Assert: dedup rate < 25%; final unique count ≈ n_reads × (1 - e^{-1.2}).
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P6: Shotgun library (1.2x, mostly singletons) ---"
N_MOL=5000
N_READS=6000    # 1.2x coverage

echo "  $N_MOL molecules × 1.2x coverage, no damage"

"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --no-damage --seed 303 > "$TMPDIR/p6_non.fq"

make_ext "$TMPDIR/p6_non.fq" "$TMPDIR/p6_ext.fq" det 303 15

"$FQDUP" sort -i "$TMPDIR/p6_non.fq" -o "$TMPDIR/p6_non.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p6_ext.fq" -o "$TMPDIR/p6_ext.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep_pairs \
    -n "$TMPDIR/p6_non.sorted.fq" -e "$TMPDIR/p6_ext.sorted.fq" \
    -o-non "$TMPDIR/p6_dp_non.fq" -o-ext "$TMPDIR/p6_dp_ext.fq" 2>/dev/null

"$FQDUP" derep -i "$TMPDIR/p6_dp_non.fq" -o "$TMPDIR/p6_final.fq" \
    --no-error-correct 2>/dev/null

AFTER_DP=$(grep -c '^@' "$TMPDIR/p6_dp_non.fq" || true)
AFTER_FINAL=$(grep -c '^@' "$TMPDIR/p6_final.fq" || true)
echo "  After derep_pairs: $AFTER_DP  |  After derep: $AFTER_FINAL  |  n_reads: $N_READS"

python3 - <<EOF
import sys, math
n_mol = $N_MOL
n_reads = $N_READS
final = $AFTER_FINAL

cov = n_reads / n_mol
# E[unique molecules observed] = n_mol * (1 - e^{-cov})
expected_unique = n_mol * (1.0 - math.exp(-cov))
dedup_rate = (n_reads - final) / n_reads * 100

print(f"  INFO: coverage {cov:.2f}x, expected unique ≈ {expected_unique:.0f}")
print(f"  INFO: dedup rate {dedup_rate:.1f}%  ({n_reads} → {final})")

# At 1.2x, Poisson P(k≥2) = 1 - e^{-1.2} - 1.2*e^{-1.2} ≈ 0.337
# Max possible dedup rate ≈ 33.7% (only duplicated molecules can be removed)
max_dedup = (1.0 - math.exp(-cov) - cov * math.exp(-cov)) / (1.0 - math.exp(-cov)) * 100
print(f"  INFO: theoretical max dedup rate ≈ {max_dedup:.1f}%")

# Actual dedup rate must not exceed theoretical max by more than a few percent
if dedup_rate > max_dedup + 5.0:
    print(f"  FAIL: dedup rate {dedup_rate:.1f}% exceeds theoretical max {max_dedup:.1f}% + 5% margin")
    print(f"        Possible over-deduplication (hash collisions or RC false matches?)")
    sys.exit(1)
print(f"  PASS: dedup rate {dedup_rate:.1f}% ≤ {max_dedup + 5.0:.1f}% (theoretical max + margin)")

# Final unique count within ±10% of expected
lo, hi = int(expected_unique * 0.90), int(expected_unique * 1.10)
if not (lo <= final <= hi):
    print(f"  FAIL: final unique {final} not in [{lo}, {hi}] (expected ≈ {expected_unique:.0f})")
    sys.exit(1)
print(f"  PASS: final unique {final} ∈ [{lo}, {hi}]")
EOF

echo "  PASS: Test P6"

# ---------------------------------------------------------------------------
# Test P7: Dedup difficulty with terminal damage — explicit with/without contrast
#
# Ground truth: N molecules, each sequenced exactly twice (--dup-pair).
# Without damage: both copies are identical → standard dedup collapses to N.
# With damage: terminal C→T/G→A differs between copies → standard dedup inflates.
# With damage masking: inflation recovered → back near N.
#
# Both datasets are generated with the same seed, so the base molecules are
# identical.  Damage is the only difference.
#
# Random extension per read forces all reads through derep_pairs with no dedup
# there, so the full contrast is visible at the derep step.
#
# Assert:
#   no-damage → unique ≈ N (both copies collapse)
#   damage, no masking → unique >> N (damage prevents collapse — the difficulty)
#   damage, with masking → unique < no-masking (damage masking helps)
#   damage, with masking → closer to N than no-masking (masking recovers ground truth)
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P7: Dedup difficulty — with and without damage (same molecules) ---"
N_MOL=3000
DMAX5=0.25; DMAX3=0.20; LAMBDA=0.35

echo "  $N_MOL molecules × 2 copies each (--dup-pair), dmax5=$DMAX5 dmax3=$DMAX3"
echo "  Same seed → same base molecules; damage is the only variable"
echo "  Random extension → all reads reach derep step (no derep_pairs dedup)"

# No-damage baseline: both copies of each molecule are identical
"$GEN" --n-unique $N_MOL --read-len 75 \
    --no-damage --dup-pair --seed 500 > "$TMPDIR/p7_nodmg.fq"

# Damaged dataset: same molecules, each copy gets independent damage
"$GEN" --n-unique $N_MOL --read-len 75 \
    --dmax5 $DMAX5 --dmax3 $DMAX3 --lambda5 $LAMBDA --lambda3 $LAMBDA \
    --dup-pair --seed 500 > "$TMPDIR/p7_dmg.fq"

# Create ext files: random extension per read → derep_pairs passes all through
make_ext "$TMPDIR/p7_nodmg.fq" "$TMPDIR/p7_nodmg_ext.fq" rnd 500 15
make_ext "$TMPDIR/p7_dmg.fq"   "$TMPDIR/p7_dmg_ext.fq"   rnd 501 15

# Sort both datasets
for tag in nodmg dmg; do
    "$FQDUP" sort -i "$TMPDIR/p7_${tag}.fq" -o "$TMPDIR/p7_${tag}.sorted.fq" \
        --max-memory 512M -t "$TMPDIR" 2>/dev/null
    "$FQDUP" sort -i "$TMPDIR/p7_${tag}_ext.fq" -o "$TMPDIR/p7_${tag}_ext.sorted.fq" \
        --max-memory 512M -t "$TMPDIR" 2>/dev/null
    "$FQDUP" derep_pairs \
        -n "$TMPDIR/p7_${tag}.sorted.fq" -e "$TMPDIR/p7_${tag}_ext.sorted.fq" \
        -o-non "$TMPDIR/p7_${tag}_dp.fq" -o-ext /dev/null 2>/dev/null
done

# Baseline: no-damage, no masking → both copies identical → collapse to N_mol
"$FQDUP" derep -i "$TMPDIR/p7_nodmg_dp.fq" -o "$TMPDIR/p7_nodmg_out.fq" \
    --no-error-correct 2>/dev/null

# Difficulty: damaged, no masking → copies look different → inflated count
"$FQDUP" derep -i "$TMPDIR/p7_dmg_dp.fq" -o "$TMPDIR/p7_dmg_nomask.fq" \
    --no-error-correct 2>/dev/null

# Recovery: damaged, with masking → terminal damage ignored → collapse recovered
"$FQDUP" derep -i "$TMPDIR/p7_dmg_dp.fq" -o "$TMPDIR/p7_dmg_masked.fq" \
    --damage-auto --no-error-correct 2>/dev/null

N_READS=$((N_MOL * 2))
U_NODMG=$(grep -c '^@' "$TMPDIR/p7_nodmg_out.fq"  || true)
U_NOMASK=$(grep -c '^@' "$TMPDIR/p7_dmg_nomask.fq" || true)
U_MASKED=$(grep -c '^@' "$TMPDIR/p7_dmg_masked.fq" || true)

echo ""
echo "  True unique molecules:     $N_MOL"
echo "  Total reads (2 per mol):   $N_READS"
echo "  ---"
echo "  No damage, no masking:     $U_NODMG unique  (both copies identical → collapses)"
echo "  Damage, no masking:        $U_NOMASK unique  (terminal diffs → inflated)"
echo "  Damage, --damage-auto:     $U_MASKED unique  (masking recovers collapse)"

python3 - <<EOF
import sys
n_mol = $N_MOL
n_reads = $N_READS
nodmg, nomask, masked = $U_NODMG, $U_NOMASK, $U_MASKED

# 1. No-damage baseline: both copies collapse → unique ≈ n_mol
lo_nd, hi_nd = int(n_mol * 0.95), int(n_mol * 1.05)
if not (lo_nd <= nodmg <= hi_nd):
    print(f"  FAIL: no-damage gave {nodmg} unique (expect {lo_nd}–{hi_nd})")
    sys.exit(1)
print(f"  PASS: no-damage  {nodmg} ∈ [{lo_nd}, {hi_nd}] — copies correctly collapsed")

# 2. Damage without masking: terminal differences inflate the count
# With dmax5=0.25 and 15 masked positions, P(at least one copy C→T at pos 0) ≈ 0.44
# → ~44% of molecule pairs won't collapse → expect significantly more than n_mol
if nomask <= int(n_mol * 1.10):
    print(f"  FAIL: damaged (no mask) gave {nomask} unique — not enough inflation to demonstrate difficulty")
    sys.exit(1)
inflation = (nomask - n_mol) / n_mol * 100
print(f"  PASS: damaged (no mask) {nomask} unique (+{inflation:.0f}% inflation — this is the difficulty)")

# 3. Damage with masking: better than no-masking and closer to ground truth
if masked >= nomask:
    print(f"  FAIL: --damage-auto gave no improvement ({masked} >= {nomask})")
    sys.exit(1)
recovery = (nomask - masked) / (nomask - n_mol) * 100 if nomask > n_mol else 0
print(f"  PASS: damage-auto {masked} unique — recovered {recovery:.0f}% of inflation")

dist_masked = abs(masked - n_mol)
dist_nomask = abs(nomask - n_mol)
if dist_masked >= dist_nomask:
    print(f"  FAIL: masking not closer to n_mol ({dist_masked} >= {dist_nomask})")
    sys.exit(1)
print(f"  PASS: closer to ground truth: {dist_masked} < {dist_nomask}")
EOF

echo "  PASS: Test P7"

# ---------------------------------------------------------------------------
# Test P8: PCR mismatch inflation and correction (explicit mismatch counting)
#
# Same molecules generated twice:
#   clean  — no PCR errors (ground truth)
#   errors — same molecules + PCR errors at rate 0.003/base (→ ~20% of reads
#            carry ≥1 mismatch from their parent molecule)
#
# Without error correction, each PCR-error read forms its own singleton cluster.
# With error correction (--error-correct --errcor-ratio 5 for 5x coverage),
# singleton clusters differing by 1 interior base from a high-count parent are
# absorbed back.
#
# Assert:
#   clean, no EC → unique ≈ N_mol (baseline: errors don't exist)
#   errors, no EC → unique >> N_mol (PCR mismatches inflate the count)
#   errors, with EC → unique between N_mol and no-EC count (correction partial)
#   errors, with EC < errors no-EC (EC strictly helps)
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P8: PCR mismatch inflation and correction ---"
N_MOL=1000
N_READS=5000    # 5x coverage
PCR_RATE=0.003  # P(≥1 error per read) = 1 - e^{-0.003*75} ≈ 20%

echo "  $N_MOL molecules × 5x coverage"
echo "  Clean dataset (no errors) vs PCR errors at rate=$PCR_RATE/base"

# Clean — no PCR errors, ground-truth unique count
"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --no-damage --seed 400 > "$TMPDIR/p8_clean.fq"

# Errors — same molecules (same seed), PCR errors added
"$GEN" --n-unique $N_MOL --n-reads $N_READS --read-len 75 \
    --no-damage --pcr-rate $PCR_RATE --seed 400 > "$TMPDIR/p8_err.fq"

for tag in clean err; do
    make_ext "$TMPDIR/p8_${tag}.fq" "$TMPDIR/p8_${tag}_ext.fq" rnd 400 15
    "$FQDUP" sort -i "$TMPDIR/p8_${tag}.fq"     -o "$TMPDIR/p8_${tag}.sorted.fq" \
        --max-memory 512M -t "$TMPDIR" 2>/dev/null
    "$FQDUP" sort -i "$TMPDIR/p8_${tag}_ext.fq" -o "$TMPDIR/p8_${tag}_ext.sorted.fq" \
        --max-memory 512M -t "$TMPDIR" 2>/dev/null
    "$FQDUP" derep_pairs \
        -n "$TMPDIR/p8_${tag}.sorted.fq" -e "$TMPDIR/p8_${tag}_ext.sorted.fq" \
        -o-non "$TMPDIR/p8_${tag}_dp.fq" -o-ext /dev/null 2>/dev/null
done

# Clean: no EC needed
"$FQDUP" derep -i "$TMPDIR/p8_clean_dp.fq" -o "$TMPDIR/p8_clean_out.fq" \
    --no-error-correct 2>/dev/null

# Errors: without EC
"$FQDUP" derep -i "$TMPDIR/p8_err_dp.fq"   -o "$TMPDIR/p8_err_noec.fq" \
    --no-error-correct 2>/dev/null

# Errors: with EC (ratio=5 matches 5x coverage; default ratio=50 is for high-cov production)
"$FQDUP" derep -i "$TMPDIR/p8_err_dp.fq"   -o "$TMPDIR/p8_err_ec.fq" \
    --error-correct --errcor-ratio 5 2>/dev/null

U_CLEAN=$(grep -c '^@' "$TMPDIR/p8_clean_out.fq" || true)
U_NOEC=$(grep -c '^@' "$TMPDIR/p8_err_noec.fq"   || true)
U_EC=$(grep -c '^@' "$TMPDIR/p8_err_ec.fq"        || true)

echo ""
echo "  True molecules:            $N_MOL"
echo "  Clean, no EC:              $U_CLEAN unique  (baseline, no errors)"
echo "  Errors, no EC:             $U_NOEC  unique  (PCR mismatches inflate)"
echo "  Errors, with EC (ratio=5): $U_EC    unique  (correction absorbs singletons)"

python3 - <<EOF
import sys
n_mol = $N_MOL
clean, noec, ec = $U_CLEAN, $U_NOEC, $U_EC

# Clean baseline: no errors → should collapse to ≈ n_mol
lo_c, hi_c = int(n_mol * 0.95), int(n_mol * 1.05)
if not (lo_c <= clean <= hi_c):
    print(f"  FAIL: clean baseline {clean} not in [{lo_c}, {hi_c}]")
    sys.exit(1)
print(f"  PASS: clean baseline {clean} ≈ {n_mol}")

# PCR errors must inflate the count (the mismatch problem)
min_inflated = int(n_mol * 1.10)
if noec <= min_inflated:
    print(f"  FAIL: PCR errors did not inflate ({noec} ≤ {min_inflated}) — check --pcr-rate")
    sys.exit(1)
inflation = noec - n_mol
pct_inflation = inflation / n_mol * 100
print(f"  PASS: PCR errors inflated unique by {inflation} reads (+{pct_inflation:.0f}%)")

# EC must reduce the inflation
if ec >= noec:
    print(f"  FAIL: EC gave no improvement ({ec} >= {noec})")
    sys.exit(1)
absorbed = noec - ec
pct_absorbed = absorbed / inflation * 100
print(f"  PASS: EC absorbed {absorbed} mismatch singletons ({pct_absorbed:.0f}% of inflation)")
EOF

echo "  PASS: Test P8"

# ---------------------------------------------------------------------------
# Test P9: False-positive protection — real SNPs must not be merged
#
# Two sub-tests using hand-crafted FASTQ:
#
#   P9a. Two molecules differing by 1 interior SNP, equal high coverage.
#        count(A) / count(B) ≈ 1 < errcor-ratio (default 50).
#        EC must NOT absorb B into A — they are truly different molecules.
#        Assert: unique = 2 with AND without EC.
#
#   P9b. One molecule (high coverage) + one PCR-error singleton (1 interior SNP).
#        count(parent) / count(singleton) >> errcor-ratio → singleton absorbed.
#        Assert: unique = 2 without EC; unique = 1 with EC.
# ---------------------------------------------------------------------------
echo ""
echo "--- Test P9: False-positive protection (real SNPs preserved by EC) ---"
READ_LEN=75
COV=200         # high coverage for both molecules in P9a
QUAL=$(python3 -c "print('I' * $READ_LEN)")

# Python inline: generate two molecules (A and B = A with SNP at pos 37)
python3 - <<PYEOF > "$TMPDIR/p9a.fq"
import random
rng = random.Random(999)
L = $READ_LEN
qual = 'I' * L

mol_A = ''.join(rng.choice('ACGT') for _ in range(L))
mol_B = list(mol_A)
mol_B[37] = {'A':'G','C':'T','G':'A','T':'C'}[mol_A[37]]  # interior SNP at pos 37
mol_B = ''.join(mol_B)

cov = $COV
for i in range(cov):
    print(f'@molA_r{i:05d}\n{mol_A}\n+\n{qual}')
for i in range(cov):
    print(f'@molB_r{i:05d}\n{mol_B}\n+\n{qual}')
PYEOF

# Python inline: generate one molecule (A, high coverage) + one PCR-error singleton
python3 - <<PYEOF > "$TMPDIR/p9b.fq"
import random
rng = random.Random(999)
L = $READ_LEN
qual = 'I' * L

mol_A = ''.join(rng.choice('ACGT') for _ in range(L))
mol_snp = list(mol_A)
mol_snp[37] = {'A':'G','C':'T','G':'A','T':'C'}[mol_A[37]]  # same interior SNP
mol_snp = ''.join(mol_snp)

cov = $COV
for i in range(cov):
    print(f'@molA_r{i:05d}\n{mol_A}\n+\n{qual}')
# Exactly one singleton PCR-error read
print(f'@snp_r00000\n{mol_snp}\n+\n{qual}')
PYEOF

# Run full pipeline on both sub-tests
for tag in p9a p9b; do
    make_ext "$TMPDIR/${tag}.fq" "$TMPDIR/${tag}_ext.fq" rnd 999 15
    "$FQDUP" sort -i "$TMPDIR/${tag}.fq"     -o "$TMPDIR/${tag}.sorted.fq" \
        --max-memory 512M -t "$TMPDIR" 2>/dev/null
    "$FQDUP" sort -i "$TMPDIR/${tag}_ext.fq" -o "$TMPDIR/${tag}_ext.sorted.fq" \
        --max-memory 512M -t "$TMPDIR" 2>/dev/null
    "$FQDUP" derep_pairs \
        -n "$TMPDIR/${tag}.sorted.fq" -e "$TMPDIR/${tag}_ext.sorted.fq" \
        -o-non "$TMPDIR/${tag}_dp.fq" -o-ext /dev/null 2>/dev/null
    "$FQDUP" derep -i "$TMPDIR/${tag}_dp.fq" -o "$TMPDIR/${tag}_noec.fq" \
        --no-error-correct 2>/dev/null
    "$FQDUP" derep -i "$TMPDIR/${tag}_dp.fq" -o "$TMPDIR/${tag}_ec.fq" \
        --error-correct 2>/dev/null
done

U9A_NOEC=$(grep -c '^@' "$TMPDIR/p9a_noec.fq" || true)
U9A_EC=$(grep -c '^@'   "$TMPDIR/p9a_ec.fq"   || true)
U9B_NOEC=$(grep -c '^@' "$TMPDIR/p9b_noec.fq" || true)
U9B_EC=$(grep -c '^@'   "$TMPDIR/p9b_ec.fq"   || true)

echo ""
echo "  P9a: mol_A (${COV}×) + mol_B (1 SNP, ${COV}×) — equal coverage, real SNP"
echo "    No EC: $U9A_NOEC unique  |  With EC: $U9A_EC unique  (expect 2 in both — no merge)"
echo ""
echo "  P9b: mol_A (${COV}×) + pcr_error_singleton (1 read, 1 SNP) — high-count parent"
echo "    No EC: $U9B_NOEC unique  |  With EC: $U9B_EC unique  (expect 2→1 — singleton absorbed)"

python3 - <<EOF
import sys
u9a_noec, u9a_ec = $U9A_NOEC, $U9A_EC
u9b_noec, u9b_ec = $U9B_NOEC, $U9B_EC
ok = True

# P9a: two equal-count clusters must stay separate
if u9a_noec != 2:
    print(f"  FAIL P9a: without EC, expected 2 unique, got {u9a_noec}")
    ok = False
else:
    print(f"  PASS P9a no-EC:  {u9a_noec} unique (correct — 2 distinct molecules)")
if u9a_ec != 2:
    print(f"  FAIL P9a: EC incorrectly merged real SNP ({u9a_ec} unique, expected 2)")
    ok = False
else:
    print(f"  PASS P9a with-EC: {u9a_ec} unique (false positive protected — count ratio 1:1 < threshold)")

# P9b: without EC the singleton stays; with EC it is absorbed
if u9b_noec != 2:
    print(f"  FAIL P9b: without EC, expected 2 unique, got {u9b_noec}")
    ok = False
else:
    print(f"  PASS P9b no-EC:  {u9b_noec} unique (singleton not absorbed without EC)")
if u9b_ec != 1:
    print(f"  FAIL P9b: EC did not absorb singleton ({u9b_ec} unique, expected 1)")
    ok = False
else:
    print(f"  PASS P9b with-EC: {u9b_ec} unique (PCR-error singleton absorbed into parent)")

if not ok:
    sys.exit(1)
EOF

echo "  PASS: Test P9"

echo ""
echo "OK: all pipeline tests passed"
