#!/usr/bin/env bash
# test_math_invariants.sh — property tests for mathematical correctness.
#
# Covers:
#   M1. gen: exact read count output
#   M2. damage: d_max estimates remain in [0, 1]
#   M3. derep: more masking → unique count ≤ less masking (monotonicity)
#   M4. derep: stricter snp_threshold → ≥ as many unique reads (less absorption)
#   M5. derep: --no-error-correct produces ≥ unique reads as --error-correct
#   M6. Phase 3: C↔T mismatch is never absorbed (damage protection invariant)
set -euo pipefail

FQDUP=${1:-fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== fqdup math/model invariant tests ==="

# ---------------------------------------------------------------------------
# M1: gen produces exactly N reads
# ---------------------------------------------------------------------------
echo ""
echo "--- M1: gen produces exactly N reads ---"
"$GEN" --n-unique 50 --n-reads 500 --read-len 60 --no-damage --seed 1 \
    > "$TMPDIR/m1.fq" 2>/dev/null
N=$(grep -c '^@' "$TMPDIR/m1.fq")
if [[ "$N" -ne 500 ]]; then
    echo "  FAIL: expected 500 reads, got $N"
    exit 1
fi
echo "  PASS: gen produced exactly 500 reads"

# ---------------------------------------------------------------------------
# M2: damage d_max estimates in [0, 1]
# ---------------------------------------------------------------------------
echo ""
echo "--- M2: damage d_max estimates in [0, 1] ---"
"$GEN" --n-unique 200 --n-reads 5000 --read-len 75 \
    --dmax5 0.30 --dmax3 0.20 --lambda5 0.40 --lambda3 0.40 \
    --seed 2 > "$TMPDIR/m2.fq" 2>/dev/null
DMAX5=$("$FQDUP" damage -i "$TMPDIR/m2.fq" 2>&1 | \
    grep "5'-end" | grep -oP 'd_max=\K[0-9.]+' | head -1)
DMAX3=$("$FQDUP" damage -i "$TMPDIR/m2.fq" 2>&1 | \
    grep "3'-end" | grep -oP 'd_max=\K[0-9.]+' | head -1)
python3 - <<EOF
import sys
d5, d3 = float("${DMAX5:-0}"), float("${DMAX3:-0}")
if not (0.0 <= d5 <= 1.0):
    print(f"  FAIL: d_max_5={d5} out of [0,1]"); sys.exit(1)
if not (0.0 <= d3 <= 1.0):
    print(f"  FAIL: d_max_3={d3} out of [0,1]"); sys.exit(1)
print(f"  PASS: d_max_5={d5:.3f} d_max_3={d3:.3f} both in [0,1]")
EOF

# ---------------------------------------------------------------------------
# M3: more masking → unique count ≤ less masking (monotonicity)
# ---------------------------------------------------------------------------
echo ""
echo "--- M3: masking monotonicity (more mask → fewer or equal unique) ---"
"$GEN" --n-unique 300 --n-reads 3000 --read-len 75 \
    --dmax5 0.25 --dmax3 0.20 --lambda5 0.35 --lambda3 0.35 \
    --seed 3 > "$TMPDIR/m3.fq" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/m3.fq" -o "$TMPDIR/m3.sorted.fq" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null

U_NODMG=$("$FQDUP" derep -i "$TMPDIR/m3.sorted.fq" -o "$TMPDIR/m3_nodmg.fq" \
    --no-damage --no-error-correct 2>/dev/null && \
    grep -c '^@' "$TMPDIR/m3_nodmg.fq")
U_DMG=$("$FQDUP" derep -i "$TMPDIR/m3.sorted.fq" -o "$TMPDIR/m3_dmg.fq" \
    --damage-auto --no-error-correct 2>/dev/null && \
    grep -c '^@' "$TMPDIR/m3_dmg.fq")

python3 - <<EOF
import sys
no_dmg, dmg = $U_NODMG, $U_DMG
if dmg > no_dmg:
    print(f"  FAIL: masking increased unique count ({no_dmg} → {dmg})")
    sys.exit(1)
print(f"  PASS: --damage-auto {no_dmg} → {dmg} unique (damage masking reduces or holds)")
EOF

# ---------------------------------------------------------------------------
# M4: stricter snp_threshold → ≥ as many unique reads (less EC absorption)
# ---------------------------------------------------------------------------
echo ""
echo "--- M4: stricter snp_threshold → ≥ unique reads (less absorption) ---"
"$GEN" --n-unique 200 --n-reads 4000 --read-len 60 \
    --no-damage --seed 4 > "$TMPDIR/m4.fq" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/m4.fq" -o "$TMPDIR/m4.sorted.fq" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null

U_LOOSE=$("$FQDUP" derep -i "$TMPDIR/m4.sorted.fq" -o "$TMPDIR/m4_loose.fq" \
    --errcor-snp-threshold 0.05 2>/dev/null && \
    grep -c '^@' "$TMPDIR/m4_loose.fq")
U_STRICT=$("$FQDUP" derep -i "$TMPDIR/m4.sorted.fq" -o "$TMPDIR/m4_strict.fq" \
    --errcor-snp-threshold 0.50 2>/dev/null && \
    grep -c '^@' "$TMPDIR/m4_strict.fq")

python3 - <<EOF
import sys
loose, strict = $U_LOOSE, $U_STRICT
if strict < loose:
    print(f"  FAIL: stricter threshold absorbed MORE ({loose} → {strict})")
    sys.exit(1)
print(f"  PASS: snp_threshold 0.05→{loose} unique, 0.50→{strict} unique (monotone)")
EOF

# ---------------------------------------------------------------------------
# M5: --no-error-correct ≥ --error-correct unique count
# ---------------------------------------------------------------------------
echo ""
echo "--- M5: EC reduces or holds unique count (no EC ≥ with EC) ---"
U_NOEC=$(grep -c '^@' "$TMPDIR/m4_strict.fq" || true)
U_EC=$("$FQDUP" derep -i "$TMPDIR/m4.sorted.fq" -o "$TMPDIR/m4_ec.fq" \
    --error-correct 2>/dev/null && grep -c '^@' "$TMPDIR/m4_ec.fq")
U_NOEC2=$("$FQDUP" derep -i "$TMPDIR/m4.sorted.fq" -o "$TMPDIR/m4_noec.fq" \
    --no-error-correct 2>/dev/null && grep -c '^@' "$TMPDIR/m4_noec.fq")

python3 - <<EOF
import sys
ec, noec = $U_EC, $U_NOEC2
if ec > noec:
    print(f"  FAIL: EC increased unique count ({noec} → {ec})")
    sys.exit(1)
print(f"  PASS: no-EC={noec} ≥ EC={ec} (EC never inflates unique count)")
EOF

# ---------------------------------------------------------------------------
# M6: C↔T mismatch never absorbed by Phase 3
# ---------------------------------------------------------------------------
echo ""
echo "--- M6: C↔T mismatch never absorbed (damage protection invariant) ---"
PARENT="ACGAACGAACGAACGAACGAACGAACGAAC"
# C→T at position 13 (interior, damage-consistent)
CT_CHILD="ACGAACGAACGAATGAACGAACGAACGAAC"
python3 - <<PYEOF > "$TMPDIR/m6.fq"
parent  = '$PARENT'
ct_sub  = '$CT_CHILD'
qual    = 'I' * 31
for i in range(100):
    print(f'@par_{i:04d}\n{parent}\n+\n{qual}')
print(f'@ct_child\n{ct_sub}\n+\n{qual}')
PYEOF
"$FQDUP" sort -i "$TMPDIR/m6.fq" -o "$TMPDIR/m6.sorted.fq" \
    --max-memory 64M -t "$TMPDIR" 2>/dev/null
U_EC=$("$FQDUP" derep -i "$TMPDIR/m6.sorted.fq" -o "$TMPDIR/m6_ec.fq" \
    --error-correct 2>/dev/null && grep -c '^@' "$TMPDIR/m6_ec.fq")
python3 - <<EOF
import sys
u = $U_EC
if u != 2:
    print(f"  FAIL: C↔T child absorbed (got {u} unique, expected 2)")
    sys.exit(1)
print(f"  PASS: C↔T child kept ({u} unique — damage protection held)")
EOF

echo ""
echo "=== All math/model invariant tests passed ==="
