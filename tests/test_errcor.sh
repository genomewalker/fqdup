#!/usr/bin/env bash
# tests/test_errcor.sh — validate Phase 3 PCR error correction
#
# Five cases:
#   1. H=1, A→T transversion (non-damage), count=1                → absorbed
#   2. H=1, C→T (deamination), count=1                            → NOT absorbed
#   3. H=2, C→G + C→T (one damage sub), count=1                   → NOT absorbed
#   4. H=2, A→T + C→G (both non-damage transversions), count=1    → absorbed (new)
#   5. H=2, C→G + A→T (both non-damage transversions), count=3    → NOT absorbed (above max_h2_count)
#
# Uses a 30-bp parent at count=100 and child counts as noted above.
#
# Case 3 uses positions 5 (C→G) and 25 (C→T) — distinct from cases 1 (pos 15)
# and 2 (pos 13) to avoid spurious H=1 edges between test sequences.
set -euo pipefail

FQDUP=${1:-fqdup}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
cd "$TMPDIR"

# 30-bp parent (repeating ACGAAC unit)
PARENT="ACGAACGAACGAACGAACGAACGAACGAAC"
# pos 15 A→T : A↔T transversion, non-damage → absorbed (H=1)
ERROR_1SUB="ACGAACGAACGAACGTACGAACGAACGAAC"
# pos 13 C→T : deamination pattern → NOT absorbed (H=1, damage sub)
DAMAGE_1SUB="ACGAACGAACGAATGAACGAACGAACGAAC"
# pos 5 C→G + pos 25 C→T : one damage sub → NOT absorbed (H=2, damage sub present)
# Positions chosen to be ≥3 away from cases 1 (pos 15) and 2 (pos 13)
TWO_SUB_DMG="ACGAAGGAACGAACGAACGAACGAATGAAC"
# pos 8 A→T + pos 21 C→G : both non-damage transversions → absorbed (H=2)
TWO_TRANS="ACGAACGATCGAACGAACGAAGGAACGAAC"

python3 -c "
parent       = '$PARENT'
error_1sub   = '$ERROR_1SUB'
damage_1sub  = '$DAMAGE_1SUB'
two_sub_dmg  = '$TWO_SUB_DMG'
two_trans    = '$TWO_TRANS'
q            = 'I' * len(parent)
for i in range(1, 101):
    print(f'@apar_{i:04d}\n{parent}\n+\n{q}')
# Case 1: H=1 A→T, count=1 → should be absorbed
print(f'@berr_0001\n{error_1sub}\n+\n{q}')
# Case 2: H=1 C→T, count=1 → should be preserved (damage sub)
print(f'@cdam_0001\n{damage_1sub}\n+\n{q}')
# Case 3: H=2 C→G+C→T, count=1 → should be preserved (damage sub present)
print(f'@ddmg_0001\n{two_sub_dmg}\n+\n{q}')
# Case 4: H=2 both transversions, count=1 → should be absorbed
print(f'@etrn_0001\n{two_trans}\n+\n{q}')
# Case 5: distinct H=2 both-transversion sequence at count=3 → should NOT be absorbed
two_trans3 = list(parent)
two_trans3[9]  = 'G'   # C→G
two_trans3[20] = 'T'   # A→T
two_trans3 = ''.join(two_trans3)
for i in range(1, 4):
    print(f'@ftrn_{i:04d}\n{two_trans3}\n+\n{q}')
" > input.fq

"$FQDUP" sort -i input.fq -o sorted.fq --max-memory 1G -t . --fast 2>/dev/null
"$FQDUP" derep -i sorted.fq -o output.fq --error-correct 2>/dev/null

# Derive the count=3 two-transversion sequence for case 5 check
TWO_TRANS3=$(python3 -c "
s = list('$PARENT'); s[9]='G'; s[20]='T'; print(''.join(s))")

# 1. Non-damage H=1 (A→T) must be absorbed
if grep -q "^${ERROR_1SUB}$" output.fq; then
    echo "FAIL: H=1 A→T child was not absorbed (expected absorption)"
    exit 1
fi
echo "OK: H=1 non-damage (A→T) absorbed"

# 2. H=1 damage sub (C→T) must be preserved
if ! grep -q "^${DAMAGE_1SUB}$" output.fq; then
    echo "FAIL: H=1 C→T child was absorbed (should be preserved as damage)"
    exit 1
fi
echo "OK: H=1 damage sub (C→T) preserved"

# 3. H=2 with one damage sub must be preserved
if ! grep -q "^${TWO_SUB_DMG}$" output.fq; then
    echo "FAIL: H=2 (C→G+C→T) child absorbed (damage sub present, should be preserved)"
    exit 1
fi
echo "OK: H=2 with damage sub (C→G+C→T) preserved"

# 4. H=2, both non-damage transversions, count=1 → must be absorbed
if grep -q "^${TWO_TRANS}$" output.fq; then
    echo "FAIL: H=2 (A→T+C→G) count=1 child was NOT absorbed (expected H=2 absorption)"
    exit 1
fi
echo "OK: H=2 both-transversion (A→T+C→G) count=1 absorbed"

# 5. H=2, both non-damage transversions, count=3 → must NOT be absorbed (above max_h2_count)
if ! grep -q "^${TWO_TRANS3}$" output.fq; then
    echo "FAIL: H=2 (C→G+A→T) count=3 child was absorbed (count > max_h2_count=2)"
    exit 1
fi
echo "OK: H=2 both-transversion count=3 preserved (above max_h2_count)"

echo "OK: test_errcor passed"
