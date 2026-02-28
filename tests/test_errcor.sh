#!/usr/bin/env bash
# tests/test_errcor.sh — validate Phase 3 PCR error correction
#
# Three cases:
#   1. 1 interior A→T substitution (transversion, no damage mechanism) → absorbed
#   2. 1 interior C→T substitution (ancient deamination pattern)       → NOT absorbed
#   3. 2 interior substitutions                                         → NOT absorbed
#
# Uses a 31-bp parent at count=100 and single-count children (ratio=100 ≥ 50).
set -euo pipefail

FQDUP=${1:-fqdup}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
cd "$TMPDIR"

# 31-bp parent
PARENT="ACGAACGAACGAACGAACGAACGAACGAAC"
# pos 15 (A) → T : A↔T transversion, non-damage → should be absorbed
ERROR_1SUB="ACGAACGAACGAACGTACGAACGAACGAAC"
# pos 13 (C) → T : C↔T, ancient deamination pattern → should NOT be absorbed
DAMAGE_1SUB="ACGAACGAACGAATGAACGAACGAACGAAC"
# pos 13 C→T + pos 15 A→T : Hamming 2 → should NOT be absorbed
TWO_SUB="ACGAACGAACGAATGTACGAACGAACGAAC"

# Build sorted FASTQ (IDs pre-ordered: apar_ < berr_ < cdam_ < dtwo_)
python3 -c "
parent  = '$PARENT'
error   = '$ERROR_1SUB'
damage  = '$DAMAGE_1SUB'
twosub  = '$TWO_SUB'
q       = 'I' * len(parent)
for i in range(1, 101):
    print(f'@apar_{i:04d}\n{parent}\n+\n{q}')
print(f'@berr_0001\n{error}\n+\n{q}')
print(f'@cdam_0001\n{damage}\n+\n{q}')
print(f'@dtwo_0001\n{twosub}\n+\n{q}')
" > input.fq

"$FQDUP" sort -i input.fq -o sorted.fq --max-memory 1G -t . --fast 2>/dev/null
"$FQDUP" derep -i sorted.fq -o output.fq --error-correct 2>/dev/null

# 1. Non-damage 1-sub must be absorbed
if grep -q "^${ERROR_1SUB}$" output.fq; then
    echo "FAIL: A→T child was not absorbed (expected absorption)"
    exit 1
fi
echo "OK: non-damage 1-sub (A→T) absorbed"

# 2. Damage 1-sub must be preserved
if ! grep -q "^${DAMAGE_1SUB}$" output.fq; then
    echo "FAIL: C→T child was absorbed (should be preserved as damage)"
    exit 1
fi
echo "OK: damage 1-sub (C→T) preserved"

# 3. Two-sub must be preserved
if ! grep -q "^${TWO_SUB}$" output.fq; then
    echo "FAIL: two-sub child was absorbed (Hamming 2 should not be absorbed)"
    exit 1
fi
echo "OK: two-sub child preserved"

echo "OK: test_errcor passed"
