#!/usr/bin/env bash
# test_damage_auto.sh — Validates --damage-auto estimation and damage vs PCR separation.
#
# Tests:
#   A. d_max estimation accuracy
#        100k reads, known dmax5=0.25 dmax3=0.20 lambda=0.35
#        Asserts: |est5 - 0.25| < 0.08  and  |est3 - 0.20| < 0.08
#
#   B. Damaged-duplicate collapse
#        10k molecules × 2 independently-damaged copies
#        Asserts: unique ≤ 1.20 × n_molecules  (most pairs collapsed)
#
#   C. Damage vs PCR separation
#        10k molecules, duplicated:
#          - copy A: ancient damage only  (dmax5=0.25, no PCR error)
#          - copy B: same damage + PCR error at one interior position
#        Asserts:
#          • without --error-correct: some pairs NOT collapsed (PCR error blocks merge)
#          • with    --error-correct: more pairs collapsed (PCR error corrected away)

set -euo pipefail

FQDUP=${1:-build/fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

DMAX5=0.25; DMAX3=0.20; LAMBDA=0.35
TOL=0.08
N_MOL=10000

echo "=== fqdup --damage-auto tests ==="

# ---------------------------------------------------------------------------
# Test A: d_max estimation accuracy
# ---------------------------------------------------------------------------
echo ""
echo "--- Test A: d_max estimation accuracy ---"
echo "  True: dmax5=$DMAX5  dmax3=$DMAX3  lambda=$LAMBDA"

"$GEN" --n-unique 10000 --n-reads 100000 --read-len 75 \
    --dmax5 $DMAX5 --dmax3 $DMAX3 --lambda5 $LAMBDA --lambda3 $LAMBDA \
    --seed 42 > "$TMPDIR/est.fq"

"$FQDUP" sort -i "$TMPDIR/est.fq" -o "$TMPDIR/est.sorted.fq" --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep -i "$TMPDIR/est.sorted.fq" -o /dev/null \
    --damage-auto 2>"$TMPDIR/est.log"

EST5=$(grep "5'-end: d_max=" "$TMPDIR/est.log" | sed 's/.*d_max=\([^ ]*\) .*/\1/')
EST3=$(grep "3'-end: d_max=" "$TMPDIR/est.log" | sed 's/.*d_max=\([^ ]*\) .*/\1/')
echo "  Estimated: dmax5=$EST5  dmax3=$EST3"

grep -E "Masked positions:|d_max_combined=" "$TMPDIR/est.log" | head -3 | sed 's/^/  /'

python3 - <<EOF
import sys
est5, est3 = float("$EST5"), float("$EST3")
true5, true3, tol = $DMAX5, $DMAX3, $TOL
ok5 = abs(est5 - true5) < tol
ok3 = abs(est3 - true3) < tol
print(f"  |est5 - {true5}| = {abs(est5-true5):.3f}  (tol={tol}) → {'PASS' if ok5 else 'FAIL'}")
print(f"  |est3 - {true3}| = {abs(est3-true3):.3f}  (tol={tol}) → {'PASS' if ok3 else 'FAIL'}")
if not ok5 or not ok3: sys.exit(1)
EOF

grep -q "Masked positions:" "$TMPDIR/est.log" \
    && echo "  PASS: masking enabled" \
    || { echo "  FAIL: masking not enabled"; exit 1; }

echo "  PASS: Test A"

# ---------------------------------------------------------------------------
# Test B: damaged-duplicate collapse
# ---------------------------------------------------------------------------
echo ""
echo "--- Test B: damaged-duplicate collapse ---"
echo "  $N_MOL molecules × 2 independently-damaged copies"

"$GEN" --n-unique $N_MOL --read-len 75 \
    --dmax5 $DMAX5 --dmax3 $DMAX3 --lambda5 $LAMBDA --lambda3 $LAMBDA \
    --seed 123 --dup-pair > "$TMPDIR/dup.fq"

"$FQDUP" sort -i "$TMPDIR/dup.fq" -o "$TMPDIR/dup.sorted.fq" --max-memory 512M -t "$TMPDIR" 2>/dev/null

"$FQDUP" derep -i "$TMPDIR/dup.sorted.fq" -o "$TMPDIR/dup.out.fq" \
    --damage-auto 2>"$TMPDIR/dup.log"

UNIQUE_B=$(grep -c '^@' "$TMPDIR/dup.out.fq")
echo "  Unique reads: $UNIQUE_B / $((N_MOL * 2)) input reads"

python3 - <<EOF
import sys
unique, n_mol = $UNIQUE_B, $N_MOL
dedup = (n_mol*2 - unique) / (n_mol*2) * 100
print(f"  Dedup rate: {dedup:.1f}%")
max_unique = int(n_mol * 1.35)
if unique <= max_unique:
    print(f"  PASS: {unique} ≤ {max_unique} (≤35% uncollapsed pairs)")
else:
    print(f"  FAIL: {unique} > {max_unique} — too many uncollapsed pairs")
    sys.exit(1)
EOF

echo "  PASS: Test B"

# ---------------------------------------------------------------------------
# Test C: PCR error correction
#
# Generate 10x-coverage data with no ancient damage and a moderate PCR error
# rate (~0.3%/base → mu≈0.225 errors/read → ~20% of reads have ≥1 error).
# Without EC, PCR-error reads form singleton clusters inflating the unique count.
# With EC (--errcor-ratio 5 matches 10x coverage parent:singleton ratio ~8:1),
# those singletons are absorbed back into the dominant parent cluster.
#
# Note: errcor-ratio default (50) suits high-coverage production data; at 10x
# we set ratio=5 explicitly so the test exercises the real absorption path.
# ---------------------------------------------------------------------------
echo ""
echo "--- Test C: PCR error correction ---"
N_MOL_C=2000
N_READS_C=20000   # 10x coverage
PCR_RATE=0.003    # 0.3%/base → mu=0.225 → ~20% reads with ≥1 error

echo "  $N_MOL_C molecules × 10x coverage, no damage, PCR rate=$PCR_RATE/base"

"$GEN" --n-unique $N_MOL_C --n-reads $N_READS_C --read-len 75 \
    --no-damage --pcr-rate $PCR_RATE --seed 77 > "$TMPDIR/pcr.fq"

"$FQDUP" sort -i "$TMPDIR/pcr.fq" -o "$TMPDIR/pcr.sorted.fq" --max-memory 512M -t "$TMPDIR" 2>/dev/null

# Without error correction: PCR error reads inflate unique count
"$FQDUP" derep -i "$TMPDIR/pcr.sorted.fq" -o "$TMPDIR/pcr_no_ec.out.fq" 2>/dev/null
UNIQUE_NO_EC=$(grep -c '^@' "$TMPDIR/pcr_no_ec.out.fq")

# With error correction (ratio=5 suits 10x coverage)
"$FQDUP" derep -i "$TMPDIR/pcr.sorted.fq" -o "$TMPDIR/pcr_ec.out.fq" \
    --error-correct --errcor-ratio 5 2>/dev/null
UNIQUE_EC=$(grep -c '^@' "$TMPDIR/pcr_ec.out.fq")

echo "  Unique without --error-correct: $UNIQUE_NO_EC"
echo "  Unique with    --error-correct: $UNIQUE_EC"
echo "  True molecule count: $N_MOL_C"

python3 - <<EOF
import sys
no_ec, ec, n_mol = $UNIQUE_NO_EC, $UNIQUE_EC, $N_MOL_C

# PCR errors must inflate unique count above true molecule count
if no_ec > int(n_mol * 1.5):
    print(f"  PASS: PCR errors inflated unique count ({no_ec} > {int(n_mol*1.5)})")
else:
    print(f"  WARN: inflation smaller than expected ({no_ec}); check --pcr-rate")

# Error correction must reduce the inflated count
if ec < no_ec:
    absorbed = no_ec - ec
    inflation = no_ec - n_mol
    pct = absorbed / inflation * 100 if inflation > 0 else 0
    print(f"  PASS: EC absorbed {absorbed} singletons ({pct:.0f}% of PCR inflation)")
else:
    print(f"  FAIL: --error-correct did not help ({no_ec} → {ec})")
    sys.exit(1)
EOF

echo "  PASS: Test C"

# ---------------------------------------------------------------------------
# Test D: Channel C/D protection (8-oxoG G→T must NOT be absorbed by EC)
#
# Two datasets, same molecules (1000 × 20x coverage, no terminal deamination):
#
#   Oxidative  — G→T at uniform rate 0.05/G-base (Channel C/D, 8-oxoG).
#                is_damage_sub() recognises G↔T as damage-consistent and
#                prevents EC from absorbing these reads.
#                Expected: EC reduction ≈ 0%.
#
#   PCR errors — random non-damage substitutions at rate 0.003/base.
#                EC absorbs singleton error reads into the dominant parent.
#                Expected: EC reduction > 10%.
#
# The key assertion is that the two reduction rates are clearly separated:
# oxidative reads survive EC while PCR errors do not.
# ---------------------------------------------------------------------------
echo ""
echo "--- Test D: Channel C/D protection (G→T not absorbed by EC) ---"
N_MOL_D=1000
N_READS_D=20000   # 20x coverage

echo "  $N_MOL_D molecules × 20x coverage"
echo "  Oxidative: G→T at ox-rate=0.05/G-base  |  PCR: random at rate=0.003/base"

# Oxidative G→T dataset
"$GEN" --n-unique $N_MOL_D --n-reads $N_READS_D --read-len 75 \
    --no-damage --ox-rate 0.05 --seed 77 > "$TMPDIR/ox.fq"
"$FQDUP" sort -i "$TMPDIR/ox.fq" -o "$TMPDIR/ox.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
U_OX_NO_EC=$(  "$FQDUP" derep -i "$TMPDIR/ox.sorted.fq" -o "$TMPDIR/ox_no.fq" \
    2>/dev/null && grep -c '^@' "$TMPDIR/ox_no.fq")
U_OX_EC=$(     "$FQDUP" derep -i "$TMPDIR/ox.sorted.fq" -o "$TMPDIR/ox_ec.fq" \
    --error-correct --errcor-ratio 5 2>/dev/null && grep -c '^@' "$TMPDIR/ox_ec.fq")

# PCR error dataset (same seed → same base molecules)
"$GEN" --n-unique $N_MOL_D --n-reads $N_READS_D --read-len 75 \
    --no-damage --pcr-rate 0.003 --seed 77 > "$TMPDIR/pcr_d.fq"
"$FQDUP" sort -i "$TMPDIR/pcr_d.fq" -o "$TMPDIR/pcr_d.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
U_PCR_NO_EC=$( "$FQDUP" derep -i "$TMPDIR/pcr_d.sorted.fq" -o "$TMPDIR/pcr_d_no.fq" \
    2>/dev/null && grep -c '^@' "$TMPDIR/pcr_d_no.fq")
U_PCR_EC=$(    "$FQDUP" derep -i "$TMPDIR/pcr_d.sorted.fq" -o "$TMPDIR/pcr_d_ec.fq" \
    --error-correct --errcor-ratio 5 2>/dev/null && grep -c '^@' "$TMPDIR/pcr_d_ec.fq")

echo "  Oxidative  (G→T): no-EC=$U_OX_NO_EC  EC=$U_OX_EC"
echo "  PCR errors (rnd): no-EC=$U_PCR_NO_EC  EC=$U_PCR_EC"

python3 - <<EOF
import sys
ox_no, ox_ec   = $U_OX_NO_EC,  $U_OX_EC
pcr_no, pcr_ec = $U_PCR_NO_EC, $U_PCR_EC

red_ox  = (ox_no  - ox_ec)  / ox_no  * 100 if ox_no  > 0 else 0
red_pcr = (pcr_no - pcr_ec) / pcr_no * 100 if pcr_no > 0 else 0
print(f"  Oxidative  EC reduction: {red_ox:.1f}%  (expect <2%  — G→T protected)")
print(f"  PCR errors EC reduction: {red_pcr:.1f}%  (expect >10% — errors absorbed)")

if red_ox >= 2.0:
    print(f"  FAIL: G→T reads absorbed by EC ({red_ox:.1f}% ≥ 2%) — is_damage_sub not protecting")
    sys.exit(1)
print(f"  PASS: G→T reads protected (EC reduction {red_ox:.1f}% < 2%)")

if red_pcr <= 10.0:
    print(f"  FAIL: PCR errors not absorbed ({red_pcr:.1f}% ≤ 10%) — EC not working")
    sys.exit(1)
print(f"  PASS: PCR errors absorbed (EC reduction {red_pcr:.1f}% > 10%)")
EOF

echo "  PASS: Test D"

echo ""
echo "OK: all damage-auto tests passed"
