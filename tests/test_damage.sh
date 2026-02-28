#!/usr/bin/env bash
# test_damage.sh — correctness tests for empirical damage-aware deduplication.
#
# Purpose of the method under test:
#   Before hashing, terminal positions with observed excess deamination are
#   replaced with a neutral byte (\x01 at C/T in the 5' zone, \x02 at G/A in
#   the 3' zone).  Two reads from the same original molecule that differ only
#   due to cytosine deamination (C→T at 5', G→A at 3') therefore produce
#   identical masked sequences and land in the same cluster.  Positions outside
#   the damage zone are untouched, so genuinely distinct molecules remain
#   separate.  The mask is applied symmetrically (same position indices from
#   both ends) so canonical_hash(seq) == canonical_hash(revcomp(seq)) holds
#   even after masking.
#
# Tests
#   1. Masked C→T must merge into 1 cluster.
#   2. Mutation outside damage zone must stay as 2 clusters.
#   3. A read and its exact reverse complement must merge (canonical hash).
#   3b. With --no-revcomp they must NOT merge.
#   4. Forward read and RC of its deaminated copy must merge (symmetry).
#
# Damage parameters used for all damage tests:
#   --damage-dmax5 0.2  --damage-lambda5 0.8  --damage-dmax3 0.0
#   populate_mask_from_model() yields:
#     pos 0: max(0.2*exp( 0), 0) = 0.200 > 0.05 → masked
#     pos 1: max(0.2*exp(-0.8), 0) = 0.090 > 0.05 → masked
#     pos 2: max(0.2*exp(-1.6), 0) = 0.040 < 0.05 → NOT masked
#   3' masking uses the same mask_pos array (symmetry), so positions 0 and 1
#   from the 3' end are also masked (for G/A bases).

set -euo pipefail

FQDUP=${1:-build/fqdup}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# 40 bp test sequence: position 0=G, 1=C (good for C→T deamination test),
# positions 3+ are well outside the damage zone.
BASE="GCAGATCGCGATCGATCGATCGATCGATCGATCGATCGAT"
QUAL=$(python3 -c "print('I' * ${#BASE})")

DMGFLAGS="--damage-dmax5 0.2 --damage-lambda5 0.8 --damage-dmax3 0.0"

# Helper: reverse complement a DNA sequence.
revcomp() { echo "$1" | tr 'ACGTacgt' 'TGCAtgca' | rev; }

# Helper: run derep and count output sequences.
count_output() {
    grep -c '^@' "$1"
}

check() {
    local label="$1" expected="$2" got="$3"
    if [ "$got" -eq "$expected" ]; then
        echo "  PASS: $label"
    else
        echo "  FAIL: $label — expected $expected unique reads, got $got"
        exit 1
    fi
}

echo "=== fqdup damage-aware deduplication tests ==="

# ---------------------------------------------------------------------------
# Test 1: C→T at a masked position must merge into 1 cluster.
#
# SEQ_A and SEQ_B differ ONLY at position 1 (C vs T).  Position 1 is in the
# damage zone (mask_pos[1] = true), so both C and T are replaced with \x01
# before hashing.  The masked sequences are identical → 1 cluster.
# ---------------------------------------------------------------------------
SEQ_A="$BASE"
SEQ_B="${BASE:0:1}T${BASE:2}"   # C→T at position 1

printf "@read_001\n%s\n+\n%s\n" "$SEQ_A" "$QUAL" >  "$TMPDIR/t1.fq"
printf "@read_002\n%s\n+\n%s\n" "$SEQ_B" "$QUAL" >> "$TMPDIR/t1.fq"

"$FQDUP" derep -i "$TMPDIR/t1.fq" -o "$TMPDIR/t1.out.fq" \
    $DMGFLAGS --no-revcomp 2>/dev/null

check "masked C→T at pos 1 merges into 1 cluster" 1 "$(count_output "$TMPDIR/t1.out.fq")"

# ---------------------------------------------------------------------------
# Test 2: Mutation at an unmasked position must stay as 2 clusters.
#
# SEQ_A and SEQ_C differ ONLY at position 5 (T vs G).  Position 5 is outside
# the damage zone (mask_pos[5] = false), so the bases are not replaced.  The
# masked sequences differ at position 5 → 2 clusters.
# ---------------------------------------------------------------------------
SEQ_C="${BASE:0:5}G${BASE:6}"   # T→G at position 5 (unmasked)

printf "@read_001\n%s\n+\n%s\n" "$SEQ_A" "$QUAL" >  "$TMPDIR/t2.fq"
printf "@read_002\n%s\n+\n%s\n" "$SEQ_C" "$QUAL" >> "$TMPDIR/t2.fq"

"$FQDUP" derep -i "$TMPDIR/t2.fq" -o "$TMPDIR/t2.out.fq" \
    $DMGFLAGS --no-revcomp 2>/dev/null

check "unmasked mutation at pos 5 stays as 2 clusters" 2 "$(count_output "$TMPDIR/t2.out.fq")"

# ---------------------------------------------------------------------------
# Test 3: A read and its exact reverse complement must merge.
#
# canonical_hash(seq) = min(h(seq), h(revcomp(seq))) = canonical_hash(revcomp(seq))
# so both reads map to the same fingerprint → 1 cluster.
# ---------------------------------------------------------------------------
RC_A=$(revcomp "$SEQ_A")

printf "@read_001\n%s\n+\n%s\n" "$SEQ_A" "$QUAL" >  "$TMPDIR/t3.fq"
printf "@read_002\n%s\n+\n%s\n" "$RC_A"  "$QUAL" >> "$TMPDIR/t3.fq"

"$FQDUP" derep -i "$TMPDIR/t3.fq" -o "$TMPDIR/t3.out.fq" 2>/dev/null

check "read and its exact revcomp merge into 1 cluster" 1 "$(count_output "$TMPDIR/t3.out.fq")"

# ---------------------------------------------------------------------------
# Test 3b: With --no-revcomp the same pair must stay as 2 clusters.
# ---------------------------------------------------------------------------
"$FQDUP" derep -i "$TMPDIR/t3.fq" -o "$TMPDIR/t3b.out.fq" \
    --no-revcomp 2>/dev/null

check "read and revcomp stay separate with --no-revcomp" 2 "$(count_output "$TMPDIR/t3b.out.fq")"

# ---------------------------------------------------------------------------
# Test 4: Symmetry invariant — forward read merges with RC of its deaminated copy.
#
# SEQ_A (C at pos 1) and RC(SEQ_B) (T at pos 1, then reverse-complemented)
# represent the same original molecule sequenced from opposite strands, where
# one copy carries a C→T deamination at the 5' end.
#
# Proof:
#   mask(SEQ_A)[1]  = \x01  (C in 5' zone)
#   mask(SEQ_B)[1]  = \x01  (T in 5' zone)
#   → h(mask(SEQ_A)) = h(mask(SEQ_B))
#
#   revcomp(SEQ_A) and revcomp(SEQ_B) differ at position L-2 (= pos 1 from 3'):
#   mask applies \x02 to G or A in the 3' zone; both G (from comp(C)) and
#   A (from comp(T)) are masked to \x02.
#   → h(mask(revcomp(SEQ_A))) = h(mask(revcomp(SEQ_B)))
#
#   canonical_hash(SEQ_A) = min(h(mask(SEQ_A)),     h(mask(revcomp(SEQ_A))))
#   canonical_hash(RC_B)  = min(h(mask(revcomp(SEQ_B))), h(mask(SEQ_B)))
#                         = min(h(mask(SEQ_A)),     h(mask(revcomp(SEQ_A))))  ✓
# ---------------------------------------------------------------------------
RC_B=$(revcomp "$SEQ_B")

printf "@read_001\n%s\n+\n%s\n" "$SEQ_A" "$QUAL" >  "$TMPDIR/t4.fq"
printf "@read_002\n%s\n+\n%s\n" "$RC_B"  "$QUAL" >> "$TMPDIR/t4.fq"

"$FQDUP" derep -i "$TMPDIR/t4.fq" -o "$TMPDIR/t4.out.fq" \
    $DMGFLAGS 2>/dev/null

check "forward + RC(deaminated) merge (symmetry invariant)" 1 "$(count_output "$TMPDIR/t4.out.fq")"

echo "OK: damage tests passed"
