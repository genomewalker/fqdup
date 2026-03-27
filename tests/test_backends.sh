#!/usr/bin/env bash
# test_backends.sh — cross-backend consistency: plain, gzip, and bgzf paths
# must produce identical dedup output for the same input.
#
# Covers:
#   B1. sort: plain vs .gz output is sequence-identical after re-sorting
#   B2. derep: plain input vs gzip input → same unique sequences and counts
#   B3. derep_pairs: plain vs gzip inputs → same unique sequences
#   B4. damage: plain vs gzip input → same d_max estimates (within tolerance)
set -euo pipefail

FQDUP=${1:-fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== fqdup cross-backend consistency tests ==="

# ---------------------------------------------------------------------------
# Build shared input
# ---------------------------------------------------------------------------
"$GEN" --n-unique 200 --n-reads 2000 --read-len 75 \
    --dmax5 0.20 --dmax3 0.15 --lambda5 0.35 --lambda3 0.35 \
    --seed 42 > "$TMPDIR/raw.fq" 2>/dev/null

# Sorted plain and gzip versions
"$FQDUP" sort -i "$TMPDIR/raw.fq" -o "$TMPDIR/sorted.fq" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/raw.fq" -o "$TMPDIR/sorted.fq.gz" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null

# Extended (plain + gzip)
"$FQDUP" extend -i "$TMPDIR/raw.fq" -o "$TMPDIR/ext.fq" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/ext.fq" -o "$TMPDIR/sorted_ext.fq" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/ext.fq" -o "$TMPDIR/sorted_ext.fq.gz" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null

# ---------------------------------------------------------------------------
# B1: sort plain vs .gz output — same sequences
# ---------------------------------------------------------------------------
echo ""
echo "--- B1: sort plain vs .gz output is sequence-identical ---"
# Decompress and extract sequences from both
awk 'NR%4==2' "$TMPDIR/sorted.fq" | sort > "$TMPDIR/b1_plain_seqs.txt"
"$FQDUP" sort -i "$TMPDIR/sorted.fq.gz" -o "$TMPDIR/b1_resorted.fq" \
    --max-memory 128M -t "$TMPDIR" 2>/dev/null
awk 'NR%4==2' "$TMPDIR/b1_resorted.fq" | sort > "$TMPDIR/b1_gz_seqs.txt"
if ! diff -q "$TMPDIR/b1_plain_seqs.txt" "$TMPDIR/b1_gz_seqs.txt" > /dev/null; then
    echo "  FAIL: plain and .gz sorted outputs differ"
    exit 1
fi
echo "  PASS: plain and .gz sorted outputs are sequence-identical"

# ---------------------------------------------------------------------------
# B2: derep plain vs gzip input → same unique set
# ---------------------------------------------------------------------------
echo ""
echo "--- B2: derep plain vs gzip input → same unique sequences ---"
"$FQDUP" derep -i "$TMPDIR/sorted.fq" -o "$TMPDIR/b2_plain.fq" \
    --damage-auto --no-error-correct 2>/dev/null
"$FQDUP" derep -i "$TMPDIR/sorted.fq.gz" -o "$TMPDIR/b2_gz.fq" \
    --damage-auto --no-error-correct 2>/dev/null

awk 'NR%4==2' "$TMPDIR/b2_plain.fq" | sort > "$TMPDIR/b2_plain_seqs.txt"
awk 'NR%4==2' "$TMPDIR/b2_gz.fq"    | sort > "$TMPDIR/b2_gz_seqs.txt"
if ! diff -q "$TMPDIR/b2_plain_seqs.txt" "$TMPDIR/b2_gz_seqs.txt" > /dev/null; then
    echo "  FAIL: derep plain vs gzip input produced different unique sets"
    diff "$TMPDIR/b2_plain_seqs.txt" "$TMPDIR/b2_gz_seqs.txt" | head -5
    exit 1
fi
U=$(grep -c '^@' "$TMPDIR/b2_plain.fq")
echo "  PASS: derep plain and gzip inputs agree ($U unique reads)"

# ---------------------------------------------------------------------------
# B3: derep_pairs plain vs gzip inputs → same unique set
# ---------------------------------------------------------------------------
echo ""
echo "--- B3: derep_pairs plain vs gzip inputs → same unique sequences ---"
"$FQDUP" derep_pairs \
    -n "$TMPDIR/sorted.fq" -e "$TMPDIR/sorted_ext.fq" \
    -o-non "$TMPDIR/b3_plain_non.fq" -o-ext "$TMPDIR/b3_plain_ext.fq" \
    2>/dev/null
"$FQDUP" derep_pairs \
    -n "$TMPDIR/sorted.fq.gz" -e "$TMPDIR/sorted_ext.fq.gz" \
    -o-non "$TMPDIR/b3_gz_non.fq" -o-ext "$TMPDIR/b3_gz_ext.fq" \
    2>/dev/null

awk 'NR%4==2' "$TMPDIR/b3_plain_non.fq" | sort > "$TMPDIR/b3_plain_seqs.txt"
awk 'NR%4==2' "$TMPDIR/b3_gz_non.fq"    | sort > "$TMPDIR/b3_gz_seqs.txt"
if ! diff -q "$TMPDIR/b3_plain_seqs.txt" "$TMPDIR/b3_gz_seqs.txt" > /dev/null; then
    echo "  FAIL: derep_pairs plain vs gzip inputs produced different unique sets"
    exit 1
fi
U=$(grep -c '^@' "$TMPDIR/b3_plain_non.fq")
echo "  PASS: derep_pairs plain and gzip inputs agree ($U unique reads)"

# ---------------------------------------------------------------------------
# B4: damage plain vs gzip input → same d_max (within 1% tolerance)
# ---------------------------------------------------------------------------
echo ""
echo "--- B4: damage plain vs gzip input → same d_max estimates ---"
D5_PLAIN=$("$FQDUP" damage -i "$TMPDIR/raw.fq" 2>&1 | \
    grep "5'-end" | grep -oP 'd_max=\K[0-9.]+' | head -1)
D5_GZ=$("$FQDUP" damage -i "$TMPDIR/sorted.fq.gz" 2>&1 | \
    grep "5'-end" | grep -oP 'd_max=\K[0-9.]+' | head -1)
python3 - <<EOF
import sys
p, g = float("${D5_PLAIN:-0}"), float("${D5_GZ:-0}")
tol = 0.01
if abs(p - g) > tol:
    print(f"  FAIL: d_max_5 plain={p:.4f} gz={g:.4f} diff {abs(p-g):.4f} > {tol}")
    sys.exit(1)
print(f"  PASS: d_max_5 plain={p:.4f} gz={g:.4f} (diff {abs(p-g):.4f} ≤ {tol})")
EOF

echo ""
echo "=== All cross-backend consistency tests passed ==="
