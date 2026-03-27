#!/usr/bin/env bash
# bench.sh — performance benchmarks for fqdup hot paths.
#
# Usage:
#   .scripts/bench.sh [fqdup_binary] [gen_synthetic_binary]
#
# Measures wall time and peak RSS for:
#   P1. derep Phase 3 on a high-depth library (stress bucket collisions)
#   P2. sort on large compressed input with many chunks
#   P3. derep_pairs on a large paired library
#   P4. extend on a large damaged library
#
# Results are printed in a table; no pass/fail assertions.
set -euo pipefail

FQDUP=${1:-build/fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# time_cmd label cmd... → prints wall time and peak RSS
time_cmd() {
    local label="$1"; shift
    local start end wall_ms rss_kb
    local out
    if command -v /usr/bin/time &>/dev/null; then
        out=$( { /usr/bin/time -v "$@" 2>&1 1>/dev/null; } 2>&1 || true )
        wall_ms=$(echo "$out" | grep "Elapsed (wall clock)" | grep -oP '\d+:\d+\.\d+' | \
            awk -F: '{printf "%d", ($1*60+$2)*1000}')
        rss_kb=$(echo "$out" | grep "Maximum resident" | grep -oP '\d+$')
        printf "  %-40s  %6s ms  %6s KB RSS\n" "$label" "${wall_ms:-?}" "${rss_kb:-?}"
    else
        start=$(date +%s%3N)
        "$@" 2>/dev/null
        end=$(date +%s%3N)
        printf "  %-40s  %6s ms  (RSS unavailable)\n" "$label" "$((end-start))"
    fi
}

echo "=== fqdup performance benchmarks ==="
echo ""

# ---------------------------------------------------------------------------
# P1: Phase 3 on high-depth library (stress bucket hashing)
# ---------------------------------------------------------------------------
echo "--- P1: derep Phase 3, high-depth (30x, 1000 molecules) ---"
"$GEN" --n-unique 1000 --n-reads 30000 --read-len 75 \
    --no-damage --seed 10 > "$TMPDIR/p1.fq" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p1.fq" -o "$TMPDIR/p1.sorted.fq" \
    --max-memory 512M -t "$TMPDIR" 2>/dev/null
time_cmd "derep --error-correct (1K mol × 30x)" \
    "$FQDUP" derep -i "$TMPDIR/p1.sorted.fq" -o "$TMPDIR/p1.out.fq"
echo ""

# ---------------------------------------------------------------------------
# P2: sort on compressed input, force many chunks via small memory budget
# ---------------------------------------------------------------------------
echo "--- P2: sort many-chunk merge (5K mol, 64M budget) ---"
"$GEN" --n-unique 5000 --n-reads 25000 --read-len 75 \
    --no-damage --seed 20 > "$TMPDIR/p2.fq" 2>/dev/null
gzip -c "$TMPDIR/p2.fq" > "$TMPDIR/p2.fq.gz"
time_cmd "sort .gz input, 64M budget" \
    "$FQDUP" sort -i "$TMPDIR/p2.fq.gz" -o "$TMPDIR/p2.sorted.fq.gz" \
        --max-memory 64M -t "$TMPDIR"
echo ""

# ---------------------------------------------------------------------------
# P3: derep_pairs on a larger paired library
# ---------------------------------------------------------------------------
echo "--- P3: derep_pairs, 5K molecules × 5x ---"
"$GEN" --n-unique 5000 --n-reads 25000 --read-len 75 \
    --dmax5 0.20 --dmax3 0.15 --lambda5 0.35 --lambda3 0.35 \
    --seed 30 > "$TMPDIR/p3.fq" 2>/dev/null
"$FQDUP" extend -i "$TMPDIR/p3.fq" -o "$TMPDIR/p3_ext.fq" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p3.fq"     -o "$TMPDIR/p3.sorted.fq"     --max-memory 512M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/p3_ext.fq" -o "$TMPDIR/p3_ext.sorted.fq" --max-memory 512M -t "$TMPDIR" 2>/dev/null
time_cmd "derep_pairs (5K mol × 5x)" \
    "$FQDUP" derep_pairs \
        -n "$TMPDIR/p3.sorted.fq" -e "$TMPDIR/p3_ext.sorted.fq" \
        -o-non "$TMPDIR/p3_dp.fq" -o-ext "$TMPDIR/p3_dp_ext.fq"
echo ""

# ---------------------------------------------------------------------------
# P4: extend on large damaged library
# ---------------------------------------------------------------------------
echo "--- P4: extend, 5K molecules × 5x, dmax5=0.25 ---"
time_cmd "extend (5K mol × 5x, damaged)" \
    "$FQDUP" extend -i "$TMPDIR/p3.fq" -o "$TMPDIR/p4_ext.fq"
echo ""

echo "=== Benchmark complete ==="
