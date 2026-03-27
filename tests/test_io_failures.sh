#!/usr/bin/env bash
# test_io_failures.sh — verify fqdup commands exit non-zero on I/O failures.
#
# Covers:
#   F1. derep: main output to /dev/full        → must fail
#   F2. derep: cluster output to /dev/full     → must fail
#   F3. derep_pairs: -o-non to /dev/full       → must fail
#   F4. sort: .gz output to /dev/full          → must fail
#   F5. derep: truncated gzip input            → must fail (not silent EOF)
#   F6. derep: truncated plain input           → must fail
#
# /dev/full is a Linux pseudo-device that always returns ENOSPC on write.
# Tests that need it are skipped on platforms without it.
set -euo pipefail

FQDUP=${1:-fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== fqdup I/O failure tests ==="

# ---------------------------------------------------------------------------
# Build a small valid sorted FASTQ for use in tests
# ---------------------------------------------------------------------------
"$GEN" --n-unique 10 --n-reads 30 --read-len 40 --no-damage --seed 1 \
    > "$TMPDIR/raw.fq" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/raw.fq" -o "$TMPDIR/sorted.fq" \
    --max-memory 64M -t "$TMPDIR" 2>/dev/null
"$FQDUP" sort -i "$TMPDIR/raw.fq" -o "$TMPDIR/sorted.fq.gz" \
    --max-memory 64M -t "$TMPDIR" 2>/dev/null

# Build a valid ext file (same as input for simplicity with short reads)
cp "$TMPDIR/sorted.fq" "$TMPDIR/sorted_ext.fq"

# ---------------------------------------------------------------------------
# Helper: assert a command exits non-zero
# ---------------------------------------------------------------------------
assert_fails() {
    local desc="$1"; shift
    if "$@" 2>/dev/null; then
        echo "  FAIL: $desc — expected non-zero exit but got 0"
        exit 1
    fi
    echo "  PASS: $desc"
}

# ---------------------------------------------------------------------------
# F1: derep main output to /dev/full
# ---------------------------------------------------------------------------
echo ""
echo "--- F1: derep main output to /dev/full ---"
if [[ -w /dev/full ]]; then
    assert_fails "derep -o /dev/full exits non-zero" \
        "$FQDUP" derep -i "$TMPDIR/sorted.fq" -o /dev/full
else
    echo "  SKIP: /dev/full not available"
fi

# ---------------------------------------------------------------------------
# F2: derep cluster output to /dev/full
# ---------------------------------------------------------------------------
echo ""
echo "--- F2: derep cluster output to /dev/full ---"
if [[ -w /dev/full ]]; then
    assert_fails "derep -c /dev/full exits non-zero" \
        "$FQDUP" derep -i "$TMPDIR/sorted.fq" -o "$TMPDIR/ok.fq" -c /dev/full
else
    echo "  SKIP: /dev/full not available"
fi

# ---------------------------------------------------------------------------
# F3: derep_pairs -o-non to /dev/full
# ---------------------------------------------------------------------------
echo ""
echo "--- F3: derep_pairs -o-non to /dev/full ---"
if [[ -w /dev/full ]]; then
    assert_fails "derep_pairs -o-non /dev/full exits non-zero" \
        "$FQDUP" derep_pairs \
            -n "$TMPDIR/sorted.fq" -e "$TMPDIR/sorted_ext.fq" \
            -o-non /dev/full -o-ext "$TMPDIR/ok_ext.fq"
else
    echo "  SKIP: /dev/full not available"
fi

# ---------------------------------------------------------------------------
# F4: sort .gz output to /dev/full
# ---------------------------------------------------------------------------
echo ""
echo "--- F4: sort .gz output to /dev/full ---"
if [[ -w /dev/full ]]; then
    assert_fails "sort -o /dev/full exits non-zero" \
        "$FQDUP" sort -i "$TMPDIR/raw.fq" -o /dev/full \
            --max-memory 64M -t "$TMPDIR"
else
    echo "  SKIP: /dev/full not available"
fi

# ---------------------------------------------------------------------------
# F5: derep with truncated gzip input
# ---------------------------------------------------------------------------
echo ""
echo "--- F5: derep truncated gzip input ---"
# Build a valid gzip, then truncate it
"$FQDUP" sort -i "$TMPDIR/raw.fq" -o "$TMPDIR/full.fq.gz" \
    --max-memory 64M -t "$TMPDIR" 2>/dev/null
FSIZE=$(wc -c < "$TMPDIR/full.fq.gz")
TRUNC=$(( FSIZE * 3 / 4 ))
dd if="$TMPDIR/full.fq.gz" of="$TMPDIR/trunc.fq.gz" \
    bs=1 count="$TRUNC" 2>/dev/null
assert_fails "derep truncated gzip exits non-zero" \
    "$FQDUP" derep -i "$TMPDIR/trunc.fq.gz" -o "$TMPDIR/trunc_out.fq"

# ---------------------------------------------------------------------------
# F6: derep with truncated plain input
# Sort a valid input, then truncate the sorted output mid-record.
# ---------------------------------------------------------------------------
echo ""
echo "--- F6: derep truncated plain input ---"
# Truncate the sorted file at 3/4 of its length (mid-record)
FSIZE=$(wc -c < "$TMPDIR/sorted.fq")
TRUNC=$(( FSIZE * 3 / 4 ))
dd if="$TMPDIR/sorted.fq" of="$TMPDIR/trunc_plain.fq" \
    bs=1 count="$TRUNC" 2>/dev/null
assert_fails "derep truncated plain FASTQ exits non-zero" \
    "$FQDUP" derep -i "$TMPDIR/trunc_plain.fq" -o "$TMPDIR/trunc_out2.fq"

echo ""
echo "=== All I/O failure tests passed ==="
