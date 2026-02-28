#!/usr/bin/env bash
# Smoke test for fqdup sort, derep_pairs, and derep
set -euo pipefail

FQDUP=${1:-fqdup}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
cd "$TMPDIR"

# ---- minimal paired FASTQ ----
cat > non.fq <<'EOF'
@read_002
CCCCCC
+
IIIIII
@read_001
AAAAAA
+
IIIIII
@read_003
TTTTTTTT
+
IIIIIIII
@read_001
AAAA
+
IIII
EOF

cat > ext.fq <<'EOF'
@read_002
GGGGGG
+
IIIIII
@read_001
AAAAAA
+
IIIIII
@read_003
TTTTTTTT
+
IIIIIIII
@read_001
AAAAAA
+
IIIIII
EOF

# ---- fqdup sort ----
"$FQDUP" sort -i non.fq -o non.sorted.fq --max-memory 1G -t . --fast
"$FQDUP" sort -i ext.fq -o ext.sorted.fq --max-memory 1G -t . --fast

if [[ ! -s non.sorted.fq ]]; then
    echo "FAIL: non.sorted.fq missing or empty"
    exit 1
fi

# Verify sorted order: read_001 should come before read_002
first=$(awk 'NR==1{print $0}' non.sorted.fq)
if [[ "$first" != "@read_001" ]]; then
    echo "FAIL: sort order wrong, first header is '$first' (expected @read_001)"
    exit 1
fi

# ---- fqdup derep_pairs ----
"$FQDUP" derep_pairs \
    -n non.sorted.fq -e ext.sorted.fq \
    -o-non out.non.fq -o-ext out.ext.fq

if [[ ! -s out.non.fq ]]; then
    echo "FAIL: out.non.fq missing or empty"
    exit 1
fi
if [[ ! -s out.ext.fq ]]; then
    echo "FAIL: out.ext.fq missing or empty"
    exit 1
fi

# Verify dedup: read_001 appears once and longest rep (AAAAAA) was chosen
count=$(grep -c '^@read_001' out.non.fq || true)
if [[ "$count" -ne 1 ]]; then
    echo "FAIL: read_001 appears $count times in output (expected 1)"
    exit 1
fi

seq=$(awk '/^@read_001/{getline; print; exit}' out.non.fq)
if [[ "$seq" != "AAAAAA" ]]; then
    echo "FAIL: longest rep for read_001 should be AAAAAA, got '$seq'"
    exit 1
fi

n_pairs=$(grep -c '^@' out.non.fq || true)
if [[ "$n_pairs" -ne 3 ]]; then
    echo "FAIL: expected 3 unique reads in out.non.fq, got $n_pairs"
    exit 1
fi

# Verify ext output is in sync: same count, read_001 ext seq is AAAAAA
n_ext=$(grep -c '^@' out.ext.fq || true)
if [[ "$n_ext" -ne 3 ]]; then
    echo "FAIL: expected 3 unique reads in out.ext.fq, got $n_ext"
    exit 1
fi
ext_seq=$(awk '/^@read_001/{getline; print; exit}' out.ext.fq)
if [[ "$ext_seq" != "AAAAAA" ]]; then
    echo "FAIL: ext mate for read_001 should be AAAAAA, got '$ext_seq'"
    exit 1
fi

# ---- fqdup derep (single-file) ----
# Run on the derep_pairs non-extended output — all 3 sequences are distinct,
# so all 3 should pass through unchanged.
"$FQDUP" derep -i out.non.fq -o non.final.fq

if [[ ! -s non.final.fq ]]; then
    echo "FAIL: non.final.fq missing or empty"
    exit 1
fi

n_final=$(grep -c '^@' non.final.fq || true)
if [[ "$n_final" -ne 3 ]]; then
    echo "FAIL: expected 3 reads in non.final.fq, got $n_final"
    exit 1
fi

# ---- RC deduplication test ----
# AAAACCCC and its RC GGGGTTTT must collapse to 1 unique (default --revcomp on),
# but stay as 2 unique with --no-revcomp.
cat > rc_test.fq <<'EOF'
@fwd_001
AAAACCCC
+
IIIIIIII
@rev_001
GGGGTTTT
+
IIIIIIII
EOF

"$FQDUP" derep -i rc_test.fq -o rc_merged.fq
rc_merged=$(grep -c '^@' rc_merged.fq || true)
if [[ "$rc_merged" -ne 1 ]]; then
    echo "FAIL: RC dedup: expected 1 unique read, got $rc_merged"
    exit 1
fi

"$FQDUP" derep -i rc_test.fq -o rc_unmerged.fq --no-revcomp
rc_unmerged=$(grep -c '^@' rc_unmerged.fq || true)
if [[ "$rc_unmerged" -ne 2 ]]; then
    echo "FAIL: RC disabled: expected 2 unique reads, got $rc_unmerged"
    exit 1
fi

# ---- error case: unequal-length paired files must fail ----
cat > short_non.fq <<'EOF'
@read_001
AAAAAA
+
IIIIII
EOF

if "$FQDUP" derep_pairs \
        -n non.sorted.fq -e short_non.fq \
        -o-non /dev/null -o-ext /dev/null 2>/dev/null; then
    echo "FAIL: derep_pairs should reject unequal-length input files"
    exit 1
fi

# ---- error case: ID mismatch in paired files must fail ----
# ext file sorted with a different ID at position 2 (read_003 vs read_002)
cat > ext_wrong_ids.fq <<'EOF'
@read_001
AAAAAA
+
IIIIII
@read_003
TTTTTTTT
+
IIIIIIII
EOF

if "$FQDUP" derep_pairs \
        -n short_non.fq -e ext_wrong_ids.fq \
        -o-non /dev/null -o-ext /dev/null 2>/dev/null; then
    echo "FAIL: derep_pairs should reject mismatched read IDs"
    exit 1
fi

echo "OK: fqdup smoke test passed"
