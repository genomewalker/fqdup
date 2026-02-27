#!/usr/bin/env bash
# Smoke test for fqdup sort and fqdup derep
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

# ---- fqdup derep ----
"$FQDUP" derep \
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

echo "OK: fqdup smoke test passed"
