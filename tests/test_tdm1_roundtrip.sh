#!/usr/bin/env bash
# fqdup profile writes a TDM1 damage summary; libtaph's reader (via tdm1_dump) must read it back as
# version 3 with the 24 fields DART consumes, and a second identical run must write the same bytes.
set -euo pipefail

FQDUP=$1; GEN=$2; DUMP=$3
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

"$GEN" --n-unique 20000 --n-reads 40000 --read-len 75 \
    --dmax5 0.25 --dmax3 0.20 --lambda5 0.35 --lambda3 0.35 --seed 7 > "$TMPDIR/reads.fq"

for run in a b; do
    "$FQDUP" profile -i "$TMPDIR/reads.fq" -p 2 --library-type ds \
        --damage-summary "$TMPDIR/$run.tdm1" > "$TMPDIR/$run.log" 2>&1 \
        || { echo "profile failed:"; cat "$TMPDIR/$run.log"; exit 1; }
done

cmp "$TMPDIR/a.tdm1" "$TMPDIR/b.tdm1" || { echo "TDM1 not deterministic"; exit 1; }
"$DUMP" "$TMPDIR/a.tdm1" > "$TMPDIR/a.tsv"
test "$(wc -l < "$TMPDIR/a.tsv")" = 24 || { echo "expected 24 fields"; cat "$TMPDIR/a.tsv"; exit 1; }
grep -qx $'version\t3' "$TMPDIR/a.tsv" || { echo "not TDM1 version 3"; exit 1; }
grep -qx $'library_type\t2' "$TMPDIR/a.tsv" || { echo "library_type is not ds"; exit 1; }
awk -F'\t' '$1=="d_max_5prime" && !($2 > 0) {exit 1}' "$TMPDIR/a.tsv" || { echo "d_max_5prime not positive"; exit 1; }
echo "TDM1_ROUNDTRIP_PASS"
