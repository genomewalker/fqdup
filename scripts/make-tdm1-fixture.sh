#!/bin/bash
# Produce the TDM1 fixture DART's reader test compares against: <out_dir>/fqdup_sample.tdm1 and
# <out_dir>/fqdup_sample.expected.tsv (every DamageSummary field as field<TAB>value, read back with
# libtaph's own reader via tdm1_dump). fqdup merges the first FQDUP_FIXTURE_PAIRS read pairs, then
# profile runs on the merged reads with length bins (raw-pair input leaves the pi fields at -1).
# Temporary files go to a private directory under <out_dir> and are removed afterwards.
# Usage: scripts/make-tdm1-fixture.sh <fqdup_bin> <out_dir>
# Environment: FQDUP_FIXTURE_R1 / FQDUP_FIXTURE_R2 (paired FASTQ, default the real Ellesmere ds
# subset), FQDUP_FIXTURE_PAIRS (default 1000000), FQDUP_FIXTURE_LIBRARY (default ds),
# FQDUP_FIXTURE_ADAPTER1 / FQDUP_FIXTURE_ADAPTER2, TDM1_DUMP (default: tdm1_dump next to fqdup).
set -euo pipefail
[ $# -eq 2 ] || { echo "usage: $0 <fqdup_bin> <out_dir>" >&2; exit 2; }
fqdup=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); out=$2
refs=/path/to/tmp/fqdup-speed/task27-refs/inputs
r1=${FQDUP_FIXTURE_R1:-$refs/ds_R1.fastq}; r2=${FQDUP_FIXTURE_R2:-$refs/ds_R2.fastq}
pairs=${FQDUP_FIXTURE_PAIRS:-1000000}; lib=${FQDUP_FIXTURE_LIBRARY:-ds}
a1=${FQDUP_FIXTURE_ADAPTER1:-AGATCGGAAGAGCACACGTCTGAACTCCAG}
a2=${FQDUP_FIXTURE_ADAPTER2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAG}
dump=${TDM1_DUMP:-$(dirname "$fqdup")/tdm1_dump}
[ -x "$fqdup" ] && [ -x "$dump" ] || { echo "need executable $fqdup and $dump" >&2; exit 2; }
mkdir -p "$out"
tmp=$(mktemp -d "$out/.fixture-XXXXXX")
trap 'rm -rf "$tmp"' EXIT
head -n $((pairs * 4)) "$r1" > "$tmp/R1.fq"
head -n $((pairs * 4)) "$r2" > "$tmp/R2.fq"
"$fqdup" --version > "$out/fqdup_version.txt"
"$fqdup" merge -1 "$tmp/R1.fq" -2 "$tmp/R2.fq" -o "$tmp/merge" --cutoff 0.316 -p "${FQDUP_JOBS:-16}" \
    --library-type "$lib" --adapter1 "$a1" --adapter2 "$a2" --outcomes none > "$tmp/merge.log" 2>&1 \
    || { tail -20 "$tmp/merge.log"; exit 1; }
for run in a b; do
    "$fqdup" profile -i "$tmp/merge/merged.fq" -p "${FQDUP_JOBS:-16}" --library-type "$lib" --length-bins auto \
        --damage-summary "$tmp/$run.tdm1" > "$tmp/profile-$run.log" 2>&1 || { tail -20 "$tmp/profile-$run.log"; exit 1; }
done
cmp -s "$tmp/a.tdm1" "$tmp/b.tdm1" || { echo "fqdup profile is not deterministic on this input" >&2; exit 1; }
cp "$tmp/a.tdm1" "$out/fqdup_sample.tdm1"
"$dump" "$out/fqdup_sample.tdm1" > "$out/fqdup_sample.expected.tsv"
(cd "$out" && sha256sum fqdup_sample.tdm1 fqdup_sample.expected.tsv > SHA256SUMS)
echo "fixture: $out/fqdup_sample.tdm1 $out/fqdup_sample.expected.tsv ($(head -1 "$out/fqdup_version.txt"))"
