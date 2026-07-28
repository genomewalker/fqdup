#!/usr/bin/env bash
# test_merge_polyg.sh — regression guard for canonical poly-G QC in `fqdup merge`.
#
# What it locks down:
#   NovaSeq/NextSeq two-color chemistry writes long 3' poly-G run-through tails.
#   On long inserts that do NOT merge, those tails land on the true 3' terminus of
#   the unmerged mate outputs. `fqdup merge` is the single QC choke point: it must
#   strip that run-through (data-derived run-length) so downstream `fqdup profile`
#   receives clean reads and its 3' composition-artifact gate never trips.
#
# Check: synthesize long-insert pairs (insert > 2*readlen, so they cannot merge)
#   with heavy 3' poly-G on both mates. Run `fqdup merge`. The unmerged R1/R2
#   outputs must come out with ~0 trailing poly-G, even though the INPUT was heavy.

set -euo pipefail
FQDUP=${1:-build/fqdup}
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

python3 - "$TMPDIR" <<'PY'
import random, gzip, sys
random.seed(7)
out=sys.argv[1]
# Long inserts (MOL=150 > 2*RL): mates never overlap -> all become unmerged, so the
# 3' poly-G run-through sits on the real terminus of the R1-out / R2-out reads.
N=40000; MOL=150; RL=60
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
bases="ACGT"
r1=gzip.open(f"{out}/R1.fq.gz","wt"); r2=gzip.open(f"{out}/R2.fq.gz","wt")
for i in range(N):
    mol="".join(random.choice(bases) for _ in range(MOL))
    a=mol[:RL]                      # R1 = molecule 5' terminus
    b=rc(mol[-RL:])                 # R2 = molecule 3' (bottom strand)
    # 3' poly-G run-through on ~70% of reads, run length 20-45 (the artifact)
    if random.random()<0.70: a=a+"G"*random.randint(20,45)
    if random.random()<0.70: b=b+"G"*random.randint(20,45)
    r1.write(f"@m{i}\n{a}\n+\n{'I'*len(a)}\n")
    r2.write(f"@m{i}\n{b}\n+\n{'I'*len(b)}\n")
r1.close(); r2.close()
PY

# 3' poly-G rate (fraction of reads ending in a run of >=6 G) for a FASTQ.gz.
polyg_rate() {
  python3 -c "
import gzip,sys
n=g=0
with gzip.open(sys.argv[1],'rt') as fh:
    for i,l in enumerate(fh):
        if i%4==1:
            s=l.strip(); n+=1
            if s.endswith('GGGGGG'): g+=1
print(f'{g/n:.4f}' if n else '0.0')
" "$1"
}

in_r1=$(polyg_rate "$TMPDIR/R1.fq.gz")
in_r2=$(polyg_rate "$TMPDIR/R2.fq.gz")

"$FQDUP" merge -1 "$TMPDIR/R1.fq.gz" -2 "$TMPDIR/R2.fq.gz" \
    -o "$TMPDIR/merged.fq.gz" \
    --r1-out "$TMPDIR/unmR1.fq.gz" --r2-out "$TMPDIR/unmR2.fq.gz" \
    --internal-panel --min-length 25 -p 4 >/dev/null 2>&1

out_r1=$(polyg_rate "$TMPDIR/unmR1.fq.gz")
out_r2=$(polyg_rate "$TMPDIR/unmR2.fq.gz")

echo "input  3' poly-G:  R1=$in_r1  R2=$in_r2"
echo "output 3' poly-G:  unmR1=$out_r1  unmR2=$out_r2"

python3 - "$in_r1" "$in_r2" "$out_r1" "$out_r2" <<'PY'
import sys
i1,i2,o1,o2=map(float,sys.argv[1:5])
assert i1>0.3 and i2>0.3, f"fixture has no poly-G to strip (in R1={i1}, R2={i2}); test invalid"
assert o1<0.01, f"FAIL: unmerged R1 still carries 3' poly-G ({o1:.4f}); merge QC not applied"
assert o2<0.01, f"FAIL: unmerged R2 still carries 3' poly-G ({o2:.4f}); merge QC not applied"
print(f"PASS: merge stripped 3' poly-G  R1 {i1:.3f}->{o1:.3f}  R2 {i2:.3f}->{o2:.3f}")
PY
