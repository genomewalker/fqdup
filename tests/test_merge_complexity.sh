#!/usr/bin/env bash
# test_merge_complexity.sh — regression guard for the data-derived low-complexity gate
# in `fqdup merge`.
#
# What it locks down:
#   Two-color chemistry emits whole-read poly-G run-through mates that carry NO insert,
#   so the 3' poly-G tail-trim cannot rescue them (nothing to keep past the run). Left in
#   the unmerged stream they corrupt the 3' (R2) damage channel and inflate a blank's d5.
#   `fqdup merge` learns a per-mate worst-window entropy/dominance reference from the
#   overlap-verified MERGED inserts (no hardcoded threshold) and drops unmerged mates more
#   extreme than that reference's 1%/99% tails.
#
# Check: a fixture with (a) enough clean mergeable short pairs to derive the gate,
#   (b) clean long-insert pairs that do NOT merge (must SURVIVE), and (c) junk pairs whose
#   R2 is a whole-read poly-G tract (must be DROPPED). After merge, the unmerged output must
#   be ~free of low-complexity reads while retaining the clean long-insert population.

set -euo pipefail
FQDUP=${1:-build/fqdup}
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

python3 - "$TMPDIR" <<'PY'
import random, gzip, sys
random.seed(11)
out=sys.argv[1]
bases="ACGT"
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
A1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"   # TruSeq R1 read-through (detection anchor)
r1=gzip.open(f"{out}/R1.fq.gz","wt"); r2=gzip.open(f"{out}/R2.fq.gz","wt")
def emit(a,b):
    r1.write(f"@r\n{a}\n+\n{'I'*len(a)}\n"); r2.write(f"@r\n{b}\n+\n{'I'*len(b)}\n")

# (a) 4000 clean SHORT-insert pairs -> overlap+adapter => MERGE => reference population.
for _ in range(4000):
    ins="".join(random.choice(bases) for _ in range(40))
    emit(ins+A1, rc(ins)+rc(A1))
# (b) 3000 clean LONG-insert pairs (insert>2*RL) -> no overlap => clean UNMERGED (survive).
for _ in range(3000):
    mol="".join(random.choice(bases) for _ in range(150)); RL=60
    emit(mol[:RL], rc(mol[-RL:]))
# (c) 3000 junk pairs: R1 clean-long, R2 = whole-read poly-G (no insert, tail-trim can't help).
for _ in range(3000):
    mol="".join(random.choice(bases) for _ in range(150)); RL=60
    emit(mol[:RL], "G"*random.randint(55,70))
# Interior-indel adapter dimer: R1 = adapter read-through with a single inserted base in the
# AGATCGGAAGAG G-run — defeats BOTH exact 5' anchors (the rigid Hamming is_adapter_dimer_5p and
# the is_adapter_fragment fast-reject 11-mer). Two deterministic streams exercise both routes;
# the adapter body must NOT survive in ANY output stream (merged/unmR1/unmR2).
A1_INDEL="AGATCGG"+"G"+A1[7:]        # 1bp insertion at pos 7, interior of the read-through prefix
# (d) 2000 UNMERGED: R2 is one FIXED clean 60-mer with no complementarity to R1, so the pair
#     never overlaps -> unmerged path. is_adapter_dimer_5p_indel must reject R1 (no unmR1 leak).
R2_FIXED="".join(random.choice(bases) for _ in range(60))   # drawn once, reused (deterministic under seed)
for _ in range(2000):
    emit(A1_INDEL, R2_FIXED)
# (e) 2000 MERGED: R2 = rc(R1) -> full overlap -> merges to ~pure adapter. is_technical_read
#     (is_adapter_fragment, Myers) must drop it; the 5' 11-mer carries the indel, so the
#     fast-reject must probe past it (middle/3' seed).
for _ in range(2000):
    emit(A1_INDEL, rc(A1_INDEL))
r1.close(); r2.close()
PY

# fraction of reads whose worst 20bp window is >90% one base (low-complexity junk)
lowcomplex_rate() {
  python3 -c "
import gzip,sys
def worst(s,W=20):
    if len(s)<=W:
        c=max(s.count(b) for b in 'ACGT'); return c/len(s) if s else 1.0
    m=0.0
    for i in range(len(s)-W+1):
        w=s[i:i+W]; c=max(w.count(b) for b in 'ACGT')
        if c/W>m: m=c/W
    return m
n=j=0
with gzip.open(sys.argv[1],'rt') as fh:
    for i,l in enumerate(fh):
        if i%4==1:
            s=l.strip(); n+=1
            if worst(s)>0.9: j+=1
print(f'{j/n:.4f}' if n else '0.0', n)
" "$1"
}

read in_lc  in_n  < <(lowcomplex_rate "$TMPDIR/R2.fq.gz")

"$FQDUP" merge -1 "$TMPDIR/R1.fq.gz" -2 "$TMPDIR/R2.fq.gz" \
    -o "$TMPDIR/merged.fq.gz" \
    --r1-out "$TMPDIR/unmR1.fq.gz" --r2-out "$TMPDIR/unmR2.fq.gz" \
    --internal-panel --min-length 25 -p 4 2>"$TMPDIR/merge.err" >/dev/null

grep -i "complexity-gate" "$TMPDIR/merge.err" || true

read out_lc out_n2 < <(lowcomplex_rate "$TMPDIR/unmR2.fq.gz")
read out_lc1 out_n1 < <(lowcomplex_rate "$TMPDIR/unmR1.fq.gz")

echo "input  R2 low-complexity: $in_lc  (n=$in_n)"
echo "output unmR2 low-complexity: $out_lc (n=$out_n2)   unmR1: $out_lc1 (n=$out_n1)"

# The gate must actually be derived (reference over floor) or the test is not exercising it.
if ! grep -qiE "complexity-gate R1\[n=[0-9]+ ent>=[0-9].* dom<=0\.[0-9]" "$TMPDIR/merge.err"; then
  echo "FAIL: complexity gate was not derived (reference below floor?)"; cat "$TMPDIR/merge.err"; exit 1
fi

python3 - "$in_lc" "$out_lc" "$out_n1" <<'PY'
import sys
in_lc, out_lc, out_n1 = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
assert in_lc >= 0.18, f"fixture R2 has too few poly-G junk reads ({in_lc}); test invalid"
assert out_lc < 0.02, f"FAIL: unmerged R2 still carries low-complexity junk ({out_lc:.4f}); gate not applied"
# The 3000 clean long-insert PAIRS emit paired unmerged reads (orphans are dropped without
# --orphan-out). The one-sided gate must retain the bulk of that clean population, not eat it.
assert out_n1 > 2000, f"FAIL: clean unmerged reads over-rejected (unmR1 n={out_n1}); gate too aggressive"
print(f"PASS: gate dropped R2 low-complexity {in_lc:.3f}->{out_lc:.3f}, kept {out_n1} clean unmerged R1")
PY

# --- interior-indel adapter-dimer leak guard (stream d) ---
# Adapter body past the AGATCGGAAGAG prefix — present only in a LEAKED dimer, never in a real
# insert (stream (a) merges and gets adapter-trimmed; (b)/(c) carry no adapter). Any occurrence
# in output means an interior-indel dimer escaped is_adapter_dimer_5p_indel (or the 3' trim).
BODY="CACACGTCTGAACTCC"
leak=0
for f in "$TMPDIR/merged.fq.gz" "$TMPDIR/unmR1.fq.gz" "$TMPDIR/unmR2.fq.gz"; do
  c=$(zcat "$f" 2>/dev/null | awk 'NR%4==2' | grep -c "$BODY" || true)
  leak=$((leak + c))
done
echo "interior-indel dimer adapter-body leak across outputs: $leak"
if [ "$leak" -ne 0 ]; then
  echo "FAIL: $leak interior-indel adapter-dimer(s) leaked into output (is_adapter_dimer_5p_indel regressed)"
  exit 1
fi
echo "PASS: no interior-indel adapter-dimer leaked into any output stream"
