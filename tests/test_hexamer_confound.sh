#!/usr/bin/env bash
# test_hexamer_confound.sh — validates CircLigase 3'-hexamer-bias confound detection.
#
# The confound fires when the 3'-terminal hexamer distribution is RC-asymmetric
# (rc_overlap_topk == 0) and the dominant 3' hexamers are not damage-consistent
# (top_damage_consistent_fraction_3prime < 0.5).  Under these conditions fqdup
# annotates d_max_3 as confounded in diagnostic_groups.adapter_position_effects.
#
# Tests:
#   SS arm:  5'-only deamination (--ss --dmax5 0.05) + 70% CircLigase CCC hexamers
#            at 3' end (--hex3-bias 0.7).
#            Expects: hexamer_end_rc_excess_jsd > 0.3,  rc_overlap_topk == 0,
#                     confounded_outputs contains "d_max_3".
#
#   DS arm:  symmetric deamination (--dmax5 0.05 --dmax3 0.04), no hexamer bias.
#            Expects: confounded_outputs empty, d_max_3prime > 0.02.
#            (Hexamer JSD / rc_overlap_topk not checked: synthetic reads lack RC-pairing.)

set -euo pipefail

FQDUP=${1:-build/fqdup}
GEN=${2:-build/gen_synthetic}

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

COMMON="--n-unique 80000 --n-reads 200000 --read-len 50 --seed 2025"

echo "=== hexamer confound: SS library with CircLigase 3'-bias ==="

"$GEN" $COMMON --ss --dmax5 0.05 --hex3-bias 0.7 --hex5-bias 0.8 > "$TMPDIR/ss.fq" 2>/dev/null
"$FQDUP" profile -i "$TMPDIR/ss.fq" --json "$TMPDIR/ss.json" --library-type ss -p 2 2>/dev/null

python3 - "$TMPDIR/ss.json" <<'PYEOF'
import json, sys
d  = json.load(open(sys.argv[1]))
dg = d['library_qc']['diagnostic_groups']
eh = dg['end_hexamer_asymmetry']
ap = dg['adapter_position_effects']
dm = d['deamination']

failures = []

jsd = eh['hexamer_end_rc_excess_jsd']
if jsd <= 0.30:
    failures.append(f"hexamer_end_rc_excess_jsd={jsd:.4f} want >0.30")

rctopk = eh['rc_overlap_topk']
if rctopk != 0:
    failures.append(f"rc_overlap_topk={rctopk} want 0")

conf = ap.get('confounded_outputs', [])
if 'd_max_3' not in conf:
    failures.append(f"confounded_outputs={conf} want 'd_max_3'")

d3 = dm['d_max_3prime']
d5 = dm['d_max_5prime']
dc = dm['d_max_combined']
if abs(dc - d5) > 0.005:
    failures.append(f"d_max_combined={dc:.4f} != d_max_5prime={d5:.4f} (delta {abs(dc-d5):.4f} > 0.005)")

if failures:
    for f in failures:
        print(f"  FAIL: {f}")
    sys.exit(1)

print(f"  PASS: jsd={jsd:.3f}  rc_overlap_topk={rctopk}  confounded=['d_max_3']"
      f"  d_max_combined={dc:.4f}=d_max_5prime={d5:.4f}")
PYEOF

echo "=== hexamer confound: DS library control (no hexamer bias) ==="
# Note: synthetic random reads have no RC-paired structure, so hexamer_end_rc_excess_jsd
# and rc_overlap_topk are not meaningful for a purely synthetic DS control.
# The important property is that the confound does NOT fire when there is no hex bias,
# because fit_offset_3prime stays at 0 (no position-0 spike at 3' end).

"$GEN" $COMMON --dmax5 0.05 --dmax3 0.04 > "$TMPDIR/ds.fq" 2>/dev/null
"$FQDUP" profile -i "$TMPDIR/ds.fq" --json "$TMPDIR/ds.json" --library-type ds -p 2 2>/dev/null

python3 - "$TMPDIR/ds.json" <<'PYEOF'
import json, sys
d  = json.load(open(sys.argv[1]))
dg = d['library_qc']['diagnostic_groups']
ap = dg['adapter_position_effects']
dm = d['deamination']

failures = []

conf = ap.get('confounded_outputs', [])
if 'd_max_3' in conf:
    failures.append(f"confounded_outputs={conf} should not contain 'd_max_3'")

d3 = dm['d_max_3prime']
if d3 < 0.02:
    failures.append(f"d_max_3prime={d3:.4f} want >0.02 (DS 3' damage should be present)")

if failures:
    for f in failures:
        print(f"  FAIL: {f}")
    sys.exit(1)

print(f"  PASS: confounded_outputs=[]  d_max_3prime={d3:.4f}")
PYEOF

echo "=== DONE ==="
