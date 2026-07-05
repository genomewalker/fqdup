#!/bin/bash
set -u
WT=/projects/caeg/scratch/kbd606/tmp/fqdup-join-wt
PERF=/maps/projects/fernandezguerra/apps/opt/conda/envs/bioinfo/bin/perf
IN=/projects/caeg/scratch/kbd606/tmp/zihao_kapk_full/sort/sorted.fq.gz
RUN=$WT/profile_run
mkdir -p "$RUN"
LOG=$RUN/derep.log
CPU=$RUN/cpu_samples.log
: > "$LOG"; : > "$CPU"

numactl --cpunodebind=1 --membind=1 "$WT/build-join/fqdup" derep \
  -i "$IN" -o "$RUN/dedup.fq.gz" -c "$RUN/clusters.tsv.gz" \
  --damage report --library-type ds -p 16 >> "$LOG" 2>&1 &
DPID=$!
echo "derep pid=$DPID"

# Continuous CPU/thread sampler (serial => ~100%, parallel => ~1600%)
( while kill -0 "$DPID" 2>/dev/null; do
    read -r cpu nlwp <<<"$(ps -o %cpu=,nlwp= -p "$DPID" 2>/dev/null)"
    tail1=$(tr '\r' '\n' < "$LOG" | grep -aE 'INFO:|κ|empirical model|Phase 3 complete' | tail -1 | cut -c1-70)
    echo "$(date +%H:%M:%S) cpu=${cpu:-NA} nlwp=${nlwp:-NA} | ${tail1}" >> "$CPU"
    sleep 3
  done ) &
SPID=$!

wait_marker() { local pat="$1"; while kill -0 "$DPID" 2>/dev/null; do
    grep -qaE "$pat" "$LOG" && return 0; sleep 2; done; return 1; }

# FIT window: begins when kappa line prints
if wait_marker 'interior κ'; then
  echo "$(date +%H:%M:%S) >> FIT window: perf record 45s" >> "$CPU"
  "$PERF" record -F 199 -g --call-graph dwarf -o "$RUN/fit.perf.data" -p "$DPID" -- sleep 45 2>>"$RUN/perf.err"
fi
# B2 window: begins when empirical-model line prints (fit done)
if wait_marker 'empirical model:'; then
  echo "$(date +%H:%M:%S) >> B2 window: perf record 45s" >> "$CPU"
  "$PERF" record -F 199 -g --call-graph dwarf -o "$RUN/b2.perf.data" -p "$DPID" -- sleep 45 2>>"$RUN/perf.err"
fi

wait "$DPID"; kill "$SPID" 2>/dev/null
echo "$(date +%H:%M:%S) >> derep exited" >> "$CPU"

for w in fit b2; do
  if [ -s "$RUN/$w.perf.data" ]; then
    echo "===== $w perf report (self-time top) =====" >> "$RUN/perf_report.txt"
    "$PERF" report -i "$RUN/$w.perf.data" --stdio --no-children -g none 2>/dev/null | \
      grep -vE '^#|^$' | head -35 >> "$RUN/perf_report.txt"
    echo >> "$RUN/perf_report.txt"
  fi
done
echo "PROFILE_DRIVER_DONE"
