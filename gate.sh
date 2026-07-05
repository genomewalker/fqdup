#!/bin/bash
set -u
WT=/projects/caeg/scratch/kbd606/tmp/fqdup-join-wt
IN=/projects/caeg/scratch/kbd606/tmp/zihao_kapk_full/sort/sorted.fq.gz
STAGE="${1:-X}"
RUN=$WT/gate_$STAGE
mkdir -p "$RUN"
LOG=$RUN/derep.log
: > "$LOG"

numactl --cpunodebind=1 --membind=1 "$WT/build-join/fqdup" derep \
  -i "$IN" -o "$RUN/dedup.fq.gz" -c "$RUN/clusters.tsv.gz" \
  --damage report --library-type ds -p 16 >> "$LOG" 2>&1 &
DPID=$!
PEAK=0
while kill -0 "$DPID" 2>/dev/null; do
  hwm=$(awk '/VmHWM/{print $2}' /proc/$DPID/status 2>/dev/null)
  [ -n "$hwm" ] && [ "$hwm" -gt "$PEAK" ] 2>/dev/null && PEAK=$hwm
  sleep 5
done
wait "$DPID"; rc=$?
echo "derep rc=$rc"
echo "PEAK_RSS_GB=$(awk -v k="$PEAK" 'BEGIN{printf "%.1f", k/1048576}')"

FA=$(zcat "$RUN/dedup.fq.gz" | md5sum | cut -d' ' -f1)
ST=$(zcat "$RUN/clusters.tsv.gz" | md5sum | cut -d' ' -f1)
echo "FASTQ_MD5=$FA  (want 0cd05a06316110e79064f152a1b47a41)"
echo "STATS_MD5=$ST  (want eab4283454b322a6bb436fa6a62cc7e7)"
[ "$FA" = "0cd05a06316110e79064f152a1b47a41" ] && echo "FASTQ MATCH" || echo "FASTQ MISMATCH"
[ "$ST" = "eab4283454b322a6bb436fa6a62cc7e7" ] && echo "STATS MATCH" || echo "STATS MISMATCH"

# window walls from log timestamps (epoch via date)
ep() { date -d "$(echo "$1" | sed -E 's/^\[([0-9-]+ [0-9:]+)\].*/\1/')" +%s 2>/dev/null; }
k=$(grep -aE 'interior κ' "$LOG" | tail -1)
m=$(grep -aE 'empirical model:' "$LOG" | tail -1)
c=$(grep -aE 'Phase 3 complete' "$LOG" | tail -1)
p3=$(grep -aE 'Phase timer: Phase 3' "$LOG" | tail -1)
[ -n "$k" ] && [ -n "$m" ] && echo "FIT window wall = $(( $(ep "$m") - $(ep "$k") )) s"
[ -n "$m" ] && [ -n "$c" ] && echo "B2 window wall  = $(( $(ep "$c") - $(ep "$m") )) s"
echo "$p3"
echo "GATE_${STAGE}_DONE rc=$rc"
