#!/bin/bash
# One-off (2026-11-13): the first 8 second planes all have more than 300 sections for C1_c04, so split_level2.py found no split.
# Waits for c1c04_after_chain.sh (PID $1) to finish its retries, then repeats the split over 24 planes with cap 1200 and runs
# the prepared resumable driver (code/run_c1c04.sh) in leased slices, as c1c04_after_chain.sh does.
PG=/c/Users/user/Documents/GitHub/erdos-triage/pg39
export PY=/c/Python314/python
export PYTHONPATH="$($PY -c 'import site; print(site.getusersitepackages())')"
LOG=$PG/logs/c1c04_after_chain.log
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }
cd "$PG" || exit 1
n=0
while tasklist //FI "PID eq $1" 2>/dev/null | grep -q " $1 "; do sleep 60; n=$((n+1)); [ $n -ge 240 ] && exit 2; done
L=/c/Users/user/.compute-broker/lease-run.sh
if [ ! -f certs/split_C1_c04.json ]; then
  say "wide split: 24 planes, cap 1200"
  bash $L --project pg39 --min-gb 2 --minutes 240 --eta-minutes 90 --preemptible --match "split_level2" \
    --why "pg39: C1_c04 second-plane split, wider search (24 planes, cap 1200)" \
    --argument "The first 8 planes gave no usable split; one python process with PySAT enumerates sections of more planes. Small memory, can wait." \
    -- bash -c "cd code && timeout 12600 $PY -u split_level2.py C1_c04 24 1200" >> logs/run_c1c04.log 2>&1
  say "wide split rc=$?"
  [ -f certs/split_C1_c04.json ] || { say "STOP: still no split plane; see results/a_slice_splits.jsonl"; exit 12; }
fi
for i in 1 2 3 4 5 6; do
  bash $L --project pg39 --min-gb 3 --minutes 600 --eta-minutes 300 --preemptible \
    --match "c1c04|cadical|cake_lpr" --why "pg39: local split of cube C1_c04 (over the WSL cap as a monolith)" \
    --argument "Last certificate not yet verified locally: the cube is split into children with a completeness certificate, one solver or checker at a time inside the 6 GB WSL cap; resumable, yields at any child boundary, can wait." \
    -- bash code/run_c1c04.sh >> logs/run_c1c04.log 2>&1
  rc=$?; say "wide slice $i rc=$rc"
  tail -1 logs/run_c1c04.log | grep -q "driver rc=0" && { say "driver finished cleanly"; break; }
  [ $rc = 2 ] || [ $rc = 3 ] || [ $rc = 4 ] && { say "STOP: setup step failed rc=$rc"; exit 12; }
  sleep 120
done
say "wide: done (run code/c1c04_finish.py after checking results/c1c04_manifest_local.jsonl)"
