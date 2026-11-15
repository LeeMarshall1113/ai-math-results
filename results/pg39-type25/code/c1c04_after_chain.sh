#!/bin/bash
# Waits for the full-regeneration chain to finish (logs/full_rerun_chain.log: "ALL STAGES DONE"), then runs the prepared
# local split of C1_c04 (code/run_c1c04.sh) under its own broker lease, in resumable slices.  Exits without running if the chain stopped on an error.
# Start once from Git Bash:  nohup bash code/c1c04_after_chain.sh > logs/c1c04_after_chain.out 2>&1 &
PG=/c/Users/user/Documents/GitHub/erdos-triage/pg39
export PY=/c/Python314/python
export PYTHONPATH="$($PY -c 'import site; print(site.getusersitepackages())')"
LOG=$PG/logs/c1c04_after_chain.log
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }
say "waiting for the regeneration chain"
while true; do
  grep -q "ALL STAGES DONE" "$PG/logs/full_rerun_chain.log" 2>/dev/null && break
  grep -q "STOP:" "$PG/logs/full_rerun_chain.log" 2>/dev/null && { say "chain stopped on an error; not running the split"; exit 12; }
  sleep 300
done
say "chain done; starting the C1_c04 split"
cd "$PG" || exit 1
for i in 1 2 3 4; do
  bash /c/Users/user/.compute-broker/lease-run.sh --project pg39 --min-gb 3 --minutes 600 --eta-minutes 300 --preemptible \
    --match "c1c04|cadical|cake_lpr" --why "pg39: local split of cube C1_c04 (over the WSL cap as a monolith)" \
    --argument "Last certificate not yet verified locally: the cube is split into children with a completeness certificate, one solver or checker at a time inside the 6 GB WSL cap; resumable, yields at any child boundary, can wait." \
    -- bash code/run_c1c04.sh >> logs/run_c1c04.log 2>&1
  rc=$?; say "slice $i rc=$rc"
  tail -1 logs/run_c1c04.log | grep -q "driver rc=0" && { say "driver finished cleanly"; break; }
  sleep 120
done
say "done (run code/c1c04_finish.py after checking results/c1c04_manifest_local.jsonl)"
