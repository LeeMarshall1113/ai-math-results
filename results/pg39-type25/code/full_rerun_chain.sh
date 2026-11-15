#!/bin/bash
# Chains the slices of the full local regeneration (see run_full_rerun.sh) until every stage list is done.
# Start once, from Git Bash:   nohup bash code/full_rerun_chain.sh > logs/full_rerun_chain.out 2>&1 &
# It first waits for an already queued/running first slice (lease id given as $1, optional) to leave the broker,
# then runs stage 1 (3 workers), stage 2 (2 workers), stage 3 (1 worker) in leased, resumable 330-minute slices.
# Stops (exit 12) on any disagreement with the cloud record or a cake_lpr failure; exit 13 after two slices with no progress.
PG=/c/Users/user/Documents/GitHub/erdos-triage/pg39
B=/c/Users/user/.compute-broker/broker.py
LEASE=/c/Users/user/.compute-broker/lease-run.sh
export PY=/c/Python314/python
export PYTHONPATH="$($PY -c 'import site; print(site.getusersitepackages())')"
LOG=$PG/logs/full_rerun_chain.log
say() { echo "$(date '+%F %T') $*" >> "$LOG"; }

state() {   # prints: remaining bad   for stage $1
  $PY - "$PG" "$1" <<'EOF'
import json, os, sys
pg, st = sys.argv[1], sys.argv[2]
ids = [j["id"] for j in json.load(open(f"{pg}/results/full_rerun_stage{st}.json"))]
done, bad = set(), 0
p = f"{pg}/results/full_rerun_dry.jsonl"
if os.path.exists(p):
    for line in open(p):
        line = line.strip()
        if not line: continue
        try: r = json.loads(line)
        except Exception: continue
        s = r.get("status", "")
        if s == "UNSAT_VERIFIED" and "flag" not in r: done.add(r.get("tag"))
        elif s in ("DISAGREE_CNF", "UNSAT_CAKE_FAILED", "ERROR", "SAT_EXACT_PASS", "SAT_EXACT_FAIL") or "flag" in r: bad += 1
        elif s == "TIMEOUT_LOCAL": done.add(r.get("tag"))      # recorded, not retried
print(len([i for i in ids if i not in done]), bad)
EOF
}

if [ -n "$1" ]; then
  say "waiting for lease $1 to leave the broker"
  while $PY $B status 2>/dev/null | grep -q "$1"; do sleep 120; done
  say "lease $1 gone"
fi

for spec in "1 3" "2 2" "3 1"; do
  set -- $spec; ST=$1; NW=$2; idle=0
  while true; do
    read REM BAD <<< "$(state $ST)"
    say "stage $ST: remaining $REM, bad records $BAD"
    [ "${BAD:-0}" -gt 0 ] && { say "STOP: disagreement or checker failure recorded; inspect results/full_rerun_dry.jsonl"; exit 12; }
    [ "${REM:-1}" -eq 0 ] && break
    bash "$LEASE" --project pg39 --min-gb 5 --minutes 480 --eta-minutes 360 --preemptible \
      --match "a_slice_driver|cadical|cake_lpr" \
      --why "pg39: full local regeneration of the remaining cube and child certificates (stage $ST)" \
      --argument "Closes the referee's last open item: every refutation regenerated and checked locally, bit-compared with the cloud record. Resumable 330-minute slices, $NW WSL workers, checkpoint per certificate, yields at any certificate boundary, can wait." \
      -- bash "$PG/code/run_full_rerun.sh" $ST $NW 330 >> "$PG/logs/full_rerun_stage$ST.log" 2>&1
    rc=$?
    read REM2 BAD2 <<< "$(state $ST)"
    say "stage $ST slice ended rc=$rc: remaining $REM -> $REM2"
    if [ "$REM2" -ge "$REM" ]; then idle=$((idle+1)); [ $idle -ge 2 ] && { say "STOP: two slices without progress"; exit 13; }; sleep 300; else idle=0; sleep 30; fi
  done
done
say "ALL STAGES DONE"
