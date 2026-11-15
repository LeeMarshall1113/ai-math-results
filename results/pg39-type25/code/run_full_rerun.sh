#!/bin/bash
# Full local regeneration of the remaining certificates (see code/full_rerun_prep.py).  Runs INSIDE the broker lease.
#   run_full_rerun.sh STAGE WORKERS SLICE_MINUTES
# Records go to results/full_rerun_dry.jsonl (--dry: certs/manifest.jsonl is NOT touched).  Resumable: completed tags are skipped.
# The driver stops launching at a broker preempt flag, after SLICE_MINUTES, when free space on C: < 40 GB, or on any disagreement / cake failure.
export MSYS_NO_PATHCONV=1
STAGE=${1:-1}; NW=${2:-3}; SL=${3:-330}
P=/mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39
PF=""
if [ -n "${BROKER_PREEMPT_FILE:-}" ]; then PF=$(wsl -d Ubuntu-24.04 -e wslpath -u "$(cygpath -w "$BROKER_PREEMPT_FILE")" | tr -d '\r'); fi
cd /c/Users/user/Documents/GitHub/erdos-triage/pg39/code
wsl -d Ubuntu-24.04 -e python3 -u $P/code/a_slice_driver.py $P/results/full_rerun_stage$STAGE.json --dry --out full_rerun --workers $NW --cap-gb 12 --min-free-gb 40 --stop-after-min $SL ${PF:+--preempt-file "$PF"}
