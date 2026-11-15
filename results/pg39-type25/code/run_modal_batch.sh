#!/bin/bash
# Launch one Modal batch under a lease (client only; the containers run on Modal).
#   bash run_modal_batch.sh NAME IDS BUDGET_USD MAXC [ETA_MIN] [TLIM_SECONDS_OVERRIDE] [BIG=1: 8 GiB containers, cake_lpr heap 6500 MB]
# IDS: comma list / class prefix / @file (see modal_cubes.py::select).  Log: ../logs/modal_NAME.log  (the run refuses when metered + budget > 48 dollars).
set -u
NAME=$1; IDS=$2; BUDGET=$3; MAXC=$4; ETA=${5:-150}; TLIM=${6:-0}; BIG=${7:-0}
cd "$(dirname "$0")"
export PY=/c/Python314/python PYTHONUTF8=1 PYTHONIOENCODING=utf-8 NO_COLOR=1 PYTHONUNBUFFERED=1 PG39_MAXC=$MAXC
unset PYTHONPATH
nohup /c/Users/user/.compute-broker/lease-run.sh --project pg39 --min-gb 1 --want-gb 1 --minutes $((ETA * 2 + 60)) --eta-minutes $ETA \
  --why "pg39 Modal batch $NAME (approved by the coordinator; budget $BUDGET dollars, guard 48): containers run on Modal, this local client only waits and writes the ledger" \
  --match "modal" -- $PY -m modal run modal_cubes.py::main --ids "$IDS" --budget $BUDGET --note "$NAME" --tlim $TLIM $([ "$BIG" = "1" ] && echo --big) > ../logs/modal_$NAME.log 2>&1 &
echo "launched $NAME: ids=$IDS budget=$BUDGET maxc=$MAXC log=logs/modal_$NAME.log"
