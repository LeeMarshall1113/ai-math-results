#!/bin/bash
# LOCAL SPLIT OF C1_c04: the single resumable command.  Runs INSIDE a broker lease, from Git Bash, e.g.
#   export PY=/c/Python314/python; export PYTHONPATH="$(/c/Python314/python -c 'import site; print(site.getusersitepackages())')"
#   /c/Users/user/.compute-broker/lease-run.sh --project pg39 --min-gb 3 --minutes 600 --eta-minutes <n> --preemptible --match "c1c04|cadical|cake_lpr" \
#       --why "pg39: local split of cube C1_c04 (over the WSL cap as a monolith)" --argument "..." -- bash code/run_c1c04.sh > logs/run_c1c04.log 2>&1
# Steps (each skipped when already done): 1 split (python+PySAT, Windows side, minutes)  2 job list  3 completeness certificate cover2_C1_c04 (CaDiCaL LRAT + cake_lpr, WSL)
# 4 children (one solver or checker at a time, WSL, memory-gated: no solver starts below 2500 MB MemAvailable).  Needs about 2.5-4.5 GB of the WSL VM when it runs: run it
# when the regeneration (3 workers) is not running.  Stops at the broker preempt flag.  Finish with: python code/c1c04_finish.py
export MSYS_NO_PATHCONV=1
PG=/c/Users/user/Documents/GitHub/erdos-triage/pg39
PY=${PY:-/c/Python314/python}
cd $PG/code || exit 1
WSLPG=/mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39
echo "start $(date '+%F %T')"
wsl -d Ubuntu-24.04 -e bash -c 'sha256sum /home/user/ext/cadical/build/cadical /home/user/lrchk/bin/cake_lpr; grep -E "MemTotal|MemAvailable|SwapFree" /proc/meminfo'
# expected: cadical f30ecc0fc740ea7a..., cake_lpr 272a75de5b8add99...
[ -f ../certs/split_C1_c04.json ] || $PY -u split_level2.py C1_c04 8 300 || exit 2
$PY c1c04_prep.py || exit 3
if ! grep -q '"tag": "cover2_C1_c04"' ../certs/manifest.jsonl; then
  $PY -u cover_cert.py level2 C1_c04 || exit 4
fi
PF=""; [ -n "${BROKER_PREEMPT_FILE:-}" ] && PF=$(wsl -d Ubuntu-24.04 -e wslpath -u "$(cygpath -w "$BROKER_PREEMPT_FILE")" | tr -d '\r')
wsl -d Ubuntu-24.04 -e python3 -u $WSLPG/code/c1c04_driver.py $WSLPG/results/c1c04_jobs.json --out c1c04_local --workers 1 --cap-gb 6 --tlim 7200 --min-free-gb 20 ${PF:+--preempt-file "$PF"}
echo "driver rc=$? $(date '+%F %T')"
