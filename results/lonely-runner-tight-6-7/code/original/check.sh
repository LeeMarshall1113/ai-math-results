#!/bin/bash
# One-command checker for implementation 1 (lr.py + tight6.py + tight7.py).
#   bash check.sh              sieves k=6 (p=11..997) and k=7 (p=11..1153), then the
#                              exact bookkeeping and a comparison with ../gates/*.json
#   bash check.sh --validate   additionally runs the brute-force cross-checks first
# One thread throughout. Set PY to choose the interpreter (default: python).
set -e
cd "$(dirname "$0")"
PY=${PY:-python}
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out
if [ "$1" = "--validate" ]; then
  $PY tight6.py validate
  $PY tight6.py validate2
  $PY validate7.py 11 13 17
  $PY kravitz_IIb.py
fi
rm -f out/k6.json out/k7.json
$PY tight6.py sieve 997 out/k6.json > out/k6.log
$PY tight7.py 11 1153 out/k7.json > out/k7.log
G=""; [ -f ../gates/gates_k6.json ] && G=../gates
$PY check_original.py out/k6.json out/k7.json $G
