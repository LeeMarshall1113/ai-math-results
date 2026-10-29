#!/bin/bash
# One-command checker for implementation 2 (ref_level1.py + ref_lift.py + ref_gates.py).
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
  $PY ref_validate.py
  $PY ref_ml.py
fi
rm -f out/k6.jsonl out/k7.jsonl
$PY ref_gates.py 6 11 997 out/k6.jsonl > out/k6.log
$PY ref_gates.py 7 11 1153 out/k7.jsonl > out/k7.log
G=""; [ -f ../gates/gates_k6.json ] && G=../gates
$PY check_referee.py out/k6.jsonl out/k7.jsonl $G
