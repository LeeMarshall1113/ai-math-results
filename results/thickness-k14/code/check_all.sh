#!/bin/bash
# One-command check. Run from anywhere:  bash code/check_all.sh
#
# Default (fast, a few minutes, < 1 GB RAM, no SAT solver binary needed):
#   certificate file hashes; upper bound (exact 4-layer drawings of K14); soundness controls for both
#   encoders; completeness of all three case trees with every leaf CNF regenerated and hash-compared.
# FULL=1: additionally re-certify every leaf: regenerate CNF, compare hash, solve with CaDiCaL (LRAT),
#   check the LRAT, re-hash it, delete it. One leaf at a time; see README for time and disk.
#
# Environment:
#   PY=python3            Python 3 with python-sat (pysat) for the controls
#   CADICAL=cadical       CaDiCaL binary (FULL=1 only). Tested with 3.0.1.
#   CADICAL_OPTS="--congruence=false"   needed for CaDiCaL 3.0.1 (rc 3 abort in congruence closure)
#   LRAT_CHECK=cake_lpr   called as: $LRAT_CHECK file.cnf file.lrat; must print "s VERIFIED UNSAT".
#                         Larger leaves need e.g. "cake_lpr --CML_HEAP_SIZE=2500 --CML_STACK_SIZE=1000"
#   TREES="sweep nosweep referee"   which trees FULL=1 replays
#   SAMPLE=k SEED=s       FULL=1 on k random leaves per tree instead of all of them
#   CONTROLS_K13=1        also run the K13/3 SAT control (about 5 min)
#   TK_WORK=../work       scratch directory (CNF and LRAT of the current leaf)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
export CADICAL=${CADICAL:-cadical}
export LRAT_CHECK=${LRAT_CHECK:-cake_lpr}
export TK_WORK=${TK_WORK:-$(cd .. && pwd)/work}
mkdir -p "$TK_WORK"

echo "== certificate files: SHA-256 against certs/CHECKSUMS.sha256"
(cd ../certs && sha256sum -c --quiet CHECKSUMS.sha256) && echo "certs OK"

echo "== upper bound: exact 4-layer drawings of K14 (referee construction + 4 drawings from the authors' models)"
$PY upper.py

echo "== soundness controls, authors' encoder: real drawings satisfy the CNF of their own hull case"
for w in k13_seamons k14_t4 k14_t4_S7 k14_t4_S5-10 k14_t4_S4-8-11 k14_t4_referee test_k10_h4 test_k10_h5 test_k10_h6; do
  printf '%s: ' "$w"; $PY check_known_prod.py ../certs/witnesses/$w.json | tail -1
done

echo "== soundness controls, referee encoder"
$PY controls.py "$( [ -n "${CONTROLS_K13:-}" ] && echo table13 || echo table )"
$PY controls.py real | grep -v "^random K"

echo "== completeness of the case trees; every leaf CNF regenerated and compared with the published SHA-256"
$PY cover_check.py 14 3 ../certs/tree_authors_sweep.jsonl --regen
$PY cover_check.py 14 3 ../certs/tree_authors_nosweep.jsonl --regen --nosweep
$PY cover_mine.py 14 3 ../certs/tree_referee.jsonl --regen

if [ -n "${FULL:-}" ]; then
  echo "== FULL: re-certify leaves (CADICAL=$CADICAL, LRAT_CHECK=$LRAT_CHECK)"
  for tr in ${TREES:-sweep nosweep referee}; do
    if [ -n "${SAMPLE:-}" ]; then
      $PY replay.py $tr --sample "$SAMPLE" --seed "${SEED:-1}" --out "$TK_WORK/replay_${tr}_sample.jsonl"
    else
      $PY replay.py $tr --out "$TK_WORK/replay_$tr.jsonl"
      case $tr in
        sweep)   $PY cover_check.py 14 3 ../certs/tree_authors_sweep.jsonl --replay "$TK_WORK/replay_sweep.jsonl" ;;
        nosweep) $PY cover_check.py 14 3 ../certs/tree_authors_nosweep.jsonl --nosweep --replay "$TK_WORK/replay_nosweep.jsonl" ;;
        referee) $PY cover_mine.py 14 3 ../certs/tree_referee.jsonl --replay "$TK_WORK/replay_referee.jsonl" ;;
      esac
    fi
  done
fi
echo "ALL CHECKS DONE"
