#!/bin/bash
# Run inside WSL. Converts Zeraoulia's DRUP traces to LRAT with drat-trim -L,
# then checks each LRAT with cake_lpr (CakeML-verified) and gate_lrat_check.
set -u
Z="$(cd "$(dirname "$0")" && pwd)/zenodo"
W=${W:-$HOME/erdos580/work}
DT=${DRAT_TRIM:-$HOME/erdos580/drat-trim/drat-trim}
CAKE=${CAKE_LPR:-~/lrchk/bin/cake_lpr}
GATE=${GATE_LRAT:-~/lrchk/bin/gate_lrat_check}
mkdir -p $W
for s in ${@:-path_support_014 path_support_0124 spider_224_core spider_233_core}; do
  cp $Z/$s.cnf $W/$s.cnf
  t0=$(date +%s)
  $DT $W/$s.cnf $Z/$s.drup -L $W/$s.lrat -t 20000 > $W/$s.dt.log 2>&1
  t1=$(date +%s)
  echo "$s drat-trim: $(grep -aE '^s ' $W/$s.dt.log | tr '\n' ' ') ($((t1-t0)) s)"
  flock $HOME/erdos580/cake.lock $CAKE --CML_HEAP_SIZE=3000 --CML_STACK_SIZE=1000 $W/$s.cnf $W/$s.lrat > $W/$s.cake.log 2>&1
  t2=$(date +%s)
  echo "$s cake_lpr: $(tr -d '\0' < $W/$s.cake.log | tail -1) ($((t2-t1)) s)"
  $GATE $W/$s.cnf $W/$s.lrat > $W/$s.gate.log 2>&1; rc=$?
  echo "$s gate_lrat_check: $(tail -1 $W/$s.gate.log) rc=$rc"
  ls -la $W/$s.lrat | awk '{print "  lrat bytes", $5}'
done
