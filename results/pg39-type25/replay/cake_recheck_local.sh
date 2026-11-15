#!/bin/bash
# re-run cake_lpr on the 19 cover and 4 planar certificates in pg39/certs
CERTS=/mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/certs
CAKE=/home/user/lrchk/bin/cake_lpr
W=$HOME/pg39work/recheck; mkdir -p $W
echo "cake_lpr sha256: $(sha256sum $CAKE | cut -d' ' -f1)"
echo "start: $(date -u +%FT%TZ)"
tags="planar_complete_k26_R0 planar_complete_k26_R1 planar_complete_k35_R0 planar_complete_k35_R1 cover1_A0 cover1_A1 cover1_A2 cover1_C0 cover1_C1 cover1_C2 cover1_C3 cover1_C4 cover1_C5 cover1_C6 cover2_A1_c04 cover2_A2_c01 cover2_A2_c04 cover2_A2_c05 cover2_A2_c17 cover2_C0_c15 cover2_C1_c04 cover2_C2_c00 cover2_C6_c08"
n=0; ok=0
for t in $tags; do
  n=$((n+1))
  cp $CERTS/$t.cnf $W/x.cnf; cp $CERTS/$t.lrat $W/x.lrat
  cs=$(sha256sum $W/x.cnf | cut -d' ' -f1); ls=$(sha256sum $W/x.lrat | cut -d' ' -f1); lb=$(stat -c %s $W/x.lrat)
  t0=$(date +%s)
  out=$($CAKE --CML_HEAP_SIZE=2000 --CML_STACK_SIZE=1000 $W/x.cnf $W/x.lrat 2>&1 | tail -1)
  if [ "$out" != "s VERIFIED UNSAT" ]; then
    first="$out"
    out=$($CAKE --CML_HEAP_SIZE=3000 --CML_STACK_SIZE=1000 $W/x.cnf $W/x.lrat 2>&1 | tail -1)
    echo "  (first try at heap 2000: $first)"
  fi
  t1=$(date +%s)
  [ "$out" = "s VERIFIED UNSAT" ] && ok=$((ok+1))
  echo "$t cnf_sha256=$cs lrat_sha256=$ls lrat_bytes=$lb cake_seconds=$((t1-t0)) result=[$out]"
  rm -f $W/x.cnf $W/x.lrat
done
echo "SUMMARY: $ok of $n verified"
echo "end: $(date -u +%FT%TZ)"
