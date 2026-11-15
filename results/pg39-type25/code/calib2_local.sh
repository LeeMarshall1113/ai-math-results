#!/bin/bash
# Second local calibration run: determinism (same LRAT sha256 as calib_local.sh) and minimum cake_lpr heap.
set -u
W=$HOME/pg39work/bench
cd $W
/usr/bin/time -v /home/user/ext/cadical/build/cadical -t 1500 --congruence=false --lrat=true --binary=true b_A1s1.cnf calib.lrat > calib.out 2> calib.time
echo "cadical rc=$?"
grep -E "c conflicts:|c total process time" calib.out
grep -E "Maximum resident|Elapsed" calib.time
ls -l calib.lrat | awk '{print "lrat bytes", $5}'
sha256sum calib.lrat
for H in 2200 1600; do
  /usr/bin/time -v /home/user/lrchk/bin/cake_lpr --CML_HEAP_SIZE=$H --CML_STACK_SIZE=1000 b_A1s1.cnf calib.lrat > calib.cake 2> calib.caketime
  echo "cake heap $H rc=$?"
  tail -2 calib.cake | cut -c1-200
  grep -E "Maximum resident|Elapsed" calib.caketime
done
rm -f calib.lrat calib.out calib.cake
