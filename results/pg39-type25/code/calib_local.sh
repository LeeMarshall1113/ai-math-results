#!/bin/bash
# Local reference run of the Modal calibration cube b_A1s1 (CaDiCaL 3.0.1 --congruence=false, binary LRAT, cake_lpr): peak RSS of both steps, sizes, sha256.
set -u
W=$HOME/pg39work/bench
cd $W
/usr/bin/time -v /home/user/ext/cadical/build/cadical -t 1500 --congruence=false --lrat=true --binary=true b_A1s1.cnf calib.lrat > calib.out 2> calib.time
echo "cadical rc=$?"
grep -E "c conflicts:|c total process time|c maximum resident" calib.out
grep -E "Maximum resident|Elapsed" calib.time
ls -l calib.lrat | awk '{print "lrat bytes", $5}'
sha256sum calib.lrat
sha256sum b_A1s1.cnf
/usr/bin/time -v /home/user/lrchk/bin/cake_lpr --CML_HEAP_SIZE=3000 --CML_STACK_SIZE=1000 b_A1s1.cnf calib.lrat > calib.cake 2> calib.caketime
echo "cake rc=$?"
tail -1 calib.cake
grep -E "Maximum resident|Elapsed" calib.caketime
rm -f calib.lrat calib.out calib.cake
