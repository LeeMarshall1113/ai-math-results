#!/bin/bash
# usage (inside WSL): cake_check.sh x.cnf x.lrat [heapMB]
# Copies both files to WSL-local disk first (cake_lpr reads Windows drives very slowly), then checks.
set -u
# one cake_lpr at a time across all our batches (WSL is capped at 6 GB)
mkdir -p $HOME/erdos580; exec 9>$HOME/erdos580/cake.lock; flock 9
W=$HOME/erdos580/cake.$$; mkdir -p $W
cp "$1" $W/f.cnf && cp "$2" $W/f.lrat
${CAKE_LPR:-~/lrchk/bin/cake_lpr} --CML_HEAP_SIZE=${3:-2500} --CML_STACK_SIZE=${CAKE_STACK_MB:-1000} $W/f.cnf $W/f.lrat
rc=$?
rm -rf $W
exit $rc
