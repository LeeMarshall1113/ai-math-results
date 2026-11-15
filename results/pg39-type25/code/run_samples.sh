#!/bin/bash
# three depth-1 cube samplers in parallel (one lease for all three); each solver process is single-threaded
cd /c/Users/user/Documents/GitHub/erdos-triage/pg39/code
PY=/c/Python314/python
$PY -u trivial_tree.py A 0 sample 3 1 3600 20000 1 > ../logs/tree_A0_R1_s1.log 2>&1 &
$PY -u trivial_tree.py A 1 sample 3 1 3600 20000 0 > ../logs/tree_A1_R1_s0.log 2>&1 &
$PY -u trivial_tree.py B 0 sample 3 1 3600 20000 0 > ../logs/tree_B0_R1_s0.log 2>&1 &
wait
