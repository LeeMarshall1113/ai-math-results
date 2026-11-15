#!/bin/bash
cd /c/Users/user/Documents/GitHub/erdos-triage/pg39/code
PY=/c/Python314/python
$PY -u trivial_tree.py C 0 sample 1 1 2700 20000 3 > ../logs/tree_C0_R1_s3.log 2>&1 &
$PY -u trivial_tree.py C 2 sample 1 1 2700 20000 0 > ../logs/tree_C2_R1_s0.log 2>&1 &
wait
