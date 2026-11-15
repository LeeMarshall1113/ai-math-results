#!/bin/bash
# level-2 split of the 5 cubes that were not closed after 2 h on Modal (stop rule), 5 python processes in parallel
cd "$(dirname "$0")"
for id in A1_c04 A2_c01 A2_c04 A2_c05 A2_c17; do
  /c/Python314/python -u split_level2.py $id 6 300 > ../logs/split_$id.log 2>&1 &
done
wait
