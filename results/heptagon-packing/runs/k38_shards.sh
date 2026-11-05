#!/bin/bash
cd /c/Users/hackathon/Documents/GitHub/erdos-triage/heptagon-packing
export PYTHONPATH=/c/Users/hackathon/AppData/Roaming/Python/Python314/site-packages OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export HEPT_MU='[-1.1684855e-12, -1.1684853e-12]' HEPT_DUAL2=1 HEPT_CHARTS=runs/charts_w030c.json
REG='{"lo": [1.8019377358048, 0, 1.2, 0.76, -1e-12, -1e-12], "hi": [4.0, 2.0, 3.5, 0.78, 0.8975979010266552, 0.8975979010266552]}'
for i in 0 1 2 3 4 5; do
  HEPT_RESUME=runs/sweep_main_a1r/shard38_$i.npy /c/Python314/python src/bnb2.py runs/f_w030.npy 3 2 1.5326754446782214 14038 "$REG" runs/sweep_main_a1r/shard38_${i}_open.npy > runs/sweep_main_a1r/shard38_$i.log 2>&1 &
done
wait
