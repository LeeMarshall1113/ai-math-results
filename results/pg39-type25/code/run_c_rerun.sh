#!/bin/bash
# Local re-run of the random sample of C cubes (results/c_rerun_sample.json, seed 20261110): WSL driver, binary LRAT + cake_lpr, records in results/c_slice.jsonl
# Runs INSIDE the lease.  Usage: bash run_c_rerun.sh [driver args]   (e.g. --workers 5 --cap-gb 17 --tlim 7200)
export MSYS_NO_PATHCONV=1
cd /c/Users/user/Documents/GitHub/erdos-triage/pg39/code
ONLY=$(/c/Python314/python -c "import json; print(','.join(json.load(open('../results/c_rerun_sample.json'))['sample']))")
wsl -d Ubuntu-24.04 -e python3 -u /mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/code/a_slice_driver.py /mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/results/c_slice_jobs.json --out c_slice --only "$ONLY" "$@"
