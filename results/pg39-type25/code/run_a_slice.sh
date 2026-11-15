#!/bin/bash
# Runs INSIDE the lease: preparation (Windows python) then the WSL driver.  Args are passed to the driver.
export MSYS_NO_PATHCONV=1
cd /c/Users/user/Documents/GitHub/erdos-triage/pg39/code
[ -f ../results/a_slice_jobs.json ] || /c/Python314/python -u a_slice_prep.py || exit 1      # lists and CNFs are deterministic: prepared once
wsl -d Ubuntu-24.04 -e python3 -u /mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/code/a_slice_driver.py /mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/results/a_slice_jobs.json "$@"
