#!/bin/bash
# Local re-run of the seeded random sample of level-2 children (results/children_rerun_sample.json, seed 20261111).  Runs INSIDE the broker lease.
# --dry only redirects the records to results/children_rerun_dry.jsonl and results/a_slice_dry_manifest.jsonl (certs/manifest.jsonl is NOT touched).
export MSYS_NO_PATHCONV=1
cd /c/Users/user/Documents/GitHub/erdos-triage/pg39/code
wsl -d Ubuntu-24.04 -e python3 -u /mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/code/a_slice_driver.py /mnt/c/Users/user/Documents/GitHub/erdos-triage/pg39/results/children_rerun_jobs.json --dry --out children_rerun --workers 3 --cap-gb 10 --tlim 7200 "$@"
