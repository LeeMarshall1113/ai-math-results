#!/bin/bash
# Runs both independent checkers, with their brute-force validations, one after the
# other on one thread, and writes the transcripts to logs/.  About 2 hours in total.
#   PY=/path/to/python bash run_all.sh
set -e
cd "$(dirname "$0")"
mkdir -p logs
export PY=${PY:-python}
{ time bash original/check.sh --validate ; } > logs/original.log 2>&1
{ time bash referee/check.sh --validate ; } > logs/referee.log 2>&1
tail -n 12 logs/original.log
tail -n 12 logs/referee.log
