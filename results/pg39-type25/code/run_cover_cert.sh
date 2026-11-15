#!/bin/bash
# level-1 and level-2 completeness certificates in parallel (2 python processes + short WSL solver/cake runs)
cd "$(dirname "$0")"
/c/Python314/python -u cover_cert.py level1 > ../logs/cover_cert_level1.log 2>&1 &
/c/Python314/python -u cover_cert.py level2 > ../logs/cover_cert_level2.log 2>&1 &
wait
