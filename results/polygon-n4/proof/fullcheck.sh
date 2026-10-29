#!/usr/bin/env bash
# Full replay of all 836 canonical faces. Fails (nonzero exit) unless every class passes:
# verify.py covers the 690 Bernstein faces and 136 trees, verify_tube.py the 10 tubes.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python}
log=$(mktemp)
trap 'rm -f "$log"' EXIT
"$PY" -u verify.py | tee "$log"
"$PY" -u verify_tube.py tube_t158.jsonl tube_tubeA.jsonl tube_tubeB.jsonl tube_tubeC.jsonl tube_tubeD.jsonl | tee -a "$log"
grep -q "^checked 826 faces .*failures: 0 \[\]" "$log"
grep -q "^ALL OK" "$log"
grep -q "^tube certificates verified: 10 all ok" "$log"
echo "FULL CHECK PASSED: 826 + 10 = 836 canonical faces certified"
