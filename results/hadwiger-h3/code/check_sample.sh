#!/bin/bash
# Quick check (about 3 min, one thread, under 1 GB). Needs Python 3 with numpy and scipy. No download needed.
#   cd results/hadwiger-h3/code && bash check_sample.sh
# It checks everything except the bulk of the leaves:
#   1. the corner counterexample to Prymak's Lemma 2.4(ii) (exact);
#   2. the counterexample showing Prymak's region D is not a fundamental domain (exact);
#   3. that D* is a fundamental domain and the 570 extension start boxes (exact, dstar_domain.py);
#   4. orbit coverage: every 1/10 box of [0,1]^6 has a cube-symmetry image among the 39,950 start boxes;
#   5. both leaf checkers (certcheck_strict.py, structure-based; replay/hcheck.py, structure-free) on
#      ../certs/sample: the centre start box (1,771 leaves, 118 chord leaves) plus 3 other start boxes
#      and 3 extension start boxes, copied verbatim from the certificate.
# The full check is check_all.sh (needs the certificate archives from Zenodo).
set -euo pipefail
PY=${PY:-python}
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CHECK_PAUSE_ON_GAMING=0
# scratch dir inside code/; pwd -W gives a C:/ path that Windows Python accepts under Git Bash
TMP=$(cd "$(mktemp -d ./_check.XXXXXX)" && (pwd -W 2>/dev/null || pwd))
trap 'rm -rf "$TMP"' EXIT
fail() { echo "FAIL: $*"; exit 1; }

echo "== 1. Lemma 2.4(ii) corner example"
$PY corner_example.py | tee "$TMP/corner.txt" | tail -1
grep -q "^RESULT: no translate" "$TMP/corner.txt" || fail "corner example"

echo "== 2. Prymak's D is not a fundamental domain"
$PY symmetry_gap.py | tee "$TMP/gap.txt" | tail -2
grep -q "images in D: 0" "$TMP/gap.txt" || fail "D counterexample"

echo "== 3. D* and the extension boxes"
(cd newmethod && $PY dstar_domain.py) | tee "$TMP/dstar.txt" | tail -2
grep -q "box orbits with no NEED box (= start boxes actually needed): 570" "$TMP/dstar.txt" || fail "dstar_domain"

echo "== 4. orbit coverage of the 1/10 grid"
$PY -c "import json; a=json.load(open('../certs/startboxes_main.json')); b=json.load(open('../certs/startboxes_ext.json')); json.dump(a+b, open('$TMP/sb.json','w'))"
(cd replay && $PY orbitboxes.py "$TMP/sb.json") | tee "$TMP/orbit.txt"
grep -q '"n_missing_orbits": 0' "$TMP/orbit.txt" || fail "orbit coverage"

echo "== 5a. certcheck_strict.py on the sample"
mkdir -p "$TMP/sample" && cp ../certs/sample/chunk_*.jsonl.gz "$TMP/sample/"
$PY newmethod/certcheck_strict.py "$TMP/sample" --only 0,1 --xcheck 0.01 | tee "$TMP/strict.txt"
[ "$(grep -c ' bad 0 ' "$TMP/strict.txt")" = 2 ] || fail "certcheck_strict"
grep -q 'disagree 0' "$TMP/strict.txt" || fail "reference cross-check"

echo "== 5b. replay/hcheck.py --strict on the sample"
$PY replay/run_all.py "$TMP/rout" "$TMP/sample" --strict --workers 1 > "$TMP/replay.log" 2>&1
(cd replay && $PY summarize.py "$TMP/rout") > "$TMP/summary.txt"
$PY -c "import json,sys; s=json.load(open('$TMP/rout/SUMMARY.json')); t=s['totals']; print('replay: leaves', t['n'], 'chord', t['chord'], 'start boxes', t['startboxes_checked'], 'tiling ok', t['tiling_ok'], 'failures', s['n_fail']); sys.exit(0 if s['n_fail']==0 and t['tiling_ok']==t['startboxes_checked']==7 else 1)" || fail "replay"

echo "SAMPLE CHECK PASSED"
