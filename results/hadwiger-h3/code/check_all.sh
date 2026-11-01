#!/bin/bash
# Full check of the H_3 <= 13 certificate.
#   DATA=/path/to/archives WORKERS=8 bash check_all.sh
# DATA must hold hadwiger-h3-fullcert.tar and hadwiger-h3-fullcert_ext.tar (Zenodo dataset record, see README).
# Needs Python 3 with numpy and scipy, and GNU tar and sha256sum. About 45 CPU-hours in total
# (certcheck_fast ~22 h, certcheck_strict ~21 h, replay ~2 h, on one 2026 desktop core); WORKERS processes run
# in parallel, each under about 1.5 GB. Every stage is resumable: rerun the script and finished chunk files are
# skipped (their verdicts are kept in ledger files keyed by the chunk's sha256 and the checker's sha256).
# Optional: LEAN=1 also builds the Lean project (needs elan; see ../lean/README.md).
set -euo pipefail
PY=${PY:-python}
WORKERS=${WORKERS:-1}
DATA=${DATA:?set DATA to the folder holding the two .tar archives}
cd "$(dirname "$0")"
CODE=$(pwd -W 2>/dev/null || pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CHECK_PAUSE_ON_GAMING=0
fail() { echo "FAIL: $*"; exit 1; }

echo "== 0. archives and chunk hashes"
(cd "$DATA" && sha256sum -c "$CODE/../certs/SHA256SUMS.archives") || fail "archive sha256"
cd newmethod
for d in fullcert fullcert_ext; do
    if [ ! -d $d ]; then
        tar -xf - < "$DATA/hadwiger-h3-$d.tar"   # via stdin: GNU tar on Windows reads "C:/..." as a remote host
        # the archives contain the ledgers of the original run; move them aside so that every chunk is re-checked
        mkdir -p $d/ledgers_original && mv $d/certcheck_*_done*.jsonl $d/ledgers_original/
    fi
done
(cd fullcert && sha256sum -c --quiet ../../../certs/SHA256SUMS.fullcert) || fail "fullcert chunk sha256"
(cd fullcert_ext && sha256sum -c --quiet ../../../certs/SHA256SUMS.fullcert_ext) || fail "fullcert_ext chunk sha256"
echo "chunk hashes OK"

shards() {   # shards N DIR NCHUNKS -> DIR/shardK.txt, chunk indices striped over N files
    for ((k = 0; k < $1; k++)); do seq -s, $k $1 $(($3 - 1)) > $2/myshard$k.txt; done
}
run_sharded() {   # run_sharded CHECKER DIR NCHUNKS
    shards $WORKERS $2 $3
    for ((k = 0; k < WORKERS; k++)); do
        $PY $1 $2 --only-file $2/myshard$k.txt --tag w$k > $2/$(basename $1 .py)_w$k.log 2>&1 &
    done
    wait
}

echo "== 1. certcheck_fast.py (exact integer arithmetic, all leaves; 1% cross-checked by the rational reference)"
run_sharded certcheck_fast.py fullcert 197
$PY certcheck_fast.py fullcert --verdict-only | tee _my_fast_verdict.txt
grep -q "^COMPLETE AND VALID" _my_fast_verdict.txt || fail "certcheck_fast fullcert"
run_sharded certcheck_fast.py fullcert_ext 57

echo "== 2. certcheck_ext.py (symmetry-correct coverage over fullcert + fullcert_ext)"
$PY certcheck_ext.py | tee _my_ext_verdict.txt
grep -q "^COMPLETE AND VALID" _my_ext_verdict.txt || fail "certcheck_ext"

echo "== 3. certcheck_strict.py (the Lean leaf predicate LeafOK: strict containment)"
run_sharded certcheck_strict.py fullcert 197
$PY certcheck_strict.py fullcert --verdict-only | tee _my_strict_verdict.txt
grep -q "^COMPLETE AND VALID" _my_strict_verdict.txt || fail "certcheck_strict fullcert"
run_sharded certcheck_strict.py fullcert_ext 57
$PY - <<'EOF' || fail "certcheck_strict fullcert_ext"
import glob, json, sys
L = {}
for f in glob.glob('fullcert_ext/certcheck_strict_done.w*.jsonl'):
    for line in open(f): v = json.loads(line); L[v['file']] = v
print('fullcert_ext strict: files', len(L), 'leaves', sum(v['leaves'] for v in L.values()),
      'bad', sum(v['bad'] for v in L.values()))
sys.exit(0 if len(L) == 57 and all(v['bad'] == 0 for v in L.values()) else 1)
EOF

echo "== 4. independent replay (replay/hcheck.py --strict, structure-free) and orbit coverage"
cd ../replay
$PY run_all.py "$CODE/replay/out_strict" "$CODE/newmethod/fullcert" "$CODE/newmethod/fullcert_ext" --strict --workers $WORKERS
$PY summarize.py "$CODE/replay/out_strict" > /dev/null
$PY - <<'EOF' || fail "replay"
import json, sys
s = json.load(open('out_strict/SUMMARY.json')); t = s['totals']; g = s['coverage']['grid_orbit_test']
print('replay: leaves', t['n'], 'chord', t['chord'], 'start boxes', t['startboxes_checked'], 'tiling ok', t['tiling_ok'],
      'failures', s['n_fail'], 'orbits with no certified box', g['n_missing_orbits'])
sys.exit(0 if (s['n_fail'] == 0 and t['n'] == 7561583 and t['tiling_ok'] == t['startboxes_checked'] == 39950
               and g['n_missing_orbits'] == 0) else 1)
EOF
cd ..

echo "== 5. the two counterexamples and D*"
bash check_sample.sh

if [ "${LEAN:-0}" = 1 ]; then
    echo "== 6. Lean"
    (cd ../lean && lake exe cache get && lake build && lake env lean Hadwiger/Audit.lean)
fi
echo "FULL CHECK PASSED"
