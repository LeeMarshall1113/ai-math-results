#!/usr/bin/env bash
# Replays everything the note refers to.  Standard library only; one process; about 3 minutes.
#   cd karp10 && bash note/run_checks.sh > note/logs/run_checks.log 2>&1
# Set PY to a Python 3 interpreter (default: python).  Do not run with python -O (the checks are asserts).
set -u
PY="${PY:-python}"
cd "$(dirname "$0")/.." || exit 70
export PYTHONDONTWRITEBYTECODE=1
echo "== python: $($PY --version 2>&1)"
echo "== 1. every certificate in certs/ + non-self-duality (replay_all.py)"
time $PY replay/replay_all.py
rc=$?
echo "== 2. Lemma Z5, details (Smith theory, Q = <(01234)(56789)>)"
$PY replay/replay_smith.py certs/lemma_z5.json
$PY replay/replay_smith.py certs/lemma_z5_euler.json
echo "== 3. non-self-duality, details"
$PY replay/replay_nonselfdual.py certs/lemma_z5_euler.json
echo "== 3b. absolute acyclicity of the Z5 patterns and the nine named graphs (both are also run by replay_all.py)"
$PY replay/check_z5_absolute.py
$PY replay/check_nine_graphs.py
echo "== 4. the groups of types 2/4/6/8 from their finite-field descriptions"
$PY note/check_explicit.py
rc2=$?
echo "== 5. type 4, certificate with verbose class table"
$PY replay/replay_cert.py certs/type4.json -v
echo "exit codes: replay_all=$rc explicit=$rc2"
