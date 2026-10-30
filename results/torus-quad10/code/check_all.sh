#!/bin/bash
# One-command replay. Run from anywhere:  bash code/check_all.sh
# Needs: python3 with numpy+scipy+pysat (only for numtest_ad.py / mus.py, optional) - the core checks use the
# standard library only; CaDiCaL (>= 1.9, LRAT output); optionally an LRAT checker.
#   CADICAL=/path/cadical        (required)
#   LRAT_CHECK="cmd"             (optional; called as: $LRAT_CHECK file.cnf file.lrat; e.g. cake_lpr)
#   LRAT_CHECK2="cmd"            (optional second checker)
set -eu
cd "$(dirname "$0")"
PY=${PY:-python3}
CADICAL=${CADICAL:?set CADICAL to a cadical binary}
W=../work; mkdir -p $W cnf
echo "== Lemma 2: map enumeration, two independent methods"
$PY enum_lattice.py | tail -2
$PY enum_backtrack.py | tail -3
echo "== exact geometric checker + encoder soundness on genuine tori (9-face, Seamons 11-face)"
$PY soundness_test.py
for w in T9_known9face T11_seamons; do $PY encode2.py cnf/e2_$w.cnf witnesses/$w.json | head -1; done
echo "== mutation tests"
$PY mutation_test.py
echo "== Lemma 5: regenerate CNFs (10 faces; and 8 faces: Q0, Q3 dropped), compare hashes, solve, check LRAT"
$PY encode.py $W/m10_core.cnf GP,PL,AD
$PY encode.py $W/m10_full.cnf GP,PL,CV,AD,INT,SYM
$PY encode2.py $W/m10_indep.cnf
$PY encode.py $W/m8_core.cnf GP,PL,AD 0,3
$PY encode2.py $W/m8_indep.cnf --drop=0,3
(cd $W && sha256sum m10_core.cnf m10_full.cnf m10_indep.cnf m8_core.cnf m8_indep.cnf) > $W/cnf.sha256
diff <(sed 's/ [ *]/  /' $W/cnf.sha256 | sort -k2) <(grep 'cnf$' ../certs/CHECKSUMS.sha256 | sort -k2) || { echo "CNF HASH MISMATCH"; exit 1; }
echo "CNF hashes match certs/CHECKSUMS.sha256"
NAMES="m10_core m10_full m10_indep m8_core m8_indep"
for n in $NAMES; do
  r=$($CADICAL -q $W/$n.cnf | head -1 || true); echo "$n: $r"
  [ "$r" = "s UNSATISFIABLE" ] || { echo "FAIL: not UNSAT"; exit 1; }
done
[ -f ../certs/m10_indep.lrat ] || { echo "regenerating m10_indep.lrat (31 MB)"; $CADICAL -q --lrat=true --binary=false $W/m10_indep.cnf ../certs/m10_indep.lrat >/dev/null || true; }
for n in $NAMES; do
  for C in "${LRAT_CHECK:-}" "${LRAT_CHECK2:-}"; do
    [ -n "$C" ] || continue
    echo "$n: $C -> $($C $W/$n.cnf ../certs/$n.lrat 2>&1 | tr -d '\000' | tail -1)"
  done
done
echo "ALL CHECKS DONE"
