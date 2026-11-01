# fullcert_ext: F1 repair of the H3<=13 certificate (2026-10-31)
- Start boxes: the 570 boxes of ../dstar_ext_min.json, one per 1/10 box orbit (48 cube symmetries, acting on p as
  signed coordinate permutations) that contains no box of the original U^D_20 need set (39,380 boxes, ../fullcert/).
  Representative = lexicographically smallest box of the orbit whose interior meets int D*. Generator and proof:
  ../dstar_domain.py (log ../dstar_domain.log). Same set as replay/orbitboxes.py.
- Certifier: fast/certify2f.py, policy L3 exactly as run_fast3.sh (CERT_NSHORT_BIG=2 CERT_NSHORT=8 CERT_ORDER=st_freq
  CERT_PRESOLVE=off CERT_CHORD_RANK_BIG=0 CERT_FALLBACK=1), CERT_SB_FILE=dstar_ext_min.json, --chunk 10 (57 chunks),
  6 workers striped (run_ext.sh 6), 2026-10-31 18:49-18:57 EDT. Pilot (pilot_ext.sh, 6 boxes) wrote parts/ only; those
  parts are not part of the certificate. No fallback was triggered, no FAIL. Every record carries "pol": "L3".
- CERT_STRICT was NOT set for this run (it was added to fast/lpfastf.py afterwards); strictness is instead verified
  directly by certcheck_strict.py: all 87,828 leaves pass with every H-rep row value > 0 (open translates).
- Checks: certcheck_fast.py (unchanged) per file, ledger certcheck_fast_done.e0.jsonl, 2% reference cross-check;
  certcheck_ext.py coverage verdict over fullcert/ + fullcert_ext/ (CERTCHECK_EXT_VERDICT.txt);
  certcheck_strict.py, ledger certcheck_strict_done.e0.jsonl (certcheck_strict_ext.log).
- Chunk files are not committed (like fullcert/); SHA256SUMS lists them.
