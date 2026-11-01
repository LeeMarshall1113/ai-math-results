# Certification policy by chunk
- Chunks completed before 2026-10-29 ~05:00 EDT: certify2.py (original policy).
- From 2026-10-29 ~05:00 EDT, worker 0 (chunks ci % 3 == 0, starting with 153) runs fast/certify2f.py with the L2 policy:
  CERT_NSHORT_BIG=2 CERT_NSHORT=8 CERT_ORDER=fast/st_freq.json CERT_PRESOLVE=off. Its log is segments_fast.log.
  This policy changes only the search order and caps; every leaf is still verified exactly, and certcheck checks each chunk independently.
- Workers 1 and 2 still run certify2.py (lease hunt-4076).
- 2026-10-29 ~22:08 EDT: all three workers switched to fast/certify2f.py policy **L3** (run_fast3.sh K, log segments_fast3.log,
  leases hadwiger-h3-4252/4253/4254). L3 = L2 + `CERT_CHORD_RANK_BIG=0` (chord rank() only on boxes of side <= 1/40)
  + `CERT_FALLBACK=1` (a start box that fails under L3 is redone with caps 40 and the old chord rule before FAIL is
  declared) + a shared queue with start-box checkpoints in fullcert/parts/ (chord-zone chunks cut into 25-start-box
  pieces). Every start-box record written from here on carries `"pol"` ("L3" or "L3+fallback40").
  Why: in the chord zone (centre within 0.15 of 1/2 in all six coordinates; chunks 151-153, 157-159, 188-191, 194-196)
  no chord leaf has ever been found, so the chord short-lists stay empty and the old rule ran three full 6000-structure
  rank() scans on every failing node: 1-2 CPU-h per start box, 32-1600 leaves each, all plain. L3 on three such boxes:
  3671 s -> 50 s, 1597 s -> 0.5 s, 5112 s -> 26 s; all plain leaves, tiling and every leaf exact under certcheck.
- Salvaged before the switch (records read from the live processes with sys.remote_exec, written as parts with "pol"):
  chunk 151 start boxes 0-68 from certify2.py worker 1 ("pol": "certify2"), chunk 153 start boxes 0-49 from the L2
  worker 0 ("pol": "L2"). Chunk files without "pol" fields are certify2.py (chunks < 151 done by workers 1/2 and
  earlier) except chunk_ci with ci % 3 == 0, ci >= 153, which are L2 or L3 as stated above.
- certcheck_chunks.py: incremental driver for certcheck.py (same checks), ledger fullcert/certcheck_done.jsonl.
