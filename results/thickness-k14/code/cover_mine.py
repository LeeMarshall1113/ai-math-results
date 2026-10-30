"""Completeness check of the referee's own cube tree (runs/<dir>/tree.jsonl).
1. the hull-pattern roots partition {0,1}^{H_3..H_{n-1}} (brute force over all 2^(n-3) patterns);
2. every root has a record; every split node has all 2^k children; every leaf is UNSAT+cake;
3. --regen: every leaf CNF regenerated from renc.py (base + cube units), SHA-256 compared.
4. --replay FILE: every leaf must have a record in FILE (written by replay.py) with matching CNF hash
   and a verified LRAT; agreement of the post-check LRAT hash with the recorded one is reported.
usage: cover_mine.py n t tree.jsonl [--regen] [--replay FILE]"""
import sys, json, itertools, hashlib
from renc import Enc
from cnc import hull_roots

n, t, d = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
enc = Enc(n, t)
roots = hull_roots(enc)
hv = [enc.H[v] for v in range(3, n)]
for bits in itertools.product((0, 1), repeat=len(hv)):
    a = {v if b else -v for v, b in zip(hv, bits)}
    m = sum(1 for c in roots if all(l in a for l in c))
    assert m == 1, (bits, m)
print("roots:", len(roots), "partition of all", 2 ** len(hv), "hull patterns: OK")
rec = {}
for line in open(d):
    r = json.loads(line)
    if r["kind"] in ("leaf", "split", "error"):
        rec[tuple(r["cube"])] = r
missing, bad, leaves, sat = [], [], [], []
st = [tuple(c) for c in roots]
while st:
    c = st.pop()
    r = rec.get(c)
    if r is None or r["kind"] == "error":
        missing.append(c); continue
    if r["kind"] == "split":
        st += [c + tuple(s * v for s, v in zip(sg, r["vars"])) for sg in itertools.product((1, -1), repeat=len(r["vars"]))]
    elif r["status"] == "UNSAT":
        leaves.append(c)
        if not r["cake"]:
            bad.append(c)
    else:
        sat.append(c)
print(f"leaves {len(leaves)}  missing {len(missing)}  cake-unverified {len(bad)}  SAT {len(sat)}  "
      f"solve CPU-s {sum(r.get('solve_s', 0) for r in rec.values()):.0f}  "
      f"LRAT GB {sum(r.get('lrat_bytes', 0) for r in rec.values()) / 1e9:.1f}")
if "--regen" in sys.argv:
    body = enc.dimacs().split("\n", 1)[1]
    mism = 0
    for c in leaves:
        s = f"p cnf {enc.nv} {len(enc.cl) + len(c)}\n" + body + "".join(f"{l} 0\n" for l in c)
        mism += hashlib.sha256(s.encode()).hexdigest() != rec[c]["cnf_sha"]
    print("regenerated leaf CNF hash mismatches:", mism)
ok = not (missing or bad or sat)
if "--replay" in sys.argv:
    rp = {}
    for line in open(sys.argv[sys.argv.index("--replay") + 1]):
        q = json.loads(line)
        rp[tuple(q["cube"])] = q
    nok = sum(1 for c in leaves if c in rp and rp[c]["cnf_match"] and rp[c]["verified"])
    same = sum(1 for c in leaves if c in rp and rp[c]["lrat_sha"] == rec[c]["lrat_sha"])
    print(f"replay: {nok}/{len(leaves)} leaves regenerated, hash-matched and verified; "
          f"LRAT hash identical to the recorded one for {same}")
    ok = ok and nok == len(leaves)
print("COVER OK" if ok else "COVER INCOMPLETE")
sys.exit(0 if ok else 1)
