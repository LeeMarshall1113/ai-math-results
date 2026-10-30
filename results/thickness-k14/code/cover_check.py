"""Completeness check of a run_all.py tree (independent of the solver driver's control flow).

For every case in case_list(n): the root cube [] has a record; every SPLIT node on variables V has a
record for each of the 2^|V| children; every leaf is UNSAT with cake_lpr verified (or SAT).
Optionally (--regen) re-generates each leaf CNF and compares its SHA-256 to the recorded one, which
ties the verified proof to the encoding in gen.py. --nosweep: the tree was made with TK_NOSWEEP=1.
--replay FILE: every reachable leaf must also have a record in FILE (written by replay.py) whose
regenerated CNF hash matches and whose LRAT was verified; the LRAT hash taken after checking is
compared with the recorded one and the agreement is reported.
usage: python cover_check.py n t tree.jsonl [--regen] [--nosweep] [--replay FILE]
"""
import hashlib, itertools, json, os, sys
if "--nosweep" in sys.argv:
    os.environ["TK_NOSWEEP"] = "1"
from run_all import case_list, case_clauses

n, t, res = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rec = {}
for line in open(res):
    r = json.loads(line)
    rec[(r["kind"], tuple(r["S"]), tuple(r["cube"]))] = r
missing, bad, leaves, sat = [], [], 0, []
for case in case_list(n):
    stack = [()]
    while stack:
        cube = stack.pop()
        key = (case[0], tuple(case[1]), cube)
        r = rec.get(key)
        if r is None:
            missing.append(key); continue
        if r["status"] == "SPLIT":
            for signs in itertools.product((1, -1), repeat=len(r["vars"])):
                stack.append(cube + tuple(s * v for s, v in zip(signs, r["vars"])))
        elif r["status"] == "UNSAT":
            leaves += 1
            if not r.get("cake"):
                bad.append(key)
        else:
            sat.append(key)
print(f"cases {len(case_list(n))}  verified leaves {leaves}  missing {len(missing)}  unverified {len(bad)}  SAT {len(sat)}")
if missing: print("first missing:", missing[:5])
if bad: print("first unverified:", bad[:5])
if sat: print("SAT:", sat[:5])
if "--regen" in sys.argv:
    mism = 0
    for case in case_list(n):
        nv, cl, _ = case_clauses(n, t, *case)
        body = f"".join(" ".join(map(str, c)) + " 0\n" for c in cl)
        for key, r in rec.items():
            if (key[0], key[1]) != (case[0], tuple(case[1])) or r["status"] != "UNSAT":
                continue
            cube = list(key[2])
            s = f"p cnf {nv} {len(cl) + len(cube)}\n" + body + "".join(f"{l} 0\n" for l in cube)
            if hashlib.sha256(s.encode()).hexdigest() != r["cnf_sha"]:
                mism += 1
    print("regenerated CNF hash mismatches:", mism)
ok = not missing and not bad and not sat
if "--replay" in sys.argv:
    rp = {}
    for line in open(sys.argv[sys.argv.index("--replay") + 1]):
        q = json.loads(line)
        rp[(q["kind"], tuple(q["S"]), tuple(q["cube"]))] = q
    reach = [k for k, r in rec.items() if r["status"] == "UNSAT"]
    nok = sum(1 for k in reach if k in rp and rp[k]["cnf_match"] and rp[k]["verified"])
    same = sum(1 for k in reach if k in rp and rp[k]["lrat_sha"] == rec[k]["lrat_sha"])
    print(f"replay: {nok}/{len(reach)} leaves regenerated, hash-matched and verified; "
          f"LRAT hash identical to the recorded one for {same}")
    ok = ok and nok == len(reach)
print("COVER OK: every case refuted with cake_lpr-verified LRAT" if ok else "COVER INCOMPLETE")
sys.exit(0 if ok else 1)
