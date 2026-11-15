#!/usr/bin/env python3
"""Independent replay of the case-split certificates (2a, 2b, 2d families).

  python replay_cases.py [NAME ...] [--cake /path/to/cake_lpr] [--all]

For each case NAME of ../certs/cases.json (written by code/certify_cases.py) this rebuilds, with the independent
geometry of replay.py (row-reduced lines, own orbit code), the CNF
    orbit CNF of the class,  |K| = k,  one unit clause per orbit touched by the case's 'ins' / 'outs' points,
checks that its SHA-256 equals the one recorded in ../certs/manifest.jsonl, and (with --cake) that cake_lpr accepts
(CNF, LRAT).  Together with the case-list arguments of NOTES.md section 4 this is the replayable certificate.
"""
import hashlib
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay as R

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, "..", "certs")


def main():
    args = sys.argv[1:]
    cake = None
    if "--cake" in args:
        i = args.index("--cake")
        cake = args[i + 1]
        del args[i:i + 2]
    cases = json.load(open(os.path.join(CERTS, "cases.json")))
    names = list(cases) if (not args or "--all" in args) else args
    classes = json.load(open(os.path.join(CERTS, "classes.json")))
    manifest = {}
    for line in open(os.path.join(CERTS, "manifest.jsonl")):
        r = json.loads(line)
        manifest[r["tag"]] = r
    points = R.build_points()
    pidx = {p: i for i, p in enumerate(points)}
    lines = R.build_lines_rref(points)
    print("PG(3,9): %d points, %d lines (row-reduced enumeration)" % (len(points), len(lines)))
    ok_all = True
    cache = {}
    for name in names:
        c = cases[name]
        cls = classes[c["cls"]]
        M = [[R.decode(e) for e in row] for row in cls["matrix"]]
        perm = R.point_perm(points, pidx, M, cls["frob"])
        orbs, oid = R.orbits_of(len(points), [perm])
        units = {}
        for p in c["ins"]:
            units[oid[p]] = 1
        for p in c["outs"]:
            assert units.get(oid[p], 0) != 1
            units[oid[p]] = 0
        cnf, comment, info = R.build_cnf_for_class(c["cls"], c["k"], points, lines, pidx, cls,
                                                   units=sorted(units.items()),
                                                   header="class %s k %d case %s amo=True" % (c["cls"], c["k"], name))
        text = cnf.dimacs([comment])
        h = hashlib.sha256(text.encode()).hexdigest()
        rec = manifest[name]
        match = h == rec["cnf_sha256"]
        ok_all &= match
        line = "%-40s k=%d  CNF %s" % (name, c["k"], "MATCH" if match else "MISMATCH")
        lr = os.path.join(CERTS, name + ".lrat")
        if cake and os.path.exists(lr):
            cn = os.path.join(CERTS, name + ".cnf")
            r = subprocess.run([cake, "--CML_HEAP_SIZE=1500", "--CML_STACK_SIZE=500", cn, lr], capture_output=True, text=True)
            line += "  cake_lpr: " + ((r.stdout + r.stderr).strip().splitlines() or ["?"])[-1]
        print(line, flush=True)
    print("ALL MATCH" if ok_all else "SOME MISMATCH")


if __name__ == "__main__":
    main()
