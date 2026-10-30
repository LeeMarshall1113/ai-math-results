"""Re-certify the leaves of a published case tree from scratch.

For each reachable leaf of the tree: regenerate the leaf CNF from the encoder (gen.py via
run_all.case_clauses for the authors' trees, renc.py for the referee's tree), compare its SHA-256
with the published one, solve it with CaDiCaL (LRAT), check the LRAT with $LRAT_CHECK, hash the
LRAT file the checker read, delete CNF and LRAT, and append a record to the output manifest.
The run is resumable: leaves already recorded as verified in the output manifest are skipped.

usage: python replay.py {sweep|nosweep|referee} [--sample K] [--seed S] [--out FILE]
  sweep    certs/tree_authors_sweep.jsonl    (gen.py, sweep lemmas on)
  nosweep  certs/tree_authors_nosweep.jsonl  (gen.py, TK_NOSWEEP=1)
  referee  certs/tree_referee.jsonl          (renc.py, independent encoder and cube partition)
Exit code 0 iff every replayed leaf had a matching CNF hash, an UNSAT answer and a verified LRAT
(and, without --sample, every reachable leaf was replayed).
"""
import argparse, hashlib, itertools, json, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("tree", choices=["sweep", "nosweep", "referee"])
ap.add_argument("--sample", type=int, default=0)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", default=None)
a = ap.parse_args()
if a.tree == "nosweep":
    os.environ["TK_NOSWEEP"] = "1"
import tools  # noqa: E402

N, T = 14, 3
WORK = os.path.abspath(os.environ.get("TK_WORK", os.path.join(HERE, "..", "work")))
os.makedirs(WORK, exist_ok=True)
out_path = a.out or os.path.join(WORK, f"replay_{a.tree}.jsonl")
tree_path = os.path.join(HERE, "..", "certs",
                         {"sweep": "tree_authors_sweep.jsonl", "nosweep": "tree_authors_nosweep.jsonl",
                          "referee": "tree_referee.jsonl"}[a.tree])

# ---- reachable leaves: (key, record, function building the leaf CNF text) ----
leaves = []
if a.tree == "referee":
    from renc import Enc
    from cnc import hull_roots
    enc = Enc(N, T)
    body = enc.dimacs().split("\n", 1)[1]
    rec = {}
    for line in open(tree_path):
        r = json.loads(line)
        if r["kind"] in ("leaf", "split"):
            rec[tuple(r["cube"])] = r
    st = [tuple(c) for c in hull_roots(enc)]
    while st:
        c = st.pop()
        r = rec[c]
        if r["kind"] == "split":
            st += [c + tuple(s * v for s, v in zip(sg, r["vars"]))
                   for sg in itertools.product((1, -1), repeat=len(r["vars"]))]
        else:
            leaves.append((("referee", (), c), r, None))
    base = lambda case: (enc.nv, len(enc.cl), body)
else:
    from run_all import case_list, case_clauses
    rec = {}
    for line in open(tree_path):
        r = json.loads(line)
        rec[(r["kind"], tuple(r["S"]), tuple(r["cube"]))] = r
    for kind, S in case_list(N):
        st = [()]
        while st:
            c = st.pop()
            r = rec[(kind, tuple(S), c)]
            if r["status"] == "SPLIT":
                st += [c + tuple(s * v for s, v in zip(sg, r["vars"]))
                       for sg in itertools.product((1, -1), repeat=len(r["vars"]))]
            else:
                leaves.append(((kind, tuple(S), c), r, None))
    _cache = {}

    def base(case):
        if case not in _cache:
            _cache.clear()
            nv, cl, _ = case_clauses(N, T, *case)
            _cache[case] = (nv, len(cl), "".join(" ".join(map(str, x)) + " 0\n" for x in cl))
        return _cache[case]

assert all(r.get("status") == "UNSAT" for _, r, _ in leaves), "tree has a non-UNSAT leaf"
print(f"{a.tree}: {len(leaves)} reachable leaves in {os.path.basename(tree_path)}", flush=True)
todo = leaves
if a.sample:
    todo = random.Random(a.seed).sample(leaves, a.sample)
todo.sort(key=lambda x: (x[0][0], x[0][1], x[0][2]))

done = {}
if os.path.exists(out_path):
    for line in open(out_path):
        q = json.loads(line)
        if q["cnf_match"] and q["verified"]:
            done[(q["kind"], tuple(q["S"]), tuple(q["cube"]))] = q
out = open(out_path, "a")
fails = 0
for i, (key, r, _) in enumerate(todo):
    if key in done:
        continue
    kind, S, cube = key
    nv, ncl, body = base((kind, S))
    text = f"p cnf {nv} {ncl + len(cube)}\n" + body + "".join(f"{l} 0\n" for l in cube)
    data = text.encode()
    cnf, lrat = os.path.join(WORK, "leaf.cnf"), os.path.join(WORK, "leaf.lrat")
    with open(cnf, "wb") as f:
        f.write(data)
    match = hashlib.sha256(data).hexdigest() == r["cnf_sha"]
    t0 = time.time()
    rc, _ = tools.solve(cnf, lrat)
    t1 = time.time()
    ok, msg, lsha, lbytes = False, "", None, None
    if rc == 20:
        ok, msg = tools.check_lrat(cnf, lrat)
        lsha, lbytes = tools.sha_file(lrat), os.path.getsize(lrat)
    t2 = time.time()
    for p in (cnf, lrat):
        if os.path.exists(p):
            os.remove(p)
    row = dict(tree=a.tree, kind=kind, S=list(S), cube=list(cube), cnf_sha=r["cnf_sha"], cnf_match=match,
               rc=rc, verified=ok, solve_s=round(t1 - t0, 1), check_s=round(t2 - t1, 1),
               lrat_sha=lsha, lrat_bytes=lbytes, lrat_sha_as_published=(lsha == r["lrat_sha"]))
    out.write(json.dumps(row) + "\n"); out.flush()
    good = match and rc == 20 and ok
    fails += not good
    print(f"[{i + 1}/{len(todo)}] {kind} {list(S)} cube={len(cube)} cnf_hash={'ok' if match else 'MISMATCH'} "
          f"rc={rc} verified={ok} solve={t1 - t0:.1f}s check={t2 - t1:.1f}s "
          f"lrat={0 if lbytes is None else lbytes >> 20}MB lrat_hash_as_published={lsha == r['lrat_sha']}",
          flush=True)
    if not ok and msg:
        print("   checker said:", msg[-200:], flush=True)
out.close()

rows = {}
for line in open(out_path):
    q = json.loads(line)
    rows[(q["kind"], tuple(q["S"]), tuple(q["cube"]))] = q
want = [k for k, _, _ in todo]
good = [k for k in want if k in rows and rows[k]["cnf_match"] and rows[k]["verified"]]
same = sum(1 for k in good if rows[k]["lrat_sha_as_published"])
print(f"{a.tree}: {len(good)}/{len(want)} leaves regenerated with matching CNF hash, UNSAT, LRAT verified; "
      f"LRAT byte-identical to the published hash for {same}/{len(good)}")
if len(good) == len(want):
    print("REPLAY OK" + ("" if a.sample else f": all {len(leaves)} reachable leaves of {a.tree}"))
    sys.exit(0)
print("REPLAY FAILED")
sys.exit(1)
