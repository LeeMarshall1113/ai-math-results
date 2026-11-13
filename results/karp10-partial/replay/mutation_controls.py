"""mutation_controls.py -- targeted mutation controls for the Karp n=10 replay: every mutation must make the replay FAIL, and the failure message is recorded.

  python replay/mutation_controls.py            layers L1 + L2 (about 10 minutes)
  python replay/mutation_controls.py --full     also layer L3 (three full replay_all runs on re-pinned trees, about 8 minutes more)

Each mutation is applied to a FRESH COPY of certs/ and replay/ in a temporary directory (the release tree is never touched).  Layers:
  L1  `replay_all.py` on the mutated copy with the released manifest: the pinned SHA-256 / exact file-set / structure checks must refuse it.
  L2  the individual replay script(s) that consume the mutated file, run DIRECTLY (no manifest involved): the mathematical checks themselves must refuse
      it (assertion, wrong result, missing file, solutions found ...).  This is the stronger test: it shows the checker, not only the manifest, rejects.
  L3  (--full) the manifest is RE-PINNED to the mutated files (what an adversary who also edits the manifest would do) and replay_all.py is run in full:
      structure, hard-coded counts, row pins and the checks themselves must still refuse.
A control is DETECTED when the replay exits non-zero (L1, L3) or when at least one consuming script exits non-zero (L2).  Mutations marked informational
are NOT expected to fail (they alter something redundant) and are reported separately.  Exit code 1 if any non-informational mutation is not detected at L2.
"""
import copy
import glob
import hashlib
import itertools
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, here)
import replay_lib as RL                       # noqa: E402   (used only to CHOOSE mutations that are guaranteed to be wrong; the verdicts come from subprocesses)

PY = sys.executable
FULL = "--full" in sys.argv
man = json.load(open(os.path.join(here, "RELEASE-MANIFEST.json")))
ROWS = {r["name"]: r for r in man["expected_rows"]}

pristine = tempfile.mkdtemp(prefix="karp_mut_pristine_")
shutil.copytree(os.path.join(root, "certs"), os.path.join(pristine, "certs"))
os.makedirs(os.path.join(pristine, "replay"))
for f in glob.glob(os.path.join(here, "*.py")) + [os.path.join(here, "RELEASE-MANIFEST.json")]:
    shutil.copy(f, os.path.join(pristine, "replay"))

records = []      # (id, description, layer, expected_fail, rc, message, detected)


def fresh():
    w = tempfile.mkdtemp(prefix="karp_mut_work_")
    shutil.rmtree(w)
    shutil.copytree(pristine, w)
    return w


def runpy(w, args, timeout=300):
    try:
        r = subprocess.run([PY, "-I", "-B"] + args, cwd=w, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return 124, "(timeout)"
    lines = [l.strip() for l in (r.stdout + "\n" + r.stderr).splitlines() if l.strip()]
    return r.returncode, (lines[-1] if lines else "(no output)")


def l1(w):
    """replay_all.py on the mutated copy with the released manifest"""
    r = subprocess.run([PY, "-I", "-B", "replay/replay_all.py"], cwd=w, capture_output=True, text=True)
    lines = [l.strip() for l in r.stdout.splitlines() if l.strip()]
    msg = next((l for l in lines if l.startswith("FAILED CHECK")), None) or next((l for l in lines if l.startswith("manifest:")), lines[-1] if lines else "(none)")
    return r.returncode, msg


def consumer_rows(cert_name):
    """the rows of the manifest that read this certificate"""
    return [r for r in man["expected_rows"] if any(a == f"certs/{cert_name}" for a in r["args"])]


def l2(w, cert_name):
    """run every row that consumes the certificate, directly; returns (detected, message of the first failing row)"""
    msgs = []
    for r in consumer_rows(cert_name):
        rc, msg = runpy(w, [os.path.join("replay", os.path.basename(r["script"]))] + r["args"])
        ok_line = r["must_contain"] in msg
        if rc != 0 or not ok_line:
            return True, f"[{r['name']}] rc={rc}: {msg}"
        msgs.append(f"[{r['name']}] rc=0")
    return False, "; ".join(msgs)


def record(mid, desc, w, cert_name, expected_fail=True, note=""):
    rc1, m1 = l1(w)
    det2, m2 = l2(w, cert_name) if cert_name else (None, "")
    records.append((mid, desc, "L1", True, rc1, m1, rc1 != 0))
    if cert_name:
        records.append((mid, desc, "L2", expected_fail, 1 if det2 else 0, m2, det2))
    tag2 = "" if cert_name is None else f" | L2 {'DETECTED' if det2 else 'NOT DETECTED'}: {m2[:150]}"
    print(f"{mid:34s} {desc[:60]:60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:90]}{tag2}{(' | ' + note) if note else ''}", flush=True)


def load(w, name):
    return json.load(open(os.path.join(w, "certs", name)))


def save(w, name, c):
    json.dump(c, open(os.path.join(w, "certs", name), "w"), indent=1)


cert_names = sorted(f for f in man["certificates"])
type_certs = [n for n in cert_names if man["certificates"][n]["kind"] == "cert"]
smith_certs = [n for n in cert_names if man["certificates"][n]["kind"] == "smith"]


def group_entries(c):
    """[(label, dict)] for every Oliver group witness inside a certificate"""
    out = [(f"group {n}", g) for n, g in c.get("groups", {}).items()]
    out += [(f"euler[{i}]", e) for i, e in enumerate(c.get("euler", []))]
    return out


def identity():
    return tuple(range(10))


def ptuple(x):
    return tuple(x)


def choose_wrong_generator(g):
    """swap two images in the first generator so that the closure has a different order than the stated one"""
    gens = [list(x) for x in g["gens"]]
    for i, j in itertools.combinations(range(10), 2):
        new = [x[:] for x in gens]
        new[0][i], new[0][j] = new[0][j], new[0][i]
        try:
            order = len(RL.closure([tuple(x) for x in new]))
        except RuntimeError:
            continue
        if order != g["order"]:
            return new, order
    return None, None


def choose_nonnormal(g):
    G = RL.closure([tuple(x) for x in g["gens"]])
    p = g["p"]
    gens = [tuple(x) for x in g["gens"]]
    for x in sorted(G, key=lambda x: len(RL.closure([x]))):
        H = RL.closure([x])
        n = len(H)
        if n == 1:
            continue
        t = n
        while t % p == 0:
            t //= p
        if t != 1:
            continue
        if any(RL.compose(RL.compose(s, h), RL.inverse(s)) not in H for s in gens for h in H):
            return [list(x)], n
    return None, None


def mutate_groups(cname):
    """group-witness mutations (all of which are mathematically wrong witnesses) for every group of the certificate"""
    c0 = load(pristine, cname)
    for label, g0 in group_entries(c0):
        def with_group(mod):
            w = fresh()
            c = load(w, cname)
            ents = dict(group_entries(c))
            mod(ents[label])
            save(w, cname, c)
            return w
        base = f"{cname}:{label}"
        G = None
        # wrong generator (different group)
        newg, order = choose_wrong_generator(g0)
        if newg is not None:
            w = with_group(lambda g: g.__setitem__("gens", newg))
            record(f"{base}:wrong-generator", f"gens[0] altered (order now {order}, stated {g0['order']})", w, cname)
        # non-normal subgroup
        nn, n = choose_nonnormal(g0)
        if nn is not None:
            w = with_group(lambda g: g.__setitem__("N_gens", nn))
            record(f"{base}:non-normal-N", f"N := a non-normal {g0['p']}-subgroup of order {n}", w, cname)
        else:
            print(f"{base}:non-normal-N  (not applicable: every {g0['p']}-subgroup of this Gamma is normal)")
        # wrong quotient generator
        Ngrp = RL.closure([tuple(x) for x in g0["N_gens"]])
        Ggrp = RL.closure([tuple(x) for x in g0["gens"]])
        if len(Ngrp) != len(Ggrp):
            w = with_group(lambda g: g.__setitem__("g", list(range(10))))
            record(f"{base}:wrong-quotient-g", "g := identity, so <g>N = N != Gamma", w, cname)
        else:
            print(f"{base}:wrong-quotient-g  (not applicable: N = Gamma)")
        # N too small
        if len(RL.closure([tuple(g0["g"])])) != len(Ggrp):
            w = with_group(lambda g: g.__setitem__("N_gens", [list(range(10))]))
            record(f"{base}:N-trivial", "N := trivial, so <g>N = <g> != Gamma", w, cname)
        # wrong order / prime / non-permutation
        w = with_group(lambda g: g.__setitem__("order", g0["order"] + 1))
        record(f"{base}:order+1", "stated order increased by one", w, cname)
        w = with_group(lambda g: g.__setitem__("p", 4))
        record(f"{base}:p-not-prime", "p := 4", w, cname)
        if len(Ngrp) > 1:
            q = next(q for q in (2, 3, 5, 7, 11) if q != g0["p"])
            w = with_group(lambda g: g.__setitem__("p", q))
            record(f"{base}:wrong-prime", f"p := {q} but |N| = {len(Ngrp)}", w, cname)

        def dup(g):
            x = list(g["gens"][0])
            x[0] = x[1]
            g["gens"][0] = x
        w = with_group(dup)
        record(f"{base}:not-a-permutation", "gens[0][0] := gens[0][1] (not a bijection)", w, cname)
        # N not inside Gamma: N := a transposition outside Gamma (p = 2)
        if g0["p"] == 2:
            t = next(((i, j) for i, j in itertools.combinations(range(10), 2) if tuple((j if x == i else i if x == j else x) for x in range(10)) not in Ggrp), None)
            if t is not None:
                tp = [t[1] if x == t[0] else t[0] if x == t[1] else x for x in range(10)]
                w = with_group(lambda g: g.__setitem__("N_gens", [tp]))
                record(f"{base}:N-outside-Gamma", f"N := the transposition {t} (not in Gamma)", w, cname)


def mutate_equations(cname):
    c0 = load(pristine, cname)
    for i, eq in enumerate(c0["equations"]):
        if eq["kind"] == "fix":
            continue
        g0 = c0["groups"][eq["group"]]
        orbs = RL.edge_orbits([tuple(x) for x in g0["gens"]])
        bad = next((RL.EDGES[k] for k in range(45) if (1 << k) not in orbs), None)
        w = fresh()
        c = load(w, cname)
        c["equations"][i]["fixed_edge"] = list(bad)
        save(w, cname, c)
        record(f"{cname}:eq{i}:fixed-edge", f"fixed_edge := {bad}, an edge that Gamma moves", w, cname)
    # remove all equations: there is nothing left to contradict
    w = fresh()
    c = load(w, cname)
    c["equations"] = []
    save(w, cname, c)
    record(f"{cname}:no-equations", "all equations removed (nothing to refute: solutions must exist)", w, cname)
    # wrong group name in an equation
    w = fresh()
    c = load(w, cname)
    c["equations"][0]["group"] = "G_does_not_exist"
    save(w, cname, c)
    record(f"{cname}:eq0:unknown-group", "equation refers to a group that is not in the certificate", w, cname)


def mutate_force(cname):
    c0 = load(pristine, cname)
    a, b = cname[len("type5_force_K"):-len(".json")].split("_K")
    # (i) assumption removed
    w = fresh()
    c = load(w, cname)
    del c["assume_out"]
    save(w, cname, c)
    record(f"{cname}:assumption-removed", "assume_out deleted (then the system is satisfiable)", w, cname)
    # (ii) wrong graph identity: K_{a'} u K_{b'} with another a'
    a2 = 5 if int(a) != 5 else 4
    w = fresh()
    c = load(w, cname)
    c["assume_out"] = [[[i, j] for i in range(a2) for j in range(i + 1, a2)] + [[i, j] for i in range(a2, 10) for j in range(i + 1, 10)]]
    save(w, cname, c)
    record(f"{cname}:wrong-identity", f"assumed graph replaced by K_{a2} u K_{10-a2} (declared identity K_{a} u K_{b})", w, cname)
    # (iii) one edge removed from the assumed graph
    w = fresh()
    c = load(w, cname)
    c["assume_out"][0] = c["assume_out"][0][1:]
    save(w, cname, c)
    record(f"{cname}:assumed-graph-1-edge-short", "one edge removed from the assumed graph", w, cname)
    # (iv) assumption turned into assume_in
    w = fresh()
    c = load(w, cname)
    c["assume_in"] = c.pop("assume_out")
    save(w, cname, c)
    record(f"{cname}:assume-in-instead", "assume_out renamed assume_in", w, cname)
    # (v) the declared identity omitted: run the row without its --assume-out-cliques flag
    w = fresh()
    rc, msg = runpy(w, ["replay/replay_cert.py", f"certs/{cname}"])
    records.append((f"{cname}:no-declaration", "run without --assume-out-cliques", "L2", True, rc, msg, rc != 0))
    print(f"{cname+':no-declaration':34s} run without the identity declaration                              | L2 {'DETECTED' if rc else 'NOT DETECTED'}: {msg[:150]}", flush=True)


def mutate_smith(cname):
    c0 = load(pristine, cname)
    # expect counts
    w = fresh()
    c = load(w, cname)
    c["expect"]["downsets"] += 1
    save(w, cname, c)
    record(f"{cname}:expect-downsets+1", "expect.downsets increased by one", w, cname)
    w = fresh()
    c = load(w, cname)
    c["expect"]["relative_survivors"] -= 1
    save(w, cname, c)
    record(f"{cname}:expect-survivors-1", "expect.relative_survivors decreased by one", w, cname)
    w = fresh()
    c = load(w, cname)
    del c["expect"]
    save(w, cname, c)
    record(f"{cname}:expect-deleted", "the pinned counts (key expect) deleted", w, cname)
    # wrong prime
    w = fresh()
    c = load(w, cname)
    c["prime"] = 3
    save(w, cname, c)
    record(f"{cname}:prime=3", "prime := 3 (the group has order 5)", w, cname)
    w = fresh()
    c = load(w, cname)
    c["prime"] = 6
    save(w, cname, c)
    record(f"{cname}:prime=6", "prime := 6 (not a prime)", w, cname)
    # wrong generator of Q: swap two images so that the order is not 5
    gens = [list(x) for x in c0["gens"]]
    for i, j in itertools.combinations(range(10), 2):
        new = [x[:] for x in gens]
        new[0][i], new[0][j] = new[0][j], new[0][i]
        if len(RL.closure([tuple(x) for x in new])) != 5:
            w = fresh()
            c = load(w, cname)
            c["gens"] = new
            save(w, cname, c)
            record(f"{cname}:wrong-Q-generator", f"gens[0] altered by swapping images {i},{j} (|<gens>| = {len(RL.closure([tuple(x) for x in new]))})", w, cname)
            break
    # a different free element of order 5 (cycle type 5^2 but a different pairing of the cycles is conjugate: use one with fixed points)
    w = fresh()
    c = load(w, cname)
    c["gens"] = [[1, 2, 3, 4, 0, 6, 5, 7, 8, 9]]
    save(w, cname, c)
    record(f"{cname}:Q-order-10", "Q := <(01234)(56)> (order 10, not a 5-group)", w, cname)
    # a claimed graph with one edge moved (a graph that is not one of the claimed ones)
    for key in ("claims_in", "claims_out"):
        for gi in (0, len(c0[key]) - 1):
            orig = [tuple(e) for e in c0[key][gi]]
            es = set(orig)
            done = False
            for e in orig:
                for f in itertools.combinations(range(10), 2):
                    if f not in es:
                        mut = [x for x in orig if x != e] + [f]
                        if not RL.isomorphic(RL.mask_of(mut), RL.mask_of(orig)):
                            w = fresh()
                            c = load(w, cname)
                            c[key][gi] = [list(x) for x in mut]
                            save(w, cname, c)
                            record(f"{cname}:{key}[{gi}]-edge-moved", f"one edge of {key}[{gi}] moved ({e} -> {f}), result not isomorphic to the original", w, cname)
                            done = True
                            break
                if done:
                    break
    # claims in/out swapped
    w = fresh()
    c = load(w, cname)
    c["claims_in"], c["claims_out"] = c["claims_out"], c["claims_in"]
    save(w, cname, c)
    record(f"{cname}:claims-swapped", "claims_in and claims_out exchanged", w, cname)
    # a claimed-in graph replaced by the complement of a claimed-out? replace first claim_in by a graph from claims_out
    w = fresh()
    c = load(w, cname)
    c["claims_in"][0] = c["claims_out"][0]
    save(w, cname, c)
    record(f"{cname}:claim-in-replaced-by-an-out-graph", "claims_in[0] := claims_out[0]", w, cname)
    # Euler group witnesses were covered by mutate_groups (euler[i] entries)


# ------------------------------------------------------------------------------------------------------------------------------------- run
t_start = time.time()
print(f"python {sys.version.split()[0]}; pristine copy {pristine}; manifest sha256 {hashlib.sha256(open(os.path.join(here, 'RELEASE-MANIFEST.json'), 'rb').read()).hexdigest()}")
print("=== group witnesses (wrong generator, non-normal N, wrong quotient, ...) and fixed edges, equations, assumptions")
for cname in type_certs:
    mutate_groups(cname)
    mutate_equations(cname)
    if cname.startswith("type5_force"):
        mutate_force(cname)
for cname in smith_certs:
    mutate_groups(cname)          # the Euler groups of lemma_z5_euler
    mutate_smith(cname)

print("=== missing input")
for cname in cert_names:
    w = fresh()
    os.remove(os.path.join(w, "certs", cname))
    record(f"{cname}:deleted", "certificate file deleted", w, cname)
for cname in type_certs[:1] + smith_certs:
    for key in (["equations", "groups", "type"] if cname in type_certs else ["gens", "claims_in", "claims_out", "prime"]):
        w = fresh()
        c = load(w, cname)
        del c[key]
        save(w, cname, c)
        record(f"{cname}:key-{key}-deleted", f"key '{key}' deleted", w, cname)
for cname in type_certs:
    w = fresh()
    c = load(w, cname)
    for k in list(c["groups"]):
        del c["groups"][k]
        break
    save(w, cname, c)
    record(f"{cname}:first-group-deleted", "the first group of the certificate deleted (equation then names a missing group)", w, cname)
for script in sorted(os.path.basename(f) for f in glob.glob(os.path.join(here, "*.py"))):
    if script in ("make_certs.py", "make_sandwich_certs.py", "make_smith_cert.py", "make_smith_euler_cert.py", "make_release_manifest.py",
                  "mutation_controls.py", "crosscheck_independent.py", "indep_graph.py", "indep_homology.py"):
        continue
    w = fresh()
    os.remove(os.path.join(w, "replay", script))
    rc1, m1 = l1(w)
    rows = [r for r in man["expected_rows"] if os.path.basename(r["script"]) == script or script == "replay_lib.py"]
    rc2, m2 = (runpy(w, [os.path.join("replay", os.path.basename(rows[0]["script"]))] + rows[0]["args"]) if rows else (None, ""))
    records.append((f"{script}:deleted", "script deleted", "L1", True, rc1, m1, rc1 != 0))
    if rows:
        records.append((f"{script}:deleted", "script deleted", "L2", True, rc2, m2, rc2 != 0))
    print(f"{script+':deleted':34s} {'script deleted':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:90]} | L2 {'DETECTED' if rc2 else 'NOT DETECTED'}: {m2[:100]}", flush=True)
w = fresh()
os.remove(os.path.join(w, "replay", "RELEASE-MANIFEST.json"))
rc1, m1 = l1(w)
records.append(("manifest:deleted", "RELEASE-MANIFEST.json deleted", "L1", True, rc1, m1, rc1 != 0))
print(f"{'manifest:deleted':34s} {'manifest file deleted':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1}", flush=True)

print("=== one-bit flips of the raw bytes of every pinned file (L1 only: the file hash must change the verdict)")
rnd = random.Random(7)
for f in [f"certs/{n}" for n in cert_names] + [f"replay/{os.path.basename(x)}" for x in sorted(glob.glob(os.path.join(here, "*.py")))]:
    w = fresh()
    path = os.path.join(w, f)
    data = bytearray(open(path, "rb").read())
    pos = rnd.randrange(len(data))
    data[pos] ^= 1 << rnd.randrange(8)
    open(path, "wb").write(bytes(data))
    rc1, m1 = l1(w)
    records.append((f"bitflip:{f}", f"byte {pos} of {f} with one bit flipped", "L1", True, rc1, m1, rc1 != 0))
    print(f"{'bitflip:'+f:34s} {'one bit flipped at byte '+str(pos):60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:110]}", flush=True)
# a bit flip inside a sha256 digest of the manifest
w = fresh()
mp = os.path.join(w, "replay", "RELEASE-MANIFEST.json")
m = json.load(open(mp))
k0 = sorted(m["certificates"])[3]
h = m["certificates"][k0]["sha256"]
m["certificates"][k0]["sha256"] = ("0" if h[0] != "0" else "1") + h[1:]
json.dump(m, open(mp, "w"), indent=1)
rc1, m1 = l1(w)
records.append(("manifest:digest-altered", f"one hex digit of the pinned digest of {k0} altered", "L1", True, rc1, m1, rc1 != 0))
print(f"{'manifest:digest-altered':34s} {'one hex digit of a pinned digest altered':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:110]}", flush=True)
# extra files
w = fresh()
shutil.copy(os.path.join(w, "certs", "type4.json"), os.path.join(w, "certs", "type4_copy.json"))
rc1, m1 = l1(w)
records.append(("extra-certificate", "an extra certificate file added to certs/", "L1", True, rc1, m1, rc1 != 0))
print(f"{'extra-certificate':34s} {'an extra certificate added':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:110]}", flush=True)
w = fresh()
open(os.path.join(w, "replay", "extra_helper.py"), "w").write("print('hi')\n")
rc1, m1 = l1(w)
records.append(("extra-script", "an extra .py file added to replay/", "L1", True, rc1, m1, rc1 != 0))
print(f"{'extra-script':34s} {'an extra script added':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:110]}", flush=True)
# wrong class: a type certificate given a 'prime' key; type key altered
w = fresh()
c = load(w, "type4.json")
c["prime"] = 5
save(w, "type4.json", c)
rc1, m1 = l1(w)
records.append(("wrong-class:prime-key", "type4.json given a 'prime' key (so a path sniffing on that key would send it to replay_smith)", "L1", True, rc1, m1, rc1 != 0))
print(f"{'wrong-class:prime-key':34s} {'type4.json given a prime key':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:110]}", flush=True)
w = fresh()
shutil.copy(os.path.join(w, "certs", "type6.json"), os.path.join(w, "certs", "type4.json"))
rc1, m1 = l1(w)
records.append(("wrong-class:type4-is-type6", "type4.json replaced by the bytes of type6.json", "L1", True, rc1, m1, rc1 != 0))
print(f"{'wrong-class:type4-is-type6':34s} {'type4.json replaced by type6.json':60s} | L1 {'DETECTED' if rc1 else 'NOT DETECTED'}: {m1[:110]}", flush=True)

# ------------------------------------------------------------------------------------------------------------------------------------- L3
def repin(w):
    """re-pin the manifest of the copy to the (mutated) files, as an adversary who also edits the manifest would"""
    mp = os.path.join(w, "replay", "RELEASE-MANIFEST.json")
    m = json.load(open(mp))
    sha = lambda p: hashlib.sha256(open(p, "rb").read().replace(b"\r\n", b"\n")).hexdigest()
    for n in list(m["certificates"]):
        p = os.path.join(w, "certs", n)
        if os.path.exists(p):
            m["certificates"][n]["sha256"] = sha(p)
    for rel in list(m["replay_scripts"]):
        p = os.path.join(w, rel)
        if os.path.exists(p):
            m["replay_scripts"][rel] = sha(p)
    json.dump(m, open(mp, "w"), indent=1)
    return m


if FULL:
    print("=== L3: manifest re-pinned to the mutated files, full replay_all.py")
    cases = []
    # a. altered order of the group of type4
    w = fresh()
    c = load(w, "type4.json")
    c["groups"]["G1279"]["order"] = 191
    save(w, "type4.json", c)
    repin(w)
    cases.append(("L3:type4-order-191", w))
    # b. force certificate with the assumption removed (and the row flag still pinned)
    w = fresh()
    c = load(w, "type5_force_K2_K8.json")
    del c["assume_out"]
    save(w, "type5_force_K2_K8.json", c)
    repin(w)
    cases.append(("L3:force-K2K8-no-assumption", w))
    # c. a claimed graph of lemma_z5_euler altered
    w = fresh()
    c = load(w, "lemma_z5_euler.json")
    c["claims_in"][0] = c["claims_in"][0][:-1] + [[0, 9]] if [0, 9] not in c["claims_in"][0] else c["claims_in"][0][:-1] + [[2, 9]]
    save(w, "lemma_z5_euler.json", c)
    repin(w)
    cases.append(("L3:z5euler-claim-altered", w))
    # d. certificate removed AND manifest/rows edited consistently (19 certificates, 22 rows): the hard-coded counts must still refuse
    w = fresh()
    os.remove(os.path.join(w, "certs", "type7_AB50.json"))
    mp = os.path.join(w, "replay", "RELEASE-MANIFEST.json")
    m = json.load(open(mp))
    del m["certificates"]["type7_AB50.json"]
    m["expected_rows"] = [r for r in m["expected_rows"] if r["name"] != "type7_AB50.json"]
    m["n_certificates"], m["n_rows"] = 19, 22
    json.dump(m, open(mp, "w"), indent=1)
    cases.append(("L3:cert-and-row-removed-consistently", w))
    # e. type key changed with the manifest re-pinned (kind 'cert' pins the type)
    w = fresh()
    c = load(w, "type4.json")
    c["type"] = 6
    save(w, "type4.json", c)
    repin(w)
    cases.append(("L3:type-key-changed", w))
    for name, w in cases:
        t0 = time.time()
        r = subprocess.run([PY, "-I", "-B", "replay/replay_all.py"], cwd=w, capture_output=True, text=True)
        lines = [l.strip() for l in r.stdout.splitlines() if l.strip()]
        fails = [l for l in lines if l.startswith("FAILED CHECK") or "rc=1" in l or l.startswith("rows that failed")]
        msg = (fails[0] if fails else lines[-1])[:200]
        records.append((name, "re-pinned manifest, full replay_all.py", "L3", True, r.returncode, msg, r.returncode != 0))
        print(f"{name:34s} {'full replay_all on a re-pinned tree':60s} | L3 {'DETECTED' if r.returncode else 'NOT DETECTED'} ({time.time()-t0:.0f} s): last={lines[-1]} ; first failure: {msg}", flush=True)

# ------------------------------------------------------------------------------------------------------------------------------------- summary
print("\n=== SUMMARY")
for layer in ("L1", "L2", "L3"):
    rs = [r for r in records if r[2] == layer]
    if not rs:
        continue
    det = sum(1 for r in rs if r[6])
    print(f"{layer}: {len(rs)} mutations, {det} detected, {len(rs) - det} NOT detected")
undetected = [r for r in records if not r[6]]
for r in undetected:
    print("NOT DETECTED:", r[0], r[1], r[2], r[5])
print(f"elapsed {time.time()-t_start:.0f} s")
shutil.rmtree(pristine, ignore_errors=True)
if any(r[3] and not r[6] for r in records):
    print("MUTATION CONTROLS: SOME NOT DETECTED")
    sys.exit(1)
print("MUTATION CONTROLS: ALL DETECTED")
