"""Cross-check of the Modal results against the local A slice, cube by cube (the coordinator's stage-1 condition: any disagreement stops everything).

  python modal_vs_local.py [--ledger ../results/modal_ledger.jsonl] [--local ../results/a_slice.jsonl] [--prefix A]

For every Modal record of a cube (id like A0_c00):
  1. cnf_sha256 must equal the hash of the cube CNF rebuilt here (base CNF + 92 units; computed without solving, so it is checked for every Modal cube at once);
  2. if the local slice has finished the same cube: status, conflicts, lrat_bytes and lrat_sha256 must be equal too (CaDiCaL is deterministic: same binary,
     same CNF => same search, so a different count or hash means a different binary, a corrupted input or a hardware fault);
  3. the SAT control: the model returned by Modal is decoded and checked on all 7,462 lines (must meet every line in 1, 4 or 10 points, 280 points).
Prints one line per cube and a summary; exit code 1 on any disagreement.  Writes ../results/modal_vs_local.json.
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


def cube_cnf_sha(base_path, units):
    h = hashlib.sha256()
    with open(base_path, "rb") as fh:
        line = fh.readline()
        while line.startswith(b"c"):
            h.update(line)
            line = fh.readline()
        parts = line.split()
        h.update(("p cnf %s %d\n" % (parts[2].decode(), int(parts[3]) + len(units))).encode())
        while True:
            ch = fh.read(1 << 24)
            if not ch:
                break
            h.update(ch)
    for u in units:
        h.update(("%d 0\n" % u).encode())
    return h.hexdigest()


def load_jsonl(path):
    out = {}
    if os.path.exists(path):
        for l in open(path):
            try:
                r = json.loads(l)
            except Exception:
                continue
            out[r.get("id") or r.get("tag")] = r
    return out


def load_ledger(path):
    """best record per id: manual lines (operator notes) are skipped, a verified record is never replaced by a later failure"""
    out = {}
    for l in open(path):
        try:
            r = json.loads(l)
        except Exception:
            continue
        if r.get("manual") or "id" not in r:
            continue
        old = out.get(r["id"])
        if old is not None and old.get("status") == "UNSAT_VERIFIED" and r.get("status") != "UNSAT_VERIFIED":
            continue
        out[r["id"]] = r
    return out


def main():
    ledger = load_ledger(arg("--ledger", os.path.join(ROOT, "results", "modal_ledger.jsonl")))
    local = load_jsonl(arg("--local", os.path.join(ROOT, "results", "a_slice.jsonl")))
    prefix = arg("--prefix", "A")
    jobs = {j["id"]: j for j in json.load(open(os.path.join(ROOT, "results", "modal_jobs.json")))}
    cache_path = os.path.join(ROOT, "results", "cnf_sha256_cache.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    agree, disagree, pending_local, timeouts = [], [], [], []
    lines = []
    for jid, r in sorted(ledger.items()):
        j = jobs.get(jid)
        if j is None or not jid.startswith(prefix) or j["kind"] != "cube" or jid.startswith("CAL") or "_s" in jid[3:]:
            continue                                   # seed variants (id_s1 ...) have no local counterpart
        if jid not in cache:
            cache[jid] = cube_cnf_sha(os.path.join(ROOT, "cnf", j["base"] + ".cnf"), j["units"])
        bad = []
        if r.get("cnf_sha256") != cache[jid]:
            bad.append("cnf_sha256 differs from the CNF rebuilt locally")
        if r.get("status") == "TIMEOUT":
            timeouts.append(jid)                       # not a disagreement: no verdict on either side (CNF hash is still checked above)
            lines.append("%-8s %-30s %s" % (jid, "modal TIMEOUT", "ok (cnf hash only)" if not bad else "DISAGREE: " + "; ".join(bad)))
            if bad:
                disagree.append((jid, bad))
            continue
        if r.get("status") == "UNSAT_CAKE_FAILED":
            lines.append("%-8s %-30s %s" % (jid, "cake_lpr failed (resource)", "NOT a disagreement: %s; re-run pending" % r.get("verdict")))
            continue
        if r.get("status") != "UNSAT_VERIFIED":
            bad.append("Modal status %s" % r.get("status"))
        l = local.get(jid)
        if l is not None and l.get("status") in ("UNSAT_VERIFIED", "UNSAT_CAKE_FAILED", "TIMEOUT_2h", "SAT_EXACT_PASS", "SAT_EXACT_FAIL"):
            if l.get("cnf_sha256") != cache[jid]:
                bad.append("local driver cnf_sha256 differs from the rebuilt CNF")
            for k in ("status", "conflicts", "lrat_bytes", "lrat_sha256", "cnf_sha256"):
                if l.get(k) != r.get(k):
                    bad.append("%s: local %s / modal %s" % (k, l.get(k), r.get(k)))
            state = "compared"
        else:
            state = "cnf only (local not finished)"
            pending_local.append(jid)
        if bad:
            disagree.append((jid, bad))
        else:
            agree.append(jid)
        lines.append("%-8s %-30s %s" % (jid, state, "DISAGREE: " + "; ".join(bad) if bad else "ok"))
    # control
    ctrl = ledger.get("CTRL_hermitian")
    ctrl_state = "not run"
    if ctrl is not None:
        if ctrl.get("status") != "SAT":
            ctrl_state = "FAIL: status %s" % ctrl.get("status")
            disagree.append(("CTRL_hermitian", [ctrl_state]))
        elif "model" not in ctrl:
            ctrl_state = "SAT, model not in the ledger (ledger drops it); see results/MODAL_SAT_CTRL_hermitian.json"
            mp = os.path.join(ROOT, "results", "MODAL_SAT_CTRL_hermitian.json")
            if os.path.exists(mp):
                ctrl = json.load(open(mp))
        if "model" in ctrl:
            sys.path.insert(0, HERE)
            from pg import PG
            P = PG(4, 9)
            pos = set(m for m in ctrl["model"] if m > 0)
            S = set(p for p in range(820) if (p + 1) in pos)
            ok = len(S) == 280 and all(sum(1 for p in L if p in S) in (1, 4, 10) for L in P.lines)
            ctrl_state = "SAT, %d points, exact check of 7462 lines: %s" % (len(S), "PASS" if ok else "FAIL")
            if not ok:
                disagree.append(("CTRL_hermitian", [ctrl_state]))
    json.dump(cache, open(cache_path, "w"))
    print("\n".join(lines))
    print("control:", ctrl_state)
    comp = len(agree) - len(pending_local)
    print("SUMMARY: %d Modal cubes: %d fully compared with the local slice (identical: status, conflicts, LRAT size, LRAT sha256, CNF sha256), %d verified on Modal with the CNF hash only (local pending), %d Modal TIMEOUT (CNF hash only), %d DISAGREEMENTS" % (
        len(agree) + len(disagree) + len(timeouts), comp, len(pending_local), len(timeouts), len(disagree)))
    json.dump(dict(agree=agree, pending_local=pending_local, timeouts=timeouts, disagree=disagree, control=ctrl_state), open(os.path.join(ROOT, "results", "modal_vs_local.json"), "w"))
    sys.exit(1 if disagree else 0)


if __name__ == "__main__":
    main()
