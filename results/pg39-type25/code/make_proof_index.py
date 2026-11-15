"""Write certs/PROOF-INDEX.tsv: one line per level-1 cube (174) and per level-2 child (1,084) with the certificate data of the chain.

  python make_proof_index.py

Columns: id, level (1 = cube, 2 = child), parent, closed_by (monolithic | children | monolithic+children), modal_status, cnf_sha256 (Modal record), lrat_sha256, lrat_bytes, conflicts,
solver_cpu_s, local_status, local_identical (yes = status, conflicts, LRAT size and LRAT sha256 equal to the Modal record; no; - = not run locally), cover_tag (the completeness certificate that covers the cube's section list).
Sources: results/modal_ledger.jsonl (best record per id: a verified record is never replaced by a later failure), results/a_slice.jsonl, results/c_slice.jsonl (local runs), results/modal_jobs.json.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def load_best(path):
    out = {}
    if not os.path.exists(path):
        return out
    for l in open(path):
        try:
            r = json.loads(l)
        except Exception:
            continue
        if r.get("manual"):
            continue
        k = r.get("id") or r.get("tag")
        if k is None:
            continue
        if k in out and out[k].get("status") == "UNSAT_VERIFIED" and r.get("status") != "UNSAT_VERIFIED":
            continue
        out[k] = r
    return out


def main():
    jobs = {j["id"]: j for j in json.load(open(os.path.join(ROOT, "results", "modal_jobs.json")))}
    modal = load_best(os.path.join(ROOT, "results", "modal_ledger.jsonl"))
    local = {}
    # level-2 children of C1_c04 exist only as LOCAL records (the cube was too big for the WSL cap as a monolith): job list results/c1c04_jobs.json, records results/c1c04_local.jsonl
    local_only = {}
    local_only_recs = {}
    jp = os.path.join(ROOT, "results", "c1c04_jobs.json")
    if os.path.exists(jp):
        for j in json.load(open(jp)):
            jobs[j["id"]] = dict(j, parent=j.get("parent", "C1_c04"))
            local_only[j["id"]] = True
        local_only_recs = load_best(os.path.join(ROOT, "results", "c1c04_local.jsonl"))
    for fn in ("a_slice.jsonl", "c_slice.jsonl", "children_rerun_dry.jsonl", "full_rerun_dry.jsonl"):
        local.update(load_best(os.path.join(ROOT, "results", fn)))
    children = {}
    for jid, j in jobs.items():
        if j.get("parent"):
            children.setdefault(j["parent"], []).append(jid)
    rows = []
    n_ident = n_cmp = 0
    for jid, j in sorted(jobs.items()):
        if j["kind"] != "cube" or jid.startswith("CAL"):
            continue
        m = modal.get(jid, {})
        l = local.get(jid)
        if jid in local_only:
            l = local_only_recs.get(jid)
            m = dict(l or {}, status=(l or {}).get("status", ""), solver_cpu=(l or {}).get("solver_cpu_seconds", ""))     # the local record is the only one: it fills the record columns
        ident = "-"
        if l is not None and l.get("status") == "UNSAT_VERIFIED" and m.get("status") == "UNSAT_VERIFIED":
            ok = jid not in local_only and all(l.get(a) == m.get(b) for a, b in (("conflicts", "conflicts"), ("lrat_bytes", "lrat_bytes"), ("lrat_sha256", "lrat_sha256"), ("cnf_sha256", "cnf_sha256")))
            ident = ("yes" if ok else "NO") if jid not in local_only else "-"
            n_cmp += jid not in local_only
            n_ident += ok
        parent = j.get("parent", "")
        level = 2 if parent else 1
        if level == 1:
            kids = children.get(jid, [])
            kids_ok = bool(kids) and all((local_only_recs.get(k, {}) if k in local_only else modal.get(k, {})).get("status") == "UNSAT_VERIFIED" for k in kids)
            mono = m.get("status") == "UNSAT_VERIFIED"
            closed = "monolithic+children" if (mono and kids_ok) else ("monolithic" if mono else ("children" if kids_ok else "OPEN"))
            cover = "cover1_%s%d" % (j["case"], j["cls"]) + ("; cover2_%s" % jid if kids else "")
        else:
            closed, cover = "child of %s" % parent, "cover2_%s" % parent
        rows.append([jid, level, parent, closed, "" if jid in local_only else m.get("status", ""), m.get("cnf_sha256", ""), m.get("lrat_sha256", ""), m.get("lrat_bytes", ""), m.get("conflicts", ""),
                     m.get("solver_cpu", ""), (l or {}).get("status", ""), ident, cover, "local only" if jid in local_only else "modal"])
    out = os.path.join(ROOT, "certs", "PROOF-INDEX.tsv")
    with open(out, "w", newline="\n") as fh:
        fh.write("\t".join(["id", "level", "parent", "closed_by", "modal_status", "cnf_sha256", "lrat_sha256", "lrat_bytes", "conflicts", "solver_cpu_s", "local_status", "local_identical", "cover_tag", "source"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    cubes = [r for r in rows if r[1] == 1]
    print("PROOF-INDEX: %d cubes (open: %d), %d children; local re-runs compared: %d, bit-identical: %d" % (
        len(cubes), sum(1 for r in cubes if r[3] == "OPEN"), len(rows) - len(cubes), n_cmp, n_ident))


if __name__ == "__main__":
    main()
