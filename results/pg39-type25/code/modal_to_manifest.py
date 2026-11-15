"""Append the verified Modal cube records (results/modal_ledger.jsonl) to certs/manifest.jsonl (via = modal), once per cube.
  python modal_to_manifest.py
The LRAT itself stays in the container (deleted with it); the certificate is (CNF sha256, LRAT sha256, size, cake_lpr verdict) and the fact that the run is
deterministic: the same binary on the same CNF gives the same LRAT (checked against the local slice cube by cube for the A cubes, modal_vs_local.py).
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def main():
    jobs = {j["id"]: j for j in json.load(open(os.path.join(ROOT, "results", "modal_jobs.json")))}
    man = os.path.join(ROOT, "certs", "manifest.jsonl")
    have = set()
    if os.path.exists(man):
        for l in open(man):
            try:
                r = json.loads(l)
                if r.get("via") == "modal":
                    have.add(r["tag"])
            except Exception:
                pass
    new = 0
    with open(man, "a") as out:
        for l in open(os.path.join(ROOT, "results", "modal_ledger.jsonl")):
            r = json.loads(l)
            j = jobs.get(r["id"])
            if j is None or j["kind"] != "cube" or r["id"].startswith("CAL") or r.get("status") != "UNSAT_VERIFIED" or r["id"] in have:
                continue
            rec = dict(tag=r["id"], kind="cube", via="modal", status="UNSAT_VERIFIED", cake_lpr=r.get("verdict"), cnf_sha256=r["cnf_sha256"], lrat_sha256=r["lrat_sha256"],
                       lrat_bytes=r["lrat_bytes"], conflicts=r["conflicts"], solver_cpu_seconds=r["solver_cpu"], cake_seconds=r["cake_wall"], container_seconds=r["total_wall"],
                       case=j["case"], cube_class=j["cls"], cube=j["cube"], orbit_size=j["orbit_size"], k=j["k"], time=r.get("time"), cost_usd=r.get("est_cost_usd"))
            out.write(json.dumps(rec) + "\n")
            have.add(r["id"])
            new += 1
    print("appended %d Modal cube records to certs/manifest.jsonl (%d Modal cubes in total)" % (new, len(have)))


if __name__ == "__main__":
    main()
