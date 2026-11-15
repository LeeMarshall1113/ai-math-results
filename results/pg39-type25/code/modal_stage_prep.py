"""Stage the local files for the Modal run (no network, no spend).  Runs under a lease (python, < 0.5 GB).

  modal_stage/bases/full_{A,C}_c<i>.cnf.gz        gzipped base CNFs of the 10 star cases (units are added per cube inside the container)
  modal_stage/bases/control_hermitian.cnf.gz      SAT control (complete CNF, no extra units)
  results/modal_jobs.json                         one job per cube (A: 52, C: 122), the control, and the calibration cube CAL_b_A1s1

A job is {id, base (no suffix), units, tlim, lrat, kind, ...}; the CAL job reproduces cnf/b_A1s1.cnf exactly (base full_A_c1 + its 92 units),
so a Modal run of it must return conflicts 2,427,498 and a binary LRAT of 592,548,314 bytes with sha256 92b303d6... (local reference run, NOTES 5.5).
"""
import gzip
import hashlib
import json
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
STAGE = os.path.join(ROOT, "modal_stage", "bases")


def gz(src, dst):
    if os.path.exists(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
        return
    with open(src, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo, 1 << 22)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 24), b""):
            h.update(ch)
    return h.hexdigest()


def main():
    os.makedirs(STAGE, exist_ok=True)
    jobs = []
    shas = {}
    for case, n in (("A", 3), ("C", 7)):
        for ci in range(n):
            name = "full_%s_c%d" % (case, ci)
            src = os.path.join(ROOT, "cnf", name + ".cnf")
            gz(src, os.path.join(STAGE, name + ".cnf.gz"))
            shas[name] = sha(src)
    src = os.path.join(ROOT, "cnf", "control_hermitian.cnf")
    gz(src, os.path.join(STAGE, "control_hermitian.cnf.gz"))
    shas["control_hermitian"] = sha(src)
    for fn in ("a_slice_jobs.json", "c_slice_jobs.json"):
        for j in json.load(open(os.path.join(ROOT, "results", fn))):
            if j["kind"] != "cube":
                continue
            jobs.append(dict(id=j["id"], kind="cube", base=j["base"][:-4], units=j["units"], tlim=7200, lrat="bin", case=j["case"], cls=j["cls"],
                             cube=j["cube"], orbit_size=j["orbit_size"], k=j["k"]))
    jobs.append(dict(id="CTRL_hermitian", kind="control", base="control_hermitian", units=[], tlim=900, lrat="none"))
    units = []
    with open(os.path.join(ROOT, "cnf", "b_A1s1.cnf")) as fh:
        tail = fh.readlines()[-92:]
    for line in tail:
        parts = line.split()
        assert len(parts) == 2 and parts[1] == "0", line
        units.append(int(parts[0]))
    jobs.append(dict(id="CAL_b_A1s1", kind="cube", base="full_A_c1", units=units, tlim=3000, lrat="bin", case="A", cls=1, cube=-1, orbit_size=0, k=287,
                     expect=dict(conflicts=2427498, lrat_bytes=592548314, lrat_sha256="92b303d6c95286d5dd79510aabf8492b4ef915fafe1eb6cecdee07a348b16699",
                                 cnf_sha256="80cb367b29be45513add4df4080832476032658da814830fa8b78d6e18a6cc92")))
    json.dump(jobs, open(os.path.join(ROOT, "results", "modal_jobs.json"), "w"))
    json.dump(shas, open(os.path.join(ROOT, "results", "modal_base_sha256.json"), "w"), indent=1)
    tot = sum(os.path.getsize(os.path.join(STAGE, f)) for f in os.listdir(STAGE))
    print("jobs: %d (cubes %d + control + calibration); staged bases: %d files, %.1f MB gz" % (len(jobs), len(jobs) - 2, len(os.listdir(STAGE)), tot / 1e6), flush=True)


if __name__ == "__main__":
    main()
