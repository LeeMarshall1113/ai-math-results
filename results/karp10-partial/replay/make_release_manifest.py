"""make_release_manifest.py -- write (or verify) replay/RELEASE-MANIFEST.json.  NOT run by replay_all.py.

  python replay/make_release_manifest.py            write the manifest from the files present now (do this only for an intended release)
  python replay/make_release_manifest.py --verify   check the files named in the manifest's "release_files" section (the note source, the explicit-group
                                                    data, the pattern table, the logs-free documents); replay_all.py checks "certificates" and "replay_scripts"

The 20 expected certificate names and the 23 expected rows are written out EXPLICITLY below, not discovered by a glob: if a certificate is added or
lost, this script stops instead of silently adapting.
"""
import glob
import hashlib
import json
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, ".."))
out = os.path.join(here, "RELEASE-MANIFEST.json")


def sha(path):
    return hashlib.sha256(open(path, "rb").read().replace(b"\r\n", b"\n")).hexdigest()


# (file name, kind, type or None)
CERTS = [
    ("lemma_z5.json", "smith", None), ("lemma_z5_euler.json", "smith", None),
    ("type1.json", "cert", 1), ("type1_AB16.json", "cert", 1),
    ("type2.json", "cert", 2),
    ("type3.json", "cert", 3), ("type3_AB50.json", "cert", 3),
    ("type4.json", "cert", 4), ("type4_scout16.json", "cert", 4),
    ("type5_force_K1_K9.json", "cert", 5), ("type5_force_K2_K8.json", "cert", 5),
    ("type5_force_K3_K7.json", "cert", 5), ("type5_force_K4_K6.json", "cert", 5),
    ("type6.json", "cert", 6), ("type6_scout16.json", "cert", 6),
    ("type7.json", "cert", 7), ("type7_AB50.json", "cert", 7),
    ("type8.json", "cert", 8),
    ("type9.json", "cert", 9), ("type9_AB16.json", "cert", 9),
]
assert len(CERTS) == 20

RELEASE_FILES = [
    "note/non_evasive_ten_vertices.tex", "note/check_explicit.py", "note/z5_table.py", "note/z5_table.tex", "note/z5_patterns.json",
    "note/crosscheck_z3.py", "note/run_checks.sh", "note/explicit/type2_explicit.json", "note/explicit/type4_explicit.json",
    "note/explicit/type6_explicit.json", "note/explicit/type8_explicit.json", "PARTIAL-THEOREM.md", "TYPE5-ATTEMPT-2.md",
]


def build():
    present = sorted(os.path.basename(f) for f in glob.glob(os.path.join(root, "certs", "*")))
    assert present == sorted(c[0] for c in CERTS), f"certs/ differs from the 20 expected names: {present}"
    man = {
        "format": 1,
        "what": "pinned release of the Karp n=10 replay; sha256 of the file with CRLF normalised to LF",
        "python_tested": "3.14.0",
        "n_certificates": 20,
        "n_rows": 23,
        "certificates": {n: {"kind": k, "type": t, "sha256": sha(os.path.join(root, "certs", n))} for n, k, t in CERTS},
        "replay_scripts": {os.path.relpath(f, root).replace(os.sep, "/"): sha(f) for f in sorted(glob.glob(os.path.join(here, "*.py")))},
        "release_files": {rel: sha(os.path.join(root, rel)) for rel in RELEASE_FILES if os.path.exists(os.path.join(root, rel))},
    }
    rows = []
    for n, k, t in CERTS:
        if k == "smith":
            rows.append({"name": n, "script": "replay/replay_smith.py", "args": [f"certs/{n}"], "must_contain": "CERTIFICATE VALID"})
        else:
            must = f"type {t} is excluded (certificate valid)" if n.startswith("type5_force") is False else "the assumed membership statement is FALSE for every non-evasive P of type 5 (certificate valid)"
            args = [f"certs/{n}"]
            if n.startswith("type5_force"):
                a_, b_ = n[len("type5_force_K"):-len(".json")].split("_K")
                args.append(f"--assume-out-cliques={a_},{b_}")          # K<a>_K<10-a>: the replay checks the identity of the assumed graph
            rows.append({"name": n, "script": "replay/replay_cert.py", "args": args, "must_contain": must})
    rows += [
        {"name": "nonselfdual(z5_euler)", "script": "replay/replay_nonselfdual.py", "args": ["certs/lemma_z5_euler.json"],
         "must_contain": "NO SURVIVING PATTERN IS SELF-DUAL"},
        {"name": "absolute_acyclicity(z5_euler)", "script": "replay/check_z5_absolute.py", "args": ["certs/lemma_z5_euler.json"],
         "must_contain": "ABSOLUTE ACYCLICITY CONFIRMS THE CERTIFICATE"},
        {"name": "nine_graphs", "script": "replay/check_nine_graphs.py", "args": [], "must_contain": "THE CLAIMED GRAPHS ARE EXACTLY THE NINE NAMED GRAPHS"},
    ]
    assert len(rows) == 23
    man["expected_rows"] = rows
    return man


if "--verify" in sys.argv:
    man = json.load(open(out))
    bad = [rel for rel, h in man["release_files"].items() if not os.path.exists(os.path.join(root, rel)) or sha(os.path.join(root, rel)) != h]
    print("release_files:", len(man["release_files"]), "pinned;", "ALL MATCH" if not bad else f"DIFFER/MISSING: {bad}")
    sys.exit(1 if bad else 0)
man = build()
with open(out, "w", newline="\n") as fh:
    json.dump(man, fh, indent=1)
    fh.write("\n")
print("wrote", out, "sha256 of the file:", hashlib.sha256(open(out, "rb").read()).hexdigest())
