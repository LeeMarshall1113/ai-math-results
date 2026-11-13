"""replay_all.py -- replay every certificate in ../certs and the finite checks behind the type-5 results (standard library only).

  python replay/replay_all.py        replays everything listed in replay/RELEASE-MANIFEST.json, ending with ALL VALID (exit code 0) or SOME FAILED (exit code 1)

What is checked, in this order, with NO way to switch any of it off (there is no flag and no environment variable that skips a check;
the former REPLAY_NO_MANIFEST=1 and --write-manifest were removed for the release; the manifest is written by replay/make_release_manifest.py,
which this script never runs):
  1. replay/RELEASE-MANIFEST.json exists and lists exactly 20 certificates and the replay scripts; the files present in certs/ and replay/
     are EXACTLY the listed ones (no extra, none missing) and each has the pinned SHA-256 (CRLF normalised to LF, so that git's line-ending
     conversion does not matter);
  2. every certificate has the key structure its manifest kind requires ("cert": keys type/groups/equations and no "prime"; "smith": key prime
     and no "groups"), and for kind "cert" its "type" equals the pinned type: the replay script is chosen by the MANIFEST, not by sniffing the file;
  3. exactly the 23 rows of the manifest are run, each in a fresh isolated interpreter (python -I -B -X pycache_prefix=<empty temp dir>:
     no PYTHON* environment variables, no user site, no stale .pyc), and each must return 0 and print its pinned final line;
  4. the number of rows run equals the manifest's expected_rows (23): an empty or short run can never print ALL VALID.
The manifest is anchored externally by the SHA-256 of RELEASE-MANIFEST.json itself, published with the release (this script cannot check its
own manifest against itself).  Do not run with python -O (replay_lib refuses to start, and this script refuses too).
"""
import glob
import hashlib
import json
import os
import subprocess
import sys
import tempfile

if not __debug__:
    raise SystemExit("replay_all: assertions are disabled (python -O); refusing to run")

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, ".."))
manifest_path = os.path.join(here, "RELEASE-MANIFEST.json")


def sha(path):
    """SHA-256 of the file with CRLF normalised to LF, so that git's line-ending conversion does not change it"""
    return hashlib.sha256(open(path, "rb").read().replace(b"\r\n", b"\n")).hexdigest()


def fail(msgs):
    for m in msgs:
        print("FAILED CHECK:", m)
    print("SOME FAILED")
    sys.exit(1)


if not os.path.exists(manifest_path):
    fail(["replay/RELEASE-MANIFEST.json is missing"])
man = json.load(open(manifest_path))
problems = []
certs = man["certificates"]
scripts = man["replay_scripts"]
rows_expected = man["expected_rows"]
EXPECTED_CERTS, EXPECTED_ROWS = 20, 23      # hard-coded here as well: a manifest edited to drop a certificate and its row is still refused
if len(certs) != man["n_certificates"] or len(rows_expected) != man["n_rows"] or (len(certs), len(rows_expected)) != (EXPECTED_CERTS, EXPECTED_ROWS):
    problems.append("manifest counts are not the expected 20 certificates and 23 rows")

# 1. exact file sets + hashes
have_c = sorted(os.path.basename(f) for f in glob.glob(os.path.join(root, "certs", "*")))
have_s = sorted(os.path.relpath(f, root).replace(os.sep, "/") for f in glob.glob(os.path.join(here, "*.py")))
for name in have_c:
    if name not in certs:
        problems.append(f"certs/{name}: present but not in the manifest")
for name in certs:
    if name not in have_c:
        problems.append(f"certs/{name}: listed in the manifest but missing")
    elif sha(os.path.join(root, "certs", name)) != certs[name]["sha256"]:
        problems.append(f"certs/{name}: SHA-256 differs from the manifest")
for rel in have_s:
    if rel not in scripts:
        problems.append(f"{rel}: present but not in the manifest")
for rel in scripts:
    if rel not in have_s:
        problems.append(f"{rel}: listed in the manifest but missing")
    elif sha(os.path.join(root, rel)) != scripts[rel]:
        problems.append(f"{rel}: SHA-256 differs from the manifest")

# 2. structure of every certificate, by manifest kind
for name, ent in certs.items():
    path = os.path.join(root, "certs", name)
    if not os.path.exists(path):
        continue
    try:
        c = json.load(open(path))
    except Exception as ex:                                     # noqa: BLE001
        problems.append(f"certs/{name}: not valid JSON ({ex})")
        continue
    if ent["kind"] == "cert":
        if not ("type" in c and "groups" in c and "equations" in c and "prime" not in c):
            problems.append(f"certs/{name}: kind 'cert' needs keys type/groups/equations and no 'prime'")
        elif c["type"] != ent["type"]:
            problems.append(f"certs/{name}: type {c['type']} differs from the pinned type {ent['type']}")
    elif ent["kind"] == "smith":
        if not ("prime" in c and "gens" in c and "claims_in" in c and "claims_out" in c and "groups" not in c):
            problems.append(f"certs/{name}: kind 'smith' needs keys prime/gens/claims_in/claims_out and no 'groups'")
    else:
        problems.append(f"certs/{name}: unknown kind {ent['kind']}")
for row in rows_expected:
    if row["script"] not in scripts:
        problems.append(f"row {row['name']}: script {row['script']} is not pinned")
    for a in row["args"]:
        if a.startswith("certs/") and a[len("certs/"):] not in certs:
            problems.append(f"row {row['name']}: argument {a} is not a pinned certificate")
print(f"manifest: {len(certs)} certificates, {len(scripts)} scripts, {len(rows_expected)} rows;",
      "OK" if not problems else "; ".join(problems))
print("python:", sys.version.split()[0], "| sha256(RELEASE-MANIFEST.json) =", sha(manifest_path))
if problems:
    fail(problems)

# 3. run exactly the manifest's rows
results = []
with tempfile.TemporaryDirectory() as pyc:
    for row in rows_expected:
        cmd = [sys.executable, "-I", "-B", "-X", f"pycache_prefix={pyc}", os.path.join(root, row["script"])] + [os.path.join(root, a) if a.startswith("certs/") else a for a in row["args"]]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=root)
        lines = [l for l in (r.stdout + r.stderr).strip().splitlines() if l.strip()]
        last = lines[-1].strip() if lines else "(no output)"
        ok = r.returncode == 0 and row["must_contain"] in last
        results.append((row["name"], r.returncode, last, ok))
for name, rc, last, ok in results:
    print(f"{name:30s} rc={rc}  {last}")
bad = [n for n, rc, last, ok in results if not ok]
# 4. count
if len(results) != man["n_rows"]:
    bad.append("row count")
if bad:
    print("rows that failed or did not print their pinned final line:", bad)
    print("SOME FAILED")
    sys.exit(1)
print(f"ALL VALID ({len(results)} rows: {man['n_certificates']} certificates + {len(results) - man['n_certificates']} further checks)")
sys.exit(0)
