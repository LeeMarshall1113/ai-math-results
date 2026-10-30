"""Solver and checker calls shared by the drivers (run_all.py, cnc.py) and the replay (replay.py).

Environment:
  CADICAL       CaDiCaL binary (default "cadical"). Must write LRAT (--lrat=true).
  CADICAL_OPTS  extra options (default "--congruence=false": CaDiCaL 3.0.1 can abort with rc 3 in
                congruence closure on some of these CNFs; the option avoids it).
  LRAT_CHECK    LRAT checker command (default "cake_lpr"), called as  $LRAT_CHECK file.cnf file.lrat
                and expected to print "s VERIFIED UNSAT". For cake_lpr on the larger leaves use
                e.g. LRAT_CHECK="cake_lpr --CML_HEAP_SIZE=2500 --CML_STACK_SIZE=1000".
"""
import hashlib
import os
import shlex
import subprocess

CADICAL = os.environ.get("CADICAL", "cadical")
CADICAL_OPTS = shlex.split(os.environ.get("CADICAL_OPTS", "--congruence=false"))
LRAT_CHECK = shlex.split(os.environ.get("LRAT_CHECK", "cake_lpr"))


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def solve(cnf, lrat, timeout=None, extra=()):
    """CaDiCaL with a textual LRAT proof. Returns (returncode or None on timeout, stdout)."""
    try:
        r = subprocess.run([CADICAL, "-q", "--lrat=true", "--binary=false", *CADICAL_OPTS, *extra, cnf, lrat],
                           capture_output=True, text=True, timeout=timeout)
        return r.returncode, r.stdout
    except subprocess.TimeoutExpired:
        return None, ""


def check_lrat(cnf, lrat):
    """Run the LRAT checker. Returns (verified, last output lines)."""
    c = subprocess.run([*LRAT_CHECK, cnf, lrat], capture_output=True)
    out = (c.stdout + c.stderr).decode("utf-8", "replace").replace("\x00", "")
    return "s VERIFIED UNSAT" in out, out.strip()[-300:]
