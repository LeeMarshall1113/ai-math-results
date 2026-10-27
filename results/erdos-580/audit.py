"""TRUSTED checker side.  For one tree (m, index) and the producer's clause list (.av):
  1. regenerate the tree list and the core of tree `index` (lks.py);
  2. check that every producer clause is a genuine rooted-core placement clause
     (lks.is_avoidance_clause) -- so each clause is implied by "H avoids T";
  3. write  CNF = aux-free base (lks.base_clauses) + those clauses;
  4. (optional, --solve) run CaDiCaL to produce a text LRAT proof.
The LRAT proof is then checked against the CNF by cake_lpr (and gate_lrat_check).

usage: python audit.py m index avfile out.cnf [--solve out.lrat]
"""
import sys, os, subprocess, time, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lks import trees, core, edge_vars, base_clauses_counter, is_avoidance_clause, write_dimacs, g6

CADICAL = os.environ.get("CADICAL",
                         "cadical")


def main():
    m, ti, avf, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    T = trees(m)[ti]
    cv, ce, mk = core(T, m)
    var = edge_vars(m)
    av = [list(map(int, l.split())) for l in open(avf) if l.strip()]
    bad = [c for c in av if not is_avoidance_clause(c, cv, ce, mk, m, var)]
    if bad:
        print("REJECT: non-placement clause", bad[0]); sys.exit(1)
    if len({tuple(sorted(c)) for c in av}) != len(av):
        print("note: duplicate clauses present")
    base, nv = base_clauses_counter(m, var, exact=True)
    write_dimacs(out, nv, base + av,
                 [f"Erdos 580 vertex form, n={2*m}, tree #{ti} g6={g6(T, m)}",
                  f"core vertices {cv} edges {ce} marked(->L) {mk}",
                  f"{len(base)} base clauses (deg=m on L via unary counter, S-degree<=m-1), {len(av)} placement clauses"])
    h = hashlib.sha256(open(out, "rb").read()).hexdigest()
    print(f"{os.path.basename(out)}: {len(av)} placement clauses verified; sha256 {h[:16]}")
    if "--solve" in sys.argv:
        lrat = sys.argv[sys.argv.index("--solve") + 1]
        t0 = time.time()
        r = subprocess.run([CADICAL, "-q", "--lrat=true", "--binary=false", out, lrat],
                           capture_output=True, text=True)
        print(f"  cadical rc={r.returncode} ({time.time()-t0:.1f}s) {r.stdout.strip()[-40:]}")
        if r.returncode != 20:
            sys.exit(2)


if __name__ == "__main__":
    main()
