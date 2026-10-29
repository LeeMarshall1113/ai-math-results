"""Soundness tests of encode.py: take KNOWN genuine polyhedral tori with convex quadrilateral faces,
compute their exact chirotope, and check that it satisfies every clause generated for their own map
(after mirroring if needed for SYM).  Also check the encoding of each known map is SAT.
  T9 : the 9-face torus {4,4}_(3,0) (coordinates as in the control of Seamons's verifier,
       independently re-verified by geomcheck.py)
  T11: Seamons's 11-face witness (zenodo 10.5281/zenodo.23750303, Table 1)
"""
import json, subprocess, sys, os
import encode
from geomcheck import check_torus

T9 = {'vertices': [[2,0,0],[-1,2,0],[-1,-2,0],[4,0,3],[-2,4,3],[-2,-4,3],[6,0,1],[-3,6,1],[-3,-6,1]],
      'faces': [[0,1,4,3],[1,2,5,4],[2,0,3,5],[3,4,7,6],[4,5,8,7],[5,3,6,8],[6,7,1,0],[7,8,2,1],[8,6,0,2]]}
T11 = {'vertices': [[0,0,0],[23,0,0],[0,138,0],[0,0,23],[-115,115,161],[115,69,-207],[23,138,0],
                    [0,138,23],[-23,23,161],[-15,135,-9],[47,3,-9]],
       'faces': [[0,1,6,2],[0,1,10,5],[0,2,7,3],[0,3,8,4],[0,4,9,5],[1,3,7,6],[1,3,8,10],[2,6,5,9],
                 [2,7,4,9],[4,7,6,8],[5,6,8,10]]}
CADICAL = os.environ.get("CADICAL", "cadical")

def main():
  ok = True
  for name, W in [("T9", T9), ("T11", T11)]:
      rep = check_torus(W['vertices'], W['faces'])
      print(name, "geomcheck:", rep)
      E, info = encode.build(faces=W['faces'])
      chi = encode.chirotope(W['vertices'])
      val = encode.assignment(E, chi)
      s = E.sym
      if not val[E.lit(s, 1)]:
          chi = {k: -v for k, v in chi.items()}; val = encode.assignment(E, chi)
      bad = encode.violated(E, val)
      print(name, "clauses", len(E.cls), "violated by true chirotope:", len(bad))
      ok &= (len(bad) == 0 and rep["valid"])
      path = f"cnf/{name}.cnf"; encode.write(E, path)
      r = subprocess.run([CADICAL, "-q", path], capture_output=True, text=True)
      print(name, "solver:", r.stdout.splitlines()[0] if r.stdout else r.returncode)
      ok &= ("s SATISFIABLE" in r.stdout)
  print("ALL SOUNDNESS TESTS PASSED" if ok else "FAILURE")
  sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
