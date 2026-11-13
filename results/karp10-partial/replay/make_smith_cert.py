"""make_smith_cert.py -- build certs/lemma_z5.json from stage2/data2/z5_forced.json (output of stage2/local_patterns analysis).
The forced-in / forced-out graphs are only CLAIMS; replay_smith.py re-derives everything and checks them."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from k10common import edges_of  # noqa: E402

d = json.load(open(os.path.join(HERE, "..", "stage2", "data2", "z5_forced.json")))
cert = {
    "type": 5,
    "prime": 5,
    "claim": ("every non-evasive nontrivial monotone P of type 5 contains C_10, the sun C_5+K_1 (pentagon with pendant edges), the Petersen "
              "graph, the pentagonal prism C_5xK_2, K_5+K_1 (K_5 with pendant edges), (K_5 u C_5) + perfect matching and K_5xK_2, and "
              "contains none of their complements (Smith theory for the cyclic group <(01234)(56789)> of order 5)"),
    "gens": [[1, 2, 3, 4, 0, 6, 7, 8, 9, 5]],
    "claims_in": [[list(e) for e in edges_of(m)] for m in d["always"]],
    "claims_out": [[list(e) for e in edges_of(m)] for m in d["never"]],
    "expect": {"downsets": 2539, "relative_survivors": 16},
}
json.dump(cert, open(os.path.join(HERE, "..", "certs", "lemma_z5.json"), "w"), indent=1)
print(len(cert["claims_in"]), "claimed-in graphs,", len(cert["claims_out"]), "claimed-out graphs written")
