"""make_certs.py -- assemble the certificate JSON files from the GAP outputs (gap/subgroups_s10.txt, gap/witness.out).
The certificates contain only: type, group generators, G^1 witnesses, and the list of equations used.
Everything else is recomputed by replay_cert.py.  (0-based permutations: g[i] = image of vertex i.)"""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def parse_perm(s):
    s = s.strip()
    if s == "[1..10]":
        return list(range(10))
    return [int(x) - 1 for x in s.strip("[]").split(",")]


def load_groups():
    gs = {}
    for line in open(os.path.join(ROOT, "gap", "subgroups_s10.txt")):
        if line.startswith("R|"):
            _, idx, order, pp, ex, cg, gens = line.rstrip("\n").split("|")
            gl = eval(gens)
            gs[int(idx)] = dict(order=int(order), gens=[[x - 1 for x in g] for g in gl])
    return gs


def load_witnesses():
    w = {}
    txt = open(os.path.join(ROOT, "gap", "witness.out")).read()
    for m in re.finditer(r"W\|(G\d+)\|(\d+)\|(\d+)\|(\d+)\|(\[\[.*?\]\]|\[\[\]\])\|(\[[0-9,.]*\])", txt):
        name, order, p, nord, ngens, g = m.groups()
        ng = re.findall(r"\[[0-9,]+\]", ngens)
        w[name] = dict(order=int(order), p=int(p), N_order=int(nord), N_gens=[parse_perm(x) for x in ng], g=parse_perm(g))
    return w


def group_entry(name, gs, ws):
    idx = int(name[1:])
    e = dict(gens=gs[idx]["gens"], order=gs[idx]["order"], p=ws[name]["p"], N_gens=ws[name]["N_gens"], g=ws[name]["g"])
    return e


def main():
    gs, ws = load_groups(), load_witnesses()
    assert set(ws) >= {"G836", "G97", "G1279", "G1132"}, ws.keys()
    os.makedirs(os.path.join(ROOT, "certs"), exist_ok=True)

    def edge_of_fixed(name):
        # the single-edge orbit of the group: recompute with the same union-find as the replay
        sys.path.insert(0, HERE)
        from replay_lib import edge_orbits, EDGES
        orbs = edge_orbits([tuple(g) for g in gs[int(name[1:])]["gens"]])
        fx = [o for o in orbs if bin(o).count("1") == 1]
        assert len(fx) == 1
        return list(EDGES[fx[0].bit_length() - 1])

    certs = {
        1: dict(groups=["G97"], eqs=[("G97", "lk")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 1"),
        9: dict(groups=["G97"], eqs=[("G97", "lk")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 9"),
        3: dict(groups=["G836"], eqs=[("G836", "fix")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 3"),
        7: dict(groups=["G836"], eqs=[("G836", "fix")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 7"),
        4: dict(groups=["G1279"], eqs=[("G1279", "lk")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 4"),
        6: dict(groups=["G1279"], eqs=[("G1279", "lk")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 6"),
        2: dict(groups=["G1132", "G1279"], eqs=[("G1132", "lk"), ("G1279", "lk")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 2"),
        8: dict(groups=["G1132", "G1279"], eqs=[("G1132", "lk"), ("G1279", "lk")], claim="no non-evasive nontrivial monotone P of Angel-Borja type 8"),
    }
    for t, c in certs.items():
        cert = dict(type=t, claim=c["claim"], groups={n: group_entry(n, gs, ws) for n in c["groups"]}, equations=[])
        for (n, kind) in c["eqs"]:
            eq = dict(group=n, kind=kind)
            if kind != "fix":
                eq["fixed_edge"] = edge_of_fixed(n)
            cert["equations"].append(eq)
        json.dump(cert, open(os.path.join(ROOT, "certs", f"type{t}.json"), "w"), indent=1)
        print("wrote certs/type%d.json" % t)


if __name__ == "__main__":
    main()
