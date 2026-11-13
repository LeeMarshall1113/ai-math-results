"""make_sandwich_certs.py -- build replayable certificates 'type 5  =>  K_a u K_{10-a} in P' (a = 1..4) from the
equation-level cores printed in logs/core_force2.out.  Needs GAP witnesses for every group used (gap/witness_auto.out,
produced by this script's --gap step).  Usage:
   python make_sandwich_certs.py --gap      # writes gap/witness_auto.g (run it in WSL) 
   python make_sandwich_certs.py --build    # reads gap/witness_auto.out and logs/core_force2.out, writes certs/
"""
import json, os, re, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)


def cores():
    out = {}
    for line in open(os.path.join(ROOT, "logs", "core_force2.out")):
        m = re.match(r"a=(\d) class \d+: minimal equation core: (\[.*?\])", line)
        if m:
            out[int(m.group(1))] = eval(m.group(2))
    return out


def groups_in(core_list):
    return sorted({int(re.match(r"G(\d+):", c).group(1)) for c in core_list})


def load_subgroups():
    gs = {}
    for line in open(os.path.join(ROOT, "gap", "subgroups_s10.txt")):
        if line.startswith("R|"):
            _, idx, order, pp, ex, cg, gens = line.rstrip("\n").split("|")
            gs[int(idx)] = dict(order=int(order), exact=int(ex), congr=cg, gens=eval(gens))
    return gs


def gap_step():
    gs = load_subgroups()
    need = sorted({i for core in cores().values() for i in groups_in(core)})
    lines = ["SizeScreen([100000,24]);",
             'lstr := function(l) return Concatenation("[", JoinStringsWithSeparator(List(l, String), ","), "]"); end;',
             "groups := ["]
    for i in need:
        lines.append(f'  ["G{i}", {json.dumps(gs[i]["gens"])}],')
    lines[-1] = lines[-1].rstrip(",")
    lines += ["];",
              "for entry in groups do",
              "  G := Group(List(entry[2], PermList));",
              "  found := false;",
              "  for p in Set(Factors(Size(G))) do",
              "    if found then continue; fi;",
              "    N := PCore(G, p);",
              "    for g in G do",
              "      if Size(ClosureGroup(N, g)) = Size(G) then",
              "        ngens := List(SmallGeneratingSet(N), x -> lstr(ListPerm(x, 10)));",
              '        Print("W|", entry[1], "|", Size(G), "|", p, "|", Size(N), "|", Concatenation("[", JoinStringsWithSeparator(ngens, ","), "]"), "|", lstr(ListPerm(g, 10)), "\\n");',
              "        found := true; break;",
              "      fi;",
              "    od;",
              "  od;",
              "  if not found then Print(\"W|\", entry[1], \"|NOT_G1\\n\"); fi;",
              "od;",
              "QUIT;"]
    open(os.path.join(ROOT, "gap", "witness_auto.g"), "w").write("\n".join(lines) + "\n")
    print("wrote gap/witness_auto.g for groups", need)


def parse_perm(s):
    s = s.strip()
    if s.startswith("[1..") or s == "[1..10]":
        return list(range(10))
    return [int(x) - 1 for x in s.strip("[]").split(",")]


def build_step():
    gs = load_subgroups()
    ws = {}
    for line in open(os.path.join(ROOT, "gap", "witness_auto.out")):
        if line.startswith("W|") and "NOT_G1" not in line:
            parts = line.rstrip("\n").split("|")
            name, order, p, nord, ngens, g = parts[1:7]
            ng = re.findall(r"\[[0-9,]+\]", ngens)
            ws[name] = dict(order=int(order), p=int(p), N_gens=[parse_perm(x) for x in ng], g=parse_perm(g))
    cs = cores()
    sys.path.insert(0, HERE)
    from replay_lib import edge_orbits, EDGES
    for a, core in cs.items():
        cert = dict(type=5, claim=f"every non-evasive nontrivial monotone P of type 5 contains K_{a} u K_{10-a}",
                    groups={}, equations=[],
                    assume_out=[[list(e) for e in itertools.chain(itertools.combinations(range(a), 2), itertools.combinations(range(a, 10), 2))]])
        for name in core:
            m = re.match(r"G(\d+):(fix|lk|del)(?:@(\d+))?", name)
            idx, kind, orb = int(m.group(1)), m.group(2), m.group(3)
            gname = f"G{idx}"
            if gname not in ws:
                print("no G^1 witness for", gname, "- skipping certificate a =", a); cert = None; break
            if gname not in cert["groups"]:
                cert["groups"][gname] = dict(gens=[[x - 1 for x in g] for g in gs[idx]["gens"]], order=gs[idx]["order"],
                                             p=ws[gname]["p"], N_gens=ws[gname]["N_gens"], g=ws[gname]["g"])
            eq = dict(group=gname, kind=kind)
            if kind != "fix":
                orbs = edge_orbits([tuple(x - 1 for x in g) for g in gs[idx]["gens"]])
                eq["fixed_edge"] = list(EDGES[orbs[int(orb)].bit_length() - 1])
            cert["equations"].append(eq)
        if cert:
            json.dump(cert, open(os.path.join(ROOT, "certs", f"type5_force_K{a}_K{10-a}.json"), "w"), indent=1)
            print(f"wrote certs/type5_force_K{a}_K{10-a}.json with {len(cert['equations'])} equations")


if __name__ == "__main__":
    if "--gap" in sys.argv:
        gap_step()
    else:
        build_step()
