"""Exploration of the star tree with adaptive plane choice.

  python tree_adaptive.py CASE CLASS N DEPTH [CAP]

For N random first-plane sections S1 (plane L1 = <v1,v2>), compute for every other quotient line L the number of sections
of the plane over L consistent with the points fixed so far; report min / median / max, the number of planes with ZERO
consistent sections (the node is refuted without any SAT call), and continue DEPTH levels choosing the plane with the
fewest (>0) sections.
"""
import random
import sys
import time

from trivial_tree import Setup


def main():
    case, ci, N, depth = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    cap = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
    S = Setup(case, ci)
    v1, v2, v3 = S.tri
    L1 = S.quotient_line_through(v1, v2)
    rng = random.Random(7 + ci)
    secs1 = S.sections(L1, {}, cap=cap)
    print("L1 sections:", len(secs1), flush=True)
    nL = len(S.P2.lines)
    for n in range(N):
        fixed = dict(rng.choice(secs1))
        used = [L1]
        trace = []
        t0 = time.time()
        for lev in range(2, depth + 1):
            counts = {}
            dead = False
            for L in range(nL):
                if L in used:
                    continue
                # only planes that contain at least one not-yet-fixed point matter
                if all(p in fixed or p in S.pt_units for p in S.line_pts[L]):
                    continue
                c = len(S.sections(L, fixed, cap=cap))
                counts[L] = c
                if c == 0:
                    dead = True
                    break
            if dead:
                trace.append("level %d: DEAD (plane %d has no consistent section) after %d planes" % (lev, L, len(counts)))
                break
            vals = sorted(counts.values())
            Lbest = min(counts, key=lambda L: counts[L])
            trace.append("level %d: min %d median %d max %d over %d planes" % (lev, vals[0], vals[len(vals) // 2], vals[-1], len(vals)))
            sec = rng.choice(S.sections(Lbest, fixed, cap=cap))
            fixed.update(sec)
            used.append(Lbest)
        print("sample %d (%.0fs): %s" % (n, time.time() - t0, " | ".join(trace)), flush=True)


if __name__ == "__main__":
    main()
