"""Predictive cost of the remaining Modal cubes (Bayesian lognormal law per class, censoring and the 2 h solver limit included).

  python modal_estimate.py [--class C] [--n 122] [--ratio 1.0] [--metered 19.5] [--cap 48] [--lrat-f 1.41] [--check-f 0.13] [--cmax 7500]

Observations of the container time (seconds) of a class:
  * Modal ledger records (results/modal_ledger.jsonl) of that class: total_wall of UNSAT_VERIFIED ones; TIMEOUT ones are right-censored at their wall time;
  * the earlier LOCAL solve times without proof (A: 352, 365, 704, 1957; C: 601, 2272 and > 3600) are converted to container seconds with
    x (lrat-f + check-f) x ratio, where ratio = Modal CPU time / local CPU time of the same search (stage 0: 325 / 292-345), used only when --prior-local is given.
The law is a lognormal; (mu, sigma) get a flat prior on a grid, so the predictive total includes the parameter uncertainty.  A cube costs at most --cmax seconds
(2 h solver limit + hashing + cake_lpr), after which it is a TIMEOUT to be split or re-run.  Prints the predictive percentiles of the cost of n cubes and the number of
expected timeouts, the structural worst case (all n cubes hit the limit) and the verdict of the coordinator's rule: metered + p99 <= cap.
"""
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RATE = 0.0473 + 4 * 0.008       # $ per container-hour (1 core + 4 GiB)


def arg(name, default, cast=float):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def posterior(obs, cens):
    """flat-prior grid posterior of a lognormal (mu, sigma) given exact obs and right-censored obs (seconds)"""
    xs = [math.log(x) for x in obs]
    cs = [math.log(x) for x in cens]
    grid, ll = [], []
    for i in range(0, 150):
        mu = 5.0 + 0.05 * i
        for j in range(0, 55):
            sg = 0.3 + 0.05 * j
            v = 0.0
            for x in xs:
                v += -math.log(sg) - 0.5 * ((x - mu) / sg) ** 2
            for c in cs:
                v += math.log(max(0.5 * math.erfc((c - mu) / (sg * math.sqrt(2))), 1e-300))
            grid.append((mu, sg))
            ll.append(v)
    m = max(ll)
    w = [math.exp(v - m) for v in ll]
    t = sum(w)
    return grid, [x / t for x in w]


def main():
    cls = arg("--class", "C", str)
    n = int(arg("--n", 122))
    ratio = arg("--ratio", 1.0)
    metered = arg("--metered", 19.5)
    cap = arg("--cap", 48.0)
    lf, cf = arg("--lrat-f", 1.41), arg("--check-f", 0.13)
    cmax = arg("--cmax", 7500.0)
    obs, cens = [], []
    path = os.path.join(ROOT, "results", "modal_ledger.jsonl")
    nmod = 0
    if os.path.exists(path):
        for l in open(path):
            r = json.loads(l)
            if r["id"].startswith(cls) and r["id"][1:2].isdigit() and r.get("status") in ("UNSAT_VERIFIED", "TIMEOUT"):
                (obs if r["status"] == "UNSAT_VERIFIED" else cens).append(r["total_wall"])
                nmod += 1
    if "--prior-local" in sys.argv or nmod == 0:
        loc = {"A": ([352, 365, 704, 1957], []), "C": ([601, 2272], [3600])}[cls]
        k = (lf + cf) * ratio
        obs += [x * k for x in loc[0]]
        cens += [x * k for x in loc[1]]
    print("class %s: %d exact, %d censored observations of container seconds (Modal records %d); median of exact %.0f s" % (
        cls, len(obs), len(cens), nmod, sorted(obs)[len(obs) // 2] if obs else float("nan")))
    grid, w = posterior(obs, cens)
    # posterior summary
    mu_m = sum(g[0] * x for g, x in zip(grid, w))
    sg_m = sum(g[1] * x for g, x in zip(grid, w))
    print("posterior mean mu %.2f (median %.0f s) sigma %.2f" % (mu_m, math.exp(mu_m), sg_m))
    rnd = random.Random(7)
    cum, acc = [], 0.0
    for x in w:
        acc += x
        cum.append(acc)
    import bisect
    tot, tmo = [], []
    for _ in range(4000):
        mu, sg = grid[min(bisect.bisect_left(cum, rnd.random()), len(grid) - 1)]
        s, k = 0.0, 0
        for _ in range(n):
            t = math.exp(rnd.gauss(mu, sg))
            if t > cmax:
                t, k = cmax, k + 1
            s += t
        tot.append(s / 3600.0 * RATE)
        tmo.append(k)
    tot.sort()
    q = lambda f: tot[int(f * (len(tot) - 1))]
    print("%d cubes: predictive cost median $%.1f  mean $%.1f  p90 $%.1f  p99 $%.1f ; expected timeouts %.1f (of %d) ; structural worst case (all hit the %d s limit) $%.1f" % (
        n, q(0.5), sum(tot) / len(tot), q(0.9), q(0.99), sum(tmo) / len(tmo), n, cmax, n * cmax / 3600.0 * RATE))
    ok = metered + q(0.99) <= cap
    print("rule: metered %.2f + p99 %.2f = %.2f %s cap %.2f -> %s" % (metered, q(0.99), metered + q(0.99), "<=" if ok else ">", cap, "FITS" if ok else "DOES NOT FIT"))
    print("structural: metered %.2f + worst case %.2f = %.2f %s cap" % (metered, n * cmax / 3600.0 * RATE, metered + n * cmax / 3600.0 * RATE,
                                                                     "<=" if metered + n * cmax / 3600.0 * RATE <= cap else ">"))


if __name__ == "__main__":
    main()
