"""Combine exact K4^- certificates into one bound on pi(K4^-) (pure Python, Fractions).

    python check_bound.py BOUND certs/cert_plain.json.gz certs/cert_X.json.gz [...]

Each certificate is checked in full by check_cert.check (catalogue completeness, PSD
Gram factors, exact coefficients for all 8157 seven-vertex K4^- -free graphs).

Logic.  Let phi be an extremal limit, pi = phi(e).  phi can be taken vertex-regular
(every vertex has degree pi almost surely; see REPORT.md), so for each interval
certificate [t1, t2] with exact value v:  pi in [t1, t2]  =>  pi <= v.
The plain certificate (no interval rows) gives pi <= P unconditionally.
If the intervals cover [t0, P] (sorted, contiguous) then
    pi <= max(t0, max_i v_i),
and an interval with v_i < t1_i contributes nothing (it is excluded outright).
The script checks the cover and that this maximum is <= BOUND.
"""
import gzip, json, sys
from fractions import Fraction
from check_cert import check


def load(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        return json.load(fh)


def main():
    bound = Fraction(sys.argv[1])
    plain = None
    itvs = []
    for p in sys.argv[2:]:
        cert = load(p)
        print(f"=== {p}")
        v = check(cert)
        it = cert.get("interval_stat")
        if it:
            itvs.append((Fraction(it["t1"]), Fraction(it["t2"]), v, p))
        else:
            plain = v if plain is None else min(plain, v)
    assert plain is not None, "need a plain certificate"
    itvs.sort()
    print(f"\nplain: pi <= {plain} = {float(plain):.10f}")
    assert itvs, "no interval certificates"
    t0 = itvs[0][0]
    reach = t0
    worst = t0
    for t1, t2, v, p in itvs:
        assert t1 <= reach, f"gap in cover before {float(t1)}"
        reach = max(reach, t2)
        tag = "excluded" if v < t1 else f"pi <= {float(v):.10f}"
        print(f"  [{float(t1):.6f}, {float(t2):.6f}]  value {float(v):.10f}  {tag}   ({p})")
        if v >= t1:
            worst = max(worst, v)
    assert reach >= plain, f"cover ends at {float(reach)} below the plain bound {float(plain)}"
    print(f"cover [{float(t0):.6f}, {float(reach):.6f}] contains [t0, plain bound]")
    print(f"=> pi(K4^-) <= {worst} = {float(worst):.12f}")
    assert worst <= bound, "BOUND NOT ESTABLISHED"
    print(f"BOUND {bound} = {float(bound):.12f}  ESTABLISHED")


if __name__ == "__main__":
    main()
