# Step 2a: group order, centralizers, class sizes and pairwise (non-)conjugacy of the
# 13 representatives of prime-order subgroup classes of PGammaL(4,9).
# Input: pgammal49_data.g (written by code/groups.py).  Usage (from this directory):
#   gap -q -o 1500m reps.g
Read("pgammal49_data.g");
G := Group(gens);;
sz := Size(G);;
q := 9;;
# |PGL(4,q)| = q^6 (q-1)^-1 * prod_{i=1}^{4} (q^i-1) / ... : |GL(4,q)| = q^6 prod (q^i-1), PGL = GL/(q-1)
szGL := q^6 * Product([1..4], i -> q^i - 1);;
szPGL := szGL / (q-1);;
Print("Size(G) = ", sz, "   expected |PGammaL(4,9)| = ", 2*szPGL, "   equal: ", sz = 2*szPGL, "\n");
Print("degree ", NrMovedPoints(G), "  primes ", PrimeDivisors(sz), "\n");
for i in [1..Length(reps)] do
  a := reps[i];
  C := Centralizer(G, a);
  Print(repnames[i], "  order ", Order(a), "  |C_G(a)| = ", Size(C), "  class size ", sz/Size(C),
        "  fixed pts ", NrMovedPoints(G) - NrMovedPoints(a), "\n");
od;
# pairwise: is <reps[i]> conjugate to <reps[j]> (same prime order)?
for i in [1..Length(reps)] do
  for j in [i+1..Length(reps)] do
    if Order(reps[i]) = Order(reps[j]) then
      r := Order(reps[i]);
      conj := false;
      for k in [1..r-1] do
        if IsConjugate(G, reps[i], reps[j]^k) then conj := true; fi;
      od;
      Print(repnames[i], " ~ ", repnames[j], " as subgroups: ", conj, "\n");
    fi;
  od;
od;
QUIT_GAP(0);
