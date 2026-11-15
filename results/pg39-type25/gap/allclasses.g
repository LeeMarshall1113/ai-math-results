# Step 2b: ALL conjugacy classes of PGammaL(4,9) (GAP's own algorithm, class sizes are
# checked by GAP to sum to |G|), then the classes of elements of prime order, grouped into
# classes of subgroups of prime order via power maps.  Compared to the 13 representatives.
# Usage (from this directory):  gap -q -o 1500m allclasses.g > ../results/gap_allclasses.txt
Read("pgammal49_data.g");
G := Group(gens);;
sz := Size(G);;
t0 := Runtime();;
cc := ConjugacyClasses(G);;
Print("number of conjugacy classes: ", Length(cc), "   time ", Runtime()-t0, " ms\n");
Print("sum of class sizes = |G| : ", Sum(cc, Size) = sz, "\n");
primeclasses := Filtered([1..Length(cc)], i -> IsPrimeInt(Order(Representative(cc[i]))));;
Print("classes of elements of prime order: ", Length(primeclasses), "\n");
# group by subgroup class: c ~ c^k classes
done := [];;
subgroupclasses := [];;
for i in primeclasses do
  if i in done then continue; fi;
  r := Order(Representative(cc[i]));
  members := [];
  for k in [1..r-1] do
    x := Representative(cc[i])^k;
    for j in primeclasses do
      if Order(Representative(cc[j])) = r and x in cc[j] then AddSet(members, j); fi;
    od;
  od;
  UniteSet(done, members);
  Add(subgroupclasses, [r, members]);
od;
Print("classes of subgroups of prime order: ", Length(subgroupclasses), "\n");
for sc in subgroupclasses do
  rep := Representative(cc[sc[2][1]]);
  Print("prime ", sc[1], "  element classes ", sc[2], "  |C(a)| = ", Size(Centralizer(G, rep)),
        "  fixed points ", 820 - NrMovedPoints(rep), "  #elements in the union ", Sum(sc[2], j -> Size(cc[j])), "\n");
od;
QUIT_GAP(0);
