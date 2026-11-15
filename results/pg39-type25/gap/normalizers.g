# Generators of the normaliser N_G(<a>) of each prime-order representative, as permutations of the
# 820 points (1-based images, one generator per line).  Output: normalizers.txt
# Usage (from this directory):  gap -q -o 1500m normalizers.g
Read("pgammal49_data.g");
G := Group(gens);;
out := OutputTextFile("normalizers.txt", false);;
SetPrintFormattingStatus(out, false);;
for i in [1..Length(reps)] do
  a := reps[i];
  H := Group(a);
  N := Normalizer(G, H);
  gl := SmallGeneratingSet(N);
  PrintTo(out, "CLASS ", repnames[i], " order ", Size(N), " gens ", Length(gl), "\n");
  for g in gl do
    AppendTo(out, JoinStringsWithSeparator(List([1..820], j -> String(j^g)), ","), "\n");
  od;
  Print(repnames[i], "  |N_G(<a>)| = ", Size(N), "  |C_G(a)| = ", Size(Centralizer(G, a)), "  #gens ", Length(gl), "\n");
od;
CloseStream(out);
QUIT_GAP(0);
