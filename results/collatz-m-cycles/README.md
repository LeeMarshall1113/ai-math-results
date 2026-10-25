# Collatz m-cycles: no nontrivial cycle with fewer than 112 local minima

**Full repository (Lean sources, method, computed tables): https://github.com/LeeMarshall1113/collatz-m-cycles**
The Lean tree is deliberately not duplicated here.

## Claim
Every nontrivial positive Collatz cycle has at least 112 local minima (no m-cycle with m <= 111).
- **m <= 111: Lean-verified**, kernel-checked, with `#print axioms` giving only `propext`, `Classical.choice`, `Quot.sound`.
  The theorem is `Collatz.no_111_cycle_2_71`, conditional on two published inputs that are parameters of the theorem and are NOT verified there:
  - `hX`: every n <= 2^71 reaches 1 (D. Barina, J. Supercomput. 81 (2025), art. 810);
  - `RhinBound`: Rhin's irrationality measure for log 3 / log 2, in the form of Simons-de Weger (Acta Arith. 117 (2005), Lemma 12).
- **m <= 127: computed and refereed, NOT formally verified.** The same method with a four-rise bound, in exact rational arithmetic, re-computed by independent AI referee sessions. The margin at m = 127 is 0.37%, and the argument is sensitive to the exactness of the rise-bound tables.

## Prior work
m <= 91 is Hercher's theorem (formalised in Lean by Alick); m <= 101 is Alick's Lean result. The Lean m <= 111 extension and the computed m <= 127 range are ours.

## Method
Chain argument on cycle shapes using rise bounds, with per-shape certificates re-checked by the Lean kernel (m <= 111). See METHOD.md in the full repository.

## How to verify
```
git clone https://github.com/LeeMarshall1113/collatz-m-cycles
cd collatz-m-cycles
lake exe cache get
lake build Collatz && lake env lean scripts/Statement.lean
```
Needs elan and at least 16 GB RAM or pagefile; the repository README describes a serial build for small machines.

## Status
Lean-verified (m <= 111, conditional on hX and RhinBound); computed and refereed (m <= 127).

## Caveats
- Conditional on the two published inputs above.
- An improvement of a known bound, not a resolution of anything central; Collatz itself is untouched.
- m = 112..127 has no formal verification; treat it as a strongly cross-checked computation.
