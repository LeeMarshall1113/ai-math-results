# Collatz m-cycles: no nontrivial cycle with fewer than 112 local minima

**Full repository (Lean sources, method, computed tables): https://github.com/LeeMarshall1113/collatz-m-cycles**
The Lean tree is deliberately not duplicated here.

## Claim
Every nontrivial positive Collatz cycle has at least 112 local minima (no m-cycle with m <= 111).
- **m <= 111: Lean-verified**, kernel-checked, with `#print axioms` giving only `propext`, `Classical.choice`, `Quot.sound`.
  The headline theorem is `Collatz.Gelfond.main_statement_111_barina_only`, conditional on one published input that is a
  parameter of the theorem and is NOT verified there:
  - `hX`: every n <= 2^71 reaches 1 (D. Barina, J. Supercomput. 81 (2025), art. 810).
  Rhin's irrationality measure for log 3 / log 2 is no longer assumed: it is replaced by an effective Gelfond-type bound
  proved inside Lean (`Collatz/Gelfond/`). The earlier form `Collatz.no_111_cycle_2_71`, which also takes `RhinBound`
  (Simons-de Weger, Acta Arith. 117 (2005), Lemma 12) as a hypothesis, is kept in the repository.
- **m <= 133: computed, NOT formally verified.** Every table value is computed by two independent codes and every m is
  checked by two independent chain implementations, in exact rational arithmetic. See the repository README for the
  margins and for what the computation assumes.

## Prior work
m <= 91 is Hercher's theorem (formalised in Lean by Alick); m <= 101 is Alick's Lean result. The Lean m <= 111 extension and the computed range up to m <= 133 are ours.

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
Lean-verified (m <= 111, conditional on hX only); computed, not formally verified (m <= 133).

## Caveats
- Conditional on the published input above (Barina's verification to 2^71).
- An improvement of a known bound, not a resolution of anything central; Collatz itself is untouched.
- m = 112..133 has no formal verification; treat it as a cross-checked computation.
