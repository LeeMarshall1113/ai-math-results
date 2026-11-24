# In plain words

What these results say, why anyone cared, and why you can believe them, without the notation. The precise statements, sources and limits are in each result's folder; where this page and a folder differ, the folder is right.

## The 30-second version

- **What these are.** Known open problems that mathematicians have worked on and published about, some of them for decades.
- **What was done.** Each result here either closes one of those problems or sets a new best on it.
- **Why believe it.** Every result comes with a file that records the proof step by step and a small program that checks it. Nothing rests on anyone's word.
- **How big.** On the main problems, the previous step forward was a paper in a research journal. None of these solves a famous problem.
- **Who has looked.** Two results have been reviewed by a specialist. None has been through a journal yet.

## How to read a computer-assisted proof

Some mathematical questions come down to a very large but finite amount of checking: millions of cases, each one simple. A person cannot do that by hand; a computer can. The danger is obvious: programs have bugs.

So the work here is split in two. A **search program** (complicated, untrusted) finds the argument and writes it down as a **certificate**: a file that records, step by step, why each case holds. A separate **checker** (small, simple, written independently) reads the certificate and confirms every step. You do not need to trust the search, or the AI that wrote it. You only need to trust the checker, and you can read it, run it, or write your own.

It is the same idea as a sudoku: finding the solution is hard, checking a filled-in grid is easy, and it does not matter who or what filled it in.

For several results the checker is one whose own correctness has been proved mathematically (`cake_lpr`), or the argument is written in Lean, a language in which the computer verifies each logical step.

**What "AI-generated" means here.** An AI system (Claude) did the search, wrote the programs and drafted the arguments, directed by Lee Marshall. A second AI system was then told to attack each result and find the hole. None of that is why you should believe a result. The reason is the certificate and the checker.

**What has not happened yet.** No journal has refereed these. Two results have been reviewed by a specialist in the area. Each folder says exactly what is still taken on trust.

## The headline results

### Covering a shape with smaller copies of itself (Hadwiger, 1957)
Take any solid convex shape: a ball, a brick, a pyramid, anything without dents. How many slightly smaller copies of it do you need to cover the original completely? A cube needs 8, one for each corner. The conjecture, open since the 1950s, says 8 is always enough in three dimensions.

Nobody can prove 8. The proven number has come down slowly: 20 (1998), 16 (1999), 14 (2023). **This repository proves 13.**

How: the problem is reduced to about 7.5 million small exact optimisation problems. Each solution is stored, and two separately written programs check every one using exact fractions, with no rounding. The reduction itself is verified in Lean.

### Packing regular heptagons (seven-sided tiles)
Squares and hexagons tile a floor with no gaps. Heptagons cannot, so the question is how much of the floor they can cover at best. A specific pattern found in 1990 covers about 89.27%. It was conjectured that nothing does better; the same question for pentagons took until 2016 to settle.

**This repository contains a computer-assisted proof that the 1990 pattern is the best possible.** Of the headline results it is the one with the lightest verification so far: the computation (about 125 million cases, using arithmetic that tracks its own rounding error) can be rerun, but it has not been formally verified, and no specialist has finished reading it. The folder says so plainly.

### A point pattern that cannot exist (PG(3,9))
Finite geometries are miniature worlds with finitely many points and lines; they underlie error-correcting codes. A natural question is whether one can choose points so that every line contains exactly 2 or exactly 5 of them. In the geometry called PG(3,9) this was singled out in the literature as the first unresolved case.

**The answer is no.** The proof splits the question into cases, a solver rules out each one, and every refutation is checked by a formally verified checker.

### Settling a 2005 conjecture on network matrices (polygon conjecture, n = 4)
A 2005 conjecture described exactly where certain numbers attached to a network (eigenvalues of its Laplacian matrix) can lie. It was known to hold for networks of up to 3 nodes and, since 2026, to fail for 5 or more. The case of 4 nodes was the last one open.

**It holds for 4**, which completes the picture. The proof is 836 exact inequalities, each with a certificate, checked by three programs; one of them uses nothing beyond basic Python. A specialist in the area reran the checks independently. (Another researcher, D. Z. Liu, posted a proof slightly earlier; this one is independent and uses different certificates.)

### The Collatz problem: how complicated a loop would have to be
Take a number; if it is even, halve it; if it is odd, triple it and add one; repeat. Every number anyone has tried ends up at 1. Whether that always happens has been open since 1937 and is nowhere near solved.

One way the conjecture could fail is a loop: a number that eventually returns to itself. **This result shows that any such loop would have to rise and fall at least 112 times**; the earlier records were 68 (2005) and 91 (2023). It does not solve the Collatz problem. It is verified in Lean and was reviewed by a specialist.

## The other eight, one line each

- **Karp's conjecture, 10 nodes:** a 1973 conjecture says certain questions about a network cannot be answered without looking at every possible connection. For 10 nodes, eight of the nine ways it could fail are now ruled out; one remains.
- **A doughnut made of 10 four-sided faces** cannot be built with flat, convex faces. With a 2026 construction for 11 faces by another author, this settles which face counts are possible.
- **Drawing the complete network on 14 points** with straight lines needs 4 layers to avoid crossings; 3 are not enough. This filled the last gap in a table open since 2000.
- **Lonely runner, extreme cases:** runners circle a track at different constant speeds; the conjecture says each is at some moment far from all the others. This classifies exactly which speed sets with 6 or 7 speeds are the worst case.
- **Erdős problem 580** verified for two new sizes (20 and 21).
- **Ulam's packing conjecture** proved for shapes extremely close to a ball.
- **Lonely runner, 8 to 10 runners:** already-known cases, now with formally checked certificates.
- **A density bound in hypergraph theory** that has since been beaten by others; kept because the certificate is exact and checkable.

## Where each result stands

| Problem | Posed | State before | State now |
|---|---|---|---|
| Hadwiger covering, 3 dimensions (conjectured answer: 8) | 1957 | 14 copies proved (2023); before that 16 (1999), 20 (1998) | **13 copies proved** |
| Densest packing of regular heptagons | pattern found 1990 | conjectured best; shown unbeatable by small changes only | **proved best among all packings** (computer-assisted; not yet formally verified) |
| Type (2,5) point sets in PG(3,9) | named the first open case in 2018 and again in 2025 | open | **settled: none exists** |
| Agaev-Chebotarev polygon conjecture | 2005 | true up to 3 nodes, false from 5; 4 nodes open until October 2026 | **4 nodes proved; the conjecture is fully settled** |
| Collatz loops | problem from 1937 | a loop needs more than 91 rises and falls (2023); 68 (2005) | **more than 111** |
| Polyhedral doughnuts with convex four-sided faces | Gruenbaum-Szilassi conjecture | 10 faces open | **10 faces impossible; with the 11-face example the possible counts are settled** |
| Layers needed to draw the complete network on n points | table started 2000 | 14 points open | **14 points needs 4; table complete to 16** |
| Karp's conjecture, 10 nodes | 1973 | nine possible types of counterexample | **eight excluded, one left** |
| Lonely runner, worst cases | problem posed by Perarnau-Serra; conjecture of Kravitz | open for 6 and 7 speeds | **fully classified** |

Three of these close a question completely (PG(3,9), the polygon conjecture, the doughnut face counts). Two are the best published result on a problem that is still open (Hadwiger, Collatz loops). One is a full solution awaiting formal verification (heptagons). None is a solution to a famous open problem, and this page does not claim one.

If you are a mathematician, the folders are written for you, and corrections are welcome as issues.
