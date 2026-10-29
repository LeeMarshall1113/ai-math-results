# Novelty check, 2026-10-29

All searches were run on 2026-10-29, between 00:42 and 00:55 EDT (04:42 to 04:55 UTC).
The raw query output is summarised below.

**Result: no competing claim found.**
* No record, preprint, repository or web page classifies the tight instances for six or
  seven speeds.
* None proves Kravitz's Conjecture 8.3 beyond n = 6 speeds.
* The most recent statements that the six-speed classification is open are Cordella
  (arXiv:2609.03444 v2, 8 Sep 2026) and Zhang (arXiv:2608.13599 v3, 18 Oct 2026). Zhang
  says "for at most five speeds the tight sets are known completely".

## Zenodo (REST API `https://zenodo.org/api/records`, sorted by most recent, all pages up to 100 hits)

| time (EDT) | query `q=` | hits |
|---|---|---|
| 00:42:41 | `"tight instances" "lonely runner"` | 53 |
| 00:42:49 | `tight instances lonely runner` | 21 621 (OR semantics; first 100 read) |
| 00:43:00 | `"lonely runner" spectrum` | 46 260 (OR; first 100 read) |
| 00:43:13 | `lonely runner spectrum` | 48 811 (OR; first 100 read) |
| 00:43:24 | `Kravitz conjecture` | 14 913 (OR; first 100 read) |
| 00:43:38 | `"Kravitz" "lonely runner"` | 230 (first 100 read) |
| 00:43:53 | `"lonely runner" sextuples` | 58 |
| 00:44:02 | `"lonely runner" "seven speeds"` | 53 |
| 00:44:10 | `"lonely runner" "tight sets"` | 55 |
| 00:44:19 | `creators.name:Cordella` | 70 |
| 00:44:27 | `creators.name:Bonfioli` | 17 |
| 00:44:30 | `creators.name:Zhang AND "lonely runner"` | 1 |
| 00:44:33 | `creators.name:Kravitz` | 100 |
| 00:44:43 | `"lonely runner"` | 53 |
| 00:44:50 | `"loneliness spectrum"` | 0 |

I read the titles, and where relevant the descriptions, of all 53 records containing
"lonely runner".

**Relevant records, none of which overlaps the result:**
* **Cordella 22345079** (code for arXiv:2609.03444). Near-tight values for 5 and 6 speeds.
* **Bonfioli 20534835** (n = 4 spectrum near 1/4).
* **Bonfioli 23445368** (2026-10-18). Six manuscripts, including "eight runners", which
  is LRC(7) via 27 prime checks plus Lean, and a Lean check of Barajas–Serra LRC(6).
  Its CLAIMS.tsv (GitHub ElVec1o/lonely-runner-covering) contains no classification of
  tight sets.
* **Zhang 21756995** (code for 2608.13599, one-entry modifications).
* **Besson 21339534 and 21393876** (converse of Goddyn–Wong for doubling and for
  multiplier 3).
* **Goddyn–Wong 8275490** (2006 paper).
* **H. Chen 22851977 and 22852990** (Lean LRC for 3–7 runners; the repair of
  Barajas–Serra §6.6).
* **Hu–Bedert 21399405** (gap bound).

**New since the 2026-10-28 check:**
* 23751011, Marshall, the ai-math-results repository itself;
* 23726067, Callens, a replication of the k = 10 sieve;
* 23679206, Baird, a bundle of computational reports.

None concerns tight instances.

**Author queries.**
* *Cordella* returns no record beyond 22345079 on this topic.
* *Kravitz* returns no Zenodo record by N. Kravitz.
* *Zhang* gives one lonely-runner record, 21756995.

## arXiv (API `export.arxiv.org/api/query`, `sortBy=submittedDate`, descending)

| time (EDT) | query | total |
|---|---|---|
| 00:46:19 | `all:"lonely runner"` | 45 |
| 00:46:23 | `abs:"lonely runner" AND abs:tight` | 4 |
| 00:46:28 | `abs:"loneliness spectrum" OR abs:"lonely runner spectrum"` | 3 |
| 00:46:32 | `abs:"lonely runner" AND abs:Kravitz` | 1 |
| 00:46:36 | `ti:"lonely runner"` | 36 |
| 00:46:40 | `abs:"view obstruction"` | 12 |
| 00:46:45 | `abs:"tight instances"` | 10 |

**Newest lonely-runner submission:** 2610.23790 (Hu–Jing, 16 Oct 2026, gap bound).
**Latest revision:** 2603.24784 v4 (Blanco–Criado–Santos, 24 Oct 2026, shifted LRC).
Neither concerns tight sets.

**Tight-set papers:**
* 2608.13599 v3 (Zhang): one-entry modifications of [n−1] only.
* 2609.06259 (Beck–Everett): integer relations satisfied by tight instances; no
  classification.

**Spectrum papers:**
* 2609.03444 (Cordella): 6 speeds up to finitely many exceptions, and says seven speeds
  are blocked by the missing classification.
* 2306.10417 (Fan–Sun), 2304.01462 (Giri–Kravitz), 2411.12684 (Jain–Kravitz).

**Kravitz's arXiv author feed** (`arxiv.org/a/kravitz_n_1`) shows no lonely-runner paper
after 2411.12684.

**Rosenfeld, 2509.14111 (LRC(7))** was read for overlap. It does not classify or discuss
tight sets: its Theorem 1 is the existence statement only. It uses MSS Theorem A
unchanged, with no strict form.

## GitHub (`gh search repos "lonely runner" --sort updated`, 00:48 EDT; `gh search code`)

Repository search returned 40 repositories, updated up to 2026-10-29. The relevant ones:
* **Qinrayn/lonely-runner-tight and -tight-II.** Zhang's code. The repository
  description says "complete classification", but the README shows that this means the
  one-entry modifications.
* **ElVec1o/lonely-runner-covering** (Bonfioli, as above).
* **JD-Jones-ASES/cl-lean.** Exact maxima for Cordella's sporadic sextuples, which are
  near-tight, not tight.
* **mia-research/lonely-runner-lean.** A Lean LRC for 3–7 runners.
* **vzsky/lonely-runners** (Sungkawichai–Trakulthongchai).
* **t-tanupat/nine-and-ten-lonely-runners.**

None classifies tight sets for six or seven speeds. Code search for
`"lonely runner" tight`, `"tight sextuples"`, `"lonely runner" "1,2,3,4,5,7,12"` and
`lonely runner kravitz` returned no results.

## OSF (API `api.osf.io/v2`, 00:48 EDT)

* `preprints?filter[title]=lonely runner`: 0.
* `preprints?filter[title]=tight instances`: 0.
* `nodes?filter[title]=lonely runner`: 1, "lonely runner conjecture proof" (5yr8s,
  2026-03-16), unrelated.

## Web (00:50 to 00:55 EDT)

Queries:
* `lonely runner tight instances six speeds classification 2026`
* `Kravitz conjecture 8.3 lonely runner "one fast runner" seven speeds`
* `"lonely runner" tight sets eight runners (1,4,5,6,7,11,13) classification`

The results were Zhang, Cordella, Barajas–Serra, Rosenfeld, Fan–Sun, Perarnau–Serra and
Wikipedia. Nothing claims the classification.

## Sanity check against known families (`code/referee/ref_ml.py`, exact)

| speeds | modifications checked | tight ones found |
|---|---|---|
| 6 | every one-entry modification of [6] with new entry w in 7..150 | none |
| 7 | every one-entry modification of [7] with w in 8..120 | exactly (1,2,3,4,5,7,12) |

This agrees with Zhang's classification of one-entry modifications, and with B1 and B2.
