# Novelty search (2026-10-29)

## Automated queries (`code/novelty.py`; full output in `code/novelty_2026-10-29.log`)

**Zenodo API** (`/api/records?q=...`):
* `Seamons torus quadrilateral`
* `"quadrilateral faces" torus`
* `toroidal polyhedron quadrilateral`
* `creators.name:Seamons`
* `Szilassi`
* `"polyhedral torus"`
* `"{4,4}" torus realization`
* `torus "ten" quadrilateral faces`
* `"equivelar" torus`

**arXiv API:**
* `abs:toroidal AND abs:quadrilateral AND abs:faces`
* `abs:Szilassi`
* `abs:"polyhedral torus"`
* `abs:equivelar AND abs:torus`
* `abs:"{4,4}" AND abs:polyhedral`
* `abs:"polyhedral realization" AND abs:torus`

**OSF** (preprints and nodes, `filter[title]`):
* `toroidal polyhedron`
* `polyhedral torus`
* `Szilassi`

**GitHub:** `gh search repos` for
* `"toroidal polyhedron"`
* `"Szilassi"`
* `"polyhedral torus quadrilateral"`
* `"torus quadrilateral faces"`

**openai/math CONTENTS.md:** grep for `torus|toroid|szilassi|quadrilateral face|polyhedral`. The only hits were
unrelated (algebraic torus).

**Result.** The only relevant record is Seamons, Zenodo 23750303 (F = 11). No F = 10 result was found.

## Web and literature
* Searches for "{4,4} (3,1) polyhedral realization torus", "Ljubic torus Simutis" and the Simutis thesis title
  found the following.
  * GS 2009, Remark 3.7: Simutis proved the Ljubić torus has no convex realisation.
  * Grünbaum 1999 (*Acoptic polyhedra*, p. 13, Fig. 13) identifies the Ljubić torus as the 9-quadrangle twisted
    triangular picture frame. Fig. 14 lists {4,4}_(3,1) as not known to be acoptically realisable.
* **The Simutis thesis was not obtained.** J. S. S. Simutis, *Geometric realizations of toroidal maps*, PhD,
  UC Davis, 1977, advisor D. W. Barnette.
  * Searches for the title and author, ProQuest/Dissertation Abstracts, and the UC Davis repository returned only
    the Mathematics Genealogy entry (id 75989) and citations: Brehm–Altshuler, Geom. Dedicata 18 (1985); a
    Mathematika article; GS 2009.
  * No full text or abstract was found. This is the main residual risk.

## Limits
* MathSciNet, the Brehm–Wills and Brehm–Schulte handbook chapters, and non-English sources were not checked.
