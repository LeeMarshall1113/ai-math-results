# Tight lonely-runner instances for six and seven speeds: full proof

This file holds the complete argument. The computational steps are stated as finite
facts. The two implementations in `code/` check them, and Section 7 says exactly which
facts they compute.

## 0. Conventions

* $k$ is the number of nonzero speeds, so there are $k+1$ runners. A speed vector is
  $v=(v_1,\dots,v_k)\in\mathbb Z_{>0}^k$. Repeated entries are allowed unless we say
  otherwise.
* $\|x\|$ is the distance from $x\in\mathbb R$ to the nearest integer.
* $\mathrm{ML}(v)=\max_{t\in\mathbb R}\min_i\|tv_i\|$. The maximum exists because the
  function is continuous and 1-periodic.
* $\mathrm{LRC}(k)$ is the statement that $\mathrm{ML}(v)\ge 1/(k+1)$ for every
  $v\in\mathbb Z_{>0}^k$.
  * For distinct speeds this is the usual Lonely Runner Conjecture with $k+1$ runners.
  * With repeated speeds it follows from the distinct case with fewer speeds, since
    $1/(j+1)\ge 1/(k+1)$ for $j\le k$.
* $v$ is *primitive* if $\gcd(v)=1$. Because $\mathrm{ML}(gv)=\mathrm{ML}(v)$ for every
  integer $g\ge1$, we may always pass to the primitive vector $v/\gcd(v)$.
* A *dilate of a permutation of* $T$ is a vector $g\,(T_{\sigma(1)},\dots,T_{\sigma(k)})$
  with $g\in\mathbb Z_{>0}$ and $\sigma$ a permutation. We simply say *dilate of* $T$.
* $\mathcal T_6=\{(1,2,3,4,5,6)\}$.
* $\mathcal T_7=\{(1,2,3,4,5,6,7),\ (1,2,3,4,5,7,12),\ (1,4,5,6,7,11,13)\}$. The second
  set is due to Wills; Goddyn and Wong list both of the last two.
* All three members of $\mathcal T_7$, and $(1,\dots,6)$, are primitive.

## 1. Results

**Theorem B1.** Let $v\in\mathbb Z_{>0}^6$ with $\mathrm{ML}(v)\le 1/7$. Then $v$ is a
dilate of $(1,2,3,4,5,6)$. In particular $\mathrm{ML}(v)=1/7$ and the entries of $v$ are
distinct.

**Theorem B2.** Let $v\in\mathbb Z_{>0}^7$ with $\mathrm{ML}(v)\le 1/8$. Then $v$ is a
dilate of $(1,\dots,7)$, of $(1,2,3,4,5,7,12)$, or of $(1,4,5,6,7,11,13)$. In particular
$\mathrm{ML}(v)=1/8$.

**Corollary 1 (LRC(6) and LRC(7)).**
* B1 together with LRC(5) gives LRC(6).
* B2 together with LRC(6) gives LRC(7).

*Proof.* Suppose $\mathrm{ML}(v)<1/7$ for some $v\in\mathbb Z_{>0}^6$. Then
$\mathrm{ML}(v)\le1/7$, so B1 makes $v$ a dilate of $(1,\dots,6)$, and hence
$\mathrm{ML}(v)=1/7$. That is a contradiction. The argument for seven speeds is the same,
using B2 and the fact that each member of $\mathcal T_7$ has $\mathrm{ML}=1/8$ (§10). ∎

**Corollary 2 (tight instances).**
* The tight 6-speed instances are exactly the dilates of $(1,\dots,6)$.
* The tight 7-speed instances are exactly the dilates of the three members of
  $\mathcal T_7$.
* This answers Problem 1 of Perarnau–Serra, the characterisation of tight instances, for
  seven and for eight runners.

**Corollary 3 (Kravitz's Conjecture 8.3 for $n=7$ and $n=8$ speeds).** See §9.

## 2. Dependencies

| Step | Used in | Published source |
|---|---|---|
| LRC(5) | B1: Lemma T (n = 6), Lemma G (k = 6) | Bohman–Holzman–Kleitman, EJC 8(2) (2001) R3; also Renault, Discrete Math. 287 (2004) 93–101 |
| LRC(6) | B2: Lemma T (n = 7), Lemma G (k = 7) | **Either** Barajas–Serra, EJC 15 (2008) R48, with H. Chen's repair of its §6.6 (Zenodo 22852990), **or** Corollary 1 of this paper (B1 + LRC(5)) |
| MSS Theorem A′ proof (Section 3 of arXiv:2411.06903, Forum Math. Sigma 13 (2025) e164) | Lemma T | Malikiosis–Santos–Schymura |
| Kravitz Thm 10.2, 10.3, Prop 10.4 (arXiv:1912.06034, Comb. Theory 1 (2021) #17) | §9 only | Kravitz |

* B2 holds whichever source of LRC(6) the reader accepts.
* If the reader uses Corollary 1, B2 depends on the B1 computation as well as its own.
* Nothing else is assumed. Every remaining step is proved below or is a finite
  computation (§7).
* An independent published proof of LRC(7) exists: Rosenfeld, arXiv:2509.14111. B2 does
  not use it. The Kravitz corollary for 8 speeds can use either Rosenfeld's proof or
  Corollary 1.

## 3. Lemma T (strict form of the MSS bound)

**Notation, following MSS.**
* Write $n$ for the number of speeds. Let $v\in\mathbb Z_{>0}^n$ be primitive.
* MSS associate to $v$ a *Lonely Runner Zonotope*
  $Z_v=\sum_{i=1}^n[0,u_i]\subset\mathbb R^{n-1}$ (MSS Def. 1.4 and (1.2)). Here
  $u_i\in\mathbb Z^{n-1}$ span $\mathbb Z^{n-1}$, every $n-1$ of them are linearly
  independent, and $\sum_i v_iu_i=0$.
* Its centre is $x_v=\frac12\sum_iu_i$.
* For $\emptyset\ne J\subseteq[n]$ put $g_J=\gcd(v_j:j\in J)$.
* MSS Cor. 2.3: $|Z_v\cap\mathbb Z^{n-1}|=\sum_{\emptyset\ne J}g_J$.
* $\lambda_1(K)$ is the first successive minimum of a convex body $K$ with respect to
  $\mathbb Z^{n-1}$.
* The *lattice length* of a segment $[a,a+sw]$, with $w$ primitive, is $|s|$.

**MSS Proposition 1.5 with Remark 1.6.** For $0<\eta<1/2$, the following are equivalent:
* there is $t\in\mathbb R$ with $\|tv_i\|\ge\eta$ for all $i$;
* $(1-2\eta)(Z_v-x_v)\cap(x_v+\mathbb Z^{n-1})\neq\emptyset$.

MSS state this for $\eta=1/(n+1)$, i.e. scale $(n-1)/(n+1)$. Remark 1.6 records that the
equivalence, which is a statement about one instance, holds for every $\eta$. The
inequality is non-strict on both sides, because the zonotope is closed.

**MSS Lemma 2.5.** Let $C\subset\mathbb R^d$ be a convex body and $\ell\in\mathbb N$. If
$|C\cap\mathbb Z^d|>\ell^d$, then $\lambda_1(C-C)\le1/\ell$.

*Proof.* Two of the lattice points agree mod $\ell\mathbb Z^d$, so their difference
divided by $\ell$ is a nonzero lattice point of $\frac1\ell(C-C)$. ∎

**Lemma T.** Let $n\ge3$ and assume LRC($n-1$). Let $m\ge1$ be an integer, and let
$v\in\mathbb Z_{>0}^n$ be primitive with
$$\sum_{\emptyset\ne J\subseteq[n]}g_J\;>\;m^{\,n-1}.$$
Then $\mathrm{ML}(v)\ge\dfrac1n-\dfrac1{2m}$.

*Proof.*

**Setting up.**
* Put $\mu=\frac{n-2}{n}$ and $\lambda=\mu+\frac1m$, and assume $\lambda<1$; otherwise
  the claimed bound is $\le0$ and there is nothing to prove.
* Put $K=\lambda(Z_v-x_v)$. This is a $0$-symmetric convex body.
* By Proposition 1.5 with $1-2\eta=\lambda$, it suffices to show
  $K\cap(x_v+\mathbb Z^{n-1})\ne\emptyset$. Then $\eta=\frac{1-\lambda}{2}=\frac1n-\frac1{2m}$.

*Step 1: a short lattice direction.*
* By MSS Cor. 2.3 and the hypothesis, $|Z_v\cap\mathbb Z^{n-1}|>m^{n-1}$.
* Lemma 2.5 with $\ell=m$, an integer (this is the only place integrality of $m$ is
  used), gives $\lambda_1(Z_v-Z_v)\le 1/m$.
* Since $K$ is $0$-symmetric, $K-K=2K$, so
  $\lambda_1(K)=2\lambda_1(K-K)=\frac{2}{\lambda}\lambda_1(Z_v-Z_v)\le\frac{2}{\lambda m}$.
* Let $w\in\mathbb Z^{n-1}$ be a primitive vector attaining $\lambda_1(K)$, and write
  $\lambda_1=\lambda_1(K)$. Then $\pm w/\lambda_1\in K$.

*Step 2: a point on a good line, using LRC($n-1$).*
This step is MSS §3.1 and the first part of the proof of Theorem A′, verbatim. The scale
$\lambda$ does not appear in it.
* Let $T:\mathbb R^{n-1}\to\mathbb R^{n-2}$ be linear with $\ker T=\mathbb Rw$ and
  $T(\mathbb Z^{n-1})=\mathbb Z^{n-2}$. Such a map exists because $w$ is primitive.
* Let $Z=T(Z_v)$ and $y_v=T(x_v)$.
* MSS Theorem 3.1 gives an LRZ $Z'\subset\mathbb R^{n-2}$, with $n-1$ generators and
  centre $z\equiv y_v\pmod{\mathbb Z^{n-2}}$, and a lattice vector $\tau$ with
  $Z'+\tau\subseteq Z$. We may take $\tau=y_v-z$. Then $Z'+\tau$ has centre exactly
  $y_v$, which is how MSS use the theorem in the first line of their proof of
  Theorem A′.
* Proposition 1.5 for $Z'$, in dimension $n-2$ with $\eta=1/n$, together with
  LRC($n-1$), gives a point of $\mu(Z'-z)\cap(z+\mathbb Z^{n-2})$.
* Translating by $\tau$ and using $Z'+\tau\subseteq Z$, we get a point
  $$b\in\mu\,(Z-y_v)\cap(y_v+\mathbb Z^{n-2}).$$
* If generators of $Z'$ coincide or are opposite, the corresponding speed vector has
  repeated entries. LRC for it follows from LRC for fewer distinct speeds, with a larger
  constant. This is the convention of §0.
* Choose $a\in\mu(Z_v-x_v)$ with $T(a)=b$.
* Since $b-y_v\in\mathbb Z^{n-2}=T(\mathbb Z^{n-1})$, the line
  $a+\mathbb Rw=T^{-1}(b)$ meets $x_v+\mathbb Z^{n-1}$ in a one-dimensional affine lattice
  of points spaced at lattice length exactly $1$.

*Step 3: the segment (MSS Lemma 3.6, with the new scale).* We show that
$K\cap(a+\mathbb Rw)$ has lattice length at least $1$. A closed segment of lattice
length $\ge1$ on that line contains a point of $x_v+\mathbb Z^{n-1}$, which proves the
lemma.

*Case $a\in\mathbb Rw$.*
* Then $K\cap(a+\mathbb Rw)=K\cap\mathbb Rw=[-w/\lambda_1,w/\lambda_1]$.
* Its lattice length is $2/\lambda_1\ge\lambda m=\mu m+1\ge1$.

*Case $a\notin\mathbb Rw$.*
* Put $a'=\frac{\lambda}{\mu}a$. Since $a\in\mu(Z_v-x_v)$, we get
  $a'\in\lambda(Z_v-x_v)=K$.
* $K$ is convex and contains $a'$ and $\pm w/\lambda_1$, so it contains the triangle
  $S=\mathrm{conv}(a',w/\lambda_1,-w/\lambda_1)$.
* For $\theta=\mu/\lambda\in(0,1)$ and either sign,
  $$\theta a'+(1-\theta)\left(\pm\frac{w}{\lambda_1}\right)=a\pm\frac{1-\theta}{\lambda_1}w\in S.$$
* So $S\cap(a+\mathbb Rw)\supseteq\left[a-\frac{1-\theta}{\lambda_1}w,\ a+\frac{1-\theta}{\lambda_1}w\right]$.
  This segment has lattice length
  $$\frac{2(1-\mu/\lambda)}{\lambda_1}\ \ge\ 2\left(1-\frac{\mu}{\lambda}\right)\frac{\lambda m}{2}\ =\ m(\lambda-\mu)\ =\ 1.$$ ∎

**Where the scale enters.** In MSS the scale $(n-1)/(n+1)$ enters in exactly two places:
* the definition of $K$, and with it the translation back to times through
  Proposition 1.5;
* the bound on $\lambda_1(K)$ in Step 1, where it meets the lattice-point threshold
  through Lemma 2.5.

Theorem 3.1, the choice of $T$, $w$, $b$ and $a$, and the convexity argument of Step 3
do not involve it. MSS remark this themselves (Remark 1.6: the arguments of §3 do not
depend on the exact value). The only constraint linking the scale to the threshold is
$m(\lambda-\mu)\ge1$, with $m$ an integer for Lemma 2.5.

* **Consistency.** $m=\binom{n+1}{2}$ gives $\lambda=\frac{n-2}{n}+\frac{2}{n(n+1)}=\frac{n-1}{n+1}$,
  which is MSS's Theorem A′ exactly.
* **Strictness.** The conclusion is $\mathrm{ML}(v)\ge\frac1n-\frac1{2m}$. This is
  strictly larger than $\frac1{n+1}$ iff $2m>n(n+1)$, i.e. $m\ge\binom{n+1}2+1$. So the
  lemma excludes tight instances as well as counterexamples.

**The two instances used.**

| $n$ | $m$ | bound $\frac1n-\frac1{2m}$ | compared with $\frac1{n+1}$ |
|---|---|---|---|
| 6 | 22 | $\frac{19}{132}$ | $> \frac17$ since $133>132$ |
| 7 | 29 | $\frac{51}{406}$ | $> \frac18$ since $408>406$ |

For $n=7$, $m=28$ gives exactly $1/8$, so $29$ is the least admissible integer.

**Corollary T.** Let $S(v)=\sum_iv_i$.
* If $v\in\mathbb Z_{>0}^6$ is primitive with $\mathrm{ML}(v)\le1/7$, then
  $S(v)\le\sum_{J}g_J\le 22^5=5\,153\,632$. This uses LRC(5).
* If $v\in\mathbb Z_{>0}^7$ is primitive with $\mathrm{ML}(v)\le1/8$, then
  $S(v)\le 29^6=594\,823\,321$. This uses LRC(6).

The singletons alone contribute $S(v)$ to $\sum_Jg_J$, and every $g_J\ge1$.

## 4. Lemma G (one speed coprime to a common factor of the others)

**Lemma G.** Let $k\ge4$ and assume LRC($k-1$). Let $v\in\mathbb Z_{>0}^k$ be primitive,
and suppose some $g\ge2$ divides $v_j$ for every $j\ne i$. Then
$\mathrm{ML}(v)\ge\min(\frac1k,\frac14)=\frac1k>\frac1{k+1}$.

*Proof.*
* By LRC($k-1$), applied to $(v_j)_{j\ne i}$ (repeated entries allowed, §0), there is
  $t_0$ with $\|t_0v_j\|\ge1/k$ for all $j\ne i$.
* For $r\in\mathbb Z$ let $t_r=t_0+r/g$. Then $t_rv_j\equiv t_0v_j\pmod1$ for $j\ne i$.
* Primitivity and $g\mid v_j$ ($j\ne i$) force $\gcd(v_i,g)=1$. So $t_rv_i$ runs, mod 1,
  over the whole coset $t_0v_i+\frac1g\mathbb Z$.
* That coset has a point in $[\frac12-\frac1{2g},\frac12+\frac1{2g}]$, at distance
  $\ge\frac12-\frac1{2g}\ge\frac14$ from $\mathbb Z$. ∎

This is also Theorem 2.1 of Fan–Sun (arXiv:2306.10417). We use it in this form:

* **$k=6$, $g=7$.** If $7$ divides at least five of the $v_i$, then $\mathrm{ML}(v)>1/7$.
  (All six is impossible for primitive $v$.)
* **$k=7$, $g=2$.** If at most one $v_i$ is odd, then $\mathrm{ML}(v)>1/8$.

## 5. Gates

Throughout this section $p$ is a prime with $p>k+1$, and $N=(k+1)p$.

### 5.1 Definitions

* **Covering.** A residue $x\in\mathbb Z/N$ *covers* the time $t\in\mathbb Z/N$ if
  $\|tx/N\|\le\frac1{k+1}$, i.e. $(k+1)\,\mathrm{dist}(tx\bmod N,\,0)\le N$. Note the
  closed inequality.
* **Closed-improper.** A tuple $x\in(\mathbb Z/N)^k$ is *closed-improper* if every
  $t\in\{1,\dots,N-1\}$ is covered by some $x_i$.
* **Lemma-G filtered.** Let $g_G=7$ for $k=6$ and $g_G=2$ for $k=7$. The tuple $x$ is
  *Lemma-G filtered* if at most one $x_i$ is nonzero mod $g_G$. Since $g_G\mid N$, this
  depends only on $x\bmod g_G$.
* **Level one.**
  * A multiset $M$ of $k$ nonzero classes mod $p$ is *level-one improper* if every
    $s\in\{1,\dots,H\}$, $H=\frac{p-1}2$, has some $c\in M$ with $\|sc/p\|\le\frac1{k+1}$.
  * Equality is impossible here, since $(k+1)\nmid p$. So the condition is
    $\mathrm{dist}(sc\bmod p,0)\le m_p:=\lfloor (p-1)/(k+1)\rfloor$.
* **Orbits.** The group generated by the following acts on such multisets:
  * multiplication by a unit $u\in(\mathbb Z/p)^\times$;
  * sign changes $c\mapsto -c$ of single elements;
  * permutations.

  Covering is invariant under all three: $s\mapsto su^{-1}$ permutes the folded times.
* **Exempt orbits.** For $T\in\mathcal T_k$ with $p\nmid\prod T$, the orbit of
  $T\bmod p$ is *exempt*. Let $\mathcal E_k(p)$ be the set of exempt orbits.
* **Gate.** $p$ is a *gate* (for $k$) if every closed-improper, not Lemma-G-filtered
  $x\in(\mathbb Z/N)^k$ with $p\nmid\prod x_i$ has $x\bmod p$ in an exempt orbit.

### 5.2 The normalisation lemma

**Lemma N.** Let $M$ be a level-one improper multiset of $k$ nonzero classes mod $p$.
Put
$$q=\left\lceil\frac{2H-km_p}{k}\right\rceil .$$
Then $q\ge1$, and some element $e$ of $M$ (counted with multiplicity) is the *unique*
element of $M$ covering at least $q$ times $s\in\{1,\dots,H\}$. Consequently, the orbit
of $M$ contains a multiset in which the class $1$ occurs, and that occurrence is the
unique coverer of at least $q$ times.

*Proof.*

**Each class covers exactly $m_p$ times.**
* For a nonzero class $c$, the map $s\mapsto\pm sc$ is a bijection of $\{1,\dots,H\}$
  onto itself, where $\pm$ means folding into $\{1,\dots,H\}$.
* So the number of times $s$ with $\mathrm{dist}(sc,0)\le m_p$ is exactly $m_p$.

**Counting.**
* Let $n_s$ be the number of elements of $M$, with multiplicity, covering $s$. Then
  $\sum_sn_s=km_p$, and improperness gives $n_s\ge1$ for all $s$.
* Let $P=\#\{s:n_s=1\}$. Then $km_p=\sum_sn_s\ge P+2(H-P)$, so $P\ge 2H-km_p$.
* A class of multiplicity $\ge2$ in $M$ covers no time with $n_s=1$. So each of the $P$
  times has a unique coverer, which is an element of multiplicity one.
* By pigeonhole, some element $e$ is the unique coverer of at least $\lceil P/k\rceil\ge q$
  times.

**$q\ge1$.** We have $km_p\le k\frac{p-1}{k+1}<p-1=2H$.

**Normalising.** Multiply $M$ by $e^{-1}$. Since $c$ covers $s$ iff $ce^{-1}$ covers $se$,
this preserves improperness and the private-time count, and it sends $e$ to $1$. ∎

**How the implementations use Lemma N.** Both implementations enumerate every orbit
representative allowed by Lemma N: multisets containing 1, with the occurrence of 1 the
unique coverer of $\ge q$ times. They do so in different ways.

* **Implementation 1** (`lr.level_one_reps`) splits into two branches, following the
  private-time argument of Allikvere (arXiv:2609.02604, §4).
  * Branch A: $k$ distinct classes forming an irredundant cover. Lemma N applies to it
    with $q_k=q$.
  * Branch B: a cover by at most $k-1$ classes, padded. This is Lemma N for $k-1$
    elements.
  * It reduces each multiset to a canonical form with a precomputed unit table.
* **Implementation 2** (`ref_level1.py`) runs a single depth-first search over covers $S$
  containing $1$.
  * It prunes on $\mathrm{private}_S(1)<q$, which is sound because private sets only
    shrink as $S$ grows.
  * It also prunes on "$>r\,m_p$ uncovered times" and on "the $r$ best gains cannot cover
    the remainder".
  * It pads leaves with arbitrary classes, filters by $\mathrm{private}(1)\ge q$, and
    canonicalises by brute force over all $p-1$ units.

Both were checked against unpruned brute force over all $k$-multisets containing $1$
(§7.3).

### 5.3 Reduction of a real vector to a residue tuple

**Lemma R.** Let $k\in\{6,7\}$. Let $v\in\mathbb Z_{>0}^k$ be primitive with
$\mathrm{ML}(v)\le\frac1{k+1}$, and let $p$ be a gate with $p\nmid\prod v_i$. Then there
are $T\in\mathcal T_k$ with $p\nmid\prod T$, a unit $u$ mod $p$, signs $\epsilon_i$ and a
permutation $\sigma$ such that $v_i\equiv\epsilon_iuT_{\sigma(i)}\pmod p$ for all $i$.

*Proof.*
* Let $x=v\bmod N$. For $t\in\{1,\dots,N-1\}$ the real time $t/N$ has some $i$ with
  $\|tv_i/N\|\le\frac1{k+1}$, because $\mathrm{ML}(v)\le\frac1{k+1}$. Since $\|tv_i/N\|$
  depends only on $v_i\bmod N$, $x$ is closed-improper.
* By Lemma G and $\mathrm{ML}(v)\le\frac1{k+1}$, $v$ has at least two entries prime to
  $g_G$. Since $g_G\mid N$, $x$ is not Lemma-G filtered.
* Also $p\nmid\prod x_i$.
* By the definition of a gate, $x\bmod p$ lies in an exempt orbit, which is the claim. ∎

**Why the computation decides gates.** This is the soundness of the sieve.

*Unit equivariance.* Let $x$ be closed-improper and not filtered, and let $u$ be a unit
mod $p$.
* By the CRT choose $U\in(\mathbb Z/N)^\times$ with $U\equiv u\pmod p$ and
  $U\equiv1\pmod{k+1}$.
* Then $Ux$ is closed-improper, because $t\mapsto tU$ is a bijection of $\mathbb Z/N$.
* $Ux$ has the same residues mod $k+1$ as $x$, so it is still not filtered.
* Signs and permutations act in the same way.

*What is enumerated.* Consequently, to decide whether $p$ is a gate it suffices to do the
following.
1. Enumerate one representative $r$ of every level-one improper orbit. Every reduction
   $x\bmod p$ of a closed-improper $x$ is level-one improper, because the times
   $s(k+1)\in\mathbb Z/N$ give the level-one times $s/p$.
2. For each non-exempt $r$, test all $(k+1)^k$ lifts $x\equiv r\pmod p$ to $\mathbb Z/N$,
   one for each digit vector mod $k+1$.
3. $p$ is a gate iff no non-exempt $r$ has a closed-improper, unfiltered lift.

*Implementation 1 restricts the times tested in the lift.*
* For $k=6$ (file `tight6.py`) it skips $t\equiv0\pmod 7$. These are level-one times,
  already covered by step 1, and at level one "$\le$" and "$<$" coincide.
* For $k=7$ (file `tight7.py`) it lifts along $p\to2p\to4p\to8p$ and at each level skips
  even $t$, which were covered at the level below.
  * At the levels $2p$ and $4p$ the closed and open criteria coincide, since
    $8\,\mathrm{dist}=Lp$ with $L\in\{2,4\}$ and $p$ odd is impossible.
  * The closed criterion is applied at $8p$.
  * The projection of a closed-improper tuple mod $8p$ to mod $4p$ and mod $2p$ is
    improper, so pruning along the chain loses nothing.

*Implementation 2 tests every time $t\in\{1,\dots,N-1\}$ directly.* It uses per-digit
bitsets and a meet-in-the-middle join.

## 6. Invariants and the bookkeeping

For $T\in\mathcal T_k$ let $E_j(T)=e_j(T_1^2,\dots,T_k^2)$, the elementary symmetric
polynomials of the squares. For $v\in\mathbb Z^k$ put $w=(v_1^2,\dots,v_k^2)$ and
$$I^T_j(v)=e_j(w)\,E_1(T)^j-e_1(w)^j\,E_j(T)\qquad(2\le j\le k).$$

**Lemma I.**
1. If $v_i\equiv\epsilon_iuT_{\sigma(i)}\pmod p$ for all $i$, then $I^T_j(v)\equiv0\pmod p$
   for all $j$.
2. If $v$ is primitive, $v\in\mathbb Z_{>0}^k$, and $I^T_j(v)=0$ for all $j=2,\dots,k$,
   then $v$ is a permutation of $T$.

*Proof of 1.* Squaring removes the signs, so $w\equiv u^2\cdot\sigma(T^2)$. Hence
$e_j(w)\equiv u^{2j}E_j$ and $e_1(w)^j\equiv u^{2j}E_1^j$, and the two terms of
$I_j^T$ agree mod $p$.

*Proof of 2.*
* Put $c=e_1(w)/E_1(T)>0$. Then $e_j(w)=c^jE_j(T)$ for every $j\ge1$.
* So $\prod_i(X-w_i)=\prod_i(X-cT_i^2)$ (Vieta), and the multisets agree:
  $w=c\,\sigma(T^2)$ for some permutation $\sigma$.
* Then $v_i=\sqrt c\,T_{\sigma(i)}$, and $\sqrt c=v_1/T_{\sigma(1)}=a/b$ is rational, in
  lowest terms.
* From $bv_i=aT_{\sigma(i)}$ we get $b\mid T_{\sigma(i)}$ for all $i$, so $b=1$ because
  $\gcd T=1$.
* Likewise $a\mid v_i$ for all $i$, so $a=1$ because $\gcd v=1$. ∎

**Proposition B (bookkeeping).** Let $k\in\{6,7\}$ and let $G$ be a set of gates for $k$.
Put $m=\binom{k+1}2+1$, $S_0=m^{k-1}$, and
$$B_k=\Big(\frac{S_0}{k}\Big)^k\prod_{T\in\mathcal T_k}\ \max_{2\le j\le k}\Big[\binom kj\Big(\frac{S_0^2}{k}\Big)^jE_1(T)^j+S_0^{2j}E_j(T)\Big].$$
If $\prod_{p\in G}p>B_k$, then every primitive $v$ with $\mathrm{ML}(v)\le\frac1{k+1}$ is a
permutation of a member of $\mathcal T_k$.

*Proof.* Let $v$ be primitive with $\mathrm{ML}(v)\le\frac1{k+1}$, and suppose $v$ is not a
permutation of any $T\in\mathcal T_k$.

*The integer $F_T$.*
* By Lemma I(2), for each $T$ some $I^T_j(v)\ne0$. Let $F_T$ be the first nonzero one.
* By Lemma R and Lemma I(1), every $p\in G$ with $p\nmid\prod v$ divides some $F_T$.
* Every other $p\in G$ divides $\prod v$.
* The elements of $G$ are distinct primes, all dividing the nonzero integer
  $\prod v\cdot\prod_TF_T$. Hence
  $$\prod_{p\in G}p\le\prod_iv_i\cdot\prod_T|F_T|.$$

*Bounding the right-hand side with $S=S(v)\le S_0$ (Corollary T).*
* AM–GM gives $\prod v_i\le(S/k)^k$.
* $e_1(w)=\sum v_i^2\le S^2$.
* Maclaurin gives $e_j(w)\le\binom kj(e_1(w)/k)^j$.
* Both terms of $I^T_j$ are nonnegative, so $|I^T_j|\le e_j(w)E_1^j+e_1(w)^jE_j$.
* Each bound increases with $S$.

Together these give $\prod_{p\in G}p\le B_k$, a contradiction. ∎

**The numbers** (natural logarithms). The computation is exact, in integers and
fractions, in both implementations.

| k | $m$ | $S_0$ | $\log(S_0/k)^k$ | $\log\max_j$ for each $T$ | $\log B_k$ |
|---|---|---|---|---|---|
| 6 | 22 | 5 153 632 | 81.981 | 201.819 | **283.800** |
| 7 | 29 | 594 823 321 | 127.805 | 303.843; 307.827; 311.468 | **1050.943** |

The sharper bookkeeping changes $\log B_k$ only in the second decimal: it uses
$S\le S_0-(2^k-k-1)$ and $|a-b|\le\max(a,b)$, giving 283.758 and 1050.917. It is not
used.

## 7. The computation

### 7.1 Facts computed

* **(F6)** Every prime $11\le p\le997$ is a gate for $k=6$ except
  $17,19,23,29,31,37,47,53,61,79$. That is 154 gates out of 164 primes.
* **(F7)** Every prime $11\le p\le1153$ is a gate for $k=7$ except
  $11,17,19,23,29,31,37,41,43,47,53,67$. That is 175 gates out of 187 primes.

The full lists are in `code/gates/gates_k6.json` and `code/gates/gates_k7.json`. For each
prime they give the verdict, the number of level-one orbits, the number of non-exempt
orbits, and the number of surviving lifts found by each implementation.

### 7.2 Margins

| k | primes | gates | $\log\prod_Gp$ | $\log B_k$ | margin (nats) |
|---|---|---|---|---|---|
| 6 | 11..997 | 154 | 915.272 | 283.800 | **631.472** |
| 7 | 11..1129 | 173 | 1056.087 | 1050.943 | 5.144 |
| 7 | 11..1153 | 175 | 1070.185 | 1050.943 | **19.242** |

Each inequality $\prod_Gp>B_k$ is checked exactly in integer arithmetic.

**How load-bearing the individual verdicts are.**
* $k=6$: the gates up to $353$ already suffice (57 primes, 288.22 nats).
* $k=7$ with primes up to $1129$: every gate $\ge173$ is individually load-bearing, since
  $\log173=5.153>5.144$. A single wrong gate verdict anywhere in $173..1129$ would break
  the bound.
* $k=7$ with primes up to $1153$: the margin $19.24$ exceeds the sum of the logs of any
  two gates ($\le\log1151+\log1153=14.10$). So any two gate verdicts could be wrong
  without affecting B2, while the three largest gates together (21.13 nats) are
  load-bearing.
* We use primes up to 1153 for B2.

### 7.3 Two implementations

The implementations share no code. Implementation 2 was written from the definitions of
§5, and it imports nothing from implementation 1.

| | Implementation 1 (`code/original`) | Implementation 2 (`code/referee`) |
|---|---|---|
| level one | `lr.level_one_reps`: two-branch DFS (Allikvere-type private-time branches), unit-table canonical form | `ref_level1.py`: single DFS with Lemma N pruning, brute-force canonical form |
| lift | `tight6.py` loop over $7^6$ digit vectors (kernel-unit symmetry) and times; `tight7.py` chain $2\to4\to8$ | `ref_lift.py`: per-digit covered-time bitsets over **all** times, meet-in-the-middle |
| bookkeeping | `check_original.py` (exact Fractions) | `ref_bounds.py` (exact Fractions) |

**Agreement.**
* Both report the same verdict and the same number of level-one orbits at every prime:
  all 164 primes for $k=6$ and all 187 primes for $k=7$.
* Implementation 2 also checks, at every prime, that every exempt orbit occurs at level
  one and has a closed-improper lift.

**Cross-checks against brute force with no pruning.**
* *Level one.* `ref_validate.py` enumerates every $k$-multiset containing $1$ and
  compares the orbit sets with the DFS: equal for $k=6$, $p\le97$, and for $k=7$,
  $p\le61$.
* *Lift.* The meet-in-the-middle count equals a naive count over every digit vector and
  every time, for every orbit at $k=6$, $p\in\{11,13,17,19\}$ and $k=7$,
  $p\in\{11,13\}$.
* *Implementation 1, $k=6$.* `tight6.py validate`/`validate2`: every multiset, every
  digit vector, every time, at $p=11,\dots,29$. The orbit sets with a closed-improper lift
  are identical (sizes 1, 1, 19, 5, 7, 6). The symmetry reduction changes no verdict
  ($p=17,19,23$).
* *Implementation 1, $k=7$.* `validate7.py` at $p=11,13,17$: orbit sets of sizes 3, 1
  and 56 are identical, including the exempt orbits.

**Observation, not used in the proof.** For $k=7$, from $p=619$ up, level one has exactly
the three exempt orbits at every prime. For $k=6$, from $p=359$ up, the only level-one
orbit is $(1,\dots,6)$.

## 8. Proofs of B1 and B2

**B1.**
* Let $v\in\mathbb Z_{>0}^6$ with $\mathrm{ML}(v)\le1/7$, and replace $v$ by
  $v/\gcd(v)$.
* Corollary T (with LRC(5)) gives $S(v)\le22^5$.
* (F6) and §7.2 give $\prod_Gp>B_6$, which is used in Proposition B.
* Proposition B (with Lemmas G, R and I, which use LRC(5) through Lemma G) gives that
  $v$ is a permutation of $(1,\dots,6)$.
* So the original vector is a dilate of $(1,\dots,6)$, and $\mathrm{ML}=1/7$ by §10. ∎

**B2.** The same argument with $k=7$, $m=29$, (F7) for primes up to 1153, and LRC(6) in
Corollary T and Lemma G. ∎

## 9. Kravitz's Conjecture 8.3 for 7 and 8 speeds

Numbering follows arXiv:1912.06034v1. The published version (Comb. Theory 1 (2021) #17)
should be checked for renumbering.

* **Conjecture 8.3** (one very fast runner). For every $n\ge4$ there is
  $f_n:\mathbb N\to\mathbb N$ such that, if $v_1<\dots<v_n$ and $v_n>f_n(v_{n-1})$, then
  either $\mathrm{ML}(v)=s/(ns+1)$ for some $s\in\mathbb N$, or $\mathrm{ML}(v)\ge1/n$.
* **Theorem 10.3.** Fix $n\ge4$, and suppose that the LRC holds for $n-1$ speeds and that
  all tight sets of $n-1$ speeds are known. Then Conjecture 8.3 holds for $n$ iff every
  primitive tight $(n-1)$-set satisfies condition (II.b) of Theorem 10.2.
* **Proposition 10.4.** $(1,\dots,n-1)$ satisfies (II.b).

**$n=7$.**
* LRC(6) holds, by Corollary 1 or by Barajas–Serra with Chen's repair.
* By B1 the only primitive tight 6-set is $(1,\dots,6)$, which satisfies (II.b) by
  Prop. 10.4.
* Hence Conjecture 8.3 holds for $n=7$.

**$n=8$.**
* LRC(7) holds, by Corollary 1 (B2 + LRC(6)) or by Rosenfeld.
* By B2 the primitive tight 7-sets are exactly the three members of $\mathcal T_7$.
  * $(1,\dots,7)$ satisfies (II.b) by Prop. 10.4.
  * For the other two, (II.b) is checked exactly by `code/original/kravitz_IIb.py`,
    which follows Kravitz's Lemma 10.1 in rational arithmetic.
* For all three sets:
  * the equality times $t$ with $\min_i\|tv_i\|=1/8$ are $\{1/8,3/8,5/8,7/8\}$;
  * $D=8$, and the only admissible residue is $Q=0$;
  * the relevant speed is $v_\rho=1$, so $D/(n\,v_\rho)=1$ is an integer and
    $Q=-nv_\rho u=0$ holds.
* Hence Conjecture 8.3 holds for $n=8$.

**Numerical sanity check, not used.** `kravitz83.py` computes $\mathrm{ML}(T\cup\{x\})$
exactly:
* for $T=(1,\dots,6)$ with $x\in[2000,2400)$;
* for each $T\in\mathcal T_7$ with $x\in[2000,2300)$.

Every value is either $\ge1/n$ or of the form $s/(ns+1)$.

## 10. Exact ML of the listed tuples

`code/referee/ref_ml.py` decides $\mathrm{ML}(v)=d$ exactly. It checks two things by
sorting rational intervals:
* the closed sets $\{t:\|tv_i\|\le d\}$ cover $[0,1]$;
* the open sets $\{t:\|tv_i\|<d\}$ do not.

**Results.**
* $\mathrm{ML}=1/7$ for $(1,\dots,6)$ and for $(3,6,\dots,18)$.
* $\mathrm{ML}=1/8$ for $(1,\dots,7)$, $(1,2,3,4,5,7,12)$, $(1,4,5,6,7,11,13)$ and
  $(2,8,10,12,14,22,26)$.
* The witnesses are $t=1/8,1/8,1/8,1/16$ respectively.
* Controls: $(1,2,3,4,5,7)$ and $(1,\dots,6,8)$ are not tight.

## 11. Limits

* **Computer-assisted.** B1 and B2 rest on the finite facts (F6) and (F7). These are
  established by two independent exhaustive enumerations that agree, and by brute-force
  cross-checks at small primes. There is no formal proof certificate.
* **LRAT certificates.** A SAT/LRAT certificate was produced and checked for a few small
  gates ($k=6$, $p=13,41,43,59$). A whole-gate certificate for every gate is infeasible
  at this size: proof size grows about tenfold from $p=41$ to $p=43$.
* **Lemma T.** It is a modification of the MSS argument, written out in full in §3. It
  has not been refereed by a journal.
* **No human peer review.**
