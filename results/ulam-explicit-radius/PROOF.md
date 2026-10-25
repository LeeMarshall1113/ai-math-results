# An explicit neighbourhood of the ball in which Ulam's packing conjecture holds

**Status: draft for hostile review (2026-10-24).** Every step below is flagged
**[P]** (proved here), **[C]** (cited, with the source given) or
**[N]** (computed: certified by `certify.py` in interval/exact arithmetic).

## Main theorem

> **Theorem.** Let $K\subset\mathbb R^3$ be an origin-symmetric convex body with
> Banach–Mazur distance $d_{BM}(K,B^3)\le 1+\delta^*$, where $\delta^* = 6\cdot 10^{-8}$.
> Then its lattice packing density satisfies
> $$\delta_L(K)\ \ge\ \delta_L(B^3)=\frac{\pi}{\sqrt{18}},$$
> with equality if and only if $K$ is an ellipsoid.
> The certificate goes through for every $\delta\le 6.6\cdot10^{-8}$ (it fails at $6.7\cdot 10^{-8}$).

Equivalently, in the Hausdorff form: the conclusion holds whenever
$(1-\epsilon)E\subseteq K\subseteq(1+\epsilon)E$ for some ellipsoid $E$ and
$\epsilon\le 3\cdot 10^{-8}$, because then $d_{BM}\le(1+\epsilon)/(1-\epsilon)<1+6.6\cdot10^{-8}$, which the certificate covers.

**Consequence for Ulam's conjecture [P+C].** For such $K$, the packing density with
congruent copies, rotations allowed, satisfies $\delta(K)\ge\delta_L(K)\ge\pi/\sqrt{18}=\delta(B^3)$.
The last equality is Hales [Ha05]. So no **origin-symmetric** counterexample to Ulam's packing conjecture
lies within Banach–Mazur distance $1+6\cdot10^{-8}$ of the ball. Non-symmetric bodies are not covered (open, as in Kallus).

This is the first explicit radius we know of for Kallus's qualitative theorem [Ka14, Thm 5].

*Numbering note (referee, 2026-10-24):* references to Kallus's results use the numbering of an earlier reading. In arXiv:1212.2551v2 ("The 3-ball is a local pessimum for packing", Adv. Math. 264 (2014)), our "Thm 3(ii)", "Lemma 4", "Lemma 6" and "Thm 4" are his Thm 4(ii), Lemma 2, Lemma 4 and Thm 5. Lemma 6 (ours) uses sign λ ∈ L^∞; measurability holds since λ is measurable. See notes/REFEREE.md.
The radius is tiny. The section *Obstacle* explains why, and what would enlarge it.

## Novelty check (citation sweep, 2026-10-24)

* Kallus, *The 3-ball is a local pessimum for packing*, Adv. Math. 264 (2014) 355–370,
  arXiv:1212.2551. It is qualitative: "there exists $\varepsilon>0$".
  Remark 3 of that paper only *speculates* that $\varepsilon$ is governed by $\min_{l\ge4}|c_l|$, with $c_{10}\approx10^{-3}$.
* Kallus, *Pessimal packing shapes*, Geom. Topol. 19 (2015) 343–363, arXiv:1305.0289.
  It lists the 3-ball result as a local minimum and gives no constant.
* Kallus, *When is the ball a local pessimum for covering?*, DCG 54 (2015); the covering analogue, also qualitative.
* OpenAlex lists 13 citing works of 1212.2551, all checked. Most are physics or self-assembly papers. The mathematical ones are
  Kallus 2015 (two papers) and Kallus–Kusner (double lattice, DCG 2016). None gives an explicit neighbourhood.
* arXiv full-text API searches ("Ulam" + packing, pessimum, pessimal, Kallus + packing,
  asphericity + lattice packing) found nothing after 2017 on this question.
* Web searches (2023–2026, "Ulam's packing conjecture" + explicit or quantitative) found nothing.
  The ProofAtlas and Wikipedia pages state the result only qualitatively.

**Verdict: no quantitative version found.** Caveat: Google Scholar was not accessible, and Semantic Scholar
rate-limited us, so coverage of non-arXiv journals relies on OpenAlex.

## Conventions

* $B=B^3$ is the unit ball, $\sigma$ the rotation-invariant *probability* measure on $S^2$, $\|\cdot\|_p$ the norms in $L^p(\sigma)$,
  and $m(f)=\int f\,d\sigma$. On $SO(3)$, $dU$ is the Haar probability measure.
* $r_K$ is the radial function and $\lambda_K=\log r_K$. For an origin-symmetric (o-symmetric) body, $\lambda_K$ is even.
* A lattice $\Lambda$ is *admissible* for $K$ if $\Lambda\cap\operatorname{int}K=\{0\}$. $d_K$ is the minimal determinant of an admissible lattice,
  and $\varphi(K)=\operatorname{vol}K/(8d_K)$.
* $\Lambda_0=D_3/\sqrt2=\{z/\sqrt2: z\in\mathbb Z^3,\ \sum z_i\ \text{even}\}$, so $d(\Lambda_0)=1/\sqrt2$.
  Its 12 minimal vectors $S$ are the unit vectors $(\pm e_a\pm e_b)/\sqrt2$, with pair representatives
  $x_1..x_6=(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)$, each divided by $\sqrt2$.
  Every non-minimal vector of $\Lambda_0$ has norm $\ge\sqrt2$.
* $H_l$ is the space of spherical harmonics of degree $l$, and $\pi_l$ the $L^2$-orthogonal projection onto it.
  $Z_0$ is the set of even continuous $f$ with $m(f)=0$ and $\pi_2f=0$.

**[C] Basic facts.**
(F1) Minkowski: $K+2\Lambda$ is a packing iff $\Lambda$ is admissible for $K$, since $K$ is o-symmetric and convex, so $(K-K)/2=K$. Hence $\delta_L(K)=\varphi(K)$, which is invariant under $GL_3$.
(F2) Gauss (1831): the critical lattices of $B$ are the rotations of $\Lambda_0$, so $d_B=1/\sqrt2$ and $\varphi(B)=\pi/\sqrt{18}$. See [CS, Ch. 1] or [Ka14, §3].
(F3) Funk–Hecke: convolution with a zonal measure $\mu$ acts on $H_l$ as multiplication by $\int P_l\,d\mu$. Also $\pi_lf(x)=(2l+1)\int P_l(\langle x,y\rangle)f(y)\,d\sigma(y)$ [Gr, Thm 3.4.1, 3.3.3].
(F4) Addition theorem: $\pi_l$ has reproducing kernel $K_l(x,y)=(2l+1)P_l(\langle x,y\rangle)$, so $\|K_l(\cdot,x)\|_2^2=2l+1$ [Gr, Thm 3.3.3].
(F5) Bernstein's inequality: $(\sin\theta)^{1/2}|P_n(\cos\theta)|<(2/(\pi n))^{1/2}$ for $n\ge1$ [Sz, Thm 7.3.3].
(F6) $H_2=\{q_A(u)=\langle u,Au\rangle: A\ \text{symmetric, traceless}\}$, and $\sup_{S^2}|q_A|=\|A\|_{op}$.
(F7) Haar disintegration: conditionally on $Ux=y$, the point $Ux'$ is uniformly distributed on the circle $\{z:\langle z,y\rangle=\langle x,x'\rangle\}$.
Write $M_t$ for the average over that circle; $M_t=M_{-t}$ on even functions.
(F8) Brouwer's fixed point theorem. Jensen's inequality.

---

## Lemma 1 (Lipschitz bound for the radial function) [P]

*If $sB\subseteq K\subseteq RB$ with $s>0$, then $\log r_K$ is $L$-Lipschitz on $S^2$ for the geodesic distance, with $L=\sqrt{R^2-s^2}/s$.*

*Proof.* Let $g$ be the gauge of $K$. It is convex, positively homogeneous and Lipschitz, and $r_K=1/g$ on $S^2$.
Let $\xi\in\partial g(u)$ with $|u|=1$. By homogeneity, $\langle\xi,u\rangle=g(u)$ and $\langle\xi,y\rangle\le g(y)\le1$ for $y\in K$.
So $s|\xi|\le h_K(\xi)\le1$. Write $\xi_\perp=\xi-\langle\xi,u\rangle u$. Then
$|\xi_\perp|^2\le s^{-2}-g(u)^2$, and
$|\xi_\perp|/g(u)\le\sqrt{r_K(u)^2/s^2-1}\le L$.
Let $c(t)$ be a unit-speed great-circle arc, so $c'\perp c$. Then $\phi=\log g\circ c$ is locally Lipschitz, hence absolutely continuous.
Its one-sided derivatives are $\max_{\xi\in\partial g(c)}\langle\xi,\pm c'\rangle/g(c)=\max\langle\xi_\perp,\pm c'\rangle/g(c)$, and each lies in $[-L,L]$.
So $|\phi'|\le L$ a.e., and $|\log r_K(c(t_1))-\log r_K(c(t_0))|\le L|t_1-t_0|$. ∎

## Lemma 2 (linear algebra of $D_3$) [N, exact]

For $v\in\mathbb R^6$ there is a unique symmetric $Q(v)$ with $\langle x_i,Q(v)x_i\rangle=v_i$. It is linear: $Q(v)=\sum v_kQ^{(k)}$. Explicitly,
$Q_{aa}=\tfrac12(s_{ab}+s_{ac}-s_{bc})$ and $Q_{ab}=\tfrac12(v_{ab+}-v_{ab-})$, where $s_{ab}=v_{ab+}+v_{ab-}$.
`certify.py` checks the following exactly:
(a) $\operatorname{tr}Q(v)=\tfrac12\sum_kv_k$ (the eutaxy coefficients of $S$ are $1/4$; cf. [Ka14, Thm 1(iv)]);
(b) with $Q_\perp x:=Qx-\langle x,Qx\rangle x$: $g_1:=\max_{k,i}|Q^{(k)}_\perp x_i|=\sqrt{3/8}$ and $G_2:=\max_k\sum_i|Q^{(k)}_\perp x_i|=\sqrt6$;
(c) $\|Q^{(k)}\|_{op}=1$, since the characteristic polynomial is $t(t-1)(t+\frac12)$.
Hence, by linearity and the triangle inequality:
$|Q(v)_\perp x_i|\le g_1\|v\|_1$, $\ \sum_i|Q(v)_\perp x_i|\le G_2\|v\|_1$, and $\|Q(v)\|_{op}\le\|v\|_1$.

## Lemma 3 (explicit version of Kallus's Theorem 3(ii)) [P; constants N]

*Let $K$ be an o-symmetric convex body with $\lambda=\lambda_K\in Z_0$, $|\lambda|\le\epsilon$ and $\operatorname{Lip}(\lambda)\le L$. Put*
$$\bar Y=12g_1(e^\epsilon-1),\quad \bar c=e^{2\epsilon}LG_2\,e^{Le^\epsilon\bar Y},\quad
\kappa=\tfrac14\epsilon e^\epsilon+\frac{e^\epsilon\bar c}{2(1-\bar c)}+72\epsilon e^{2\epsilon}+288\epsilon^2e^{3\epsilon},$$
*and assume $\bar c\le\frac12$ and $(1-12\epsilon e^\epsilon)\sqrt2\ge e^\epsilon$. Then for every $U\in SO(3)$*
$$\frac{\varphi(B)}{\varphi(K)}-1\ \le\ A\lambda(U)+\kappa\,G\lambda(U),\qquad
A\lambda(U):=\tfrac12\sum_{i=1}^6\lambda(Ux_i),\quad G\lambda(U):=\sum_{i=1}^6|\lambda(Ux_i)|.$$

*Proof.* Replacing $K$ by $U^{-1}K$, we may take $U=I$. Write $\lambda_i=\lambda(x_i)$, $r_i=e^{\lambda_i}$, $\rho_i=r_i-1$,
$s_\lambda=\sum|\lambda_i|\le6\epsilon$ and $s_\rho=\sum|\rho_i|\le e^\epsilon s_\lambda$.
For $\tau\in\mathbb R^6_{\ge0}$ put $\Psi_i(\tau)=e^\epsilon\big(\exp(Le^\epsilon|Q(\rho+\tau)_\perp x_i|)-1\big)$.
Set $T=\bar c s_\rho/(1-\bar c)\le s_\rho$ and $D=\{\tau\ge0:\sum\tau_i\le T\}$.
For $\tau\in D$, Lemma 2 gives $|Q(\rho+\tau)_\perp x_i|\le g_1(s_\rho+T)\le 2g_1s_\rho\le\bar Y$.
Using $e^z-1\le ze^z$, we get $\sum_i\Psi_i(\tau)\le e^{2\epsilon}Le^{Le^\epsilon\bar Y}G_2(s_\rho+T)\le\bar c(s_\rho+T)=T$.
So $\Psi(D)\subseteq D$, and Brouwer gives a fixed point $\tau^*$. Let $Q'=Q(\rho+\tau^*)$ and $p_i=(I+Q')x_i$.

*Admissibility.* We have $\langle p_i,x_i\rangle=r_i+\tau_i^*$. Let $\gamma_i=\angle(p_i,x_i)$; then $\tan\gamma_i=|Q'_\perp x_i|/(r_i+\tau^*_i)$, so $\gamma_i\le e^\epsilon|Q'_\perp x_i|$.
Lemma 1, along the arc from $x_i$ to $p_i/|p_i|$, gives $r_K(p_i/|p_i|)\le r_ie^{L\gamma_i}$. Hence
$|p_i|\ge r_i+\Psi_i(\tau^*)\ge r_i+r_i(e^{L\gamma_i}-1)\ge r_K(p_i/|p_i|)$, using $r_i\le e^\epsilon$.
So $\pm p_i\notin\operatorname{int}K$.
Every other $v\in\Lambda_0$ has $|v|\ge\sqrt2$, and $\|Q'\|_{op}\le s_\rho+T\le 2e^\epsilon s_\lambda\le12\epsilon e^\epsilon$.
Hence $|(I+Q')v|\ge(1-12\epsilon e^\epsilon)\sqrt2\ge e^\epsilon\ge\max r_K$.
So $(I+Q')\Lambda_0$ is admissible, and $d_K\le\det(I+Q')\,d(\Lambda_0)$.

*Determinant.* Let $\mu_j$ be the eigenvalues of $Q'$, with $m:=\|Q'\|_{op}\le2e^\epsilon s_\lambda$.
Then $\det(I+Q')\le1+\operatorname{tr}Q'+3m^2+m^3$, and $\operatorname{tr}Q'=\frac12\sum\rho_i+\frac12\sum\tau^*_i$ by Lemma 2(a).
Since $e^x-1-x\le\frac12x^2e^{|x|}$, we have $\frac12\sum\rho_i\le\frac12\sum\lambda_i+\frac14\epsilon e^\epsilon s_\lambda$.
Also $\frac12\sum\tau^*_i\le\frac T2\le\frac{e^\epsilon\bar c}{2(1-\bar c)}s_\lambda$, $3m^2\le72\epsilon e^{2\epsilon}s_\lambda$ and $m^3\le288\epsilon^2e^{3\epsilon}s_\lambda$.
So $\det(I+Q')\le1+A\lambda(I)+\kappa G\lambda(I)$.

*Volume.* By Lemma 5, $\operatorname{vol}K\ge\operatorname{vol}B$. With (F1) and (F2),
$\varphi(B)/\varphi(K)=\frac{d_K}{d_B}\frac{\operatorname{vol}B}{\operatorname{vol}K}\le\det(I+Q')$. ∎

## Lemma 4 (explicit normalisation; explicit version of Kallus's Lemma 4) [P; constants N]

*Let $B\subseteq K\subseteq(1+\delta)B$. Put $\eta=\frac12\log(1+\delta)$, $L_0=\sqrt{(1+\delta)^2-1}$ and $P_2=10/(3\sqrt3)$.
Put $E_1(a)=\frac{R_a}2+\frac{z_a^2}{4(1-z_a)}$, with $R_a=(1-a)^{-2}-1-2a$ and $z_a=2a+R_a$, and put $\theta_a=\arcsin\frac a{1-a}$.
Suppose $a<1/100$ satisfies $P_2(\eta+E_1(a)+L_0\theta_a)\le a$.
Then some linear image $K_2$ of $K$ has $\lambda_{K_2}\in Z_0$, with*
$$\sup|\lambda_{K_2}|\le\epsilon:=\log(1+\delta)+2(a+E_1(a)),\qquad
\operatorname{Lip}\lambda_{K_2}\le L:=L_0\tfrac{1+a}{1-a}+\tfrac12\Big(\tfrac{(1+a)^2}{(1-a)^2}-1\Big).$$

*Proof.* $\lambda_0:=\lambda_K\in[0,2\eta]$, and $\operatorname{Lip}\lambda_0\le L_0$ by Lemma 1 (with $s=1$, $R=1+\delta$).
For traceless symmetric $A$ with $\|A\|_{op}\le a$, put $M=(I+A)^{-1}$, $v_A(u)=Mu/|Mu|$ and $\ell_A(u)=-\log|Mu|$.
Then $\lambda_{(I+A)K}=\lambda_0\circ v_A+\ell_A$.
Write $\ell_A=q_A+e_1$ and $\lambda_0\circ v_A=\lambda_0+e_2$.

*$e_1$:* with weights $w_j=\langle u,f_j\rangle^2$ in an eigenbasis of $A$ (eigenvalues $\alpha_j$),
$\ell_A=-\frac12\log(1-2q_A+R)$, where $R=\sum w_j[(1+\alpha_j)^{-2}-1+2\alpha_j]\in[0,R_a]$ by convexity.
So $e_1=-\frac R2+[-\frac12\log(1+z)+\frac z2]$ with $|z|\le z_a$, and $|e_1|\le E_1(a)$ because $|\log(1+z)-z|\le z^2/(2(1-|z|))$.

*$e_2$:* for symmetric $C$ with spectrum in $[b_1,b_2]$, $|Cu-\langle Cu,u\rangle u|\le(b_2-b_1)/2$.
So $\sin\angle(u,Mu)\le\frac{a/(1-a^2)}{1/(1+a)}=\frac a{1-a}$, and $|e_2|\le L_0\theta_a$.

Let $\mathcal Q^{-1}:H_2\to\mathrm{Sym}_0$ invert $A\mapsto q_A$; by (F6), $\|\mathcal Q^{-1}h\|_{op}=\|h\|_\infty$.
By (F3) and Archimedes, $\|\pi_2\|_{\infty\to\infty}=5\int_0^1|P_2|=P_2$.
The map $\Gamma(A)=-\mathcal Q^{-1}\pi_2[(\lambda_0-\eta)+e_1+e_2]$ is continuous on the ball $\{\|A\|_{op}\le a\}$,
because $(A,u)\mapsto\lambda_{(I+A)K}(u)$ is uniformly continuous. It maps the ball into itself, since $\|\Gamma(A)\|\le P_2(\eta+E_1+L_0\theta_a)\le a$.
Brouwer gives a fixed point $A^*$, which means $\pi_2\lambda_{(I+A^*)K}=0$; recall $\pi_2$ kills constants.
Let $K_2=e^{-m}(I+A^*)K$ with $m$ the mean of $\lambda_{(I+A^*)K}$. Then $\lambda_{K_2}\in Z_0$, and
$\sup|\lambda_{K_2}|\le\operatorname{osc}\le 2\eta+2(a+E_1)$.
For the Lipschitz bound: $u\mapsto v_A(u)$ has differential norm $\le\|M\|/\min|Mu|=\frac{1+a}{1-a}$, and
$|\nabla_{S^2}\ell_A|=|(M^2u)_\perp|/|Mu|^2\le\frac12[(1-a)^{-2}-(1+a)^{-2}](1+a)^2$. ∎

## Lemma 5 (volume) [P]

If $\lambda_K\in Z_0$, then $\operatorname{vol}K=\operatorname{vol}B\int e^{3\lambda_K}d\sigma\ge\operatorname{vol}B\,e^{3m(\lambda_K)}=\operatorname{vol}B$ by Jensen. ∎

## Lemma 6 (the averaging constant; replaces Kallus's Lemma 6) [P; constants N]

*For every $\lambda\in Z_0$: $\ \|\lambda\|_1\le N_{\rm eff}\,\|A\lambda\|_{L^1(SO(3))}$, with $N_{\rm eff}\le 2\cdot31+\sqrt{\mathcal S}\le180.733$, where*
$$\mathcal S=\sum_{l\ge4\ \rm even}\frac{(2l+1)(1-c_l)^6}{1.5\,c_l},\qquad c_l=1+4P_l(\tfrac12)+P_l(0).$$

*Facts about $c_l$ [N, exact, plus C for the tail].*
$c_2=0$, $c_4=7/32$, $c_{10}=65/65536$.
$c_l>0$ for all even $4\le l\le20000$ (computed in exact rationals); the minimum is at $l=10$, and the next smallest are $c_4=0.219$ and $c_{22}=0.319$.
For $l>20000$, (F5) at $\theta=\pi/3,\pi/2$ gives $|1-c_l|\le B_0/\sqrt l$, with $B_0=\sqrt{2/\pi}\,(4(2/\sqrt3)^{1/2}+1)<4.2275$, so $c_l>0.97$.
Hence $\mathcal S\le 14097.3099+0.1962$; the $l=10$ term alone is $14031.65$ (14097.0 is the whole finite sum).

*Proof.* (i) *Adjoint.* Define $A^*:L^\infty(SO(3))\to L^\infty(S^2)$ by $\int A\lambda\cdot W\,dU=\int\lambda\cdot A^*W\,d\sigma$.
By (F7), $A^*W(y)=\frac14\sum_{x\in S}\mathbb E[W(U)\mid Ux=y]$, where $A\lambda=\frac14\sum_{x\in S}\lambda(Ux)$ for even $\lambda$.

(ii) *Pull-backs.* For even $\omega$ and $W(U)=\omega(Ux_1)$: the inner products $\langle x,x_1\rangle$, $x\in S$, take the values $1,\frac12,0,-\frac12,-1$ with multiplicities $1,4,2,4,1$.
So $A^*W=\frac14(2+8M_{1/2}+2M_0)\omega=\frac12\Phi\omega$, where $\Phi=I+\nu$ and $\nu=4M_{1/2}+M_0$ (a positive zonal measure of mass 5).
By (F3), $\Phi$ has multipliers $c_l$ (this is Kallus's operator).

(iii) *Harmonic pieces.* For $G\in H_l$ ($l$ even), $A^*AG=\frac1{16}\sum_{x,x'\in S}M_{\langle x,x'\rangle}G=1.5\,c_l\,G$.
This uses $\sum_{x'}P_l(\langle x,x'\rangle)=2c_l$ for each $x$, a sum of 12 terms.
By (F4) and Cauchy–Schwarz, $|AG(U)|=\frac14|\langle G,\sum_xK_l(\cdot,Ux)\rangle|\le\frac14\|G\|_2\big(\sum_{x,x'}K_l(x,x')\big)^{1/2}=\frac14\|G\|_2\sqrt{24(2l+1)c_l}$.
In particular $c_l\ge0$ for all even $l$.

(iv) *Duality.* Let $g=\operatorname{sign}\lambda$, so $\|g\|_\infty\le1$ and $\|g\|_2\le1$.
Put $\omega=2(I-\nu+\nu^2)g$; then $\|\omega\|_\infty\le2(1+5+25)=62$.
Let $W_f(U)=\omega(Ux_1)$, so $A^*W_f=\Phi(I-\nu+\nu^2)g=(I+\nu^3)g$, because $(1+u)(1-u+u^2)=1+u^3$.
The remainder $-\nu^3g$ has components $(1-c_l)^3g_l$, with $g_l=\pi_lg$.
For $l\ge4$ put $G_l=(1-c_l)^3g_l/(1.5c_l)$ and $W=W_f+\sum_{l\ge4\ \rm even}AG_l$.
By (iii) and Cauchy–Schwarz over $l$,
$\|W\|_\infty\le62+\sum_l|1-c_l|^3\|g_l\|_2\sqrt{(2l+1)/(1.5c_l)}\le62+\sqrt{\mathcal S}$.
The series converges uniformly, since $\mathcal S<\infty$.
Then $A^*W-g\in H_0\oplus H_2$, which is orthogonal to $\lambda$, so
$\|\lambda\|_1=\langle\lambda,g\rangle=\langle\lambda,A^*W\rangle=\langle A\lambda,W\rangle\le\|A\lambda\|_1\|W\|_\infty$. ∎

*Comparison.* Kallus averages only over the fibres $\{U:Ux_1=y\}$, which amounts to taking $W=W_f$ alone.
Any such argument has constant $\ge2/c_{10}=2016.5$. The direct-on-$SO(3)$ step (iii) pays only $c_{10}^{-1/2}$ instead of $c_{10}^{-1}$ for the bad degree.
This gains a factor of about 11 in $N$ and about 125 in $\delta$.

## Proof of the theorem [P; numbers N]

Take a linear image with $B\subseteq K\subseteq(1+\delta)B$, $\delta\le\delta^*$; $\varphi$ is $GL_3$-invariant (F1).
`certify.py` verifies the following for $\delta=6\cdot10^{-8}$:
* Lemma 4 applies with $a=5.7774\cdot10^{-8}$, giving $\epsilon\le1.7555\cdot10^{-7}$ and $L\le3.46526\cdot10^{-4}$.
* Lemma 3 applies, with $\bar c=8.4881\cdot10^{-4}\le\frac12$ and $\kappa\le4.37450\cdot10^{-4}$.
* Lemma 6 gives $N_{\rm eff}\le180.733$, and $12\kappa N_{\rm eff}\le0.94874<1$.

Let $\lambda=\lambda_{K_2}\in Z_0$. Lemma 3 holds for every $U$, and $A\lambda$ has Haar mean 0.
Let $\Omega=\{A\lambda<0\}$. Averaging the bound of Lemma 3 over $\Omega$ gives
$$\frac{\varphi(B)}{\varphi(K)}-1\le\frac{1}{|\Omega|}\Big(-\tfrac12\|A\lambda\|_1+\kappa\|G\lambda\|_1\Big)
\le\frac1{|\Omega|}\Big(-\frac1{2N_{\rm eff}}+6\kappa\Big)\|\lambda\|_1,$$
using $\|G\lambda\|_{L^1(SO(3))}=6\|\lambda\|_1$ and Lemma 6.
If $\lambda\not\equiv0$, the right side is $<0$, so $\varphi(K)>\varphi(B)$.
If $\lambda\equiv0$, then $K_2=B$ and $K$ is an ellipsoid; conversely, every ellipsoid has $\varphi=\varphi(B)$. ∎

No monotonicity in δ is needed: if $B\subseteq K\subseteq(1+\delta)B$ with $\delta\le\delta^*$, then $B\subseteq K\subseteq(1+\delta^*)B$. So the single certified run at $\delta^*$ covers every smaller δ.
For the record, `python certify.py 33/500000000` also passes ($12\kappa N_{\rm eff}=0.9965$), while $67/10^9$ fails.

---

## Where Kallus's proof is non-quantitative

Nothing in [Ka14] is a genuine compactness or contradiction argument. Every gap is an unspecified "for ε small enough" constant:

| Kallus [Ka14] | What is implicit | Here |
|---|---|---|
| Lemma 4 (normalisation via Brouwer: "fix a small closed ball") | the ball radius and ε′(ε) | Lemma 4: $a\approx P_2\eta$, $\epsilon\approx2.93\delta$; the Lipschitz constant is *not* degraded ($L\approx\sqrt{2\delta}$) |
| Thm 3(ii) (radial-function construction): constants $C$, $\varepsilon'$ with $\beta\le2\sqrt\varepsilon$ | the size of the dilation error, **of order $\sqrt\varepsilon\,\sum|\rho_x|$** | Lemma 3, with a targeted correction τ instead of a uniform dilation: $\kappa\approx\frac{\sqrt6}2L$ |
| Lemma 6 ($\Phi^{-1}$ bounded on $Z$ via a non-evaluated sum) | $\|\Phi^{-1}\|_{L^1\to L^1}$, which is $\ge 1/c_{10}=1008.2$ | Lemma 6 (SO(3) duality): $N_{\rm eff}\le180.7$ |
| Thm 4 ("quadratic and higher order terms absorbed") | second-order constants | explicit in κ |
| (implicit) non-minimal lattice vectors stay admissible; det > 0 | — | checked in Lemma 3 |

## Obstacle: why δ* is only 6·10⁻⁸

The binding condition is $\frac{\sqrt6}{2}L\lesssim\frac1{12N_{\rm eff}}$ with $L\approx\sqrt{2\delta}$, so
$\delta\approx\frac{1}{2}\big(\frac{1}{6\sqrt6\,N_{\rm eff}}\big)^2$. Two effects multiply.

1. **The near-degenerate degree, $c_{10}=65/65536$.**
   Perturbations in $H_{10}$ are almost invisible to the 12 contact points of $D_3$, at first order.
   This contributes $\sqrt{14/c_{10}}=118.8$ of $N_{\rm eff}=180.7$. It is intrinsic to the first-order method; Kallus's Remark 3 already points at it.
   Sanity check: Monte Carlo of $\|A\lambda\|_1/\|\lambda\|_1$ for random $\lambda\in H_{10}$ gives about 0.038 (`sanity.py`), against the 1/180.7 = 0.0055 we can certify.
   So a sharp $L^\infty$ bound for the degree-10 part (an interval branch-and-bound of $\sup_{SO(3)}|AG|$) could gain up to about 7× in $N$, which is about 50× in δ.
2. **The $\sqrt\delta$ Lipschitz loss.**
   A convex body within $1+\delta$ of the ball can have radial slopes of order $\sqrt{2\delta}$ (for example, small facets).
   Kallus's construction moves the contact points by angle $O(|\lambda|)$, so it pays $L\cdot|\lambda|$. That is first order in $\lambda$, with a factor $\sqrt\delta$ rather than δ.
   This turns the natural scale $\delta\sim1/N$ into $\delta\sim1/N^2$.
   *Route:* place the six contact pairs **exactly** on $\partial K$ (a fixed point, as in Lemma 3), and compare $\int\lambda(y_i(U))\,w(U)\,dU$ with $\int\lambda(Ux_i)w(U)\,dU$ by a change of variables on $SO(3)$.
   The error is then $\int\lambda\,\nabla\lambda$-type, which is second order after integration by parts. This could plausibly give $\delta\sim10^{-5}$–$10^{-4}$. It is not done here; it is a research step, not a computation.

**Kill-test verdict (PORTFOLIO lane G):** δ* = 6·10⁻⁸ lies between "park" (< 1e-10) and "pass" (≥ 1e-3).
It is a rigorous, explicit, apparently first constant. Level 4 seems generous at this size; route 2 is what would make it a real level-4 result.

## Certificate

* `certify.py` uses python-flint 0.9 (Arb ball arithmetic, 256-bit, outward rounding) and exact `fmpq` and `Fraction`. It runs in about 1 s, single thread.
  It checks Lemma 2 exactly; all $c_l$ for $l\le20000$ in exact rationals; the tail via (F5); and the constants of Lemmas 3, 4 and 6 and the final inequality.
  Exit code 0 iff all checks pass. Output for the default δ: `ALL CHECKS PASSED`, `12*kappa*N_eff = [0.94873876 +/- 2.55e-9]`.
* `sanity.py` (floating point, not part of the proof) runs the Lemma 3 construction on a degree-10 body.
  Over 300 random rotations, the admissibility margin is ≥ 3e-11 and the det bound has slack ≥ 1.5e-8. It also estimates the Lemma 6 ratio.

## References

* [Ka14] Y. Kallus, The 3-ball is a local pessimum for packing, Adv. Math. 264 (2014) 355–370, arXiv:1212.2551.
* [Ka15] Y. Kallus, Pessimal packing shapes, Geom. Topol. 19 (2015) 343–363, arXiv:1305.0289.
* [Ha05] T. C. Hales, A proof of the Kepler conjecture, Ann. Math. 162 (2005) 1065–1185.
* [CS] J. H. Conway, N. J. A. Sloane, Sphere Packings, Lattices and Groups, 3rd ed., Springer 1999.
* [Gr] H. Groemer, Geometric Applications of Fourier Series and Spherical Harmonics, CUP 1996 (Funk–Hecke Thm 3.4.1; addition theorem 3.3.3).
* [Sz] G. Szegő, Orthogonal Polynomials, AMS Colloq. Publ. 23, 4th ed. 1975, Thm 7.3.3 (Bernstein's inequality).
