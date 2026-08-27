# Provenance of the Hessian conjecture and its bridges — verbatim extracts from primary sources

Compiled 2026-08-08 by literature-verification agent. All quotes below are verbatim from the
primary sources stored in this directory unless marked otherwise.

Local primary sources:
- `C:\Users\mouni\hc4\lit\meng\HessianConjecture.tex` — arXiv e-print math-ph/0308035 (gzip
  was "HessianConjecture.TEX", dated Jan 31 2005 = final arXiv version; published as
  Appl. Math. Lett. 19 (2006) 503-510).
- `C:\Users\mouni\hc4\lit\debondt-vandenessen\ams_proc.pdf` (+ `ams_proc.txt`, pypdf text dump) —
  published PAMS paper S 0002-9939(05)07570-2, Proc. Amer. Math. Soc. 133 (2005), no. 8, 2201-2205.
- `C:\Users\mouni\hc4\lit\zhao\HNS-JC.tex` — arXiv e-print math/0409534 (gzip was "HNS-JC.tex",
  Nov 2004; published as Trans. Amer. Math. Soc. 359 (2007) 249-274).
- `C:\Users\mouni\hc4\lit\watanabe-debondt\watanabe_debondt.tex` — arXiv:1703.07624, J. Watanabe,
  "On the theory of Gordan-Noether on homogeneous forms with zero Hessian" (vanishing Hessian,
  tangential only).

---

## (a) G. Meng, "Legendre Transform, Hessian Conjecture and Tree Formula" (math-ph/0308035)

Title page: author Guowu Meng, HKUST; date "January 31, 2005". Original arXiv posting Aug 2003
(paper number 0308035).

### Abstract (verbatim)

> Let $\varphi$ be a polynomial over $K$ (a field of characteristic $0$) such that the Hessian of
> $\varphi$ is a nonzero constant. Let $\bar\varphi$ be the formal Legendre Transform of $\varphi$.
> Then $\bar\varphi$ is well-defined as a formal power series over $K$. The Hessian Conjecture
> introduced here claims that $\bar\varphi$ is actually a polynomial. This conjecture is shown to be
> true when $K=\bb{R}$ and the Hessian matrix of $\varphi$ is either positive or negative definite
> somewhere. It is also shown to be equivalent to the famous Jacobian Conjecture. Finally, a tree
> formula for $\bar\varphi$ is derived; as a consequence, the tree inversion formula of Gurja and
> Abyankar is obtained.

### Setup and definition of the formal Legendre transform (verbatim, Sec. 1.1)

> Let $K$ be a field of characteristic zero, $\varphi$ a polynomial in $n$ variables with
> coefficients in $K$, i.e., $\varphi\in K[x_1,\cdots, x_n]$. The Hessian matrix $H_\varphi(x)$ is a
> symmetric matrix whose $(i,j)$-entry is $\partial_i\partial_j\varphi(x)$. By definition, the
> determinant of $H_\varphi(x)$ is called the Hessian of $\varphi$ at $x$, denoted by
> $h_\varphi(x)$.
>
> Suppose that $h_\varphi\neq 0$ at $x=0$, then $y=\nabla\varphi(x):=(\partial_1
> \varphi(x),\cdots,\partial_n\varphi(x))$ has a formal inverse $x=g(y)$ --- a formal power series
> in $y$. Let $\bar \varphi(y)$ be the (formal) Legendre transform of $\varphi$, i.e., $\bar
> \varphi(y)$ is a formal power series in $y$ defined by equation
> $$\bar\varphi(y)=[xy-\varphi(x)]|_{x=g(y)}. \tag{LT}$$
> It is clear that $x=\nabla\bar\varphi(y)$, so $\bar\varphi$ is a potential function for $g$.

### The Hessian Conjecture (verbatim, Conjecture 1.2 in the TeX numbering)

> **Conjecture (Hessian Conjecture).** Let $\varphi$ be a polynomial over $K$ whose Hessian is a
> nonzero constant, $\bar\varphi$ the formal Legendre transform of $\varphi$. Then $\bar\varphi$ is
> also a polynomial.

So the exact formulation is: constant nonzero Hessian **determinant** + formal Legendre transform
is a polynomial; over an arbitrary field of characteristic zero. (No dimension index appears in
Meng's statement; the dimension-indexed family HC_n is a later convention, cf. Meng-Yang 2026.)

### Real positive-definite case (verbatim, Theorem 1.3)

> **Theorem.** The Hessian Conjecture is true when $K=\bb{R}$ and the Hessian matrix is definite
> (either positive or negative) somewhere. Therefore, if
> $\varphi(x)={1\over 2} x^2+\mbox{higher order terms}$
> is a real polynomial with $h_{\varphi}=1$ everywhere, then $\bar\varphi$ is also a polynomial.

(The proof uses a spectral-flow argument for constant signature, a mean-value argument for
injectivity of $\nabla\varphi$, and then **Theorem 2.1 of Bass-Connell-Wright [BCW82]** —
injective polynomial maps have polynomial inverse.)

### Equivalence with JC and the doubling HC_{2n} => JC_n (verbatim, Proposition 1.4 + proof)

> **Proposition.** The Hessian Conjecture is equivalent to the Jacobian Conjecture.
>
> *Proof.* If the Jacobian Conjecture is true, then equation (LT) implies that the Hessian
> Conjecture is also true. On the other hand, assume the Hessian Conjecture is true, then the
> Jacobian Conjecture is also true, and this can be proved by the following trick: Let $f$:
> $K^n\rightarrow K^n$ be a polynomial map whose Jacobian is $1$ everywhere. Let
> $\varphi (v,x)=v\cdot f(x)$, then $\varphi$ is a polynomial function on $K^{2n}$ whose Hessian is
> $(-1)^n$ everywhere. Then $\bar\varphi$ is also a polynomial function by the assumption. Now
> $\bar\varphi(w,y)=w\cdot f^{-1}(y)$ where $f^{-1}(y)$ is the formal inverse of $f$, so
> $f^{-1}(y)$ is also a polynomial. $\square$

**Verification note:** the doubling construction IS in this paper, exactly as claimed: the
potential $\varphi(v,x)=v\cdot f(x)$ lives on $K^{2n}$ and has constant Hessian $(-1)^n$; a
polynomial $\bar\varphi$ yields the polynomial inverse of $f$. Meng states the equivalence at the
level of the un-indexed conjectures; the per-dimension bookkeeping "HC in $2n$ variables implies
JC in $n$ variables" is exactly what the displayed proof establishes (and is stated with explicit
indices in de Bondt-van den Essen Theorem 1.1 below, and as Proposition 2.6 of Meng-Yang 2026).

### Meng's reduction theorem (verbatim, Sec. 1.2 "A Reduction Theorem")

> **Theorem.** The Hessian Conjecture is true $\Leftrightarrow$ for each integer $n\ge 1$ and for
> each polynomial map $\varphi$ : $\bb{C}^{2n}\to\bb{C}$ of the form
> $$\varphi(x)={1\over 2}x^2+\hbox{a homogeneous quartic polynomial in $x$},$$
> if the hessian of $\varphi$ is constant, then $\bar\varphi$ is a polynomial.

### Priority acknowledgment inside Meng's paper (verbatim, from the Introduction)

> Significant work has appeared pertaining to polynomial maps with symmetric Jacobian matrices
> [JC4SymmetricJ, NilpotentSymmetric, SymmetricReduction]. In [SymmetricReduction] M. de Bondt and
> A. van den Essen describe a Hessian conjecture virtually identical to the one we formulate, but
> which does not involve the Legendre transform. They also show its equivalence to the Jacobian
> conjecture and prove the reduction theorem of section [1.2] in this paper. We learned of that
> work only after this paper was originally written. Although the Hessian conjecture has been
> articulated as such only recently, the first result in this area - the case $n=2$ - was proved in
> 1991 [SymmetricJacobian].

where [SymmetricJacobian] = Franki Dillen, *Polynomials with constant Hessian determinant*, J. Pure
Appl. Algebra 71 (1991), 13-18, and [SymmetricReduction] = de Bondt & van den Essen, Report 0308,
Univ. of Nijmegen, June 2003, "to appear in Proc. AMS".

---

## (b) M. de Bondt, A. van den Essen, "A reduction of the Jacobian Conjecture to the symmetric case", Proc. AMS 133 (2005) 2201-2205

Primary source: published AMS PDF (downloaded from ams.org, file `ams_proc.pdf`). Received by
editors June 30, 2003; published electronically March 4, 2005. Everything is over **C** (the
construction uses $i=\sqrt{-1}$ essentially; see $f_H$ below).

### Abstract (verbatim)

> The main result of this paper asserts that it suffices to prove the Jacobian Conjecture for all
> polynomial maps of the form $x + H$, where $H$ is homogeneous (of degree 3) and $JH$ is nilpotent
> and symmetric. Also a 6-dimensional counterexample is given to a dependence problem posed by de
> Bondt and van den Essen (2003).

### Their Hessian Conjecture — the *nilpotent* form (verbatim)

> **Hessian Conjecture HC(n).** Let $f \in C[x]$. If $h(f)$ is nilpotent, then
> $F := (x_1 + f_{x_1},\ldots,x_n + f_{x_n})$ is invertible.

(Here $h(f)$ is the Hessian matrix of $f$. Note the contrast with Meng's HC, which is about
constant nonzero Hessian determinant and the Legendre transform; see their own comparison below.)

### Main theorem with dimension bookkeeping (verbatim)

> **Theorem 1.1.** The Jacobian Conjecture is equivalent to the Hessian Conjecture. More precisely,
> if HC(2n) holds, then $x + H$ is invertible for every $H : C^n \to C^n$ with $JH$ nilpotent.

**Dimension bookkeeping: the reduction doubles the dimension — JC_n (nilpotent Keller maps in
dimension n) follows from the symmetric/Hessian statement in dimension 2n.** The key construction
(their Eq. (3)) is
> $$f_H := (-i)H_1(x_1+iy_1,\ldots,x_n+iy_n)\,y_1 + \ldots + (-i)H_n(x_1+iy_1,\ldots,x_n+iy_n)\,y_n,$$
a polynomial in the $2n$ variables $(x,y)$, together with the invertible linear map
$S := (x_1 - iy_1,\ldots,x_n - iy_n,\,y_1,\ldots,y_n)$, so that $g_H := f_H\circ S =
(-i)H_1(x)y_1+\ldots+(-i)H_n(x)y_n$ and (their Eq. (4))
> $$h(g_H) = \begin{pmatrix} * & (-i)(JH)^t \\ (-i)JH & 0\end{pmatrix}.$$

> **Lemma 1.2.** Let $H = (H_1,\ldots,H_n): C^n \to C^n$ and let $f_H \in C[x,y]$ be as defined in
> (3). Then $JH$ is nilpotent iff $h(f_H)$ is nilpotent.

Proof of Theorem 1.1 (structure): $JH$ nilpotent => $h(f_H)$ nilpotent (Lemma 1.2) => by HC(2n),
$F = (x+ f_{H,x},\, y + f_{H,y})$ invertible in dimension $2n$ => $S^{-1}\circ F\circ S =
(x_1+H_1(x),\ldots,x_n+H_n(x),*,\ldots,*)$ invertible => $x+H$ invertible in dimension $n$.

### Homogeneity/nilpotency structure (verbatim)

> **Corollary 1.3.** It suffices to prove the Jacobian Conjecture for all $n \ge 2$ and all $F$ of
> the form $F = (x_1 + f_{x_1},\ldots,x_n + f_{x_n})$, where $h(f)$ is nilpotent and $f$ is
> homogeneous of degree 4 (or equivalently for all $n \ge 2$ and all $F$ of the form $F = x + H$
> with $JH$ nilpotent and symmetric and $H$ homogeneous of degree 3).
>
> *Proof.* Follows immediately from Theorem 1.1 and Corollary 2.2 of [1 = Bass-Connell-Wright 1982].

### Dependence problems and the SDP(6) counterexample

> **(Homogeneous) Dependence Problem (H)DP(n).** Let $H := (H_1,\ldots,H_n)$ with $H(0)=0$ be
> (homogeneous of degree $d\ge1$) such that $JH$ is nilpotent. Are the $H_i$ linearly dependent
> over $C$?

> **(Homogeneous) Symmetric Dependence Problem (H)SDP(n).** Let $H$ with $H(0) = 0$ be (homogeneous
> of degree $d \ge 1$) such that $JH$ is nilpotent and symmetric. Are the $H_i$ linearly dependent
> over $C$?

> **Theorem 2.1 ([2, Theorem 2.1]).** i) If SDP(p) has an affirmative answer for all $p \le n$,
> then HC(n) holds.
> ii) If SDP(p) has an affirmative answer for all $p \le n-2$ and HSDP(p) for $p = n-1$ and
> $p = n$, then HC(n) holds for all homogeneous $f \in C[x]$.

(Here [2] = de Bondt-van den Essen, *Nilpotent symmetric Jacobian matrices and the Jacobian
Conjecture*, J. Pure Appl. Algebra 193 (2004), no. 1-3, 61-70 — so the SDP => HC(n) machinery is
in that 2004 paper.)

> **Proposition 2.2.** If $n$ is minimal such that (H)DP(n) does not hold, then (H)SDP(2n) does not
> hold either.

Applied to the classical 3-dimensional DP counterexample $H_1 = x_2 - x_1^2$, $H_2 = x_3 +
2x_1(x_2-x_1^2)$, $H_3 = -(x_2-x_1^2)^2$ (their Eq. (8); DP(2) holds), this gives: **$h(f_H)$ is a
counterexample to SDP(6)**.

> **Added in proof.** In a recent paper the authors gave an affirmative answer to HDP(3). Also, the
> first author found counterexamples to HDP(n) for all $n \ge 5$.

### Dimension-4/5 facts recorded in the Introduction (verbatim)

> In [12] and [7] the cubic homogeneous cases in dimension 3 (resp. 4) were treated by Wright
> (resp. Hubbers).
>
> Recently, in [6] Washburn and the second author treated one more special case, namely they showed
> that if $n \le 4$, then the Jacobian Conjecture holds for all polynomial maps of the form
> $F = x + H$, where $JH$ is homogeneous, nilpotent and symmetric.
>
> ... Finally we would like to mention that in [3] the authors have obtained the following
> extensions of the results from [6]: the Jacobian Conjecture holds for all $F$ of the form
> $x + H$, where $JH$ is nilpotent and symmetric in the case $n \le 4$ ($H$ need not be
> homogeneous) and in the case $n = 5$ when $H$ is homogeneous.

with [6] = A. van den Essen and S. Washburn, *The Jacobian Conjecture for symmetric matrices*
[sic; journal title: "... for symmetric Jacobian matrices"], J. Pure Appl. Algebra 189 (2004),
no. 1-3, 123-133; [3] = M. de Bondt and A. van den Essen, *Nilpotent symmetric Jacobian matrices
and the Jacobian Conjecture II*, J. Pure Appl. Algebra 196 (2005), 135-148; [7] = E. Hubbers,
*The Jacobian Conjecture: Cubic homogeneous maps in Dimension Four*, Master's thesis, Univ. of
Nijmegen, 1994; [12] = D. Wright, Linear and Multilinear Algebra 34 (1993), 85-97.

### Their comparison with Meng's conjecture (verbatim, Sec. 3 "Final remarks")

> Almost three months after this paper was submitted, the authors were notified by David Wright
> that the paper [9] by Guowu Meng had appeared on the internet, in which he obtained a result
> similar to ours. He also formulates a Hessian Conjecture and shows that the Jacobian Conjecture
> is equivalent to his Hessian Conjecture. Meng's Hessian Conjecture states that the Jacobian
> Conjecture holds for all gradient maps $\nabla f := (f_{x_1},\ldots,f_{x_n})$. The difference
> between our Hessian Conjecture and the one formulated by Meng is that he considers all polynomial
> maps of the form $\nabla f$ with $\det h(f) \in C^*$, where we only need to consider all
> polynomial maps of the form $x + \nabla f$, with $h(f)$ nilpotent. So our reduction is more
> refined in the sense that it preserves the nilpotency as formulated in the classical reduction
> theorems of [1] and [13].

---

## (c) W. Zhao, "Hessian Nilpotent Polynomials and the Jacobian Conjecture" (math/0409534; TAMS 359 (2007) 249-274)

Primary source: arXiv e-print TeX (Nov 2004 version). Everything over C.

### Abstract (verbatim)

> Let $z=(z_1, \cdots, z_n)$ and $\Delta=\sum_{i=1}^n \frac{\partial^2}{\partial z^2_i}$ the
> Laplace operator. The main goal of the paper is to show that the well-known Jacobian conjecture
> without any additional conditions is equivalent to the following what we call *vanishing
> conjecture*: for any homogeneous polynomial $P(z)$ of degree $d=4$, if $\Delta^m P^m(z)=0$ for
> all $m \geq 1$, then $\Delta^m P^{m+1}(z)=0$ when $m>>0$, or equivalently, $\Delta^m
> P^{m+1}(z)=0$ when $m> \frac 32 (3^{n-2}-1)$. It is also shown in this paper that the condition
> $\Delta^m P^m(z)=0$ ($m \geq 1$) above is equivalent to the condition that $P(z)$ is Hessian
> nilpotent, i.e. the Hessian matrix $\text{Hes } P(z)=(\frac {\partial^2 P}{\partial z_i\partial
> z_j})$ is nilpotent. The goal is achieved by using the recent breakthrough work of M. de Bondt,
> A. van den Essen [BE1] and various results obtained in this paper on Hessian nilpotent
> polynomials.

### Hessian-nilpotency criterion (verbatim, Theorem `Crit-1`)

> **Theorem.** For any $P(z)\in \bC[[z]]$ with $o(P(z))\geq 2$, the following statements are
> equivalent.
> (1) $P(z)$ is HN [Hessian nilpotent].
> (2) $\Delta^m P^m=0$ for any $m\geq 1$.
> (3) $\Delta^m P^m=0$ for any $1\leq m\leq n$.

### Heat-equation formulation (verbatim, Theorem `Heat`)

> **Theorem.** Let $P(z)\in \bC [[z]]$ be HN with $o(P(z)) \geq 2$ and $Q_t(z)$ its deformed
> inversion pair. For any non-zero $s\in \bC$, set $U_{t,s}(z)=\exp(sQ_t(z))$. Then, $U_{t,s}(z)$
> is the unique formal power series solution of the following Cauchy problem of the Heat equation.
> $$\frac{\partial U_{t,s}}{\partial t}(z) = \frac 1{2s} \Delta U_{t,s}(z), \qquad
> U_{t=0,s}(z) = \exp(s P(z)).$$

(The *deformed inversion pair* $Q_t$ of $P$ is defined by: $G(z)=z+t\nabla Q_t(z)$ is the inverse
of $F(z)=z-t\nabla P(z)$.) The companion heat/Burgers form of JC (their Conjecture 3.1): for any
homogeneous HNP of degree $d\ge2$, the solutions of these Cauchy problems must be a polynomial in
$(z,t)$ and the exponential of a polynomial in $(z,t)$, respectively.

### Powers formula linking $Q_t$ and $\Delta^m P^{m+k}$ (verbatim, Theorem `T3.4`)

> **Theorem.** Suppose $P(z)\in \bC[[z]]$ with $o(P(z))\geq 2$ is HN. Then, for any $k \geq 1$,
> $$Q_t^k(z)=k! \sum_{m=0}^\infty \frac {t^m} {2^{m}m!(m+k)!} \Delta^{m} P^{m+k}(z).$$

### The exact JC equivalence (verbatim, Conjecture `VC` and Theorem `T6.2`)

> **Conjecture (Vanishing Conjecture).** For any HN (not necessarily homogeneous) polynomial $P(z)$
> of degree $d\geq 2$, its deformed inversion pair $Q_t(z)$ is a polynomial in both $t$ and $z$.
> More precisely, $\Delta^k P^{k+1}=0$ when $k>>0$.

> **Theorem.** The following statements are equivalent.
> (1) The vanishing conjecture for homogeneous HNP of degree $d=4$.
> (2) The vanishing conjecture for homogeneous HNP of degree $d\geq 2$.
> (3) The vanishing conjecture.
> (4) The Jacobian conjecture.

The proof uses "the gradient reduction in [BE1] and the homogeneous reduction in [BCW], [Y]": JC
holds iff it holds for maps $F(z)=z-\nabla P(z)$ with $P$ homogeneous HNP of degree 4.

> **Conjecture (Homogeneous Vanishing Conjecture).** For any homogeneous HNP $P(z)$ of degree
> $d\geq 2$, we have
> (1) $\Delta^m P^{m+1}=0$ for any $m> \alpha_{[n,d]}:=\frac 1{d-2}((d-1)^{n-1}-(d-1))$.
> (2) For any $k\geq 1$, $\Delta^m P^{m+k}=0$ for any $m> k\alpha_{[n,d]}$.

(Proposition `P6.4`: (1)<=>(2), and HVC for $d\ge2$ and HVC for $d=4$ are both equivalent to JC.
For $d=4$: $\alpha_{[n,4]} = \frac32(3^{n-2}-1)$, matching the abstract.)

### Known cases recorded by Zhao (verbatim, end of Section 7) — directly relevant to (d)

> - S. Wang [Wa] proved that the Jacobian conjecture holds for any polynomial map $F(z)$ of $\deg
>   F(z)\leq 2$. Hence Conjecture VC and HVC hold for any HNP $P(z)$ of degree $d\leq 3$.
> - For any symmetric polynomial map $F(z)=z-H(z)$ with $o(H(z))\geq 2$ and $JH(z)$ nilpotent,
>   A. van den Essen and S. Washburn [EW] showed that the Jacobian conjecture holds when $n\leq 4$
>   and $H(z)$ is homogeneous. Later, M. de Bondt and A. van den Essen [BE2]-[BE4] further proved
>   that the Jacobian conjecture holds either $n\leq 4$ without $H(z)$ being homogeneous, or $n=5$
>   with $H(z)$ being homogeneous. (For an exposition discussion on these results, see [BE5].)
>   From the results above, we see that, Conjecture VC has an affirmative answer when $n\leq 4$,
>   and Conjecture HVC is true when $n\leq 5$.
> - By Theorem 4.1 in [EW] and similar arguments there, it is easy to show that the only HNP's
>   (not necessarily homogeneous) $P(z)\in {\mathbb R}[z]$ $(o(P(z))\geq 2)$ with real
>   coefficients are $P(z)=0$. Hence, Conjecture VC and HVC hold trivially in this case.
> - Recently, D. Wright [Wr1] showed that the Jacobian conjecture holds for any symmetric
>   polynomial map $F(z)=z-H(z)$ with $H(z)$ homogeneous and $JH^3(z)=0$. Hence Conjecture HVC
>   holds for any homogeneous HNP $P(z)$ with $\text{Hes}^3(P(z))=0$.

---

## (d) Known partial results touching HC_4 / symmetric Keller maps in dimension 4

### d.1 The *nilpotent-symmetric* (de Bondt-van den Essen) side — dimension 4 IS settled there

- **van den Essen & Washburn**, *The Jacobian Conjecture for symmetric Jacobian matrices*, J. Pure
  Appl. Algebra 189 (2004), no. 1-3, 123-133. Statement (as recorded verbatim in the Proc AMS
  paper's introduction, primary source): "if $n \le 4$, then the Jacobian Conjecture holds for all
  polynomial maps of the form $F = x + H$, where $JH$ is homogeneous, nilpotent and symmetric."
  Also (via Zhao): its Theorem 4.1 implies the only real-coefficient Hessian-nilpotent polynomials
  are $0$ — the nilpotent-symmetric problem is genuinely complex.
- **de Bondt & van den Essen**, *Nilpotent symmetric Jacobian matrices and the Jacobian
  Conjecture*, J. Pure Appl. Algebra 193 (2004), 61-70 (= Nijmegen Report 0307, June 2003).
  Contains Theorem 2.1 quoted verbatim in (b) above: SDP(p) affirmative for all $p\le n$ =>
  HC(n) (nilpotent form); posed the symmetric dependence problems.
- **de Bondt & van den Essen**, *Nilpotent symmetric Jacobian matrices and the Jacobian Conjecture
  II*, J. Pure Appl. Algebra 196 (2005), 135-148. Statement (from the Proc AMS introduction,
  primary, and the ScienceDirect abstract, secondary): **JC holds for all $F = x + H$ with $JH$
  nilpotent and symmetric when $n \le 4$ ($H$ need not be homogeneous), and when $n = 5$ with $H$
  homogeneous.**
- **Hubbers** (Master's thesis, Nijmegen 1994): the cubic homogeneous case of JC in dimension 4
  (classification of cubic homogeneous nilpotent maps in dim 4) — cited in Proc AMS intro.
- **Drużkowski**, *The Jacobian Conjecture: symmetric reduction and solution in the symmetric
  cubic linear case*, Ann. Polon. Math. 87 (2005), 83-92: JC holds for cubic-linear maps
  $x + (Ax)^{*3}$ with symmetric Jacobian, in **all** dimensions, over char-0 fields (secondary:
  search-verified, EUDML/ResearchGate metadata; independent of the de Bondt-van den Essen line).
- **Wright** (cited by Zhao as [Wr1], "Ideal Membership Questions Relating to the Jacobian
  Conjecture", then to appear): JC holds for symmetric $F=z-H$, $H$ homogeneous, $JH^3=0$, any $n$.
- **de Bondt**, *Symmetric Jacobians*, Cent. Eur. J. Math. 12 (2014) (arXiv:1206.2865): "it
  suffices to prove the Jacobian conjecture for polynomial maps $x + H$ over $C$ such that $JH$
  satisfies all symmetries of the square, where $H$ is homogeneous of arbitrary degree $d \ge 3$"
  (abstract, secondary via arXiv abs page).

**Consequence recorded by Zhao (verbatim above): the Vanishing Conjecture holds for $n \le 4$ and
the Homogeneous Vanishing Conjecture for $n \le 5$.** I.e., in the nilpotent-symmetric formulation
the 4-variable case was PROVED (2004/2005). Caution for the HC_4 project: this does *not* settle
Meng's HC_4 — the nilpotent-symmetric dimension-$n$ results feed the *doubled* reduction (HC(2n)
=> JC_n needs the symmetric statement in dimension $2n$, i.e., JC_2 would need the nilpotent
symmetric case in dimension 4 *for all degrees with no nilpotency-free analogue*; note the proved
$n\le4$ nilpotent-symmetric case does NOT imply JC_2 because Theorem 1.1's input HC(2n) is the
full Hessian conjecture HC(4) = "h(f) nilpotent => x+∇f invertible" for arbitrary $f$ in 4
variables — and that IS exactly what BE-II proved for $n\le4$!... **BUT only for
$F = x + \nabla f$ maps arising with $JH = h(f)$ nilpotent symmetric; combined with Theorem 1.1
this yields precisely the classical result that JC_2 holds for Keller maps of the form
$x + H$, $JH$ nilpotent** — which is consistent: in dimension 2, nilpotent Keller maps
$x+H$ with $JH$ nilpotent are known invertible anyway (see van den Essen's book, e.g. via
Jung-van der Kulk); the BCW cubic reduction changes dimension, so no contradiction with JC_2
remaining open. The chain JC_2 <= HC(4)-nilpotent applies only to *nilpotent* 2-dim Keller maps,
not to all of JC_2.)

### d.2 The *constant-Hessian-determinant* (Meng) side — dimension 4 remains open

- **Dillen**, *Polynomials with constant Hessian determinant*, J. Pure Appl. Algebra 71 (1991),
  13-18: the case $n=2$ (cited as such in Meng's introduction, primary).
- **de Bondt**, *Polynomials with constant Hessian determinants in dimension three*, J. Pure Appl.
  Algebra 219 (2015), 3743-3754 (arXiv:1203.6605): abstract states "the Jacobian conjecture holds
  for gradient maps in dimension $n \le 3$ over a field $K$ of characteristic zero" and the paper
  notes its structural lemmas "do not hold for larger $n$" (arXiv abs page, secondary but
  direct). Via Meng's Prop. (LT substitution), this yields HC_n true for $n \le 3$ — the
  de Bondt anchor used by Meng-Yang 2026.
- **Meng's real-definite theorem** (primary, quoted in (a)): HC true over $\mathbb R$ whenever the
  Hessian matrix is definite somewhere — valid in every dimension, including 4. This is the only
  unconditional HC_4-relevant positive result found that is not a degree restriction.
- **Degree restriction** via Wang's theorem (S. Wang, J. Algebra 65 (1980), 453-494: JC holds for
  $\deg F \le 2$): if $\deg\varphi \le 3$ then $\nabla\varphi$ is a Keller map of degree $\le 2$,
  hence invertible with polynomial inverse, hence $\bar\varphi$ polynomial. So HC_n (Meng form)
  holds for $\deg\varphi \le 3$, all $n$. [Inference — Wang's theorem itself is verified only via
  Zhao's citation; the two-line bridge to HC is ours, flagged as such.]
- **Meng-Yang 2026** (local, `meng-yang\jc_hc_status_note.tex`, primary): "no result known to us
  settles $\HC_4$"; their Remark 4.2: the Schur descent self-terminates at 5 variables and
  "Reaching $\HC_4$ therefore appears to require genuinely new ideas."
- **Searches performed** (Aug 8, 2026; WebSearch): "Hessian conjecture dimension four HC_4 2026",
  "symmetric Jacobian conjecture dimension four", "Hessian nilpotent dimension four", "Xiaosong
  Sun symmetric cubic", de Bondt/van den Essen nilpotent symmetric. **No paper was found claiming
  any progress on Meng's HC_4 beyond the above** (the only 2026 items surfaced were Meng-Yang
  2607.22198 itself, Gao 2608.00222, and T. Shaska, *Graded Keller maps and the Jacobian
  Conjecture*, arXiv:2607.20210 — the latter checked: no Hessian/symmetric/dimension-4 content).
  Xiaosong Sun's related output found: quadratic linear Keller maps of nilpotency index 3 (Liu-Du-
  Sun, LAA 2008) and Guo-de Bondt-Du-Sun (LAA 2012) on invertible sums of Jacobian matrices —
  neither is specifically about symmetric dimension-4 maps.

### d.3 Summary table (all statements over char-0 fields; symmetric side over C)

| Result | Formulation | Dimensions covered | Source (status) |
|---|---|---|---|
| Dillen 1991 | const. Hessian det (Meng HC) | n = 2 | JPAA 71:13-18 (cited in Meng, primary) |
| de Bondt 2015 | JC for gradient maps / HC | n <= 3 | JPAA 219:3743-3754, arXiv:1203.6605 (abstract fetched) |
| Meng 2003/2006 | HC over R, definite Hessian | all n | math-ph/0308035 Thm (primary TeX) |
| Wang 1980 + LT | HC for deg phi <= 3 | all n | inference from Zhao citation |
| vdEssen-Washburn 2004 | x+H, JH homog. nilpotent symmetric | n <= 4 | JPAA 189:123-133 (quoted in PAMS, primary) |
| de Bondt-vdEssen 2005 (II) | x+H, JH nilpotent symmetric | n <= 4 (any H); n = 5 (homog.) | JPAA 196:135-148 (quoted in PAMS, primary) |
| Wright | JH^3 = 0, homog. symmetric | all n | quoted in Zhao (secondary) |
| Druzkowski 2005 | symmetric cubic linear | all n | Ann. Polon. Math. 87:83-92 (secondary) |
| Hubbers 1994 | cubic homogeneous JC | n = 4 | Master's thesis (cited in PAMS, primary) |
| de Bondt-vdEssen 2005 (PAMS) | SDP(6) counterexample | n = 6 | PAMS 133:2201-2205 (primary PDF) |
| Meng-Yang 2026 | HC_5 false; HC_4 open | n = 5 (and >= 5) | arXiv:2607.22198 v2 (local TeX, primary) |
