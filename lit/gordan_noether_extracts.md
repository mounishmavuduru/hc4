# Gordan-Noether theory: verified extracts from primary sources

Compiled 2026-08-08 for the HC_4 project (input to leading-form arguments in <= 4 variables).
All quotes below are verbatim from TeX source files downloaded from arXiv into this
directory; file paths and line numbers refer to those local copies.

## 0. Sources on disk (primary, TeX from arXiv e-prints)

| Key | Paper | Local file |
|---|---|---|
| [WdB] | J. Watanabe, M. de Bondt, *On the theory of Gordan-Noether on homogeneous forms with zero Hessian (improved version)*, arXiv:1703.07624; published in Springer Proc. Math. Stat. 319 (2020), Polynomial Rings and Affine Algebraic Geometry; earlier version Proc. Sch. Sci. Tokai Univ. 49 (2014) | `lit\watanabe-debondt\watanabe_debondt.tex` |
| [GR] | A. Garbagnati, F. Repetto, *A geometrical approach to Gordan-Noether's and Franchetta's contributions to a question posed by Hesse*, Collect. Math. 60 (2009), no. 1, 27-41; arXiv:0802.0959 | `lit\garbagnati-repetto\gr.tex` |
| [CRS] | C. Ciliberto, F. Russo, A. Simis, *Homaloidal hypersurfaces and hypersurfaces with vanishing Hessian*, Adv. Math. 218 (2008), no. 6, 1759-1805; arXiv:math/0701596 (revised 2022 TeX) | `lit\crs\crs_homaloidal.tex` |

Additional sources checked online (not downloaded): Bricalli-Favale-Pirola arXiv:2201.07550
(Selecta Math. 2023); Gondim-Russo arXiv:1312.1618 (J. Algebra); Russo's book chapter
(Springer, paywalled - see Section 6).

**Variable-count dictionary (critical to avoid off-by-one errors).** [WdB] works with
forms in `n` variables (affine/algebraic convention). [GR] and [CRS] work with
hypersurfaces in `P^n` resp. `P^r`, i.e. `n+1` (resp. `r+1`) variables. So:

- "Hesse's claim true for n <= 3" in [GR]/[CRS] (projective) = forms in **<= 4 variables** = [WdB] n <= 4.
- First counterexamples in P^4 = forms in **5 variables**.
- HC_4's "four variables" corresponds to P^3, where Gordan-Noether says vanishing Hessian => cone.

---

## 1. The Gordan-Noether theorem (forms in <= 4 variables => cone)

### 1.1 Ground-field convention of [WdB]

> "Throughout the paper we denote by $K$ an algebraically closed field of characteristic
> zero." ([WdB] tex line 378, opening of Section 2 "Notation and preliminaries")

[GR] line 81-83 similarly: "Let $f=f(x_0, \ldots, x_n)\in k[x_0, \ldots, x_n]$ be a
non-zero irreducible homogeneous polynomial over an algebraically closed field $k$ of
characteristic zero."

### 1.2 Binary forms (n = 2) - Hesse's claim, sharpest form

> "\begin{theorem}
> Assume that $n=2$ and let $f \in K[x_1, x_2]$ be a form of degree $d$ with zero Hessian.
> Then $f=(a_1x_1+a_2x_2)^d$ for some $a_1, a_2 \in K$.
> \end{theorem}"
> ([WdB] tex lines 1021-1025, Section 5 "Binary and ternary forms with zero Hessian";
> proof: the ideal I(f) of relations among the partials is prime, hence principal generated
> by a linear form.)

So: **a binary form with identically vanishing Hessian is the d-th power of a linear form.**
Historical attribution ([WdB] lines 206-207): "it should have been easy to see that Hesse's
claim is true for binary forms as well as quadrics"; the claim originates with O. Hesse
(J. reine angew. Math. 42 (1851), 117-124 and 56 (1859), 263-269 - citation as in [GR]
bibliography, lines 925-928; **note** the [WdB] bibliography (line 2267) misprints
"52 (1851)" for volume 42, while [WdB]'s own introduction says Bd. 42).

Projective version ([GR] lines 538-543):

> "Let $X=V(f)\subset\mathbb{P}^2$ be a reduced hypersurface of degree $d\geq 2$. Then
> $X=V(f)$ has vanishing Hessian if and only if $X$ is a cone, i.e. if and only if $X$
> consists of $d$ distinct lines through a point."

### 1.3 Ternary forms (n = 3)

> "\begin{theorem}
> Suppose that $f=f(x) \in K[x_1, x_2, x_3]$ is a form with zero Hessian. Then
> a variable can be eliminated from $f$ by means of a linear transformation of variables.
> \end{theorem}"
> ([WdB] tex lines 1075-1078, Theorem `ternary_form_with_zero_hessian_is_trivial`)

Projective version ([GR] lines 554-562, Prop. `case n=3`, stated in P^3 = 4 variables; see 1.4).

### 1.4 Quaternary forms (n = 4) - the case needed for HC_4

> "\begin{theorem}
> Let $f=f(x_1, x_2, x_3, x_4) \in R = K[x_1, x_2, x_3, x_4]$ be a form with zero Hessian.
> Then $f$ can be transformed into a form with three variables via a linear transformation
> of the variables $x_1, x_2, x_3, x_4$.
> \end{theorem}"
> ([WdB] tex lines 1438-1441, Theorem `quaternary_form_with_zero_hessian`, Section 7
> "Quaternary and quinary forms with zero Hessian")

Projective version, with the precise cone structure ([GR] lines 554-562):

> "\begin{prop}
> Let $X=V(f)\subset\mathbb{P}^3$ be a reduced hypersurface of degree $d\geq 3$. Then
> $X=V(f)$ has vanishing Hessian if and only if $X$ is a cone. More precisely, $X=V(f)$
> has vanishing Hessian if and only if either $X$ is a cone over a curve of vertex a point
> or $X$ consists of $d$ distinct planes through a line. In the first case $Z(f)$ is a
> plane in $\mathbb{P}^{3*}$ while in the second case it is a line in $\mathbb{P}^{3*}$.
> \end{prop}"

(The d >= 3 hypothesis in [GR] is only because d <= 2 is classical: a quadric with vanishing
Hessian is a singular quadric, i.e. a cone - [GR] Remark, lines 313-319. The [WdB]
polynomial-level statement has **no degree restriction and no reducedness/irreducibility
restriction** - important for leading-form arguments, where the leading form need not be
reduced.)

Modern one-line statement (Bricalli-Favale-Pirola, *A theorem of Gordan and Noether via
Gorenstein rings*, Selecta Math. 2023, arXiv:2201.07550, from the abstract, fetched
2026-08-08): "an hypersurface $X=V(F)\subseteq \mathbb{P}^n$ with $n\leq 3$ is a cone if
and only if $F$ has vanishing hessian", and "the statement is false if $n\geq 4$".

### 1.5 Hypotheses audit

- **Characteristic 0: required.** All sources assume it. (Standard caveat, *not* from
  these sources: in characteristic p the binary form $f = x_1^p x_2$ has identically zero
  Hessian matrix but depends on both variables in an essential way, so the theorem fails
  badly in char p.)
- **Algebraic closure: assumed in all three sources; harmless for the statement over
  subfields of C.** Justification (our observation, standard linear algebra, *not* a quote):
  being a cone is equivalent to linear dependence of the first partials
  ([GR] lines 284-289: "the hypersurface $X=V(f)$ is a cone if and only if the partial
  derivative of $f$ are linearly dependent", Prop. `cone iff degenerate`, citing Ein).
  Linear dependence over $\bar{k}$ of polynomials with coefficients in $k$ is a linear
  system over $k$, hence implies linear dependence over $k$. So for $f \in k[x_1,..,x_4]$,
  $k \subseteq \mathbb{C}$, vanishing Hessian implies $f$ is degenerate already after a
  $k$-linear change of variables. (The finer *normal forms*, e.g. $f = \ell^d$ in the
  binary case, need the closure of $k$ for the roots of $\ell$-factorizations, but the
  cone statement itself descends.)
- **Degree: no restriction** in the [WdB] statements (trivial for d <= 2).
- **Reducedness/irreducibility: not needed** at the polynomial level ([WdB]); the
  [GR]/[CRS] projective statements assume reduced (their Hessian-vanishing notion depends
  only on the support, cf. [GR] Theorem `polar degree` quoting Dimca-Papadima).

---

## 2. Gordan-Noether's general structure theorems (all n)

Zero Hessian = algebraically dependent gradient ([WdB] lines 852-857):

> "\begin{proposition}
> A homogeneous polynomial $f \in R$ is a form with zero Hessian if and only if the
> partial derivatives $\frac{\pa f}{\pa x_1},\frac{\pa f}{\pa x_2},\ldots, \frac{\pa f}{\pa x_n}$
> are algebraically dependent.
> \end{proposition}"

The correct general form of Hesse's claim, valid for every n ([WdB] lines 911-915):

> "\begin{theorem}[Gordan-Noether]
> Suppose that $f(x) \in K[x_1, \ldots, x_n]$ is a form with zero Hessian.
> Then a variable can be eliminated from $f$ and its partial derivatives
> simultaneously by means of a birational transformation of the variables.
> \end{theorem}"

Mechanism (from [WdB] introduction, lines 264-279): if hess(f) = 0, pick a polynomial
null vector $(h_1,...,h_n)$ of the Hessian matrix; then $f$ and each $\partial f/\partial x_j$
satisfy the PDE $h_1 F_{x_1} + ... + h_n F_{x_n} = 0$, and Gordan-Noether's key discovery is
that **each coefficient $h_j$ itself satisfies the same PDE** ("self-vanishing system",
Yamada's terminology). For n >= 3 one has moreover
$\mathrm{rk}(\partial h_i/\partial x_j) = \mathrm{trdeg}_K K(\boldsymbol{h}) \le n-2$
([WdB] Lemma `rank_of_hessian_matrix_corrected`, lines 1034-1044), which is what makes
n <= 4 rigid and n = 5 classifiable.

---

## 3. Five variables (P^4): counterexamples and complete classification

### 3.1 The [WdB] classification of quinary forms with zero Hessian

> "\begin{theorem}
> Let $R=K[x_1, x_2, x_3, x_4, x_5]$, and let $\Delta$ be a homogeneous polynomial of the form
> $$\Delta = p_3(x_1,x_2) x_3 + p_4(x_1,x_2) x_4 + p_5(x_1,x_2) x_5$$
> Then any element in the algebra $K[x_1, x_2][\Delta]$ is a polynomial with zero Hessian.
>
> Conversely, let $f$ be a homogeneous form in five variables with zero Hessian and assume
> that $f$ properly involves five variables. Then we can choose $\Delta$ such that $f$ can
> be transformed into a homogeneous polynomial in the algebra $K[x_1, x_2][\Delta]$ by
> means of a linear change of variables.
> \end{theorem}"
> ([WdB] tex lines 1498-1511, Theorem `thm45`)

Here $p_3,p_4,p_5$ are (homogeneous, same degree) binary forms in $x_1,x_2$. Supporting
structure result ([WdB] lines 1465-1472, Prop. `svs_for_forms_with_zero_hessian_in_5_variables`):
if there is no *linear* relation among the partials of $f$, then after a linear change
the reduced self-vanishing system satisfies $h_1=h_2=0$ with $h_3,h_4,h_5$ polynomials in
$x_1,x_2$ only.

### 3.2 The cubic case: the Perazzo cubic, essentially unique

> "It is easy to see that any degree three homogeneous polynomial $F \in K[x_1,x_2,\Delta]$
> in Theorem [thm45] which properly involves five variables can be transformed into the
> canonical form $x_1^2\,x_3+ x_1x_2\,x_4+ x_2^2\,x_5$ by means of a linear transformation
> of the variables. This form is also known as the Macaulay dual of the trivial extension
> of the algebra $K[x_1, x_2]/ (x_1,x_2)^3$ by the canonical module."
> ([WdB] tex lines 1664-1670, Remark `rem_by_H-Nasu`; cf. also the claim in the
> introduction, lines 345-347: "In this paper it is proved that a cubic form in five
> variables is essentially unique".)

Same cubic in [GR] (lines 310-311): "An easy example for $n=4$ is the following cubic
polynomial $f(x_0,x_1,x_2,x_3,x_4)=x_0x_3^2+2x_1x_3x_4+x_2x_4^2$." (Equivalent to the
[WdB] canonical form up to relabeling and scaling.) Gondim-Russo arXiv:1312.1618 use
$V(x_0x_3^2 + x_1x_3x_4 + x_2x_4^2) \subset \mathbb{P}^4$ and call it the Perazzo cubic;
per their abstract, "special Perazzo cubic hypersurfaces" *exhaust* cubic hypersurfaces
with vanishing hessian that are not cones for $\mathbb{P}^N$, $N \le 6$ (wording via
WebFetch of the abs page, 2026-08-08). Geometry: $(\Delta=0)\subset\mathbb{P}^4$ for the
canonical cubic is the projection of the Segre 3-fold
$\mathbb{P}^1\times\mathbb{P}^2 \subset \mathbb{P}^5$ from a general external point
([WdB] Remark, lines 1651-1662).

Original source of the cubic classification: U. Perazzo, *Sulle varieta cubiche la cui
hessiana svanisce identicamente*, Giornale di matematiche (Battaglini) **38** (1900),
337-354 (citation from [GR] bibliography lines 933-935; [GR] intro lines 107-111:
Perazzo "considered the case of cubic hypersurfaces with vanishing hessian and obtained
the classification of these cubics in $\PP^4$, $\PP^5$ and $\PP^6$"). We did not fetch
Perazzo's 1900 paper itself.

### 3.3 Geometric classification of vanishing-Hessian hypersurfaces in P^4

[GR] main theorem (tex lines 824-840, Theorem `HGNF`; attributed to Gordan-Noether and,
independently, Franchetta):

> "\begin{theorem}
> Let $X=V(f)\subset\mathbb{P}^4$ be an irreducible and reduced hypersurface of degree
> $d\geq 3$, not a cone. The following conditions are equivalent:
> i) $X=V(f)$ has vanishing Hessian.
> ii) $X=V(f)$ is a Franchetta hypersurface.
> iii) $X^*=V(f)^*$ is a scroll surface of degree $d$, having a line directrix $L$ of
> multiplicity $e$, sitting in a $3$-dimensional rational cone $W(f)$ with vertex $L$,
> and the general plane ruling of the cone cuts $V(f)^*$ off $L$ along $\mu \leq e$ lines
> of the scroll, all passing through the same point of $L$.
> iv) $X=V(f)$ is a general GN-hypersurface of type $(4,2,1,s)$, with $\mu =[\frac{d}{s}]$,
> which has a plane of multiplicity $d-\mu$.
> In particular, $X^*=V(f)^*$ is smooth if and only if $d=3$, $X^*=V(f)^*$ is a rational
> normal scroll of degree $3$ and $X=V(f)$ contains a plane, the orthogonal of the line
> directrix of $X^*=V(f)^*$, with multiplicity $2$.
> \end{theorem}"

[CRS] version (tex lines 1824-1841, Theorem `frank`, attributed to Franchetta [Fr2], with
the remark that "according to [Lossen] this result is contained in [GN]"): same three-way
equivalence for a *reduced* hypersurface of degree d in P^4 (their wording does not repeat
the "not a cone" hypothesis inside the theorem; the scroll description in (iii) applies to
the non-cone case). A Franchetta hypersurface ([GR] Def. lines 706-717) is a reduced
hypersurface of degree d in P^4 swept by a 1-dimensional family of planes all tangent to a
plane rational curve C on it, with the tangency/incidence condition of the definition.
Permutti proved Franchetta hypersurfaces = GN-hypersurfaces of type (4,2,1,s)
([GR] Remark lines 719-728).

Summary sentence ([GR] lines 737-741): "the hypersurfaces in $\mathbb{P}^4$ with vanishing
Hessian are either cones or Franchetta hypersurfaces and ... there are no other
possibilities. A similar result is not known in higher dimension."

### 3.4 The Gordan-Noether counterexample machine (all n >= 4 projective)

[GR] lines 629-673 / [CRS] Section (label `GNP`): for $r\ge 4$, integers $t\ge m+1$,
$2\le t\le r-2$, $1\le m\le r-t-1$, binary data $h_i(y_0,...,y_m)$ and
$\psi_j(x_{t+1},...,x_n)$, one forms the determinants $Q_\ell$ ([GR] display at lines
636-649) and $f := \sum_{k=0}^{\mu} P_k(Q_1,\ldots,Q_{t-m},x_{t+1},\ldots,x_n)$
([GR] eq. (labelled equGN), line 658-660): a *GN-polynomial of type (n,t,m,s)*. Then
([GR] Prop. `GNpol`, lines 671-673): "Every GN-polynomial has vanishing Hessian."
(Proof in [CRS] Prop. 2.9.) General GN-hypersurfaces with $\mu > n-t-2$ are not cones
([GR] Prop. 2.11(iv) = [CRS] Prop. 2.11, lines 692-703). This is the historical refutation
of Hesse's claim in >= 5 variables (Gordan-Noether 1876).

---

## 4. Original literature (citations as printed in the sources on disk)

- O. Hesse, *Ueber die Bedingung, unter welcher eine homogene ganze Function von n
  unabhaengigen Variabeln durch lineaere Substitutionen ... auf eine homogene Function
  sich zurueckfuehren laesst, die eine Variable weniger enthaelt*, J. reine angew. Math.
  **42** (1851), 117-124. ([GR] bib line 925; [WdB] bib line 2260-2267 misprints vol. as 52.)
- O. Hesse, *Zur Theorie der ganzen homogenen Functionen*, J. reine angew. Math. **56**
  (1859), 263-269.
- M. Pasch, *Zur Theorie der Hesseschen Determinante*, J. reine angew. Math. **80** (1875),
  169-176. (Proved Hesse's claim for ternary and quaternary **cubics**; [WdB] lines 208-209, bib line 2301-2303.)
- P. Gordan, M. Noether, *Ueber die algebraischen Formen, deren Hesse'sche Determinante
  identisch verschwindet*, Math. Ann. **10** (1876), 547-568. ([WdB] bib lines 2275-2277.)
- U. Perazzo, *Sulle varieta cubiche la cui hessiana svanisce identicamente*, Giornale di
  matematiche (Battaglini) **38** (1900), 337-354.
- A. Franchetta, classification of vanishing-Hessian hypersurfaces in P^4 - [CRS] cite it
  as [Fr2]; we did not fetch the original (1954, Ricerche Mat., per the literature).
- R. Permutti, Ricerche di Mat. **6** (1957), 3-10 and **13** (1964), 97-105.
- C. Lossen, *When does the Hessian determinant vanish identically? (On Gordan and
  Noether's Proof of Hesse's Claim)*, Bull. Braz. Math. Soc. **35** (2004), 71-82.
- F. Russo, *On the Geometry of Some Special Projective Varieties*, Lecture Notes of the
  Unione Matematica Italiana **18**, Springer, 2016. Chapter 7 "Hypersurfaces with
  Vanishing Hessian", pp. 177-220. (Existence and metadata verified on the Springer
  landing page 2026-08-08; **full text paywalled**, so no theorem numbers/verbatim
  statements from the book are quoted here.)
- J. Watanabe, M. de Bondt, arXiv:1703.07624 (see Section 0) - the fully-proved modern
  treatment; note their warning that the 2014 Tokai Proceedings version contains a serious
  error (its Lemma 5.2), corrected in this improved version ([WdB] lines 172-179, 350-361).
  Use the arXiv v-latest / Springer 2020 version, not the 2014 one.
- D. Bricalli, F. F. Favale, G. P. Pirola, *A theorem of Gordan and Noether via Gorenstein
  rings*, Selecta Math. (2023), arXiv:2201.07550 - independent modern proof of GN via
  Lefschetz properties of standard Artinian Gorenstein algebras.

## 5. Not verified / caveats

- Russo book Chapter 7 theorem numbers: NOT obtained (paywall). The content is covered by
  [CRS]+[GR] (Russo is a coauthor of [CRS] and the [GR] project was suggested by Russo).
- Perazzo 1900, Franchetta original, Gordan-Noether 1876 originals: not fetched; their
  statements are used here only as quoted/attributed by [WdB], [GR], [CRS].
- Gondim-Russo wording was obtained through a summarizing fetch of the arXiv abstract page;
  treat the N <= 6 exhaustion claim as reliably sourced but re-quote from the paper before
  citing verbatim.
- The char-p failure example and the descent-to-non-closed-fields remark in 1.5 are
  standard facts supplied by us, not quotes from these sources.

## 6. Use in HC_4 leading-form arguments (what is now safe to cite)

1. **(GN, n<=4)** For $f$ homogeneous in $\le 4$ variables over a field of characteristic 0
   (algebraically closed WLOG, and the cone conclusion descends): $\det \mathrm{Hess}(f)
   \equiv 0$ iff the partials of $f$ are linearly dependent iff, after a linear change of
   variables, $f$ involves at most 3 variables. No degree, irreducibility, or reducedness
   hypotheses. Source: [WdB] Thm (lines 1438-1441) + [GR] Prop (lines 554-562) +
   Bricalli-Favale-Pirola.
2. **(n=5 structure)** Every quinary form with zero Hessian properly involving 5 variables
   is, after a linear change, an element of $K[x_1,x_2][\Delta]$ with
   $\Delta = p_3(x_1,x_2)x_3+p_4(x_1,x_2)x_4+p_5(x_1,x_2)x_5$; conversely all such forms
   have zero Hessian. Source: [WdB] Thm `thm45` (lines 1498-1511).
3. **(Perazzo cubic)** The unique (up to linear equivalence) quinary cubic with zero
   Hessian properly involving 5 variables is $x_1^2x_3+x_1x_2x_4+x_2^2x_5$. Source: [WdB]
   Remark (lines 1664-1670); classical source Perazzo 1900.
4. Any leading-form argument in $\le 4$ variables may conclude from
   $\det\mathrm{Hess} \equiv 0$ that the leading form is degenerate (a "cone"), and the
   degenerating linear change is available over the ground field (char 0).
