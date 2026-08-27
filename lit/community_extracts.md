# Community extracts: status of JC_2 / HC_4 and the counterexample mechanism (as of Aug 8, 2026)

Compiled by literature-verification agent, 2026-08-08. All items below were fetched from the
primary URLs cited; quotes are verbatim as extracted by WebFetch unless marked otherwise.
Local primary source already on disk: `C:\Users\mouni\hc4\lit\meng-yang\jc_hc_status_note.tex`
(Meng-Yang, arXiv:2607.22198 v2, "Version 4 (arXiv v2) — July 26, 2026").

---

## 0. Bottom line (cross-source synthesis)

- **JC_2: open.** Every source that addresses dimension 2 says it is untouched. Tao: "The
  conjecture remains open in two dimensions, and is easy to establish in one dimension."
  jacobianfun.org: "These three-variable constructions do not settle the separate
  two-variable case." No 2026 arXiv paper on the plane Jacobian conjecture exists (newest
  in the feed is Lee-Li, Magnus formula IV, Aug 2024).
- **HC_4: open, and NOBODY has claimed it as of Aug 8, 2026.** An arXiv full-text search for
  "Hessian conjecture" (sorted newest first) returns exactly two papers: Meng-Yang
  arXiv:2607.22198 and Meng's 2003 paper math-ph/0308035. No preprint proves or refutes
  HC_4, and no one has published an obstruction theorem for Schur descent to dimension 4;
  the only statement on record is Meng-Yang's informal Remark ("Why the descent stops at
  five variables"), which is heuristic, not a theorem. **Track C territory is unclaimed.**
- The community's mechanism consensus: constant-Jacobian non-injectivity is a
  **non-properness / escape-to-infinity** phenomenon (étale but not proper covers of C^n),
  realized concretely by **multiplication/factorization of polynomials** (Tao, Speyer) and
  generalized by **tangent-sweep duality** (Gao). None of the mechanism accounts extend to
  gradient maps / Hessian potentials except through Meng-Yang's doubling + Schur descent.

---

## 1. Terence Tao, "A digestion of the Jacobian conjecture counterexample"
**URL:** https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/
**Date:** July 21, 2026.

**Mechanism as Tao presents it.** The counterexample is rebuilt from the multiplication map
on symmetric powers: F : Sym^1(C^2) x Sym^2(C^2) -> Sym^3(C^2), F(L,Q) := LQ; in coordinates
"F( (a,b), (c,d,e) ) = (ac, ad + bc, ae + bd, be)". A generic cubic factors as C = L1 L2 L3,
so the pairs "(L1, L2L3), (L2, L1L3), (L3, L1L2)" all map to the same cubic — the map is
generically **three-to-one**. Scaling symmetry is removed by imposing the resultant
normalization Res(L,Q) = 1 (SL_2-invariant; Res(lambda L, lambda^{-1} Q) = lambda Res(L,Q)).
The key step producing an honest map C^3 -> C^3 is that when an auxiliary third-order
differential operator D (dual to a hyperplane V) has a **repeated root**, the resulting
variety becomes isomorphic to C^3, via explicit polynomial coordinates in both directions
(fibering as a C^2-bundle over C^1 in the variable a, with smooth gluing at a = 0). Tao:
"I do not have a completely satisfactory geometric explanation for this miracle."

**Dimension 2.** Verbatim: "The conjecture remains open in two dimensions, and is easy to
establish in one dimension." No obstruction analysis for JC_2 is given; the 3-to-1
cubic-factorization mechanism has no stated 2-variable analogue.

**Hessian / gradient / Legendre content: none.** Tao does not mention gradient maps,
Hessians, symmetric Keller maps, or Legendre transforms anywhere in the post.

**Degree bounds / Moh: not mentioned.**

---

## 2. David Speyer, Secret Blogging Seminar, "The new counterexample to the Jacobian conjecture"
**URL:** https://sbseminar.wordpress.com/2026/07/20/the-new-counterexample-to-the-jacobian-conjecture/
**Date:** July 20, 2026.

**Geometric explanation.** The map (with Jacobian determinant -2) is "generically three to
one, not bijective". Speyer's picture: a map pi : P^1 x P^2 -> P^3 built from symmetric
functions (same factorization multiplication map as Tao's), with divisors H and R removed
from the source; the complement is an affine three-space despite the map being
three-sheeted. Will Sawin (comments) explains the C^3-structure as an **iterated
A^1-bundle**: "X as an iterated A^1 bundle... if we can figure out some sort of iterated
bundle structure on X, we probably win" — the space factors as an A^1-bundle over an
A^1-bundle over A^1, forcing X ~ A^3.

**Dimension 2.** The post itself contains no discussion of JC_2. A late comment by
Meng-Yang notes "JC_2, the classical two-dimensional Jacobian conjecture" remains
undecided, linked to HC_4.

**Hessian versions.** Only via the Meng-Yang comment (HC_5 false, HC_4 => JC_2, both open).
Another comment (wmayner): "the Dixmier conjecture fails for A_n, n >= 3."

---

## 3. Zihan Zhang, "Direct Consequences of the Three-Dimensional Counterexample to the Jacobian Conjecture"
**URL:** https://zzhang-iu.github.io/papers/direct-consequences-jacobian/
**Date on page:** July 20, 2026.

Consequence cascade derived from the explicit 3D Keller counterexample (credited to Claude
Fable 5; identity-coordinate stabilization gives all n >= 3):
- "the Mathieu Conjecture for SU(3) is false";
- Gaussian Moments Conjecture: "GMC(N) is false for at least one finite N" (existential, N unspecified);
- "the Vanishing Conjecture [Zhao] is false in some finite dimension" (existential);
- "the Image Conjecture is false in some finite dimension" (existential).

**Hessian conjecture content: none** (only "Hessian nilpotency of H" as a hypothesis inside
Zhao's framework). **No JC_2 discussion, no degree bounds, no August 2026 references.**
Relevance to us: the vanishing-conjecture failure is the same circle of ideas as
Zhao's TAMS 2007 Hessian-nilpotent formulation, but Zhang nowhere touches HC_n proper.

---

## 4. jacobianfun.org/jacobian-explained
**URL:** https://jacobianfun.org/jacobian-explained
**Maintainer:** "Created and curated by Alexis Gallagher, relying mainly on GPT-5.6-sol."

**Mechanism (analytic viewpoint): non-properness.** Verbatim quotes:
- "The number of (real) preimages can change only when a preimage runs away to infinity
  while its target stays bounded."
- "local invertibility controls every finite neighborhood, but it does not control what
  happens along a sequence that leaves every bounded region."
- "Finite merging would create a critical point. A fiber can instead lose points through a
  nonproper sequence."
- "local invertible maps can overlap themselves globally"; "the sheets of this map can
  never touch, anywhere."
The three preimages of (-1/4, 0, 0) lie on sheets (one flat, one curved, diverging like
1/x^2 as x -> 0) that never meet at any finite point.

**Dimension 2.** Verbatim: "These three-variable constructions do not settle the separate
two-variable case."

**Hessian / HC_4 content: none.**

---

## 5. arXiv scan, August 2026: has anyone settled or claimed HC_4? — NO

arXiv API full-text query all:"Hessian conjecture", newest first (fetched 2026-08-08):
returns ONLY (1) Meng-Yang arXiv:2607.22198v2 (Jul 24/27, 2026) and (2) Meng
math-ph/0308035 (2003). **No paper claims HC_4 either way; no Schur-descent obstruction
theorem for dimension 4 has been published.** The nearest thing on record is Meng-Yang's
Remark 4.2 (rem:selfterm, local TeX lines 359-370): "The Schur descent that produces HC_5
from HC_6 cannot be iterated to reach HC_4: the mechanism self-terminates... A single
descent already spends that linearity... Reaching HC_4 therefore appears to require
genuinely new ideas." — explicitly a remark, not a theorem.

Jacobian-conjecture papers, July-Aug 2026 (arXiv API all:"Jacobian conjecture", newest first):

| arXiv id | date | authors | content | JC_2/HC_4 relevance |
|---|---|---|---|---|
| 2608.05392 | Aug 5 | Castañeda, Honorato, Valenzuela-Henríquez | weak Markus-Yamabe conjecture fails for n >= 14 (chain realization + JC counterexamples) | none |
| 2608.00222 | Jul 31 | Shuhong Gao | "explicit examples of étale coverings C^n -> C^n that are not proper... fail to be injective only through points escaping to infinity"; tangent-sweep mechanism ("sweeps the tangent lines of a plane curve — a map that classical duality forces to hit most points several times"); counterexamples in every dim > 2 with arbitrarily large geometric degree; explicit degrees 4, 5, 10 (and 6, 12) in dims 3-5 | silent on JC_2; **no Hessian/gradient content** |
| 2608.02634 | Jul 29 | Romy Mondello | dimension-TWO counterexample to the **separable** Jacobian conjecture in **characteristic 2** | NOT char-0 JC_2; do not confuse |
| 2607.22198v2 | Jul 24/27 | Meng, Yang | HC_5 false; status table; HC_4 => JC_2 | the primary source on disk |
| 2607.21572v2 | Jul 23/30 | Piotr Migus | generic degrees of real Keller maps with non-dense image are exactly even d >= 4, n >= 3 | none |
| 2607.20968 | Jul 23 | Irit Huq-Kuruvilla | explicit char-2 separable JC counterexample | char 2 only |
| 2607.20597 | Jul 22 | Zbigniew Jelonek | polynomial automorphisms with Jac = 1 form a Zariski-closed set; generic elements of certain components are counterexamples for n >= 3, d >= 6 | none |
| 2607.20210v2 | Jul 22 | T. Shaska | graded/equivariant Keller maps: "an equivariant Keller map is always an automorphism" when all weights positive | positive-grading rigidity — possibly useful as a structural constraint for Track B/C |

Also found (non-arXiv): **Zenodo record 21504303** (Felipe Santibañez-Leal, CAOS
open-research programme, Santiago; v0.03, July 23, 2026), "The consequence cascade of the
Jacobian counterexample, with an explicit dimension-48 witness against the Hessian
conjecture": from Thompson's cubic-homogeneous form in dim 24, the de Bondt-van den Essen
construction "yields a homogeneous Hessian-nilpotent quartic potential in 48 variables (382
monomials over Q(i)) whose gradient map has an exact two-point collision"; "The Hessian
conjecture therefore fails in dimension 48 with an explicit witness." Predates/parallels
Meng-Yang's dim-5 result and is entirely superseded by it; makes **no claim about HC_4 or
JC_2**.

Plane-JC feed check (all:"two-dimensional Jacobian conjecture" OR "plane Jacobian
conjecture", newest first, fetched 2026-08-08): **no 2026 entries at all.** Most recent:
Lee-Li arXiv:2408.01279 (Magnus' formula revisited IV, Aug 2024); Moskowicz
arXiv:2407.13795 (no Keller maps with prime-degree field extension, Jul 2024).

---

## 6. Moh's theorem (JC_2 degree bound) — verified references

- **Moh's theorem (1983).** T. T. Moh, "On the Jacobian conjecture and the configurations
  of roots", *Journal für die reine und angewandte Mathematik* **340** (1983), 140-212.
  DOI: 10.1515/crll.1983.340.140 (verified at degruyterbrill.com and eudml.org/doc/152524).
  Statement (secondary-source phrasing, consistent across sources): the two-variable
  Jacobian conjecture holds for all Keller maps F = (P,Q) with deg F <= 100, where deg F is
  the maximum of deg P, deg Q. The proof is a (partly computer-assisted) analysis of
  possible degree pairs / configurations of roots at infinity. I did not obtain the 1983
  paper text itself; the theorem statement above is corroborated by multiple independent
  secondary sources including the paper below, which explicitly "confirm[s] the lower bound
  of 100 obtained by Moh".
- **Best known bound as of 2026: 108, not 100.** J. A. Guccione, J. J. Guccione,
  R. Horruitiner, C. Valqui, "Increasing the degree of a possible counterexample to the
  Jacobian Conjecture from 100 to 108", arXiv:2204.14178 (Apr 2022). Verbatim from
  abstract: they "discard them all, except the pair (72,108) (and the symmetric pair
  (108,72)), thus we confirm the lower bound of 100 obtained by Moh and raise it up to
  108." I.e., any counterexample to JC_2 must have max(deg P, deg Q) >= 108, and the only
  surviving degree pair at that level is (108,72)/(72,108). No improvement beyond 108
  surfaced in any 2026 search.
- Consequence for us: HC_4 => JC_2 means any HC_4 proof strategy needs no degree
  hypothesis, but conversely a JC_2 counterexample (which would kill HC_4 via the
  contrapositive of Meng's doubling only in dim 4 = 2x2) must have degree >= 108 —
  structured low-degree searches for JC_2-derived HC_4 counterexamples below degree ~108
  are provably futile.

---

## 7. Gaps / verification caveats

1. Moh's exact theorem wording is from secondary sources (the 1983 Crelle paper is
   paywalled); the degree-<=100 formulation with deg = max degree is uniformly attested.
2. WebFetch summaries of Tao's and Speyer's posts are model-condensed; the quoted fragments
   were returned as verbatim by the extraction pass, but for load-bearing use (e.g. citing
   Tao's exact coordinates), re-fetch the post and copy the formulas directly.
3. arXiv full-text search has indexing lag of a few days; a paper submitted Aug 6-8, 2026
   might not yet be indexed. As of the Aug 5 indexing horizon, no HC_4 claim exists.
4. Gao's paper was checked at the abstract level only; its TeX source has not been
   downloaded. If Track C wants the tangent-sweep mechanism in detail, pull
   arxiv.org/e-print/2608.00222.
