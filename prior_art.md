# Prior art — honest attribution ledger

Result of an adversarial novelty audit (2026-08-10) whose instruction was to
FIND prior art, not to confirm novelty. Several project "results" turned out
to be known. They are demoted here and cited properly. This file overrides
any stronger claim made elsewhere in the project.

## KNOWN — must be cited, not claimed

**Identity E (Euler contraction).** ∇eᵀadj(Hess e)∇e = (d/(d−1))·e·det Hess e
for e homogeneous of degree d ≥ 2, any n.
- D. J. F. Fox, *Equiaffine geometry of level sets and ruled hypersurfaces
  with equiaffine mean curvature zero*, Math. Nachr. 290 (2017) 293–320
  (arXiv:1503.09108), §3 — verbatim, with the same proof (Euler +
  H·adj(H)·H = det(H)H).
- A. Del Pia, R. Hildebrand, R. Weismantel, K. Zemmer, arXiv:1408.4711,
  Lemma 5.2 — same identity for any C² homogeneous function (more general
  than our polynomial statement); attributed further to D. Hemmer,
  *Limiting curvature near singular points of algebraic curves*,
  manuscript, July 1995 (n = 2 case).
- Also the cone-specialization of Dolgachev's (1.20) below.
Verified equivalent to our statement in `identityE_is_known.py`.

**Identity H (homogenization identity).**
det Hess_{n+1}(E)|_{x₀=1} = d(d−1)·e·det Hess_n(e) − (d−1)²·B(e).
- I. V. Dolgachev, *Classical Algebraic Geometry: A Modern View*, CUP 2012,
  §1.1.4, equations (1.20)–(1.21), pp. 17–18 — "the affine equation of the
  Hessian", general n, **with our exact proof** (row reduction by Euler's
  relation, then the mirrored column reduction). Dolgachev uses the ternary
  case explicitly on p. 20.
- The lemma ledger previously called this "derived and used here". That was
  **not defensible** and is retracted.

**Theorem D (developability).** B(e) ≡ 0 ⟺ every level hypersurface has
vanishing Gauss–Kronecker curvature.
- R. C. Reilly, *Affine geometry and the form of the equation of a
  hypersurface*, Rocky Mountain J. Math. 16 (1986) 553–565, Prop. 4.
- Fox, op. cit., §2: names exactly our invariant, U(F) = U^{ij}F_iF_j, proves
  U(F) = 𝒦·|dF|^{n+2}, and credits Reilly.
- The underlying determinant identity det[[A,b],[cᵀ,d]] = d·det A − cᵀadj(A)b
  is **Cauchy's**.
- Classical differential geometry: Spivak vol. 3 p. 204; Knoblauch 1913
  pp. 89–94; collected as eq. (4.1) of Goldman, CAGD 22 (2005) 632–658.
- Algebraic-geometry side: Dolgachev, Remark 1.1.18 + Theorem 1.1.20.

**Lemma S3 (isotropic vector for all-singular symmetric subspaces).** Known
in strictly greater generality (all n, bounded rank).
- R. Gow, arXiv:1502.05547, Theorem 1 — our statement for general n over any
  field with |K| ≥ r+1, **with our proof line for line**; attributed there to
  Fillmore–Laurie–Radjavi, *On matrix spaces with zero determinant*, Linear
  and Multilinear Algebra 18 (1985) 255–266, proof of Lemma 1.
- Equivalently the classical **Bertini theorem for a linear system of
  quadrics**; stated verbatim with the maximal-rank refinement in
  J. M. Landsberg, arXiv:math/0609507.

**Lemma S1 (2×2 symmetric of linear forms, det ≡ 0 ⟹ ℓ·uuᵀ).** The r = 1
case of Meshulam / Loewy–Radwan on symmetric matrix spaces of bounded rank
(R. Meshulam, LAA 114–115 (1989) 261–271; R. Loewy, N. Radwan, *Spaces of
symmetric matrices of bounded rank*).

**Theorem F.** Our homogenization proof was unnecessary: it is a three-line
corollary of a published classification we had not cited —
M. de Bondt, A. van den Essen, *Singular Hessians*, J. Algebra 282 (2004)
195–204, Theorem 3.3 (n = 3, det Hess h = 0 ⟹ h is linearly equivalent to
a₁(x₁,x₂) + λx₃ or to a₁(x₁) + a₂(x₁)x₂ + a₃(x₁)x₃).

## FOLKLORE / ELEMENTARY — state, don't claim

- **Theorem C** (doubling structure, det Hess = (Jac(b,e))²): the block
  determinant is elementary; the doubling and its Jacobian relation are
  Meng, Appl. Math. Lett. 19 (2006) 503–510, Prop. 1.4, and Meng–Yang
  arXiv:2607.22198 Prop. 2.3. Our contribution is only the observation that
  the class is *exactly* this shape and the resulting JC₂ equivalence.
- **T0, T1** (nondegenerate quadratic part; leading form has vanishing
  Hessian determinant): implicit in Meng–Yang §2 and in de Bondt
  arXiv:1203.6605, proof of Thm 1.1(i).
- **T3** and the branch graded identities: no prior statement located, but
  the method (weight/leading-part bookkeeping on det Hess) is de Bondt's.
- **T4** (midpoint identity): the gradient case of Oda's proof of Wang's
  theorem, as printed in BCW Bull. AMS 7 (1982) p. 298.

## UNDETERMINED — no prior art located

- **Theorem G** (B(e) ≡ 0 ⟹ det Hess e ≡ 0, any n). Reilly and Fox supply
  the framework and the identity B = det Hess e · ⟨∇e, (Hess e)⁻¹∇e⟩, but
  **both work under the standing hypothesis that this quantity does not
  vanish** — Fox's paper is about affine spheres and assumes U(F) ≠ 0 and
  H(F) ≠ 0 throughout. The content of Theorem G is precisely that the
  complementary case ⟨∇e,(Hess e)⁻¹∇e⟩ ≡ 0 with det Hess e ≢ 0 is
  impossible over C for polynomials. Over R with definite Hessian this is
  vacuous, which may be why it is not in the affine-geometry literature.
  Status: no prior art found; needs expert review before being claimed.
- **The pivot dichotomy** (HC₄-with-a-pivot ⟺ JC₂) and **Theorem A**'s
  σ = 0 classification: no prior art located, but the σ ≠ 0 half is a
  natural completing-the-square argument that may be folklore.
- **D1/D1⁺** and **G1**: the objects are weeks old; novelty is near-certain
  but untested by referees.

## Audit gaps

Three novelty agents (Monge–Ampère literature, reduction techniques, and the
HC₄/recent-arXiv sweep) died on quota or a stalled stream before reporting.
What is known from partial work: Jörgens–Calabi–Pogorelov (convex entire
solutions of det D²u = 1 are quadratic) does **not** apply — it needs
convexity/reality, and over C there are many non-quadratic polynomial
solutions (every Meng doubling). de Bondt's Theorem 4.2 is the algebraic
analogue and holds only for n ≤ 3 or K = R.
