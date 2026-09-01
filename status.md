# Status of the Jacobian and Hessian conjectures — verified as of 2026-08-08/09

Every claim below was checked against a primary source during this project
(downloads under `lit/`; extract files cited). Machine verifications are in
`certificates/` (exact rational arithmetic, sympy 1.14, fail-closed asserts).

## Definitions

- **Keller map**: polynomial F: k^n → k^n with Jac F ∈ k^× (char k = 0 here).
- **JC_n** (Keller 1939): every Keller map in dimension n is invertible with
  polynomial inverse.
- **HC_n** (Meng, Appl. Math. Lett. 19 (2006) 503–510, arXiv:math-ph/0308035;
  independently de Bondt–van den Essen in nilpotent form): if
  f ∈ k[x_1..x_n] has det Hess f ∈ k^×, the formal Legendre transform f^L is
  a polynomial — equivalently (Cynk–Rusek/BBR), ∇f is injective, equivalently
  a polynomial automorphism.
  CAUTION: de Bondt–van den Essen's "HC(n)" (Proc. AMS 133 (2005) 2201–2205)
  is the different NILPOTENT statement: Hess f nilpotent ⟹ x + ∇f invertible.
  Their Theorem 1.1: nilpotent-HC(2n) ⟹ JC_n for nilpotent Keller maps (the
  construction f_H = −i·Σ H_k(x+iy)y_k uses i essentially). The
  nilpotent-symmetric case is SETTLED for n ≤ 4 (any H; de Bondt–van den
  Essen, JPAA 196 (2005) 135–148) and n = 5 homogeneous. This does NOT settle
  Meng's HC_4, which is the open problem attacked here.

## The status table (all entries source- and machine-verified)

| n        | 1    | 2        | 3     | 4        | ≥ 5   |
|----------|------|----------|-------|----------|-------|
| JC_n     | true | **open** | false | false    | false |
| HC_n     | true | true     | true  | **open** | false |

- **JC_3 false** (hence JC_n false ∀ n ≥ 3 by padding): Alpöge's explicit map
  (announced 2026-07-19, credited to the AI system Fable; recorded in
  arXiv:2607.22198 and arXiv:2608.00222), components of degrees 7,6,4,
  Jac ≡ −2, with the 3-point fiber F(0,0,−1/4) = F(1,−3/2,13/2) =
  F(−1,3/2,13/2) = (−1/4,0,0).
  **Machine-verified here**: `certificates/verify_alpoge_3d.py`.
  Community status: multiple independent machine verifications reported;
  formally "resolution reported — pending peer review" (openconjectures.org).
  Mechanism (Tao 2026-07-21; Speyer 2026-07-20; Gao arXiv:2608.00222): the
  map sweeps tangent lines of a plane curve — e.g. Tao's model
  F(L,Q) = LQ on binary forms, generically 3:1; non-injectivity happens
  entirely through non-properness (étale, no finite merging).
- **Generalization** (Gao, arXiv:2608.00222, 2026-07-31): counterexamples in
  every dimension > 2 and arbitrarily large geometric degree; explicit maps:
  G (3-dim, fiber 4), F_4, F_5 (4-dim, fibers 5 and 10, Jac −44/9 and
  160/29 — **both independently re-verified here** by direct symbolic
  determinants: `certificates/gradient_equivalence_test.py`), F_6, F_7
  (5-dim, fibers 6 and 12). Gallagher (2026-07-20, Zenodo): infinite family,
  every geometric degree ≥ 3.
- **HC_5 false** (hence HC_n false ∀ n ≥ 5 by ⊕½t²): Meng–Yang
  arXiv:2607.22198 (v2, 2026-07-27): Ψ = A² + 13A + 2B ∈ Z[x1,x2,y1,y2,y3],
  degree 14, 42 monomials, det Hess Ψ ≡ 128, gradient collision at
  (±1,∓3/2,0,0,0). Obtained from the 6-variable doubling of Alpöge's map by
  a one-variable **Schur descent** (partial Legendre transform); the full
  two-parameter family B + (λ/2)A² + μA has det Hess ≡ 4λ.
  **Machine-verified here** (two determinant algorithms):
  `certificates/verify_meng_yang_hc5.py`, `verify_mengyang_fast.py`.
  NOTE: extensive search (SBS comments, Tao's blog, arXiv, Zenodo, GitHub)
  found NO third-party verification of this counterexample prior to ours;
  nothing disputes it either. (A parallel Zenodo record of 2026-07-23,
  Santibañez-Leal, gives a dimension-48 Hessian witness independently.)
- **HC_2 true**: Dillen, JPAA 71 (1991) 13–18 (attribution per Meng and
  de Bondt). **HC_3 true**: de Bondt, "Polynomials with constant Hessian
  determinants in dimension three", JPAA 219 (2015) 3743–3754
  (arXiv:1203.6605v3) — Theorem 2.1: for char-0 K and n ≤ 3, ∇f satisfies
  the Jacobian conjecture, with explicitly polynomial inverse (Lemma
  "antitri": x_i ∈ A[F_{n+1−i},…,F_n]). Verified from the full TeX,
  independently twice (agents debondt-hc3 and hc3-check).
- **JC_2 open**: untouched by the 2026 counterexamples (Gao, Tao,
  jacobianfun.org all explicit on this). Degree bounds: Moh (Crelle 340
  (1983) 140–212, verified from page scans): true for deg ≤ 100; improved:
  Guccione–Guccione–Horruitiner–Valqui (arXiv:2204.14178): the only
  surviving counterexample degree pair up to 108 is (72,108)/(108,72).
- **HC_4 open**: "no result known to us settles HC_4" (Meng–Yang); confirmed
  by arXiv full-text search 2026-08-08: only Meng math-ph/0308035 and
  Meng–Yang 2607.22198 mention the Hessian conjecture; no partial results,
  no obstruction theorems, no claims.

## The two bridges (proofs re-verified line-by-line here)

- JC_n ⟹ HC_n (restriction to gradient maps).
- **HC_2n ⟹ JC_n** (Meng 2006, doubling φ(x,y) = ⟨y,F(x)⟩,
  det Hess φ = (−1)^n(Jac F)²; verbatim-verified in Meng's TeX). With n = 2:
  **HC_4 ⟹ JC_2** — also stated explicitly in de Bondt's paper
  ("an affirmative answer … in dimension four would imply the planar
  Jacobian conjecture"). So HC_4 is the STRONGER of the two survivors:
  proving HC_4 settles JC_2 affirmatively; a JC_2 counterexample would
  refute HC_4; and any JC_2 counterexample has max component degree ≥ 108.

## What is known about HC_4 specifically (complete list, as far as verifiable)

1. deg f ≤ 3: true (Wang's theorem, deg ≤ 2 Keller maps invertible — BCW
   Thm 2.4, with Oda's midpoint proof; verified from BCW).
2. Real definite case: true over R when Hess f is definite somewhere (Meng;
   via Jörgens–Calabi–Pogorelov). Also de Bondt Thm 4.2: anisotropic
   quadratic part at a point + det Hess ∈ K^× ⟹ deg f = 2 (n ≤ 3 or K = R;
   open for general K, n > 3 — de Bondt's Problem). Over C vacuous.
3. de Bondt's structural machinery FAILS at n = 4: both the anti-triangular
   normal form and the weight theorem have explicit 4-variable
   counterexamples (f = (x1+x2²)x3 + (x2+(x1+x2²)²)x4 — an invertible
   doubled planar Keller map — and Example 2.3/weightcounter). These defeat
   the METHOD, not the conjecture.
4. Meng–Yang's Schur descent "self-terminates" before dimension 4
   (their Remark 4.2 — a heuristic, not a theorem, prior to this project).

## New results of this project (see lemma_ledger.md, main_result.tex)

- **D1/D1+**: no member ψ_{λ,μ} of the Meng–Yang family admits ANY direction
  v ≠ 0 with D_v²ψ ≡ const (zero or not), in any affine coordinates: the
  known HC_5 counterexample family has no Schur pivot and no quadratic
  pivot; the known descent mechanism provably cannot be iterated to reach
  HC_4. (Rigorous form of Meng–Yang Remark 4.2 for the known family.)
- **G1**: neither Gao F_4 nor F_5 is equivalent to a gradient map under
  two-sided affine composition (symmetrizer space = 0); the known dim-4
  Jacobian counterexamples cannot refute HC_4 by affine massage. Consistency
  control: Alpöge's F also has symmetrizer space 0, as HC_3 demands.
- **Theorem A (pivot theorem)**: HC_4 holds for every f that is (in some
  affine coordinates) quadratic in a variable with NONZERO constant leading
  coefficient — with no hypothesis on the pivot coefficient — via a forced
  pencil identity + de Bondt's HC_3; and, when that leading coefficient is
  zero, for pivot coefficients of degree ≤ 2 (fiberwise Dillen HC_2 for
  affine e1; a rank collapse + PDE classification for quadratic e1).
  **Corollary: an HC_4 counterexample admits no quadratic pivot in any
  affine coordinates.**
- **Theorem B**: HC_4 holds for every f of degree ≤ 4.
- **Theorem C (doubling structure)**: a(x1,x2) + x3b + x4e has
  det Hess = (Jac(b,e))² independently of a, and ∇f is injective iff (b,e)
  is, for EVERY inert a; HC_4 restricted to this class is EQUIVALENT to
  JC_2. (Caution: the affine-pivot class is NOT exhausted by these
  doublings — see the trichotomy below.)
- **Theorems D/D′**: the affine-pivot constraint is exactly vanishing of
  the pivot coefficient's bordered Hessian (developable level surfaces);
  its leading form is therefore a Gordan–Noether cone.
- **Theorem G (rigidity, dimension-free)**: for e ∈ C[x_1..x_n], if the
  bordered Hessian ∇eᵀadj(Hess e)∇e vanishes identically then
  det Hess e ≡ 0. Geometrically: *if every level hypersurface of a
  polynomial is developable, its Hessian determinant vanishes identically.*
  Proved via the Legendre transform: the field V = (Hess e)⁻¹∇e is
  gradient-related to the Euler field, so a vanishing bordered Hessian
  makes the Legendre transform satisfy (E² − E)L = 0, forcing the inverse
  gradient to be bounded along rays — impossible for a polynomial gradient.
- **Corollary (pivot planarity, formerly Conjecture E)**: every
  3-variable polynomial with vanishing bordered Hessian is affinely
  2-variable.
- **THEOREM (the pivot trichotomy)**: a 4-variable potential with constant
  nonzero Hessian determinant admitting a pivot v (γ = D_v²f, e1 = D_vf)
  falls into exactly one of: (i) γ ≠ 0 — settled outright; (ii) γ = 0 and
  deg e1 = 1 — settled outright, but NOT always a doubling (witness
  w = x1²/2 + x1x2(1+x3) + x2²(2x3+x3²)/2 + x3x4, det Hess = 1, whose
  isotropic cone is a single line — certificate W13); (iii) γ = 0 and
  deg e1 ≥ 2 — exactly a Meng doubling of a planar Keller map plus an
  inert planar potential. **Hence HC_4 restricted to potentials admitting
  a pivot is EQUIVALENT to JC_2.** The remaining open case is potentials
  with no pivot at all.
- **Theorem G now has FOUR proofs** — the eigensplit proof (paper), the
  flow proof (advG), a Puiseux-at-infinity proof (pzG) and a
  derivation-eigenvector proof (atkG, shortest, purely algebraic) — three
  of them written by independent auditors instructed to refute it. The
  n = 2 case is machine-DECIDED for every polynomial of degree ≤ 5,
  linear terms included (`atkG_n2_decision.py`).
- **Theorem G, second proof (independent, advG)**: elementary flow proof —
  B ≡ 0 forces Jac(V)V = −V for V = (Hess e)⁻¹∇e, so the local flow is
  x(s) = ae^{−s} + b with ∇e(x(s)) = e^s p0; comparing Laurent
  polynomials kills p0. Valid over any characteristic-0 field (formal
  power series + transcendence of e^{−s}). The original write-up's
  analytic-continuation paragraph was found defective by the independent
  audit and has been replaced (in theorem_G_paper.tex) by an
  eigensplit-and-extend argument; both proofs are certified
  (advG_flow_proof.py, paperG_writeup_checks.py).

## New reformulation (2026-08-27): HC_4 as "f is a polynomial in its partials"

- **Lemma R**: for f with det Hess f ∈ C^×, Meng's HC holds for f
  ⟺ f ∈ C[∂₁f,…,∂_nf] (f is a polynomial in its own first partials)
  ⟺ Λ = ⟨x,∇f⟩ − f ∈ C[∇f]
  ⟺ the n commuting polynomial derivations W_i := adj(Hess f)·e_i are all
  locally nilpotent.
  So **HC_4 ⟺ every 4-variable constant-Hessian potential is a polynomial
  in its four partial derivatives.**
- Structural identities behind it (all machine-proved, fully generic):
  div(adj(H)∇f) = n·det H; W_i(∂_j f) = det H·δ_ij; the pullback along ∇f
  of {∂_i, Euler} is a Lie-algebra homomorphism ([W_i,W_j] = 0,
  [V,W_i] = −W_i), giving a commuting polynomial framing off {det H = 0};
  and the Legendre relations W_iΛ = x_i, Vf = B/det H.
- Necessary condition (fibre constancy): ∇f(p) = ∇f(q) ⟹ f(p) = f(q). The
  Meng–Yang counterexample recovers HC-failure via non-injectivity of ∇Ψ
  (two distinct points, equal gradients); NOTE the Ψ-values there coincide
  (both 0), so fibre-constancy alone does not witness it — necessary, not
  sufficient (self-corrected, caught by a fail-closed assert).
- Cert: `euler_pullback_reformulation.py`. Ledger entry: §R.
- Degree-5 attack (rank-3, no-isotropic-direction, pivot-free branch),
  weighted leading form F = a5(y) + x4·b3(y) + ½x4²·y1, five closed-form
  graded equations E4…E0:
  - **Generic-b3 emptiness — PROVED and independently certified.** E4 = 0
    splits the cubic b3 into two rational families. Treating the cubic
    coefficients c as PARAMETERS over ℚ(c), the 21-unknown a5-system
    J = ⟨coeff_y E0…E3⟩ has std(J) = (1) in each branch — the unit ideal
    over ℚ(c) (instant in Singular; out of sympy Gröbner reach with c as
    unknowns). A Nullstellensatz certificate Σgᵢ·Jᵢ = 1 was extracted
    (Singular `lift`) and re-verified exactly in sympy
    (`verify_paramcert.py`). So both branches are EMPTY off the exceptional
    c-loci Z0 = V(972 c8⁷c9⁴(3c0c9−c1c8)⁶), Z1 = V(3 c5(2c0c5−c1c4)⁵).
    Ledger: R-D5-GEN.
  - **Tail (the loci Z0, Z1) — closed by a finite recursive stratification**
    (`_d5_close.py`): each c-substratum is decided by std(J) (= (1) ⇒ empty)
    or a cone test; the recursion is Noetherian and every branch reduces to a
    single degenerate family, to which every survivor leaf collapses:
    **b3 = y1²·(c2y1 + c0y2 + c1y3)** (det B ≡ 0, Hess b3 rank 2, kernel
    v* = (0,c1,−c0)). This family is NOT a rank-collapse cone — the E-system
    genuinely has rank-3 solutions (det Hess₃a5 ∉ √J, unbiased full-variety
    Rabinowitsch, 12/12 across 3 primes). BUT every solution's a5 carries the
    isotropic/linear direction v*: **D²_{v*}a5 = v*ᵀ(Hess₃a5)v* ∈ √J**, with
    the explicit power certificate (v*ᵀAv*)² ∈ ⟨J⟩. Such a5 has a linear
    direction, so it FAILS the branch's "no-isotropic-direction" hypothesis —
    it lies OUTSIDE this branch (which assumes no linear direction). Hence the
    rank-3, no-isotropic-direction, pivot-free branch is **empty, conditional on
    the char-0 lift of the modular fact D²_{v*}a5 ∈ √J** (which alone suffices;
    a cone would also be excluded). LABEL: the structural core S1–S4
    (`d5_survivor_family.py`) is fully deductive/char-0, but the two decisive
    radical facts are certified only MODULARLY (multi-prime, sampled c and λ;
    the single-prime F_p power reduction (v*ᵀAv*)² ≡ 0 mod std(J) is NOT char-0),
    because the char-0 parametric lift over ℚ(c) is compute-bound. Per the
    project's rule (no char-0 theorem from modular samples) the emptiness is
    [machine-modular], not a hand proof; also the solve()-based stratification is
    not independently proved exhaustive, and the c0=c1=0 sub-locus (b3=c2y1³,
    v*=0) is closed separately as a rank-collapse cone, det Hess₃a5 ∈ √J
    (`_d5_surv_c0c1.py`, 3/3 primes). Ledger: R-D5-TAIL. NOTE: closing this one branch
    does not settle HC_4 — the ultimate blocker remains JC_2 because HC_4 ⟹ JC_2;
    the WITH-pivot (affine-pivot) class is JC_2-equivalent by the trichotomy,
    while the pivot-free residual is precisely the part NOT reduced to JC_2 (the
    genuine open heart). The rank-3 f5 that HAS a linear direction is NOT a separate
    leading-form branch: its normal form f5 = P(x2,x3)+x1Q(x2,x3) pushed into
    the tower gives J = (1) at generic b3 (char 0; `d5_lindir_verify.py`),
    subsumed by R-D5-GEN — so its leading forms ARE the survivors above. The
    genuine open frontier is one level up — LIFTING a survivor leading form to
    a full quintic f with det Hess f ∈ C^× and deciding pivot-freeness (NOT
    covered by Theorem A: a linear direction of f5 gives D²_v f = D²_v(f4+f3+f2),
    degree 2, not a pivot).
