# Certificates — reproducibility

Environment: Windows 11, Python 3.12 (invoke as `py`, not `python` — the
Store stub shadows PATH), sympy 1.14.0 (`py -m pip install --user sympy`).
All scripts are fail-closed: any assertion failure exits nonzero. All
arithmetic is exact (rational/symbolic); no floating point anywhere.

Proof-status convention: a symbolic identity check whose expression is
LINEAR in generic symbolic coefficients is a complete proof for all
polynomials of that shape (multilinearity); such checks are marked
"fully generic" / "proof". Random-instance checks are consistency checks
only; the corresponding hand proofs live in ../lemma_ledger.md.

## Landscape verification (2026 counterexamples)

- verify_alpoge_3d.py — Alpöge's JC_3 counterexample: det J ≡ −2, three
  distinct rational points with a common image. (Berkowitz determinant.)
- verify_meng_yang_hc5.py — Meng–Yang: HC_6 doubling φ (det Hess ≡ −4,
  3-point gradient collision), HC_5 Ψ (deg 14, 42 monomials,
  det Hess ≡ 128, 2-point collision), Schur family det ≡ 4λ (symbolic
  λ, μ). (Berkowitz.)
- verify_mengyang_fast.py — same claims, independent determinant
  algorithm (domain-ge). Two-algorithm redundancy.

## Tower and pivot machinery (this project's theorems)

- ../src/hc4lib.py (`py src\hc4lib.py`) — T0/T1/T3 spot checks + fully
  symbolic proof of the midpoint identity T4 (all quartics).
- t3_generic.py — T3 fully generic (proof for all degree-4 f).
- verify_q2_core.py — Theorem Q2: Schur identity (3 instance classes incl.
  a nontrivial pencil), symbolic collision-transfer step.
- ap4_rank_reduction.py — AP3/AP2/AP1: K·adjK·K = detK·K (generic
  symmetric 3×3), adj-decomposition, the two contradiction evaluations.
- ap4_pivot_closure.py — AP1 collapse: adj(P+tE11) linear in t (generic),
  cofactor identity g^T adjP g = −(x1b2−b1)² (generic jet), PDE solution
  family, det Hess ≡ κ² for FULLY GENERIC a (deg 4) + generic cubic β +
  symbolic κ, and a saturated collision-Gröbner emptiness instance.
- deg4_branch_closures.py — S1 (2×2 lemma key branch), S2 ([det]_6 and
  [det]_5 fully generic), S3 (isotropic-vector identities: adj M0 =
  d1d2·v0v0ᵀ and the t-coefficient formula, generic).
- dillen4_ap4_control.py — known-good control: de Bondt's dillen4 example
  satisfies the AP identities with c = 1 and has empty saturated
  collision system (injective), as required.

## Structure theorems (affine-pivot class, all degrees)

- doubling_structure.py — Theorem C: det[[A,B],[Bᵀ,0]] = det(B)² at the jet
  level (proof), the polynomial-level check with fully generic inert
  potential a, the AP1 family's Jacobian and determinant, completeness of
  the planar Keller classification with e = x2 + x1²/2, and a collision
  Gröbner instance with a ≠ 0.
- developability.py — Theorem D: the bordered-Hessian identity
  det[[K,g],[gᵀ,0]] = −gᵀadj(K)g (fully symbolic ⇒ proof); the doubling
  class satisfies it; nondegenerate quadric pivots do not; a homogeneous
  cubic search over the (c2) coefficient ideal.
- cone_leading_form.py — Identity E (Euler contraction) verified fully
  symbolically for ternary forms of degrees 2,3,4; the top-degree
  extraction on a generic cubic with arbitrary lower-order terms;
  Theorem D′ (leading form of the pivot coefficient is a cone) checked
  against six ternary-cubic representatives.
- homogenization_identity.py — Identity H:
  det Hess₄(E)|_{x0=1} = d(d−1)·e·det Hess₃(e) − (d−1)²·B(e), verified with
  fully generic e of degrees 2 and 3.
- theorem_F.py — Theorem F (Conjecture E holds when det Hess₃(e) ≡ 0, in
  particular for all homogeneous pivot coefficients): translation
  invariance of B (fully symbolic), the converse-of-Euler step, and
  end-to-end examples.
- conjE_degree3.py — Conjecture E in degree 3: radical-membership
  (Rabinowitsch) certificates testing whether every cubic pivot
  coefficient with vanishing bordered Hessian is affinely 2-variable.
  NOTE: this is a heavy Gröbner computation (14 variables per minor) and
  may run for hours; conjE_search.py is the fast complement.
- conjE_search.py — fast exact search for a counterexample to Conjecture E
  over cubic and quartic pivot coefficients with cone leading form, plus
  structured developable ansätze. Reports any exact solution that is not
  affinely 2-variable (none found).

## Obstruction certificates

- pivot_obstruction_d1.py — D1/D1+: 95 homogeneous v-quadratics from the
  degree-≥6 coefficients of the Meng–Yang family; Gröbner leading-term
  ideal contains a pure power of each v_i ⇒ common zero cone = {0} ⇒ no
  affine-coordinates pivot (D_v²ψ ≡ const, any const) for any λ ≠ 0, μ.
- gradient_equivalence_test.py — G1: parses Gao's F4, F5 from the LaTeX
  source (../lit/gao/Jacobian_CE.tex), checksums by re-deriving
  det J ≡ −44/9 and 160/29, then computes the linear symmetrizer space
  {W : W·J_F symmetric} for Alpöge F, F4, F5 — all {0}.

## Sharp form of Theorem A

- theoremA_sharp.py — the σ ≠ 0 branch of the pivot theorem holds for pivot
  coefficients of ANY degree: the Schur/product-rule identity and the
  collision-transfer identity are verified fully symbolically for generic
  e0, e1 of degree 4 (so K = Hess₃e1 is non-constant), plus four explicit
  instances with non-constant K and empty saturated collision systems.
  Corollary: an HC_4 counterexample admits no quadratic pivot in any
  affine coordinates.

## Theorem G and the pivot dichotomy

- theorem_G.py — Theorem G: the two structural identities fully symbolic for
  n = 2,3,4; the Legendre identities EL − L = e∘ψ on exact examples; the
  homogeneous eigenvalue d/(d−1) ∉ {0,1} sanity check; and the decisive
  radical-membership test in THREE variables (every coefficient of det Hess
  vanishes on the whole variety B ≡ 0). PASSED.
- theorem_G_n2_elementary.py — the n = 2 case: elementary proof in the
  header (adj H = JᵀHJ, so B = Tᵀ(Hess e)T with T tangent to the level
  curves), plus radical membership. Degrees 2, 3 PASSED; the run was
  killed mid-degree-4 by a mass process kill. CAVEAT (found by the pzG
  audit): its scan DROPS LINEAR TERMS, on which B depends, so it is NOT
  exhaustive as its header claims; and the hand proof's "reduced fibre"
  appeal should be "smooth generic fibre". The genuinely full decision
  is pzG_n2_decision.py.
- theorem_G_n2_deg45.py — degrees 4–6, fast by design: (1) 180 exact
  rational instances with det Hess ≠ 0, all with B ≠ 0 (a single B ≡ 0
  hit would refute Theorem G); (2) an exact decision on a slice in
  degrees 4 and 5; (3) positive control e = φ(L). PASSED. CAVEAT (found
  by the atkG audit): the TEST-2 slice parametrises only the TOP FORM,
  i.e. it is HOMOGENEOUS, where Theorem G is a one-liner by Identity E —
  so its "complete decision" does NOT cover inhomogeneous e. The
  genuinely complete degree ≤ 5 decision (linear terms included, via a
  certified normalisation e = x2^d + lower terms) is atkG_n2_decision.py.

NOTE ON RUNTIME: this environment kills background jobs after roughly ten
minutes with exit code 255. Keep certificates short; prefer normalised
slices and exact sampling over full-coefficient-space Gröbner bases.
- pivot_dichotomy.py — the corrected TRICHOTOMY: the fully symbolic identity
  [x4] det Hess f = −(e0)_{x3x3}·B₂(e), branch (b)'s determinant formula,
  and an instance of each branch. PASSED. (The first write-up asserted only
  branch (a); this is the correction.)
- identityE_is_known.py — ATTRIBUTION: proves the project's Identity E is the
  s-linear coefficient of a published proposition. Cases n = 3 (d = 2,3,4)
  and n = 4 (d = 2) PASSED before the mass kill — sufficient, since the
  equivalence is a formal algebraic fact.
- why_dimension_four.py — a REFUTED hypothesis, recorded: the Meng–Yang
  leading form IS a cone (x1⁶x2⁶y1²), so Gordan–Noether failure does NOT
  explain why that counterexample has no pivot; the obstruction lives in
  lower graded pieces (degrees 12 and 8).

Superseded and not re-run: conjE_search.py and conjE_degree3.py hunted for
counterexamples to what was then Conjecture E. It is now a theorem
(Theorem G + Theorem F), and theorem_G.py's three-variable radical test
covers the same ground far more cheaply.

## Red-team artifacts

Three red-team workflow runs were terminated by model/session quota before
any agent returned a structured result — but the agents had already written
their test scripts to disk. Those were recovered and run:

- rt_novelty_q2_anydegree.py — **found a real defect**: the deg e1 ≤ 2
  hypothesis is unnecessary when σ ≠ 0. Superseded by theoremA_sharp.py.
- rt_instances_*.py, rt_instances_x2_*.py — hostile instance attacks on
  Theorems A and B: random sparse quartic ansätze with det Hess f forced
  constant by Gröbner elimination, then pivot-cone emptiness (Rabinowitsch
  saturation) and injectivity tests on every solution. Fail-closed: a
  non-injective Keller solution or an empty pivot cone would be a serious
  gap (or an HC_4 counterexample).
- rt_obstructions_parser.py — audit of the F4/F5 LaTeX parser used in G1.
- rt_ref_landscape.py — independent re-derivation of the 2026 landscape
  from the primary TeX only, not from this project's own scripts.
- deg4_q2v_*.py — engine-agent certificates for the Q2 chain (bordered
  Schur form, pencil forcing, collision transfer, AP identities).

Some of these are long-running (minutes to hours). conjE_degree3.py in
particular computes a 14-variable Gröbner basis per minor and is better
run in Singular or Macaulay2 than sympy.

## Six-agent push (agent keys advG, d5, pe, nax, tgc, paperG)

All recovered from disk after the Aug 10 quota kill and re-run to
completion on Aug 18 (all exit 0 unless noted).

- advG_flow_proof.py — INDEPENDENT adversarial verification of Theorem G.
  Verdict: theorem TRUE, but the original analytic-continuation paragraph
  was a non-sequitur; the header contains a complete replacement proof
  (local flow x'' + x' = 0, Laurent-polynomial contradiction) needing no
  continuation, valid over any characteristic-0 field. Certifies the two
  structural identities J1/J2 fully generically (n = 2,3,4; random exact
  jets n = 5,6), the Step-3 conclusion Jac(V)V = −V, non-vacuity on
  rational witnesses of Euler-degree 0, Identity E cross-checks, and the
  converse-failure control. ALL PASSED.
- d5_graded_tower.py — degree-5 machinery: universal closed forms of the
  graded pieces C12–C9 of det Hess f for f = f5+f4+f3+f2 with f5 x4-free
  (proof via universal determinant expansion). C11 = (f4)_x4x4·det3 A5;
  C10|_{d4=0} = (f3)_x4x4·det A5 − qᵀadj(A5)q. ALL PASS.
- d5_pivotfree_normalform.py — the rank-3, no-pivot branch in degree 5:
  weighted normal form f = a(y) + x4 b(y) + x4² y1/2, the weighted-T1
  analogue, and the five x4-graded components of det Hess of the weighted
  leading form.
- pe_pivot_cone.py / pe_lowdeg_pivot.py / pe_perturbation.py — pivot
  existence: exact Gröbner pivot-cone decision tool + baseline (every
  catalogued 4-variable potential has a pivot; Meng–Yang HC_5 does not);
  PE-1 (deg ≤ 4 pivot existence), PE-2 (dim span{N_α} ≤ 3 forces a
  pivot), PE-3 (pivot-free needs ≥ 2 homogeneous parts of degree ≥ 3);
  PE-4 (adj-isotropic linear-form perturbations preserve det Hess,
  fully generic proof) + catalogue scan: the perturbation route produces
  NO pivot-free potential from the current catalogue. ALL PASS.
  (pe_perturbation.py needs PYTHONIOENCODING=utf-8 on Windows consoles —
  one banner contains "∩"; the mathematics is encoding-independent.)
- nax_d1_nonaffine.py — non-affine scope of D1. PROVED: the naive
  strengthening ("no polynomial change of coordinates gives Meng–Yang a
  pivot") is FALSE — an explicit non-affine triangular automorphism tau
  makes every Psi_{lam,mu} affine-linear in a coordinate; but tau
  destroys the constant Hessian (det Hess(Psi∘tau) = lam·R², R
  non-constant), so the corrected statement ("no automorphism gives a
  pivot WHILE keeping det Hess constant") survives this attack. Plus the
  MASTER FORMULA det Hess_5 Theta = lam·(Jac F)²|_{x3=lam⟨a,z⟩+mu}
  (fully generic proof) and the fiber-linear reduction to a Keller
  condition. ALL PASS.
- nax_pivot_transfer_5var.py — honest GAP record: the proposed
  Gordan–Noether route to a full non-affine D1 breaks at exactly the
  Perazzo phenomenon (A0 = w2+w1w3+w1²w4 has B_4 ≡ 0, det Hess_4 ≡ 0,
  and is NOT affinely 3-variable). The corrected non-affine D1 remains
  OPEN. ALL stated identities PASS.
- tgc_debondt_chain.py — independent second proof of Corollary E via
  de Bondt–van den Essen's Singular-Hessians classification (branch
  identities B = lam²·det Hess_2 a and B = −W(a2,a3)², jet-level
  proofs). Corrects the task's false premise: the dBvdE exceptional
  family does NOT always have B ≠ 0 — B vanishes exactly on its
  affinely-2-variable members. ALL PASS.
- tgc_geometry_and_sharpness.py — the curvature dictionary
  (B = Gauss–Kronecker curvature × |dF|^{n+2}, re-derived; prior art
  Reilly 1986 Prop. 4, Fox Math. Nachr. 290 (2017)); SHARPNESS: the
  algebraic non-polynomial F = x1 − sqrt(x1²−2x2) has B ≡ 0,
  det Hess ≠ 0, and its local inverse has EXACTLY the form
  psi(tp) = t⁻¹a + b predicted by the proof — only the polynomiality
  step fails. Theorem G is genuinely algebraic. ALL PASSED.
- tgc_beyond_hc4.py — placement of Theorem G: R2 (B ≡ 0 ⟺ det Hess_n e
  ≡ 0 AND det Hess_{n+1} E ≡ 0, both classical conditions); R3
  (Corollary E FALSE for n ≥ 4 — dehomogenized Perazzo witnesses); R4
  (4-variable classification via Watanabe–de Bondt); R5 (homogeneous
  case is one line via Identity E — the content of Theorem G is the
  non-conical case); R6 (B ≡ 0 propagates to the whole pencil of
  projective closures). ALL PASSED.
- paperG_structure.py — RED at check (P4), by design of the correction:
  it tests ideal membership where radical membership is required
  (Gröbner basis {v1²+2v3v4, v1v2, v2², v1v3, v2v3, v3²} contains no
  pure v1 power). P1–P3 and the first half of P4 pass; the corrected
  test is (W13) in paperG_writeup_checks.py. Do not "fix" this script;
  it documents the defect.
- paperG_writeup_checks.py — the full W1–W15 checklist of
  ../theorem_G_paper.tex, including (W12) doubling collision transfer
  for arbitrary inert a, (W13) the not-a-doubling witness w by
  Rabinowitsch radical membership, and (W15) sharpness witnesses.
  ALL PASSED (~4 min).

## Second independent pass (Aug 19–20; agent keys pzG, wG, gcv)

- pzG_audit_identities.py — SECOND independent adversarial audit of
  Theorem G. Verdict: TRUE; the original continuation step is repairable
  exactly as in theorem_G_paper.tex, and the header contains a THIRD
  proof — purely algebraic, embedding C(x) into the Puiseux field at
  t = infinity, valid over every algebraically closed characteristic-0
  field, with no analytic function theory. Machine pillars: H·adj(H) =
  det(H)·I (generic, n = 2,3,4), the Legendre identities fully symbolic,
  the Puiseux mechanism on an exact branch, and the final
  no-positive-powers step proved generically. ALL PASS.
- pzG_n2_decision.py — repairs a REAL DEFECT found by pzG: the
  "exhaustive" scan in theorem_G_n2_elementary.py dropped linear terms,
  but B depends on them (B(x1²/2 + x2) = 1, B(x1²/2) = 0), so it was not
  exhaustive as claimed. This script decides Theorem G in two variables
  over the FULL coefficient space (all monomials of degree 1..d),
  degrees 2–4, by radical membership over the whole coefficient ring.
  Also corrects the hand proof: the generic fibre is smooth (avoids
  critical values), hence its line components are parallel —
  "reduced" alone would not suffice.
- wG_writeup_checks.py — the writeup agent's own certificate for the
  revised theorem_G_paper.tex: bordered-determinant formula, the
  [x4]-coefficient identity, the doubling determinant, the division-free
  Schur identity for quadratic pivots, the case-(ii) formulas, the
  collision-transfer identity — all fully generic — plus the cautionary
  example h = x1 + x2² showing "det Hess_2 e ≡ 0 ⟹ e = φ(L)" is FALSE
  as stated (the correct hypothesis is B_2(e) ≡ 0); main_result.tex's
  trichotomy proof has been corrected accordingly.
- gcv_family_wronskian.py — the de Bondt–van den Essen exceptional
  family h = g1(l)x1 + g2(l)x2 + g3(l)x3 has B(h) = −det[a | g | g′]²
  (adapted: −W(g2,g3)²), proved fully generically for deg ≤ 4 plus a
  matrix-level weight-2 covariance proof; so B(h) ≡ 0 ⟺ h affinely
  2-variable. Corrects the task's literal claim ("always B ≠ 0"). PASS.
- atkG_algebraic_proof.py — THIRD independent adversarial audit of Theorem
  G (agent key atkG). Verdict: TRUE; the paper's eigensplit proof
  re-derived line by line and found complete; the original defect
  classified as expositional. Header: a FOURTH proof — purely algebraic
  via derivation eigenvectors (δ = D_V, δg = g, δV = −V, finite Taylor
  expansion of the polynomial g at the δ-constant point b = x + V splits
  into eigencomponents with distinct eigenvalues −k, forcing g = 0). No
  analysis, no Puiseux series, any characteristic-0 field. Pillars fully
  generic, ALL PASS (44 s).
- atkG_n2_decision.py — Theorem G in TWO variables DECIDED for every
  polynomial of degree ≤ 5, linear terms included: certified
  normalisation (B affine-covariant; top form has det Hess ≡ 0 by
  Identity E; a binary form with zero Hessian is a d-th power, by
  Hankel-minor radical membership) reduces to the complete slice
  e = x2^d + lower terms, on which radical membership proves e ∈ C[x2].
  Degrees 2–5 ALL PASS in 23 s. Supersedes theorem_G_n2_elementary.py
  (linear terms dropped) and theorem_G_n2_deg45.py TEST 2 (homogeneous
  slice only).
- atkG_n3_search.py — computational REFUTATION ATTEMPTS on Theorem G in
  three variables over complete coefficient slices (modes: deg3,
  deg4:p4/p31/p22/p211/p1111, pivot, sample). Results: `sample` (exact
  contrapositive sampling, degrees 3–5) ALL PASS; `deg4:p211` complete
  slice ALL PASS (52 s); the 17–25-unknown Gröbner bases of `deg3`,
  `deg4:p4`, `deg4:p31`, `deg4:p22`, `deg4:p1111`, `pivot` do not finish
  within the environment's runtime kill (no DONE line) — NOT refutations;
  the deg-3 complete decision is the same computation already PASSED in
  theorem_G.py. No counterexample found anywhere.
- pzG_n2_decision.py NOTE: degree 4's 14-variable Gröbner does not finish
  under the runtime kill (degrees 2–3 PASS); superseded by
  atkG_n2_decision.py.
- refG_paper_integrity.py — referee-pass certificate for
  ../theorem_G_paper.tex: structural compile-check (no undefined refs /
  cites, no duplicate labels, nested environments, balanced braces,
  single document skeleton, no macro clashes) plus an independent exact
  re-derivation of the manuscript's six load-bearing new identities
  (I)–(VI) with fully generic coefficients. ALL CHECKS PASSED.
- gcv_ny_homogeneous_converse.py — the homogeneous CONVERSE of Theorem
  G: B(e) = (d/(d−1))·e·det Hess e makes B ≡ 0 ⟺ det Hess ≡ 0 on cones;
  Identity E and Nagaoka–Yazawa Prop. 2.3 [J. Algebra 577 (2021)
  175–202] proved formally equivalent in BOTH directions. Off cones the
  converse fails (x1²/2 + x2). PASS.

## New mathematics (session 2026-08-27)

- euler_pullback_reformulation.py — **Lemma R** and its structural
  identities. For f with det Hess f in C^*, Meng's HC for f is equivalent
  to "f is a polynomial in its own first partials", f in C[grad f]
  (equivalently Lambda = <x,grad f> - f in C[grad f]; equivalently the n
  commuting polynomial derivations W_i = adj(Hess f) e_i are all locally
  nilpotent). Identities certified fully generically (n = 2,3,4):
  (I1) div(adj(H) grad f) = n det H; (I2) W_i(d_j f) = det H * delta_ij and
  W_i(x_j) symmetric; (I3) [W_i,W_j] = 0, [V,W_i] = -W_i (the pullback of
  {d_i, Euler} along grad f is a Lie-algebra homomorphism); (I4) the
  Legendre relations W_i Lambda = x_i, V f = B/det H. B1: the automorphism
  x1x2 + x1^3/3 is a polynomial in its gradient and its W_i are LND. B2:
  the Meng-Yang HC_5 counterexample recovers HC-failure via
  NON-INJECTIVITY of grad Psi. SELF-CORRECTION recorded: an earlier draft
  claimed Psi violates fibre-constancy at its collision; it does NOT (both
  Psi-values are 0 there) -- fibre-constancy is necessary, not sufficient,
  and this collision does not witness the failure. ALL CHECKS PASSED.
  (NOTE: B2 imports verify_mengyang_fast, which recomputes a degree-14
  determinant -- ~12 min; the mathematics of B2 needs only Psi and its
  gradient at two points.)
- d5_rank3_pivotfree_decision.py — reduces the degree-5, rank-3,
  no-isotropic-direction, PIVOT-FREE branch to a finite decision. Verifies
  the closed forms E4..E0 of det Hess_4 of the weighted leading form
  F = a5 + x4 b3 + x4^2 y1/2; solves E4 == 0 (exactly two rational
  branches; the b3 == 0 sub-branch is EMPTY by an elementary domain
  argument); solves E3 (linear in a5); exports each residual E2=E1=E0
  radical-membership decision to d5_branchK_famM.sing (Singular) since
  sympy's Groebner does not terminate on ~17 unknowns here; and runs a
  fail-closed exact sampling of the residual variety -- 400 exact points
  per branch, NO rank-3 survivor found. A YES from the exported ideal
  computation closes degree-5 HC_4 in this branch; a NO produces the first
  pivot-free 4-variable constant-Hessian potential. TERMINATES cleanly.
- d5_branch1_normalized.py — completes the E4-branch-1 decision by first
  applying the GL_2 normalisation c4 = 0 (the (y2,y3)-part of b3 is a
  perfect square, so its linear form can be set to y3 without loss),
  shrinking the parameter count enough to attempt the radical-membership
  test directly in sympy.
- NOTE: Singular / Macaulay2 / Magma are NOT installed on this machine;
  installing one is the single change that closes the exported degree-5
  decisions (and the other runtime-killed Groebner certificates).
