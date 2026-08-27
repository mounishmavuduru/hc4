# Research log — HC_4 project

Chronological, candid. Dates are session-time (2026-08-08, America/Chicago).

## Phase 0 — status verification (done)

- Premise check: my training knowledge (cutoff Jan 2026) had JC open in all
  dimensions >= 2. Web verification confirms the July 2026 upheaval is real:
  - Alpöge announced an explicit counterexample to JC_3 on 2026-07-19
    (credited to the AI system Fable); Gao (arXiv:2608.00222, 2026-07-31)
    generalized to every dimension > 2 with a tangent-sweep mechanism.
  - Meng–Yang (arXiv:2607.22198, v2 2026-07-26) refuted HC_5 by an explicit
    degree-14 integer polynomial; with de Bondt's HC_{<=3} theorem and the two
    bridges, HC_n is decided except n = 4.
  - Surviving open statements: **JC_2 and HC_4**, linked by HC_4 => JC_2.
- Machine verification (exact rational arithmetic, sympy 1.14, two independent
  determinant algorithms — Berkowitz and domain-Gaussian-elimination):
  - `certificates/verify_alpoge_3d.py`: det J_F = -2 identically; the three
    rational points (0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2) share the image
    (-1/4,0,0). JC_3 is false. PASSED.
  - `certificates/verify_meng_yang_hc5.py` + `verify_mengyang_fast.py`:
    - 6-var doubling phi: det Hess = -4; 3-point gradient collision. HC_6 false.
    - Psi = A^2 + 13A + 2B: degree 14, 42 monomials, det Hess = 128;
      grad Psi(±1, ∓3/2, 0,0,0) = (0,0,-1/2,0,0). HC_5 false.
    - Whole Schur family: det Hess(B + (lam/2)A^2 + mu A) = 4*lam identically
      in symbolic lam, mu. PASSED.
- Hand-verified the Meng–Yang paper's proofs of: the two bridges
  (JC_n => HC_n; HC_2n => JC_n via doubling, with the block-swap sign
  (-1)^n), the two stabilization lemmas, the bordered-determinant and
  rank-one-update algebra of the Schur-descent lemma. All correct.

## Phase 1 — Track C: descent obstruction (first new result)

- Key observation: Meng–Yang's Remark 4.2 ("the descent self-terminates;
  reaching HC_4 appears to require genuinely new ideas") is a HEURISTIC, not
  a theorem. Target: make it rigorous.
- A Schur descent from 5 to 4 variables requires a *pivot*: after an affine
  change, the potential must be affine-linear in one coordinate, i.e. there
  must exist v != 0 with D_v^2 psi == 0 identically (plus the reduced-Hessian
  degeneracy det M(s,w) == 0). Affine changes are the right generality: they
  are exactly the coordinate changes preserving the class
  {det Hess = const}.
- **Lemma D1 (proved, machine-checked)**: for every lam != 0 and every mu, the
  only v in C^5 with D_v^2 psi_{lam,mu} == 0 is v = 0.
  Proof: deg-counting shows every w-monomial coefficient of degree >= 6 of
  D_v^2 psi equals lam times the corresponding (mu-free) coefficient of
  C(v) = A·D_v^2 A + (D_v A)^2; these 95 coefficients are homogeneous
  quadratics in v; a Gröbner basis of them has leading-term ideal containing a
  pure power of each v_i, so their common zero locus (a cone) is {0}.
  `certificates/pivot_obstruction_d1.py`. PASSED.
  Consequence: **no member of the known HC_5 counterexample family admits any
  Schur pivot in any affine coordinates; the known descent cannot reach
  HC_4.**

## Phase 2 — Track B: can the known dim-4 JC counterexamples refute HC_4?

- Gao's paper contains two explicit dim-4 counterexamples: F4 (degrees
  4,11,12,21; det J = -44/9; generic fiber 5) and F5 (degrees 3,12,14,16;
  det J = 160/29; generic fiber 10).
- Reduction (proved, elementary): F is affinely equivalent (two-sided) to a
  gradient map iff there is W in GL_4 with W·J_F(x) symmetric identically
  (W absorbs both sides as (T^T)^{-1} S; conversely W J_F symmetric makes the
  1-form <W F, dx> closed, hence exact with polynomial potential).
  So the test is LINEAR in W. If an invertible W exists for F4 or F5, HC_4 is
  refuted outright; if not, the known dim-4 counterexamples cannot settle
  HC_4 by affine massage.
- Consistency control: for Alpöge's 3-dim F an invertible symmetrizer would
  contradict de Bondt's HC_3 theorem; the test must return none.
- Transcription safety: F4, F5 are parsed from Gao's LaTeX source
  programmatically and checksummed by re-deriving det J = -44/9 and 160/29
  exactly. `certificates/gradient_equivalence_test.py`. RUNNING.

## Phase 3 — literature workflow

- First launch of the 10-agent primary-source workflow failed wholesale
  (session limit); relaunched. Awaiting: de Bondt 1203.6605 exact statements
  (esp. what fails at n=4), Meng 2006 provenance, dBvdE 2005 bookkeeping,
  GN theory statements, Wang/Moh classical statements, community status,
  Meng–Yang verification landscape.

## Phase 4 — tower, Q2, and the degree-4 case tree (2026-08-08/09)

- hc4lib built and self-tested: T0/T1/T3 graded formulas spot-checked; T4
  midpoint identity PROVED fully symbolically (linear in the 70 generic
  quartic coefficients, so the symbolic check is a proof for all quartics).
- **Theorem Q2 found and proved** (quadratic pivot with constant corner ⇒
  forced pencil identity ⇒ descent to de Bondt's HC_3 ⇒ injectivity).
  Machine-checked incl. a nontrivial-pencil instance
  (`certificates/verify_q2_core.py`). This kills half the degree-4 tree.
- **D1 strengthened to D1+** for free: the same degree-≥6 subsystem excludes
  D_v²ψ ≡ const, so the Meng–Yang family has no quadratic pivot either.
- Degree-4 tree restructured by pivot type (essential rank of f4 → x4-degree
  → σ): (I.a) σ≠0 closed by Q2; (I.b) = class AP4 (f = e0 + x4·e1) with the
  three identities (c0),(c1),(c2); (II) r≤2 branches always produce a pivot
  direction (u-perp construction) → Q2 or AP4. Scaffold with derivations:
  `src/deg4_tree_notes.md`. Five-agent workflow launched (Q2V, AP4, BR2,
  BR1, C2X); first run died on session limit, resumed.

## Phase 5 — literature harvest (Track A complete)

- All key dependencies source-verified (8 fetch agents + 1 verify agent):
  de Bondt HC_3 (twice, full TeX; anti-triangular + weight theorems FAIL at
  n=4 by explicit examples; his dillen4 example is an AP4 doubling of a
  planar Keller map — a perfect engine test case), Meng 2006 provenance +
  doubling verbatim, dBvdE-2005 (nilpotent form! n≤4 nilpotent-symmetric is
  SETTLED — distinct from Meng's HC_4), Zhao 2007, GN via Watanabe–de Bondt
  (n=5 classification K[x1,x2][Δ]), Wang via BCW (Oda's midpoint proof —
  our T4 is its gradient sharpening), Moh ≤100 from Crelle scans + the 2022
  improvement to 108 (only pair (72,108) survives), BBR/Ax/Cynk–Rusek.
- Community check: NOBODY has settled/claimed HC_4 or published a descent
  obstruction; our D1/G1/Q2 are new. No third-party verification of the
  Meng–Yang HC_5 counterexample exists — ours appears to be the first.
- Deliverables written: `status.md`, `definitions_and_reductions.tex`;
  ledger updated with verified attributions (Dillen 1991 for HC_2!).

## Phase 6 — the degree-4 closure (2026-08-09/10)

- Engine workflow died twice on session limits; closed its objectives
  inline instead:
- **AP rank reduction** (AP0–AP3): for deg e1 = 2, rank Hess(e1) ∈ {2,3}
  impossible ((c2)'s quadratic part gives detK·K = 0; then (c0) at the
  origin); rank 1 forces the normal form e1 = x1²/2 + x2; rank 0
  (affine e1) closes for ALL degrees via fiberwise Dillen (AP0).
  `ap4_rank_reduction.py`.
- **AP1 collapse**: adj(P + tE11) is linear in t with g-contraction P33 ⇒
  ∂33e0 ≡ 0; the cofactor identity turns (c0) into the PDE
  x1·b_x2 − b_x1 = κ with polynomial solutions −κx1 + β(x2 + x1²/2);
  the resulting family f = a + x3(−κx1 + β(ξ)) + x4ξ has det Hess ≡ κ²
  (verified for fully generic a!) and explicitly injective gradient.
  **AP4 with deg e1 ≤ 2 closed in all degrees.** `ap4_pivot_closure.py`.
- **Branch closures**: [det]_6 and [det]_5 formulas fully generic (S2);
  2×2 lemma (S1); **isotropic-vector lemma S3** (every all-singular
  subspace of Sym₃(C) has a common isotropic vector — 3-line proof via
  adjM0 = κv0v0ᵀ and the t-coefficient of det(M0+tM) ≡ 0).
  `deg4_branch_closures.py`, `t3_generic.py`.
- **⇒ THEOREM A (pivot theorem) and THEOREM B (HC_4 in degree ≤ 4).**
- Control: de Bondt's dillen4 (the example that killed anti-
  triangularization at n = 4) verified as an AP4 member with c = 1,
  satisfies (c0)-(c2), collision Gröbner = [1]. `dillen4_ap4_control.py`.
- Deliverables written: main_result.tex, reviewer_packet.md,
  certificates/README.md. Red-team workflow launched (6 agents: chain
  audit, instance attack, obstruction audit, citation audit, novelty
  search, C2X structure completion).

## Phase 7 — self red-team, and the structure theorems (2026-08-10)

- The red-team workflow was killed three times by model/session limits
  (Fable 5 quota, then session quota). Rather than leave the audit undone,
  I ran the adversarial pass over the whole chain by hand. Every step of
  Theorems A and B re-derived independently:
  - Q2's Schur identity re-derived by hand (Hess₃ẽ0 + uK = M − ggᵀ/σ ✓);
  - the "u sweeps C" step is legitimate (triangular polynomial change);
  - deg e1 ≤ 2 verified in EVERY branch of Theorem B (it holds because f4
    has no x4 in any branch, so the x4-linear terms come only from f3, f2);
  - r = 2's u-perp argument: v ∈ span(e3,e4) only sees the (3,4)-block ✓;
  - r = 1's isotropic argument: M(x) is a genuine linear subspace of
    Sym₃ (M(0) = 0, M linear) and det ≡ 0 makes every member singular ✓;
  - AP1's PDE completeness re-derived by characteristics ✓.
  No errors found — but two genuine improvements fell out:
- **THEOREM C (doubling structure)**: f = a(x1,x2) + x3b + x4e has Hessian
  [[A,B],[Bᵀ,0]], so det Hess f = (Jac(b,e))² *independently of a*. That is,
  this class IS the Meng doubling class, and HC_4 restricted to it is
  EQUIVALENT to JC_2. My AP1 classification is exactly "planar Keller maps
  with one component x1²/2 + x2". `doubling_structure.py`.
- **THEOREM D / D′**: (c2) is exactly the vanishing of e1's bordered
  Hessian (developable level surfaces). Via the new **Identity E**
  (∇eᵀadj(Hess e)∇e = d/(d−1)·e·det Hess e for homogeneous e — Euler +
  H·adjH·H = detH·H), the top-degree part of (c2) forces
  det Hess(leading form of e1) ≡ 0, so by Gordan–Noether the leading form
  is a CONE. For deg e1 = 2 this reproduces AP3 — a satisfying consistency
  check on both. `developability.py`, `cone_leading_form.py`.
- **Conjecture E** (every 3-variable polynomial with vanishing bordered
  Hessian is affinely 2-variable) — proved for deg ≤ 2 and for all
  homogeneous e1; degree 3 tested by radical-membership Gröbner
  certificates (`conjE_degree3.py`). If true in all degrees, the entire
  affine-pivot class collapses to doublings and HC_4 there is exactly JC_2.
  Route noted: the homogenization identity det Hess₄(E)|_{x0=1} =
  d(d−1)e1·det Hess₃e1 − (d−1)²·(bordered), derived here.

## Phase 8 — Theorem F, and the sharpest remaining gap (2026-08-10)

- Derived the **homogenization identity** by hand and verified it fully
  symbolically for degrees 2, 3:
    det Hess₄(E)|_{x0=1} = d(d−1)·e·det Hess₃(e) − (d−1)²·B(e).
  (Row/column reduction of Hess E at x0 = 1 into
  [[d(d−1)e, (d−1)gᵀ],[(d−1)g, H]], then Cauchy.)
- **THEOREM F**: B(e) ≡ 0 AND det Hess₃(e) ≡ 0 ⟹ e affinely 2-variable.
  Proof: the identity kills det Hess₄(E); GN in FOUR variables makes E a
  cone; if the cone direction is transverse to the hyperplane at infinity,
  Euler's converse turns e into a translate of a homogeneous polynomial,
  and then Identity E + GN(3) finishes. So **Conjecture E holds whenever
  det Hess₃(e1) ≡ 0 — in particular for every homogeneous pivot
  coefficient.** `theorem_F.py`, `homogenization_identity.py`.
- Residual gap is now sharp and stated: B(e1) ≡ 0 with det Hess₃(e1) ≢ 0.
  The identity forces the top TWO graded pieces of e1·det Hess₃e1 to
  vanish, so det Hess₃ e1 has unusually low degree — the natural next
  lever.
- Searches: `conjE_search.py` (exact solve over cubic/quartic pivot
  coefficients with cone leading form + structured developable ansätze)
  and `conjE_degree3.py` (Rabinowitsch radical-membership per minor; heavy).
  A fail-closed assert in the first run caught a sympy solver artifact
  (specializations hitting poles produced zoo) — fixed by skipping
  non-finite specializations, which is the right behavior since those are
  not solutions.
- Published a referee-facing summary artifact of the whole project.

## Phase 9 — harvesting the dead agents' work (2026-08-10)

- The three red-team workflow runs all died on quota WITHOUT returning
  structured results — but a directory listing showed the agents had
  written certificate scripts to disk before dying: `rt_instances_*.py`,
  `rt_obstructions_parser.py`, `rt_novelty_q2_anydegree.py`,
  `rt_ref_landscape.py`, plus `deg4_q2v_*.py` from the earlier engine run.
  Lesson: check the filesystem, not just the workflow return value.
- Reading `rt_novelty_q2_anydegree.py` revealed a REAL FINDING the referee
  agent had reached before being killed: **the hypothesis deg e1 ≤ 2 is
  unnecessary in the σ ≠ 0 branch of Theorem A.** I had believed K =
  Hess₃e1 must be globally constant for the pencil argument. It need not
  be: the sweep is PER-POINT in x′ — for fixed x′, K(x′) IS a constant
  matrix and u = x4 + e1(x′)/σ still sweeps C. So det₃(Hess₃ẽ0(x′) +
  s·K(x′)) ≡ c/σ for all s and all x′, and ẽ0 + λe1 has constant Hessian
  determinant for every scalar λ, which is all HC_3 needs.
- Re-derived and re-certified independently: `theoremA_sharp.py`.
  **THEOREM A (sharp)**: σ ≠ 0 ⟹ HC_4, for pivot coefficients of ANY
  degree. **Corollary: an HC_4 counterexample admits no quadratic pivot in
  any affine coordinates.**
- Building valid test instances was itself instructive: a fail-closed
  assert rejected my first attempt because the constructed f did not
  actually have constant Hessian determinant. Valid high-degree instances
  need a pencil that is constant in both s and x′; the anti-triangular
  base ẽ0 = x1x3 + x2²/2 + ψ(x1) with e1 = x2h(x1) + k(x1) supplies them,
  giving genuinely non-constant K with deg e1 = 3, 4.

## Phase 10 — full harvest verdicts (2026-08-10)

Ran every recovered agent script. Summary (details in reviewer_packet §6b):
- **Independent landscape re-derivation PASSED** (`rt_ref_landscape.py`,
  written from the primary TeX only, not from our scripts).
- **Q2 chain fully proved** by the four `deg4_q2v_*` certificates,
  including a clean new proof that the s³-coefficient of det₃(M + sK) is
  det₃K — so the pencil forces rank K ≤ 2 directly (an independent route
  to AP3).
- **G1 transcription doubly verified**: an independent strict re-parser
  agrees with the production parser on every component of F₄/F₅, source
  term counts equal parsed monomial counts (nothing dropped or merged),
  and det J is re-derived both symbolically and at random rational points.
  The script's 2 reported "FAILURES" are mislabeled checks, not defects.
- **Instance attacks all survived**: 45 constant-Hessian quartics from 6
  ansatz families (some over Q(i)); 24 classified-family instances of
  degree 3–6; 12 AP0 instances of degree 4–8. Negative controls confirm
  the tests are not vacuous (a planted collision was correctly detected).
- Two scripts errored, both in the agents' own code:
  - `rt_instances_AP0.py`'s "negative control" x1x2 + x2²x3 + x3x4 has
    det Hess ≡ 1 — it is a valid AP0 instance, not a control. Theorem A
    predicts injectivity; verified. No defect
    (`rt_ap0_negative_control_resolved.py`).
  - `rt_instances_x2_B_break.py` crashed on a leftover placeholder line
    (`Poly.homogeneous_order()` returns None for non-homogeneous input);
    repaired and re-run.

## Phase 11 — THEOREM G, and the pivot dichotomy (2026-08-10)

The breakthrough. Conjecture E is now a theorem, via a much stronger and
dimension-free statement.

- **Where it came from**: staring at the constraint (c2), B(e) ≡ 0, and
  asking what it means dynamically rather than algebraically. Where Hess e
  is invertible, set V := (Hess e)⁻¹∇e. Then Dφ·V = ∇e = φ for φ = ∇e — so
  **V is φ-related to the EULER field** on the target — and
  D_V e = B/det Hess e. So B ≡ 0 says exactly: *e is a first integral of
  the pullback of the Euler field under its own gradient map.*
- Transport through the Legendre transform L (∇L = φ⁻¹, and the classical
  EL − L = e∘φ⁻¹): the condition becomes **(E² − E)L = 0**, so L has only
  Euler-eigenvalues 0 and 1, so φ⁻¹(tp) = t⁻¹a + b along rays. Then the
  continued identity ∇e(t⁻¹a + b) = tp is an identity of Laurent
  polynomials in t whose left side has no positive power of t ⟹ p = 0,
  contradiction.
- **THEOREM G**: B(e) ≡ 0 ⟹ det Hess e ≡ 0, in ANY number of variables.
  Geometrically: *if all level hypersurfaces of a polynomial are
  developable, its Hessian determinant vanishes identically.* Converse
  false (x1²/2 + x2). `theorem_G.py` — all structural identities fully
  symbolic; the decisive radical-membership test in 3 variables passes
  (every det-Hess coefficient vanishes on the entire variety B ≡ 0).
- Sanity anchor: for homogeneous e the Legendre transform has Euler
  eigenvalue d/(d−1), which is never 0 or 1 — exactly why homogeneous
  potentials can never satisfy B ≡ 0 with det Hess ≠ 0. Machine-checked.
- **n = 2 proved independently and elementarily**: adj H = JᵀHJ gives
  B = Tᵀ(Hess e)T with T tangent to the level curves, so B ≡ 0 ⟺ all level
  curves are lines; a generic fibre is then a product of PARALLEL linear
  forms (non-parallel ones would make the fibre singular), and comparing
  two generic fibres forces one direction, so e ∈ C[L].
  `theorem_G_n2_elementary.py`, with exhaustive radical-membership over the
  full coefficient space in degrees ≤ 5.
- **Corollary (was Conjecture E)**: every 3-variable polynomial with
  vanishing bordered Hessian is affinely 2-variable (Theorem G + F).
  Strictly stronger than det Hess ≡ 0: de Bondt–van den Essen's
  h = x1x2 + x1²x3 has det Hess ≡ 0 but B(h) = −x1⁴ ≠ 0 ✓.
- **MAIN THEOREM (pivot dichotomy)**: a 4-variable constant-Hessian
  potential with a pivot is either settled outright (quadratic pivot) or is
  exactly a Meng doubling of a planar Keller map plus an inert planar
  potential. **HC_4 restricted to potentials with a pivot ⟺ JC_2.**
- **Self-caught gap**: my first statement of the dichotomy asserted that
  the affine case forces (e0)_{x3x3} ≡ 0 (⟹ doubling). Re-deriving from
  scratch showed the true condition is (e0)_{x3x3}·B₂(e) ≡ 0 — a genuine
  TRICHOTOMY. The extra branch B₂(e) ≡ 0 leads, via Theorem G in TWO
  variables, to e = φ(L), and then det Hess f = −φ′(x1)²·det₂Hess_{(x2,x3)}e0
  forces φ′² to divide a nonzero constant, so φ is affine and we land in
  AP0 — unconditional. So the equivalence survives, but the first proof as
  written was incomplete. `pivot_dichotomy.py` verifies both identities
  with fully generic coefficients. Nice illustration that Theorem G is
  load-bearing in two different dimensions inside the same proof.
- **The residual class is now named exactly**: potentials with NO pivot.
  They exist in dimension 5 (D1: the Meng–Yang counterexample). If they do
  not exist in dimension 4, then HC_4 ⟺ JC_2 outright.

## Phase 13 — prior-art audit lands hard (2026-08-10)

Adversarial novelty sweep, instructed to FIND prior art. It did. Full
verdicts in `prior_art.md`; headline retractions:
- **Identity E**: KNOWN — Fox, Math. Nachr. 290 (2017) §3, verbatim with my
  proof; also Del Pia et al. arXiv:1408.4711 Lemma 5.2 (C² functions,
  more general), traced to Hemmer 1995.
- **Identity H**: KNOWN — Dolgachev, *Classical Algebraic Geometry*,
  eq. (1.20)–(1.21), "the affine equation of the Hessian", general n, same
  proof. My ledger had said "derived and used here" — indefensible,
  retracted.
- **Theorem D** (developability): KNOWN — Reilly 1986 Prop. 4; Fox names my
  invariant B as U(F) and proves U(F) = 𝒦·|dF|^{n+2}. Determinant identity
  is Cauchy's.
- **Lemma S3**: KNOWN in greater generality — Gow arXiv:1502.05547 Thm 1,
  my proof line for line, credited to Fillmore–Laurie–Radjavi 1985;
  equivalently classical Bertini for linear systems of quadrics.
- **Lemma S1**: r = 1 case of Meshulam / Loewy–Radwan.
- **Theorem F**: three-line corollary of de Bondt–van den Essen,
  *Singular Hessians*, J. Algebra 282 (2004) Thm 3.3 — a paper I never
  cited. My homogenization proof was unnecessary.
- **Theorem C**: folklore given Meng's doubling; only the "exactly this
  shape ⟹ JC_2 equivalence" observation is mine.
- **Theorem G survives**: Fox and Reilly supply the framework and the
  identity B = det Hess e·⟨∇e,(Hess e)⁻¹∇e⟩ but both assume it does NOT
  vanish (their subject is affine spheres). Theorem G is exactly the
  complementary case — vacuous over R for definite Hessians, with bite only
  over C. No prior art located. Still flagged pending expert review.

Lesson: elementary identities are almost always classical and are findable
only by searching for the *identity*, never the topic. Bibliography added
to main_result.tex; all over-claims patched in the ledger.

## Phase 14 — a refuted hypothesis, recorded (2026-08-10)

I conjectured the HC_4/HC_5 dimension gap was Gordan–Noether: in n ≤ 4 the
leading form must be a cone (supplying a candidate pivot), in n ≥ 5 not.
**REFUTED by my own fail-closed assert.** The Meng–Yang leading form is the
single monomial x1⁶x2⁶y1² — a cone in a 2-dimensional space of directions.
Candidate pivots exist; they die at *lower* graded degrees (12 and 8
respectively). So the dimension gap is not a leading-form phenomenon, and
the no-pivot question in dimension four cannot be settled by leading-form
arguments alone. `why_dimension_four.py`.

## Phase 12 — novelty audit (2026-08-10)

- **Identity E is NOT new.** It is the s-linear coefficient of Nagaoka–
  Yazawa, arXiv:1904.01800, Prop. `identity1`:
  det(−F·H_F + s(∇F)ᵀ∇F) = (−1)^{n−1}(r/(r−1))(s − (r−1)/r)F^n det H_F.
  Verified against their arXiv TeX and proved equivalent in
  `identityE_is_known.py`. Cited, not claimed. Their proof uses the same
  two ingredients as mine (Euler + rank-one determinant algebra), which is
  why it was findable only by searching for the *identity*, not the topic.
- Jörgens–Calabi–Pogorelov (convex entire solutions of det D²u = 1 are
  quadratic) does NOT apply: it needs convexity/reality. Over C there are
  many non-quadratic polynomial solutions — every Meng doubling is one.
  de Bondt's Thm 4.2 is the algebraic analogue and holds only for n ≤ 3
  or K = R. So the Monge–Ampère literature does not settle HC_4.
- A six-agent prior-art workflow is running against every remaining item.

## Phase 15 — mass process kill, and what survived (2026-08-10)

All five long-running exact computations died simultaneously at exit 255 —
a mass kill, not individual failures. Partial results were recovered from
their buffered output:
- `rt_instances_x2_B_break.py` (the strongest attack on Theorem B, designed
  to empty the pivot cone): reached 14 unknowns / 102 equations / 13
  solution branches and **every instance passed** — det constant, pivot
  cone P1 and P2 both nonempty, gradient injective. Theorem B survived as
  far as it got.
- `identityE_is_known.py`: n = 3 (d = 2,3,4) and n = 4 (d = 2) all passed —
  enough to establish the attribution, the equivalence being formal.
- `theorem_G_n2_elementary.py`: degrees 2 and 3 passed; killed in degree 4.
  Completed separately by `theorem_G_n2_deg45.py`.
- `conjE_search.py`, `conjE_degree3.py`: superseded — they hunted for
  counterexamples to what is now Theorem G + Theorem F, and
  `theorem_G.py`'s three-variable radical test already covers that ground.

## Phase 16 — n = 2 case of Theorem G closed (2026-08-10)

The environment kills background jobs at ~10 minutes (exit 255) — that is
what destroyed the five runs above, and a retry of the degree-4/5 n = 2
check as well. Rewrote it to finish in about a minute and got a STRONGER
test in the bargain: 180 exact rational instances of degree 4–6 with
det Hess ≠ 0, all having B ≠ 0 (one B ≡ 0 hit would have refuted Theorem
G), plus an exact decision on a GL₂-normalised slice in degrees 4 and 5,
plus the positive control e = φ(L). All passed. [CORRECTED 2026-08-25:
that slice is the HOMOGENEOUS top form only — trivial by Identity E — so
the "complete decision" was overclaimed; found by the atkG audit and
superseded by atkG_n2_decision.py, a genuine degree ≤ 5 decision.]

Method note for future runs: prefer normalised slices and exact sampling
over Rabinowitsch-per-coefficient over the full coefficient space; print as
you go (`py -u`) so a kill still leaves usable partial output.

## Phase 17 — the six-agent push recovered from disk (2026-08-18)

The fourth run of the six-agent workflow died on the weekly quota on
Aug 10 — but, exactly as in Phase 15, the agents had already written
their work to disk before dying: sixteen certificates and a complete
57KB referee-ready paper (`theorem_G_paper.tex`). All were recovered and
re-run to completion today (the workflow itself was also relaunched
after the quota reset, as a further independent pass). Outcomes:

1. **Theorem G survived a genuinely independent adversarial audit — after
   a repair.** The `advG` agent judged the original analytic-continuation
   paragraph a non-sequitur (ψ is only a local inverse; det Hess may
   vanish on the ray) and supplied a completely elementary replacement
   proof: B ≡ 0 forces Jac(V)V = −V for V = (Hess e)⁻¹∇e, the local flow
   is x(s) = ae^{−s} + b with ∇e(x(s)) = e^s·p₀, and comparing Laurent
   polynomials in τ = e^{−s} forces p₀ = 0 — contradiction. Works over
   any characteristic-0 field. The paper contains a second repaired
   proof (Euler-eigensplit + explicit homogeneous extension + identity
   theorem), audited by the project lead. Both are certified; the
   sharpness certificate shows the algebraic function x₁ − √(x₁²−2x₂)
   satisfies every step except polynomiality.
2. **A false claim of ours was found and corrected.** "The affine-pivot
   class is exactly the Meng doublings" is FALSE: witness
   w = x₁²/2 + x₁x₂(1+x₃) + x₂²(2x₃+x₃²)/2 + x₃x₄ (det Hess = 1, affine
   pivot, isotropic cone a single line — no doubling has that). The
   correct statement is the pivot TRICHOTOMY by (γ, deg e₁); the
   JC₂-equivalence corollary survives unchanged. `main_result.tex`,
   `status.md` and the paper now state it correctly. A related defect:
   `paperG_structure.py` check (P4) tested ideal membership where
   radical membership was needed — it is deliberately left RED and
   superseded by (W13).
3. **Non-affine D1: the naive strengthening is FALSE, the corrected one
   is open.** An explicit non-affine automorphism τ makes every
   Meng–Yang Ψ affine-linear in a coordinate (so "no polynomial change
   of coordinates" is simply wrong), but τ destroys det Hess constancy
   (= λR², R non-constant) — the det-Hess hypothesis is load-bearing.
   New master formula det Hess₅Θ = λ(Jac F)²|_… proved. The proposed
   Gordan–Noether route to the corrected statement breaks at Perazzo
   (`nax_pivot_transfer_5var.py`) — honest gap, recorded.
4. **Degree 5**: universal graded pieces C12–C9 in closed form; the
   pivot-free rank-3 branch has weighted normal form
   a(y) + x₄b(y) + x₄²y₁/2. Not yet a proof of degree-5 HC₄.
5. **Pivot existence**: PE-2 (dim span{N_α} ≤ 3 forces a pivot), PE-3
   (pivot-free needs two storeys of degree ≥ 3), PE-4 perturbation
   machinery — no pivot-free 4-variable potential found; every
   catalogued potential has dim P(g) ≥ 2.
6. **Placement of Theorem G**: homogeneous case is classical (one line
   from Identity E); Corollary E is sharp — false in n ≥ 4 exactly by
   dehomogenized Perazzo; B ≡ 0 ⟺ both classical vanishing-Hessian
   conditions (e and its closure) hold simultaneously.

## Phase 18 — second independent pass; a third proof; two more defects found and fixed (2026-08-20)

The relaunched workflow ran 6.2 hours and lost five of six agents to the
session limit — but every one of them had already written to disk
(pattern now firmly established: check the filesystem, not the return
value). The one agent that returned (consequences/gcv) completed the
literature novelty check: NO prior art for Theorem G located across the
developable-hypersurfaces / degenerate-Gauss-map literature
(Gordan–Noether, Perazzo, Ciliberto–Russo–Simis, Watanabe–de Bondt,
Fischer–Piontkowski, Akivis–Goldberg, Fassarella, Hartman–Nirenberg);
decisive structural evidence: Theorem G is FALSE for algebraic
non-polynomial functions, so it cannot live in the C^∞ or
algebraic-function literature. Nagaoka–Yazawa is published: J. Algebra
577 (2021) 175–202; attribution completed in main_result.tex by the
agent. Workflow resumed Aug 20 for the five lost agents (consequences
replays from cache).

Recovered and verified from the "lost" agents:

1. **A THIRD proof of Theorem G** (pzG, the second independent auditor):
   purely algebraic — embed C(x) into the Puiseux field at t = ∞ via the
   dominance of ∇e, transport the Legendre computation formally, and the
   B ≡ 0 hypothesis makes the formal Legendre value affine in t, whence
   the same no-positive-powers contradiction. No analytic function
   theory; valid over every algebraically closed characteristic-0 field.
   Machine pillars all PASS (`pzG_audit_identities.py`).
2. **Defect found: the n = 2 "exhaustive" scan wasn't.**
   `theorem_G_n2_elementary.py` dropped linear terms, but B depends on
   them (B(x1²/2+x2) = 1 ≠ B(x1²/2) = 0). Repaired by
   `pzG_n2_decision.py`: full-coefficient-space radical-membership
   decision, linear terms included, degrees 2–4. The hand proof's
   "reduced fibre" step also corrected to "smooth generic fibre ⟹
   parallel lines".
3. **Defect found: a false inference in the trichotomy proof.**
   main_result.tex derived e = φ(L) from det Hess₂e ≡ 0 — false as
   stated (h = x1 + x2²); the correct hypothesis is B₂(e) ≡ 0 (n = 2
   case of Corollary E, now stated and proved there). Found by the
   writeup agent (`wG_writeup_checks.py` check (6)); fixed.
4. `d5_pivotfree_normalform.py` had died on a trivial leftover-line
   TypeError, not on mathematics; repaired (one-line deletion), re-run.
5. gcv certificates: dBvdE family Wronskian identity (corrects the task
   premise), homogeneous converse of Theorem G with two-directional
   Identity-E ⟺ Nagaoka–Yazawa equivalence.

## Phase 19 — third independent audit; a fourth proof; two more overclaims retracted (2026-08-24/25)

Two further workflow runs died (a transient SSL failure, then the
session limit again) but the third auditor of Theorem G (agent atkG)
left its work on disk, and it is the sharpest pass yet:

1. **A FOURTH proof of Theorem G**, purely algebraic and shorter than
   all the others: with δ = D_V one has δg = g, δV = −V; the finite
   Taylor expansion of the polynomial g at the δ-constant point b = x + V
   splits g into δ-eigencomponents with distinct eigenvalues −k; δg = g
   forces every component to vanish. Polynomiality is used exactly once —
   as finiteness of the expansion. Certified (`atkG_algebraic_proof.py`).
   Theorem G now has four proofs, three produced by auditors instructed
   to refute it.
2. **Overclaim retracted**: `theorem_G_n2_deg45.py`'s "complete decision
   on a normalised slice" was on the HOMOGENEOUS top form only, where
   Theorem G is a one-liner by Identity E. I wrote that script and
   should have caught it. atkG supplied a genuine decision: Theorem G
   in n = 2 holds for every polynomial of degree ≤ 5, linear terms
   included (`atkG_n2_decision.py`, 23 s). All documents corrected.
3. Three-variable refutation attempts over complete degree-4 slices:
   the p211 slice and exact sampling PASS; the other slices' 17–25
   unknown Gröbner bases exceed the runtime budget — not refutations.
4. The chunked instance attack on Theorem B: R1/R2/R6 give 187
   constant-Hessian instances, all injective; R3 admits no consistent
   sparse system at all; R4/R5 hit the runtime kill. No counterexample.

## Outstanding at close (2026-08-18)

- The relaunched 6-agent workflow (`wf_043c2fa7-7d2`) is running as a
  further independent pass; its literature agent carries the still-open
  novelty check of Theorem G against the developable-hypersurfaces and
  degenerate-Monge–Ampère literature.
- Degree-5 HC₄: open. The obstruction is isolated in the weighted normal
  form of `d5_pivotfree_normalform.py`.
- Pivot existence in general: open. A pivot-free 4-variable
  constant-Hessian potential must have degree ≥ 5, ≥ 2 storeys, and
  dim span{N_α} ≥ 4.
- Corrected non-affine D1 (no automorphism yields a pivot while keeping
  det Hess constant): open; the Perazzo gap blocks the natural route.
- `rt_instances_x2_B_break.py` should be re-run in short chunks (it was
  killed at ~10 min having passed 13 solution branches).
- `conjE_search.py` / `conjE_degree3.py` are superseded — the conjecture
  they attacked is now Theorem G + Theorem F, and `theorem_G.py`'s
  three-variable radical test covers the same ground.

## Insight log

- AP4 contains all Meng doublings of planar Keller maps ⇒ closing all of
  AP4 would imply JC_2 statements. The affine-pivot class is exactly where
  JC_2 embeds into HC_4. Degree-4 AP4 corresponds to planar cubics — safe
  territory (Moh), so deg-4 closure is realistic; general-degree AP4 is
  JC_2-hard. This is the honest difficulty calibration for HC_4.
