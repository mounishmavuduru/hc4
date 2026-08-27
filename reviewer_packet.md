# Reviewer packet — HC_4 project (August 2026)

For expert human referees. Everything referenced lives in
`C:\Users\mouni\hc4\` (papers under `lit/`, exact certificates under
`certificates/`, run each with `py <script>`; all are fail-closed asserts,
exact rational arithmetic, sympy 1.14).

## 1. One-paragraph summary

After the July 2026 refutations (JC_n false for n ≥ 3; HC_n false for
n ≥ 5), the surviving open statements are JC_2 and HC_4, with
HC_4 ⇒ JC_2. This project (i) independently machine-verified all four
underlying 2026 counterexamples — for the Meng–Yang HC_5 polynomial
apparently the first third-party verification; (ii) proved a **pivot
theorem**: HC_4 holds for every potential that, in some affine
coordinates, is quadratic in a variable with constant leading coefficient
and pivot coefficient e1 of degree ≤ 2 — with a complete classification
in the affine-pivot case; (iii) proved **HC_4 for all f of degree ≤ 4**;
(iv) proved the Meng–Yang counterexample family admits **no pivot of
either kind in any affine coordinates** (their Remark 4.2, made rigorous
for the known family), so the Schur descent self-terminates on known
data; (v) proved Gao's two dimension-4 Jacobian counterexamples are
**not affinely equivalent to gradient maps** (linear symmetrizer test).

## 1a. Correction found by the red-team (recorded for transparency)

The first version of Theorem A carried the hypothesis deg e1 ≤ 2 in BOTH
branches. A recovered red-team script showed it is unnecessary when σ ≠ 0:
the pencil sweep is **per-point in x′**, so K = Hess₃e1 need only be
constant for fixed x′ (automatic), not globally constant. Re-derived and
re-certified in `theoremA_sharp.py`. The sharp statement is:
σ ≠ 0 ⟹ HC_4 for pivot coefficients of ANY degree; σ = 0 needs deg e1 ≤ 2.
**Corollary: an HC_4 counterexample admits no quadratic pivot in any affine
coordinates.** Referees should check this per-point argument carefully —
it is the one place where the first derivation was over-cautious, and the
same style of error could hide elsewhere.

## 1c. THE HEADLINE RESULT (Theorem G and the pivot dichotomy)

**Theorem G (rigidity; dimension-free).** For e ∈ C[x_1..x_n], if
B(e) := ∇eᵀadj(Hess e)∇e ≡ 0 then det Hess e ≡ 0. Geometrically: *if every
level hypersurface of a polynomial is developable, its Hessian determinant
vanishes identically.* The converse is false.

Proof idea: where Hess e is invertible, V := (Hess e)⁻¹∇e satisfies
Dφ·V = φ for φ = ∇e, so V is φ-related to the Euler field, and
D_V e = B/det Hess e. So B ≡ 0 makes e a first integral of V. Passing to
the local Legendre transform L (∇L = φ⁻¹, EL − L = e∘φ⁻¹) this becomes
(E² − E)L = 0, so L has only Euler-eigenvalues 0 and 1 and
φ⁻¹(tp) = t⁻¹a + b. The continued identity ∇e(t⁻¹a + b) = tp compares
coefficients of t¹: the left side is a polynomial in t⁻¹, so p = 0 —
contradiction.

**Corollary (pivot planarity).** Every 3-variable polynomial with vanishing
bordered Hessian is affinely 2-variable. (Theorem G + Theorem F.) This was
Conjecture E; it is now proved.

**Theorem (pivot trichotomy).** A 4-variable potential with constant nonzero
Hessian determinant admitting a pivot v (γ = D_v²f, e1 = D_vf) falls into
exactly one of three cases:
(i) γ ≠ 0 (quadratic pivot) — settled outright (Theorem A);
(ii) γ = 0 and deg e1 = 1 — settled outright (AP0);
(iii) γ = 0 and deg e1 ≥ 2 — after an affine change,
a(x1,x2) + x3b(x1,x2) + x4e(x1,x2): a Meng doubling of the planar Keller
map (b,e) plus an inert planar potential, injective iff (b,e) is.
**Hence HC_4 restricted to potentials admitting a pivot is EQUIVALENT to
JC_2** (case (iii) is exactly JC_2; (i) and (ii) are unconditional).
The split comes from the identity
[x4] det Hess f = −(e0)_{x3x3}·B₂(e). **Two real corrections along the
way, both recorded**: case (ii) was omitted in the first draft; and a
later draft claimed the affine-pivot class "is exactly the doublings" —
false: the witness w = x1²/2 + x1x2(1+x3) + x2²(2x3+x3²)/2 + x3x4
(det Hess = 1) has an affine pivot but is not affinely equivalent to any
doubling (its isotropic cone is a single line; a doubling's contains a
2-plane). Certs: `pivot_dichotomy.py`, `paperG_writeup_checks.py` (W13).

WHAT THIS LEAVES: HC_4 is still open, and the residual class is now named
exactly — potentials admitting **no pivot at all**. (In dimension 5 such
potentials exist: D1 shows the Meng–Yang counterexample has empty pivot
cone. Whether they exist in dimension 4 with constant nonzero Hessian
determinant is the precise remaining question; if they do not, then
HC_4 ⟺ JC_2.)

**Referee focus (updated Aug 18).** The original write-up's load-bearing
step — "φ⁻¹ continues analytically along the punctured ray" — was judged
a NON-SEQUITUR by the independent audit (`advG_flow_proof.py`): φ⁻¹ is
only a local inverse and det Hess e may vanish on the ray, so no inverse
need exist there. TWO independent repairs now exist, and both avoid
continuation entirely:
1. (paper, `theorem_G_paper.tex` §2) split L = α + β into its
   Euler-eigencomponents near p0 (available exactly because B ≡ 0) and
   extend each over the punctured cone Γ = C^×·W by its OWN homogeneity,
   via an explicit formula; ∇e∘Ψ = id then propagates by the identity
   theorem, and comparing Laurent coefficients of t¹ kills p0.
2. (advG) purely local flow argument: B ≡ 0 gives MV = 2g hence
   Jac(V)V = −V, so the flow of V is x(s) = ae^{−s} + b with
   ∇e(x(s)) = e^s·p0; then τ·∇e(aτ + b) − p0 is a polynomial vanishing at
   infinitely many τ but not at τ = 0. Valid over any characteristic-0
   field (formal power series).
A THIRD proof was produced by the second independent audit (pzG, Aug 20):
purely algebraic, via an embedding of C(x) into the Puiseux field at
t = ∞ — no analytic function theory at all, valid over every
algebraically closed field of characteristic 0
(`pzG_audit_identities.py`, header; machine pillars all PASS).
Sharpness is certified: the non-polynomial algebraic function
x1 − sqrt(x1²−2x2) satisfies every step except the final polynomiality
comparison (`tgc_geometry_and_sharpness.py`). The n = 2 case is proved
independently and elementarily; CAUTION: the first "exhaustive degrees
≤ 5" scan in `theorem_G_n2_elementary.py` dropped linear terms (on which
B depends) and was not exhaustive as claimed — the genuinely full
decision (linear terms included, degrees ≤ 4, radical membership over
the whole coefficient ring) is `pzG_n2_decision.py`, plus the complete
slice decision in `theorem_G_n2_deg45.py` — CAVEAT: that slice is
homogeneous (top form only), where the theorem is trivial by Identity E;
the genuinely complete degree ≤ 5 decision over ALL polynomials, linear
terms included, is `atkG_n2_decision.py` (third independent auditor).
Certificates: `theorem_G.py`, `theorem_G_n2_elementary.py`,
`theorem_G_n2_deg45.py`, `pzG_n2_decision.py`, `advG_flow_proof.py`,
`pzG_audit_identities.py`, `paperG_writeup_checks.py`.

## 1b. Structure theorems (added after the self red-team pass)

- **Theorem C (doubling structure).** f = a(x1,x2) + x3·b + x4·e has
  Hessian [[A,B],[Bᵀ,0]] with B = ∂(b,e)/∂(x1,x2), so
  **det Hess f = (Jac(b,e))², independent of a**. This class *is* the Meng
  doubling class, and **HC_4 restricted to it is equivalent to JC_2**.
  Theorem A's last case is therefore exactly: planar Keller maps with one
  component (affinely) x1²/2 + x2 are b = −κx1 + β(x2+x1²/2), invertible.
- **Theorem D (developability).** det[[K,g],[gᵀ,0]] = −gᵀadj(K)g, so the
  affine-pivot constraint (c2) says the pivot coefficient has identically
  vanishing bordered Hessian: its level surfaces are developable.
  Nondegenerate quadric pivots are impossible.
- **Identity E.** For homogeneous e of degree d ≥ 2 (any n):
  ∇eᵀ adj(Hess e) ∇e = (d/(d−1))·e·det Hess e (Euler + H·adjH·H = detH·H).
- **Theorem D′ (cone leading form).** Hence the leading form of any pivot
  coefficient has vanishing Hessian determinant and is, by Gordan–Noether
  (n = 3), a cone. For deg e1 = 2 this reproduces the rank collapse used in
  Theorem A — a mutual consistency check.
- **Conjecture E.** Every 3-variable polynomial with vanishing bordered
  Hessian is affinely 2-variable. Proved for deg ≤ 2 and for all
  homogeneous pivot coefficients. If true in general, the whole
  affine-pivot class collapses to doublings, i.e. **HC_4 on affine-pivot
  potentials is exactly JC_2**.

## 2. Statement dependency graph (bottom-up)

- External, source-verified: Gordan–Noether (n ≤ 4); Dillen 1991 (HC_2);
  de Bondt 2015 (HC_3, polynomial inverse); Wang/Oda (deg ≤ 2);
  BBR/Ax/Cynk–Rusek (injective ⇒ automorphism); Meng 2006 (bridges);
  Meng–Yang 2026 + Alpöge/Gao 2026 (landscape).
- Project-internal, machine-certificated:
  T0–T4 (graded tower + midpoint identity) → T5 (branch structure)
  → Q2 (quadratic pivot, uses HC_3) + AP0 (affine pivot, affine e1, uses
  HC_2 fiberwise) + AP1–AP3 (rank collapse + classification, uses nothing
  external) = **Theorem A**;
  S1 (2×2 lemma) + S2 (branch graded formulas, fully generic) + S3
  (isotropic-vector lemma) close the three branches = **Theorem B**;
  D1/D1+ (Gröbner pure-power certificate) and G1 (symmetrizer space)
  are self-contained.

## 2b. PRIOR ART — read `prior_art.md` before assessing novelty

An adversarial audit found that six items we had presented as our own are
in the literature. They are retracted and cited: Identity E (Fox 2017; Del
Pia et al.; traced to Hemmer 1995), Identity H (Dolgachev, *Classical
Algebraic Geometry*, (1.20)–(1.21) — we had wrongly written "derived and
used here"), Theorem D (Reilly 1986; Fox names our invariant U(F)), Lemma
S3 (Gow, Thm 1, in greater generality; Fillmore–Laurie–Radjavi 1985;
= Bertini for linear systems of quadrics), Lemma S1 (Meshulam /
Loewy–Radwan), and Theorem F (a corollary of de Bondt–van den Essen,
*Singular Hessians*, J. Algebra 282 (2004), Thm 3.3, which we had not
cited). Theorem C is folklore given Meng's doubling.

**Theorem G survives the audit.** Reilly and Fox give exactly our framework
and the identity B = det Hess e·⟨∇e,(Hess e)⁻¹∇e⟩, but both assume that
quantity is nonvanishing (their subject is affine spheres). Theorem G is
the complementary case, vacuous over R for a definite Hessian and with bite
only over C. No prior art located — but not claimed as new pending review.

## 3. Novelty claims (for referees to check)

1. **Theorem A** (pivot theorem, incl. the classification
   f ≅ a(x1,x2) + x3(−κx1 + β(x2+x1²/2)) + x4(x2+x1²/2), det Hess = κ²):
   we found no prior statement. The σ ≠ 0 half (Q2) is a natural
   completing-the-square + HC_3 bootstrap that experts may know
   informally; the σ = 0 rank collapse and PDE classification we believe
   new. CHECK against: de Bondt's JPAA 2015 paper (his weight machinery is
   different and provably fails at n = 4); unpublished folklore.
2. **Theorem B** (HC_4, deg ≤ 4): no prior claim found; the natural prior
   art would be symmetric-cubic-Keller classifications (Hubbers-type,
   dBvdE nilpotent-symmetric n ≤ 4) — those concern the NILPOTENT
   x + ∇f setting, not Meng's HC_4; we found no overlap.
3. **Lemma S3** (all-singular subspaces of Sym_3(C) have a common
   isotropic vector): possibly classical (spaces of singular symmetric
   matrices literature: Meshulam, Loewy, bounded-rank subspaces). Our
   3-line proof is included; novelty immaterial to correctness.
4. **D1/D1+, G1**: new (the objects are three weeks old).

## 4. Likely weak points (be hostile here)

1. **Case exhaustiveness of Theorem B**: the essential-rank trichotomy of
   f4 and the three rotations-to-pivot arguments. The r = 2 case needs
   Hess_{(x3,x4)}f3 = ℓ·uuᵀ (S1) — check the degenerate subcases
   (u = 0, ℓ = 0) which we handle by "any v works".
2. **Q2's sweep argument**: det Hess f = σ·det₃(pencil at u(x', x4)) as a
   polynomial identity, then substituting arbitrary s for u — this uses
   that (x', x4) ↦ (x', u) is a polynomial automorphism of C⁴ (triangular),
   so the identity in (x', x4) is equivalent to the identity in (x', s).
   Verify you agree.
3. **AP0's fiberwise HC_2**: for each fixed a ∈ C, e0(·,·,a) has constant
   Hessian determinant −c; Dillen's theorem is applied over C for each a.
   No uniformity in a is needed since injectivity is tested pointwise.
   Verify you agree.
4. **AP1 collapse**: the step "adj(P + tE11) linear in t" is a 3×3 fact
   (no 2×2 minor uses the (1,1) entry twice) — machine-verified
   symbolically. The PDE completeness argument (characteristic change of
   variables ξ = x2 + x1²/2, polynomial in (x1, ξ)) — check the claim
   H ∈ C[ξ] (H polynomial in (x1, ξ) with ∂H/∂x1 ≡ 0).
5. **Scope of D1**: "no pivot in any affine coordinates" — non-affine
   automorphisms are excluded from the claim (they do not preserve the
   constant-Hessian class in general); a referee should confirm the
   framing is honest (we state it explicitly). UPDATE (Aug 18): the
   framing is now PROVED necessary — `nax_d1_nonaffine.py` exhibits an
   explicit non-affine automorphism τ making every Meng–Yang member
   affine-linear in a coordinate (so the naive strengthening is false),
   while τ destroys det Hess constancy (= λR², R non-constant). The
   corrected strengthening (no automorphism yields a pivot while
   PRESERVING constant Hessian) remains open; the natural
   Gordan–Noether route breaks at Perazzo (`nax_pivot_transfer_5var.py`).
6. **Scope of G1**: two-sided AFFINE equivalence only. Composition with
   nonlinear automorphisms is open and stated as such.
7. **Transcription risk for F4/F5**: mitigated by re-deriving Gao's exact
   Jacobian constants from the parsed polynomials (a 1-in-astronomically-
   many coincidence if mis-parsed); still, referees can re-parse.

## 5. Exact checks a referee can run (minutes each)

    py certificates\verify_alpoge_3d.py        # JC_3 counterexample
    py certificates\verify_mengyang_fast.py    # HC_5/HC_6 + family
    py certificates\pivot_obstruction_d1.py    # D1/D1+
    py certificates\gradient_equivalence_test.py  # G1 + Gao determinants
    py src\hc4lib.py                           # tower + midpoint proofs
    py certificates\verify_q2_core.py          # Q2 chain
    py certificates\ap4_rank_reduction.py      # AP3/AP2/AP1 identities
    py certificates\ap4_pivot_closure.py       # AP1 collapse + closure
    py certificates\deg4_branch_closures.py    # S1,S2,S3 + branch formulas
    py certificates\t3_generic.py              # T3 fully generic
    py certificates\dillen4_ap4_control.py     # known-good control

## 5a. Status of the machine verification (be precise about this)

A mass process kill (exit 255 on five concurrent jobs) truncated several
runs. What is fully verified versus partial:

| Certificate | Status |
|---|---|
| `theorem_G.py` | COMPLETE — incl. the 3-variable radical-membership test |
| `pivot_dichotomy.py` | COMPLETE — both branches, fully symbolic |
| `theoremA_sharp.py`, `doubling_structure.py`, `cone_leading_form.py`, `theorem_F.py`, `homogenization_identity.py`, `deg4_branch_closures.py`, `t3_generic.py`, `ap4_*.py`, `verify_q2_core.py`, `dillen4_ap4_control.py` | COMPLETE |
| `verify_alpoge_3d.py`, `verify_mengyang_fast.py`, `verify_meng_yang_hc5.py`, `pivot_obstruction_d1.py`, `gradient_equivalence_test.py` | COMPLETE |
| `theorem_G_n2_elementary.py` | degrees 2–3 complete (exhaustive); killed in degree 4 |
| `theorem_G_n2_deg45.py` | 180 exact instances degrees 4–6 PASS; its "complete slice decision" is on the HOMOGENEOUS slice only (overclaimed; found by atkG) — superseded by `atkG_n2_decision.py` |
| `identityE_is_known.py` | 4 of 6 cases before the kill — sufficient (formal equivalence) |
| `rt_instances_x2_B_break.py` | re-run per family (chunked CLI, Aug 24–25): R1 10 instances, R2 6 (16 draws), R6 171 instances — every constant-Hessian solution injective with nonempty pivot cone; R3 admits NO consistent system in 20 draws (the r = 3 ansatz is over-constrained; attack vacuous there); R4 passes 39 instances then a single branch's injectivity Gröbner exceeds the 10-minute kill; R5's first system never finishes. No counterexample, no gap. |
| `conjE_search.py`, `conjE_degree3.py` | superseded by `theorem_G.py` |
| `advG_flow_proof.py` | COMPLETE (Aug 18) — independent audit of Theorem G: found the continuation gap, certified a replacement elementary proof |
| `pzG_audit_identities.py` | COMPLETE (Aug 20) — second independent audit; third proof (Puiseux at infinity); all machine pillars PASS |
| `pzG_n2_decision.py` | degrees 2–3 COMPLETE over the full coefficient space (linear terms included); degree 4's 14-variable Gröbner does not finish under the 10-minute kill — superseded by `atkG_n2_decision.py` |
| `atkG_algebraic_proof.py` | COMPLETE (Aug 25) — third independent audit; FOURTH proof (derivation eigenvectors, purely algebraic); pillars fully generic, ALL PASS |
| `atkG_n2_decision.py` | COMPLETE (Aug 25) — Theorem G in n = 2 DECIDED for all polynomials of degree ≤ 5, linear terms included, via certified normalisation; 23 s |
| `atkG_n3_search.py` | mode deg4:p211 ALL PASS; modes deg3 (13-unknown complete decision — the same computation already PASSED in `theorem_G.py` on Aug 10) and the other deg-4 slices: see §5a addendum below as runs complete |
| `wG_writeup_checks.py`, `refG_paper_integrity.py` | COMPLETE — the rewritten `theorem_G_paper.tex` is referentially sound and its displayed identities re-certified |
| `gcv_family_wronskian.py`, `gcv_ny_homogeneous_converse.py` | COMPLETE (Aug 20) — dBvdE family Wronskian identity; homogeneous converse; NY ⟺ Identity E both directions |
| `paperG_writeup_checks.py` | COMPLETE (Aug 18) — W1–W15, the full checklist of `theorem_G_paper.tex` |
| `paperG_structure.py` | RED at (P4) BY DESIGN — documents an ideal-vs-radical membership defect; superseded by (W13) |
| `nax_d1_nonaffine.py`, `nax_pivot_transfer_5var.py` | COMPLETE (Aug 18) — non-affine D1: naive strengthening refuted, corrected version open |
| `d5_graded_tower.py`, `d5_pivotfree_normalform.py` | graded pieces C12–C9 closed-form (COMPLETE); pivot-free weighted normal form (re-run in progress) |
| `pe_pivot_cone.py`, `pe_lowdeg_pivot.py`, `pe_perturbation.py` | COMPLETE (Aug 18) — pivot-existence tool, PE-1/2/3/4 |
| `tgc_debondt_chain.py`, `tgc_geometry_and_sharpness.py`, `tgc_beyond_hc4.py` | COMPLETE (Aug 18) — independent 2nd proof of Corollary E; sharpness; placement (R1–R6) |

## 5b. Later additions (structure theorems)

    py certificates\doubling_structure.py      # Theorem C
    py certificates\developability.py          # Theorem D
    py certificates\cone_leading_form.py       # Identity E + Theorem D'
    py certificates\homogenization_identity.py # Identity H
    py certificates\theorem_F.py               # Theorem F
    py certificates\conjE_search.py            # Conjecture E: counterexample search
    py certificates\conjE_degree3.py           # Conjecture E, degree 3 (SLOW)

**Theorem F** (the strongest partial result on Conjecture E): if
e ∈ C[x1,x2,x3] has vanishing bordered Hessian AND vanishing Hessian
determinant, it is affinely 2-variable. Hence Conjecture E holds for all
homogeneous pivot coefficients. The residual case is exactly
B(e) ≡ 0 with det Hess₃(e) ≢ 0; there Identity H forces the top two
graded pieces of e·det Hess₃(e) to vanish, so det Hess₃(e) sits two
degrees below generic — the natural lever for finishing the conjecture.

## 6. Questions for specialist review

1. Is Lemma S3 in the literature (singular symmetric subspaces)? If so,
   where — and does the standard version cover dim W up to 5?
2. Does de Bondt's unpublished work or the Nijmegen school folklore
   contain Theorem A's σ ≠ 0 half?
3. Conjecture C2 (descent rigidity, main_result.tex §5): is the
   Kronecker-pencil route over C(w) viable, or is there a slicker
   argument that admissible 5-variable data forces doubling structure?
4. The degree-5 frontier: the tower's next layer meets pivot coefficients
   of degree 3 (containing doubled planar quartic Keller maps — Moh
   territory). Is a fiberwise-JC_2-bounded-degree version of AP0
   formalizable to push Theorem B to degree 5?

## 6b. Red-team results actually obtained

The workflow runs died on quota, but the agents' scripts were recovered from
disk and executed. Verdicts:

| Script | Result |
|---|---|
| `rt_ref_landscape.py` | **PASS** — independent re-derivation of the whole 2026 landscape from the primary TeX only (not from this project's scripts): Jac F ≡ −2, degrees 7/6/4, triple collision; det Hess φ ≡ −4; deg Ψ = 14, 42 monomials, det Hess Ψ ≡ 128, collision; family det ≡ 4λ. |
| `deg4_q2v_bordered_schur.py` | **PASS** — bordered form and Schur identity; explicitly notes Q2(i) needs no constancy of K (independently anticipating the sharpening). |
| `deg4_q2v_collision_transfer.py` | **PASS** — collision transfer, fully symbolic. |
| `deg4_q2v_pencil_forcing.py` | **PASS** — proves the s³-coefficient of det₃(M+sK) is det₃K, so the pencil forces rank K ≤ 2; rank 0/1/2 instances exhaust the cases. |
| `deg4_q2v_ap_identities.py` | **PASS** — adj decomposition (18 symbols), bordered-determinant identity, (AP) ⟺ (c0)&(c1)&(c2), all fully symbolic. |
| `rt_obstructions_parser.py` | **substantively PASS** — an independent strict re-parse of Gao's F₄/F₅ agrees with the production parser on every component; source term counts match parsed monomial counts exactly (no dropped/merged terms); degrees match; det J ≡ −44/9 and 160/29 by the independent parse and at random rational points. Its 2 reported "FAILURES" are mislabeled checks (they assert only `\tfrac` occurs, which is the desired outcome). |
| `rt_instances_B_ansatz.py` | **PASS** — 6 ansatz families, 15 solution branches, 45 concrete constant-Hessian quartics (some over Q(i)), all injective. Theorem B survives. |
| `rt_instances_A_classified.py`, `rt_instances_x2_A_classified.py` | **PASS** — 24 instances of the classified family, degrees 3–6 (beyond Theorem B's range), all det ≡ κ² and all injective; negative controls confirm the tests are not vacuous (a planted collision was correctly detected). |
| `rt_instances_x2_AP0.py` | **PASS** — 12 AP0 instances, degrees 4–8, all pass det and injectivity; negative control gives non-constant det as required. |
| `rt_novelty_q2_anydegree.py` | **FOUND A REAL DEFECT** — the deg e1 ≤ 2 hypothesis is unnecessary when σ ≠ 0. See §1a. |
| `rt_instances_AP0.py` | **assertion failed — resolved, no defect.** Its "negative control" x₁x₂ + x₂²x₃ + x₃x₄ has det Hess ≡ 1, i.e. it is a valid AP0 instance rather than a control; Theorem A predicts injectivity and a saturated Gröbner confirms it (`rt_ap0_negative_control_resolved.py`). |
| `rt_instances_x2_B_break.py` | crashed on a leftover placeholder line in the agent's own code (`Poly.homogeneous_order()` returns None for non-homogeneous input); repaired and re-run — this is the strongest attack, designed to empty the pivot cone. |

## 7. Candid limitations

- HC_4 itself remains OPEN; nothing here settles it in general degree.
- Theorem B's proof relies on machine-verified polynomial identities
  (each labelled proof vs consistency-check in lemma_ledger.md);
  a referee should audit the labels.
- **Independence status of the adversarial pass (updated Aug 18).** Four
  runs of a six-agent verification workflow were killed by quota limits,
  but the agents wrote their work to disk before dying; all of it was
  recovered and re-run to completion. The audit of Theorem G
  (`advG_flow_proof.py`) was performed by an agent instructed to refute,
  working from the statement and the write-up, and it (a) judged the
  original analytic-continuation paragraph a non-sequitur and (b)
  certified a complete replacement proof. A separate write-up agent
  independently repaired the same gap a second way (`theorem_G_paper.tex`,
  Remark "where the care is needed") and found a false structural claim
  of the project lead's (the doubling exhaustion — see §1a and
  Remark rem:notdoubling of `main_result.tex`). A fifth workflow run is
  in flight as a further pass. This is materially more independent than
  the Aug 10 state, but it is still all machine-side: no human expert
  has reviewed anything. Treat §4 accordingly.
- The novelty search for Theorem G is still not complete: the equiaffine
  sources (Reilly, Fox) were checked and contain only the pointwise
  dictionary and the opposite-direction examples, but the
  developable-hypersurfaces and degenerate-Monge–Ampère literature has
  not been swept. The paper flags this explicitly for referees. Theorem
  A's σ ≠ 0 half and Lemma S3 remain the items most likely to be known.
- Former Conjecture E is now a theorem (Theorem G + Theorem F, plus an
  independent second proof via de Bondt–van den Essen in
  `tgc_debondt_chain.py`); it is SHARP — false in every dimension ≥ 4,
  by dehomogenized Perazzo cubics (`tgc_beyond_hc4.py`).
