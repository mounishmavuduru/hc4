# Lemma ledger — HC_4 project

Status labels: **[known]** (literature, source-verified), **[proved here]**
(full proof in this project), **[machine]** (exact machine certificate in
`certificates/`), **[conjectured]**, **[open]**.
Every entry: statement / dependencies / proof or source / adversarial check.

Conventions. k a field of char 0 (results over C unless stated); for
f ∈ k[x_1..x_n], Hess f = (∂²f/∂x_i∂x_j); f_d = degree-d homogeneous part.
HC_n: det Hess f ∈ k^× ⟹ the formal Legendre transform f^L is a polynomial
(equivalently, by L-B below, ∇f has a polynomial inverse).

---

## B0. Status inputs (verified counterexamples)

**B0.1 [known] [machine]** Alpöge's map F (components deg 7,6,4) has
Jac F ≡ −2 and F(0,0,−1/4) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−1/4,0,0).
Hence JC_3 false; with stabilization (L-S1) JC_n false for all n ≥ 3.
Cert: `verify_alpoge_3d.py` (Berkowitz), re-run in `gradient_equivalence_test.py`
(domain-GE). Source: Alpöge announcement 2026-07-19 as recorded in
arXiv:2607.22198 and arXiv:2608.00222.
Adversarial: points distinct ✓, char-0 only ✓ (integer/rational data).

**B0.2 [known] [machine]** Meng–Yang Ψ = A² + 13A + 2B ∈ Z[x1,x2,y1,y2,y3]:
deg 14, 42 monomials, det Hess Ψ ≡ 128, ∇Ψ(±1,∓3/2,0,0,0) = (0,0,−1/2,0,0).
Hence HC_5 false; with L-S2, HC_n false for all n ≥ 5.
Cert: `verify_meng_yang_hc5.py`, `verify_mengyang_fast.py` (two det algorithms).
Also verified: the 6-var doubling φ has det Hess ≡ −4 with a 3-point gradient
collision (HC_6 false directly), and the family identity
det Hess(B + (λ/2)A² + μA) ≡ 4λ for symbolic λ, μ.
Adversarial: "det Hess ∈ k^×" ✓ (128 ≠ 0); collision points distinct ✓.

---

## L. Core lemmas

**L-B (bridges) [known, proofs re-verified here]**
(i) JC_n ⟹ HC_n. (ii) HC_2n ⟹ JC_n via the doubling
φ(x,y) = ⟨y, F(x)⟩, det Hess φ = (−1)^n (Jac F)², and triangularity of ∇φ.
(iii) If det Hess f ∈ k^× and ∇f not injective, f^L is not a polynomial.
Source: Meng 2006 (statement; provenance check in progress), full proofs in
arXiv:2607.22198 §2, each step re-verified by hand here (block-swap sign
(−1)^{n²} = (−1)^n; lower-block-triangular determinant; formal-identity
globalization). Dependencies: none.
Adversarial: (ii) needs the y-block linear-invertibility from Jac F ∈ k^× ✓;
(iii) needs char 0 for formal inverse ✓.

**L-S1, L-S2 (stabilization) [known, re-verified]** JC_n false ⟹ JC_{n+1}
false (pad by identity); HC_n false ⟹ HC_{n+1} false (add ½t²; Legendre
transform separates). Source: arXiv:2607.22198 §3; proofs re-checked.
Adversarial: HC padding needs the +½t² block to keep det Hess in k^× ✓ and
separated Legendre transform ✓.

**L-INJ [known — source-verified]** An injective polynomial map C^n → C^n is
bijective with polynomial inverse. Sources (primary, fetched):
Białynicki-Birula–Rosenlicht, Proc. AMS 13 (1962) 200–203 (surjectivity);
Cynk–Rusek, Ann. Polon. Math. 56 (1991), Thm 2.2 (injective ⟺ automorphism,
alg. closed char 0); for Keller maps also BCW Bull. AMS 7 (1982), Thm 2.1
(char 0: injective ⟺ invertible, set-theoretic inverse automatically
polynomial). (So: HC_n over C ⟺ every gradient Keller map is injective; a
counterexample needs only two exact points with equal gradient.)

**L-HC3 [known — source-verified twice]** For any char-0 field K and
f ∈ K[x_1..x_n] with n ≤ 3: ∇f satisfies the Jacobian conjecture; when
det Hess f ∈ K^× the inverse is explicitly polynomial (Lemma "antitri":
x_i ∈ A[F_{n+1−i},…,F_n]). de Bondt, JPAA 219 (2015) 3743–3754 =
arXiv:1203.6605v3, Theorem 2.1; full TeX read by two independent agents;
verbatim quotes in `lit/debondt_extracts.md`. n = 2 case originally Dillen,
JPAA 71 (1991) 13–18. NOTE for n = 4: the same paper shows its own two main
tools FAIL at n ≥ 4 — the anti-triangular normal form (explicit
counterexample f = (x1+x2²)x3 + (x2+(x1+x2²)²)x4, an invertible doubled
planar Keller map, i.e. an AP4-class polynomial in our terminology) and the
weight theorem (Example 2.3, robust failure) — and explicitly notes
HC_4 ⟹ planar JC. Used for: Q2's base case; consistency control in G1.

**L-WANG [known — source-verified]** Every Keller map of degree ≤ 2 (any n,
char ≠ 2, not necessarily homogeneous) is invertible. S. S.-S. Wang,
J. Algebra 65 (1980) 453–494; statement + Oda's midpoint proof verified from
BCW, Bull. AMS 7 (1982), Thm 2.4 and p. 298 (`lit/bcw/bcw1982.pdf`). Our T4
midpoint identity is the gradient-case sharpening of exactly Oda's argument.

---

## T. Graded tower for det Hess f ∈ k^× (n = 4 unless stated)

**T0 [proved here]** For any n and f with det Hess f = c: the degree-0 graded
piece of det Hess f is det Hess f_2. Hence **the quadratic part of f is
nondegenerate**, det Hess f_2 = c; in particular f_2 ≠ 0 and, after a linear
change and dropping affine terms, f = q + f_3 + ... + f_d with
q = ½Σx_i², c = 1.
Proof: every entry of Hess f_k is homogeneous of degree k−2 ≥ 1 for k ≥ 3;
a permutation product has degree 0 iff all four entries come from Hess f_2. ∎
Adversarial: uses only char 0 and the grading; normalization of q to ½Σx²
uses alg. closedness (or any field where the quadric is split) — over
non-closed k keep S general.

**T1 [proved here]** If d = deg f ≥ 3 then det Hess(f_d) ≡ 0 (any n with
n(d−2) > 0).
Proof: the degree-n(d−2) piece of det Hess f is det Hess f_d (any product
using an entry of some Hess f_k, k < d, has degree < n(d−2)); a constant has
no positive-degree part. ∎

**T2 [known — source-verified] (Gordan–Noether, n ≤ 4)** A form in ≤ 4
variables over an algebraically closed char-0 field with identically
vanishing Hessian determinant is a cone: after a LINEAR change it involves
at most n−1 variables. Verified verbatim from Watanabe–de Bondt,
arXiv:1703.07624 (improved version; Springer PMS 319, 2020) — use this
version, the 2014 proceedings version has a known error in its Lemma 5.2.
False for n ≥ 5: Perazzo cubic x1²x3 + x1x2x4 + x2²x5 (canonical form,
WdB Remark); full n = 5 classification: after linear change, homogeneous
elements of K[x1,x2][Δ], Δ = p3(x1,x2)x3 + p4(x1,x2)x4 + p5(x1,x2)x5
(WdB Thm 4.5). Binary Hesse case: f = (a1x1 + a2x2)^d. Extracts:
`lit/gordan_noether_extracts.md`. Consequence with T1: **the leading form of
any HC_4 candidate of degree ≥ 3 is, WLOG, in C[x1,x2,x3]**.
Adversarial: char 0 essential (x1^p x2 in char p); alg. closedness assumed
in the sources, and the cone conclusion descends to non-closed char-0 fields
via linear dependence of partials (remark, not quoted).

**T3 [proved here]** With f_d ∈ C[x1,x2,x3] (n = 4, d ≥ 3): the degree-(4d−9)
piece of det Hess f equals det₃(Hess₃ f_d) · ∂²f_{d−1}/∂x4², which therefore
vanishes: **either det₃ Hess₃ f_d ≡ 0 (then, by GN₃, WLOG f_d ∈ C[x1,x2]) or
f_{d−1} has x4-degree ≤ 1.**
Proof: in a permutation product of degree 3(d−2)+(d−3), exactly three entries
come from Hess f_d; since row 4 and column 4 of Hess f_d vanish, the fourth
entry must occupy position (4,4); summing over the 3-permutations gives the
3×3 determinant as stated. ∎ (Machine spot-check planned in hc4lib tests.)

**T4 [proved here] (midpoint identity, deg ≤ 4)** For deg f ≤ 4 and any
m, v: ∇f(m+v) − ∇f(m−v) = 2[Hess f(m)·v + ∇f_4(v)].
Hence ∇f is non-injective iff ∃ v ≠ 0, m with Hess f(m)v = −∇f_4(v).
In particular (Wang's case deg f ≤ 3, f_4 = 0): Hess f(m)v = 0 with
det Hess f(m) = c ≠ 0 forces v = 0 — **HC_n holds for deg f ≤ 3, all n**.
Proof: gradient components have degree ≤ 3; odd part of a cubic map at m ± v
is JG(m)v + (cubic-in-v term), and the pure v³-term of ∇f(m+v) is ∇f_4(v)
(compare weights; machine-verified with fully symbolic coefficients —
planned in hc4lib test, linear in coefficients hence a proof for all f). ∎
Adversarial: the identity is specific to deg ≤ 4 ✓; injectivity ⟹ HC uses
L-INJ ✓.

---

## D. Descent obstruction (Track C — new results)

**D0 [proved here]** The coordinate changes of C^N preserving the class
{f : det Hess f ∈ C^×} include exactly the affine ones among polynomial
automorphisms acting by f ↦ f∘σ in general position
(Hess(f∘σ) = Jσ^T (Hess f)∘σ Jσ + Σ_k (∂_k f)∘σ · Hess σ_k; for affine σ the
second term vanishes and det multiplies by (det Jσ)²).
Scope note: this does NOT claim no exotic non-affine σ can preserve constancy
for a specific f; D1 below is a statement about the affine orbit, which is
the natural equivalence of the Hessian setting.

**D1 [proved here] [machine]** For every λ ≠ 0 and every μ, the Meng–Yang
potential ψ_{λ,μ} = B + (λ/2)A² + μA admits no v ∈ C^5 \ {0} with
D_v²ψ_{λ,μ} ≡ 0. Consequently no affine change of coordinates makes any
member of the known HC_5 counterexample family affine-linear in a variable,
so **the Schur descent of arXiv:2607.22198 cannot be applied to its own
output: the known 5-variable counterexample family admits no pivot, in any
affine coordinates**. (Rigorous form, for the known family, of their
heuristic Remark 4.2.)
Proof: D_v²ψ = D_v²B + μD_v²A + λ(A·D_v²A + (D_vA)²) with
deg_w D_v²B, deg_w D_v²A ≤ 5; so every w-monomial coefficient of degree ≥ 6
equals λ·[same coefficient of C(v) := A·D_v²A + (D_vA)²] — μ-free. The 95
such coefficients are homogeneous quadratics in v; a Gröbner basis has
leading-term ideal containing a pure power of each v_i, so their common zero
cone is {0}. Cert: `pivot_obstruction_d1.py`.
Adversarial: pure-power-LT criterion proves the affine variety is finite,
and a finite cone is {0} ✓; homogeneity of all generators asserted ✓;
λ ≠ 0 used exactly once ✓; μ eliminated by degree ≥ 6 restriction ✓;
translations don't change D_v² ✓.

**D1+ [proved here] [machine] (strengthening of D1)** The same certificate
proves more: the degree-≥6 coefficient subsystem is implied not only by
D_v²ψ ≡ 0 but by D_v²ψ ≡ γ for ANY constant γ (a constant has no
positive-degree part). Hence **no member of the Meng–Yang family admits a
quadratic pivot (D_v²ψ ≡ const ≠ 0) either**: the family is immune both to
Schur descent (affine pivot) and to the quadratic-pivot reduction Q2 below.
No new computation needed; pure logic over the same 95-generator system.

**Q2-SHARP [proved here] [machine] — SUPERSEDES the degree hypothesis below**
The hypothesis deg e1 ≤ 2 is NOT needed when σ ≠ 0. The pencil sweep is
PER-POINT in x′: for fixed x′, K(x′) is a constant matrix and
u = x4 + e1(x′)/σ still sweeps C, so det₃(Hess₃ẽ0(x′) + s·K(x′)) ≡ c/σ for
all s ∈ C and all x′ — hence for each SCALAR λ the 3-variable polynomial
ẽ0 + λe1 has constant nonzero Hessian determinant, and HC_3 applies. The
Schur identity Hess₃(ẽ0) + uK = Hess₃(e0 + x4e1) − ggᵀ/σ and the
collision-transfer identity both hold for generic e1 of degree 4 (verified
fully symbolically), so constancy of K is never used.
**Therefore: HC_4 holds for every f with a quadratic pivot (σ ≠ 0), in any
degree of e1.** Corollary: an HC_4 counterexample admits NO quadratic pivot
in any affine coordinates.
Cert: `theoremA_sharp.py` (fully symbolic identities + four instances with
non-constant K, incl. deg e1 = 3, 4, each with an empty saturated collision
system). Provenance: the unnecessary hypothesis was caught by the hostile-
referee pass (agent-authored `rt_novelty_q2_anydegree.py`), then
re-derived and re-certified independently.

**Q2 [proved here; machine-checked chain] (quadratic-pivot reduction)**
Let f = e0(x') + x4·e1(x') + (σ/2)x4² with x' = (x1,x2,x3), σ ∈ C^×, any
degrees, det Hess f = c ∈ C^×. Then:
(i) with ẽ0 := e0 − e1²/(2σ) and K := K(x') := Hess₃ e1 (constant only when
deg e1 ≤ 2), the Schur identity
det Hess f = σ·det₃(Hess₃ ẽ0 + (x4 + e1(x')/σ)·K) holds;
(ii) [CORRECTED — see Q2-SHARP above] since u = x4 + e1(x')/σ sweeps C for
each FIXED x', at which K(x') is a constant matrix, the PENCIL IDENTITY
det₃(Hess₃ ẽ0(x') + s·K(x')) ≡ c/σ holds for all s ∈ C and all x' — with NO
degree hypothesis on e1. (The original entry required deg e1 ≤ 2 here; that
was over-cautious.)
(iii) any collision ∇f(p) = ∇f(q) forces (4th component) a common
λ = p4 + e1(p')/σ = q4 + e1(q')/σ and then ∇'(ẽ0 + λe1)(p') =
∇'(ẽ0 + λe1)(q'); by (ii), ẽ0 + λe1 is a 3-variable polynomial with
constant Hessian determinant, so de Bondt's HC_3 (L-HC3) gives p' = q',
hence p = q. **Conclusion (sharp): HC_4 holds for every f that is, in some
affine coordinates, quadratic in a variable with NONZERO constant leading
coefficient — with no hypothesis on deg e1.**
For deg f = 4 the tower (T3) automatically puts every essential-ternary-
leading-form candidate in this shape with deg e1 ≤ 2 (see T5 below).
Dependencies: L-HC3 (de Bondt), L-INJ. Certificate:
`verify_q2_core.py` (Schur identity on 3 instance classes incl. nontrivial
pencil; symbolic collision-transfer step). Hand proof: `src/deg4_tree_notes.md`
§1a.
Adversarial: "u sweeps C for fixed x'" is where σ ≠ 0 and CONSTANCY of the
x4²-coefficient are both used ✓; deg e1 ≤ 2 needed so that K is constant and
ẽ0 + λe1 has constant Hessian for a SCALAR λ ✓; over non-closed char-0
fields the argument survives since HC_3 is stated over any char-0 field ✓.

**T5 [proved here, pending agent re-verification] (deg-4 affine-pivot
reduction)** For deg f = 4 with det Hess f ∈ C^×:
(a) if the leading form f4 has essential rank 3, then T3 forces
(f3)_{x4x4} = 0 and the x4²-coefficient of f is the constant S44/2; hence f
is covered by Q2 (S44 ≠ 0) or is affine in x4 (S44 = 0);
(b) if f4 has essential rank 2 (WLOG f4 ∈ C[x1,x2] binary, det₂Hess₂f4 ≠ 0),
the degree-6 graded piece of det Hess f equals
det₂(Hess₂f4)·det₂(Hess_{(x3,x4)}f3), forcing Hess_{(x3,x4)}f3 = ℓ(x)·uu^T
with u constant; any 0 ≠ v ∈ span(e3,e4) with u·v = 0 then has
D_v²f = v^TSv constant — again Q2 (if ≠ 0) or affine pivot (if = 0);
(c) essential rank 1 (f4 = x1⁴): richer; [det]_5 forces
det₃(Hess_{(x2,x3,x4)}f3) ≡ 0; constant-pivot existence under analysis
(branch BR1).
**Consequence (modulo (c) and branch closures): every degree-4
counterexample to HC_4 must be, after an affine change, of the affine-pivot
form f = e0(x1,x2,x3) + x4·e1(x1,x2,x3) (class AP4).**
Status: (a) proved here (T3 + Q2); (b) proved modulo the 2×2 lemma
(agent BR2 verifying); (c) in progress (agent BR1).

**G1 [proved here] [machine]** For F ∈ {Alpöge F, Gao F₄, Gao F₅}:
the space {W ∈ C^{n×n} : W·J_F(x) symmetric for all x} is {0}. Hence
**neither known dimension-4 Jacobian counterexample is equivalent to a
gradient map under two-sided affine composition** (and the dim-3 result is
the consistency check demanded by L-HC3: an invertible symmetrizer would
have contradicted it).
Reduction [proved here]: L1∘F∘L2 is a gradient map for some affine L1, L2
⟺ ∃W ∈ GL_n with W·J_F symmetric identically (W = (T^T)^{-1}S; converse by
the polynomial Poincaré lemma). Cert: `gradient_equivalence_test.py`, which
also re-derives det J(F₄) ≡ −44/9, det J(F₅) ≡ 160/29 by direct symbolic
determinant (independent verification of Gao Thms F4/F5, determinant part;
transcription checksum).
Adversarial: two-sided AFFINE only — polynomial (non-affine) equivalence is
NOT excluded and remains open; the parse of F₄/F₅ is trusted only through
the determinant checksum (a single wrong coefficient would a.s. destroy
constancy of det J — checked to be exactly the stated constants).

---

## AP. The affine-pivot class AP4 (f = e0(x') + x4·e1(x'), det Hess f = c ≠ 0)

Setting: x' = (x1,x2,x3); the determinant identity is equivalent to
(c0) g^T adj(P) g ≡ −c, (c1) g^T L(P,K) g ≡ 0, (c2) g^T adj(K) g ≡ 0,
where g = ∇'e1, P = Hess₃e0, K = Hess₃e1 (constant iff deg e1 ≤ 2),
adj(P+tK) = adjP + tL + t²adjK. Cert for the identities and all matrix
facts below: `ap4_rank_reduction.py`, `dillen4_ap4_control.py` (control:
de Bondt's dillen4 example satisfies (c0),(c1),(c2) with c = 1 and has
empty saturated collision system).

**AP3 [proved here] [machine]** deg e1 = 2 ⟹ rank K ≤ 2: the x'-quadratic
part of (c2) is det K · x'^T K x' (identity K·adjK·K = detK·K), so
det K = 0.

**AP2 [proved here] [machine]** rank K = 2 impossible: normalize
K = diag(1,1,0); translations kill the Im K-part of the affine part of g;
(c2) = γ3² kills the rest, so g(0) = 0; evaluating (c0) at x' = 0 gives
0 = −c, contradiction.

**AP1 [proved here] [machine]** rank K = 1 ⟹ after affine change and
scaling, **e1 = x1²/2 + x2** (γ-part off Im K forced nonzero by (c0) on
{x1 = 0}; rotate (x2,x3), rescale). Hence: **the only surviving quadratic
pivot shape in AP4 is f = e0(x1,x2,x3) + x4(x1²/2 + x2)**, subject to (c0)
and (c1) as identities on e0 alone (with g = (x1,1,0)).

**AP0 [proved here] (all degrees of e0!)** deg e1 ≤ 1 (rank 0): g = const
≠ 0 (g = 0 kills the determinant); normalize e1 = x3 (rotation, x4-scaling,
x3-translation). Then the determinant identity reads
det₂ Hess_{(x1,x2)} e0 ≡ −c identically — so every x3-fiber e0(·,·,a) is a
2-variable polynomial with constant Hessian determinant. A collision of
∇f = (∂1e0, ∂2e0, ∂3e0 + x4, x3) forces equal x3 (comp 4), then equal
(x1,x2) by Dillen's HC_2 applied fiberwise, then equal x4 (comp 3).
**HC_4 holds for every AP4 potential with affine e1, of any degree.**
Dependencies: L-HC3 (n = 2 case = Dillen), L-INJ.
Adversarial: fiberwise use of HC_2 is legitimate — for each fixed a ∈ C,
e0(·,·,a) ∈ C[x1,x2] with det Hess₂ = −c ≠ 0 ✓; the normalization chain is
affine ✓; no division anywhere ✓.

**Consequence [modulo T5 branch closures]: a degree-4 counterexample to
HC_4 must be, after an affine change, f = e0(x1,x2,x3) + x4(x1²/2 + x2)
with deg e0 ≤ 4 and (c0),(c1) — the RESIDUAL NORMAL FORM. Its collision
analysis: non-injectivity ⟺ ∃p' ≠ q' with e1(p') = e1(q'),
∂3e0(p') = ∂3e0(q'), and (p1 ≠ q1, always x4-solvable) or (p1 = q1 and the
map x3 ↦ (∂3e0, ∂1e0 − x1∂2e0) non-injective on the (x1,x2)-fiber).**

## S. Branch-closing lemmas and the degree-4 theorem

**S1 [proved here] [machine] (2×2 lemma)** A symmetric 2×2 matrix of linear
forms over C with det ≡ 0 equals ℓ(x)·uu^T for a linear form ℓ and a
CONSTANT vector u. Proof: l11l22 = l12² in the UFD C[x]; if l12 = 0 one
diagonal entry vanishes; if l12 ≠ 0, the factorizations of l12² force
{l11,l22} = {γl12, γ⁻¹l12}, so M = (l12/γ)·(γ,1)(γ,1)^T.
Cert: `deg4_branch_closures.py` (ii).

**S2 [proved here] [machine] (branch graded formulas)** For f = f4+f3+f2:
if f4 ∈ C[x1,x2]: [det Hess f]_6 = det₂(Hess₂f4)·det₂(Hess_{(x3,x4)}f3);
if f4 = x1⁴: [det Hess f]_5 = 12x1²·det₃(Hess_{(x2,x3,x4)}f3).
Verified with FULLY GENERIC symbolic coefficients (multilinear in graded
pieces ⟹ the symbolic check is a proof). Cert: `deg4_branch_closures.py`
(i),(iv).

**S3 [KNOWN — cite, do not claim]** Gow arXiv:1502.05547 Thm 1 proves this
for all n over any field with |K| ≥ r+1, with our proof line for line, and
attributes it to Fillmore–Laurie–Radjavi (1985), proof of Lemma 1;
equivalently the classical Bertini theorem for a linear system of quadrics
(Landsberg arXiv:math/0609507). S1 is the r = 1 case of Meshulam /
Loewy–Radwan. See `prior_art.md`. Statement as used here: Every linear
subspace W ⊆ Sym₃(C) consisting of singular matrices has a common isotropic
vector v ≠ 0 (v^TMv = 0 ∀M ∈ W). Proof: if some M0 ∈ W has rank 2, then
adj M0 = κv0v0^T with κ ≠ 0 and M0v0 = 0; det(M0 + tM) ≡ 0 in t (all
elements singular), and its t-coefficient is tr(adjM0·M) = κ·v0^TMv0, so
v0 works. If all elements have rank ≤ 1, W = C·vv^T and any u ⊥ v works.
Cert: `deg4_branch_closures.py` (iii).
Application: in branch r = 1, M(x) := Hess_{(x2,x3,x4)}f3 is a linear
family with det ≡ 0 (by S2/[det]_5), so {M(x)} is an all-singular subspace
and a constant isotropic direction v ∈ span(e2,e3,e4) exists ⟹ D_v²f =
v^TSv = const ⟹ pivot.

**THEOREM A [proved here] (pivot theorem, all degrees)** Let f ∈ C[x1..x4],
det Hess f = c ∈ C^×, and suppose some 0 ≠ v ∈ C⁴ has D_v²f ≡ γ constant
with the adapted pivot coefficient e1 of degree ≤ 2 (in coordinates with
v = e4: f = e0(x') + x4e1(x') + (γ/2)x4²). Then ∇f is a polynomial
automorphism:
- γ ≠ 0: Q2 (pencil identity + HC_3);
- γ = 0, deg e1 ≤ 1: AP0 (fiberwise HC_2);
- γ = 0, deg e1 = 2: AP1–AP3 rank collapse + full classification
  f ≅ a(x1,x2) + x3(−κx1 + β(x2 + x1²/2)) + x4(x2 + x1²/2), κ² = c,
  with explicit injectivity. Certs: `verify_q2_core.py`,
  `ap4_rank_reduction.py`, `ap4_pivot_closure.py`.

**THEOREM B [proved here] (HC_4 in degree ≤ 4)** Every f ∈ C[x1..x4] with
deg f ≤ 4 and det Hess f ∈ C^× has ∇f a polynomial automorphism.
Proof: deg ≤ 3 by T4/Wang. deg = 4: T1+T2 put f4 ∈ C[x1,x2,x3]; by
essential rank r of f4:
r = 3: T3 forces (f3)_{x4x4} = 0, so D_{e4}²f = S44 const and e1 deg ≤ 2 —
Theorem A applies;
r = 2: S2+S1 give Hess_{(x3,x4)}f3 = ℓuu^T; the u-perp direction
v ∈ span(e3,e4) has D_v²f = v^TSv const, pivot coefficient deg ≤ 2 —
Theorem A;
r = 1: S2+S3 give a constant isotropic v ∈ span(e2,e3,e4); same conclusion.
Dependencies: T2 (Gordan–Noether, known), L-HC3 (Dillen/de Bondt, known),
L-INJ (known), + this project's T0–T5, Q2, AP0–AP3, S1–S3, all
machine-certificated. Cert of the closing identities:
`deg4_branch_closures.py`.
Adversarial notes: every "WLOG" is an affine change (class-preserving);
(f3)_{44}·det₃ ≡ 0 splits correctly in the domain C[x]; e1 ≡ const cases
are impossible (zero row in Hess ⟹ det = 0); all Q2/AP applications have
deg e1 ≤ 2 by degree counting in each branch.

## C/D. Structure of the affine-pivot class (all degrees) — hand red-team pass

**THEOREM C [FOLKLORE/ELEMENTARY, with one new observation]** The block
determinant is elementary and the doubling is Meng 2006 Prop. 1.4 /
Meng–Yang Prop. 2.3. What is ours is the resulting JC_2 equivalence for
this class. Statement: For any
a, b, e ∈ C[x1,x2], the potential f = a(x1,x2) + x3·b + x4·e has
Hess f = [[A, B],[Bᵀ, 0]] with B = ∂(b,e)/∂(x1,x2), hence
**det Hess f = (Jac(b,e))², independent of a**. So f is precisely the Meng
doubling ⟨(x3,x4),(b,e)⟩ plus an inert planar potential, and det Hess f is
a nonzero constant iff (b,e) is a planar KELLER map. Moreover ∇f is
injective IF AND ONLY IF (b,e) is injective, for EVERY inert a — the
converse direction solves the linear system Bqᵀw = ∇a(p) − ∇a(q) + Bpᵀu
(cert (W12) in `paperG_writeup_checks.py`). **Hence HC_4
restricted to this class is EQUIVALENT to JC_2.**
Cert: `doubling_structure.py` (jet-level determinant identity = proof;
plus fully generic polynomial check and a collision-Gröbner instance).
**RETRACTED SIDE-CLAIM (Aug 18):** an earlier draft asserted the
affine-pivot class is *exhausted* by these doublings. FALSE — witness
w = x1²/2 + x1x2(1+x3) + x2²(2x3+x3²)/2 + x3x4, det Hess w = 1, affine
pivot e4 with pivot coefficient x3 of degree 1, isotropic cone
{v : D_v²w ≡ 0} = C·e4 a single line (radical-membership cert (W13)),
whereas every doubling's isotropic cone contains span(e3,e4). The correct
statement is the pivot TRICHOTOMY by (γ, deg e1): only the γ = 0,
deg e1 ≥ 2 case is a doubling. Adversarial check: the invariant used
(isotropic cone) is affine-invariant, so the non-equivalence is proved,
not numerical; independently re-verified by the project lead (Gröbner
basis of the coefficient ideal recomputed from scratch).
Consequence: the AP1 classification of Theorem A is exactly the statement
that planar Keller maps with one component of the form x1²/2 + x2 are
b = −κx1 + β(x2 + x1²/2), all invertible — which is why that case closes
unconditionally.

**THEOREM D [KNOWN — cite, do not claim]** The equivalence "B(e) ≡ 0 ⟺ all
level hypersurfaces have vanishing Gauss–Kronecker curvature" is Reilly,
Rocky Mountain J. Math. 16 (1986) 553–565, Prop. 4; Fox (op. cit.) §2 names
exactly our invariant U(F) = ∇Fᵀadj(Hess F)∇F and proves U(F) = 𝒦·|dF|^{n+2}.
The determinant identity is Cauchy's. See `prior_art.md`.
The bordered-Hessian
identity det[[K, g],[gᵀ, 0]] = −gᵀadj(K)g (fully symbolic ⟹ proof) shows
constraint (c2) says exactly: **the pivot coefficient e1 has identically
vanishing bordered Hessian, i.e. all its level surfaces are developable**
(zero Gaussian curvature). This has real content: nondegenerate quadric
pivots (e.g. e1 = x1x2 + x3, bordered Hessian ≡ 1) are impossible.
Cert: `developability.py`.

**IDENTITY E [KNOWN — cite, do not claim]** Fox, Math. Nachr. 290 (2017)
293–320 (arXiv:1503.09108) §3, verbatim and with our proof; also Del Pia–
Hildebrand–Weismantel–Zemmer arXiv:1408.4711 Lemma 5.2 (any C² homogeneous
function), attributed there to Hemmer 1995. See `prior_art.md`.
For e homogeneous of degree d ≥ 2 in
any number of variables, ∇eᵀ adj(Hess e) ∇e = (d/(d−1))·e·det Hess e.
Proof: Euler ((Hess e)x = (d−1)∇e, ∇e·x = de) plus H·adj(H)·H = det(H)·H.
Verified fully symbolically for ternary forms of degrees 2, 3, 4.
Cert: `cone_leading_form.py`.

**THEOREM D′ [proved here] [machine] (cone leading form)** In every
affine-pivot potential with det Hess f ∈ C^×, the LEADING FORM e_d of the
pivot coefficient satisfies det Hess₃(e_d) ≡ 0, hence by Gordan–Noether
(n = 3) **e_d is a cone: after a linear change it involves at most two
variables**. Proof: the top-degree part of (c2) is
∇e_dᵀadj(Hess e_d)∇e_d, which by Identity E equals
(d/(d−1))·e_d·det Hess e_d; e_d ≠ 0 in the domain C[x′].
For d = 2 this recovers AP3 (rank ≤ 2) as a special case.
Cert: `cone_leading_form.py` (incl. the top-degree extraction on a generic
cubic with arbitrary lower-order terms).

**IDENTITY H [KNOWN — retraction; see prior_art.md]** This is Dolgachev,
*Classical Algebraic Geometry: A Modern View*, CUP 2012, §1.1.4,
eq. (1.20)–(1.21), pp. 17–18 ("the affine equation of the Hessian"), for
general n and **with our exact proof**. The label "derived and used here"
below was NOT defensible and is retracted; cite Dolgachev.
For e ∈ C[x1,x2,x3]
of degree d with homogenization E (degree d in four variables):
det Hess₄(E)|_{x0=1} = d(d−1)·e·det Hess₃(e) − (d−1)²·B(e), where
B(e) := ∇eᵀadj(Hess e)∇e is the bordered Hessian. Proof: row/column
reduction of Hess E at x0 = 1 using E_{0i} = (d−1)e_i − (Hx)_i and
E_{00} = d(d−1)e − 2(d−1)x·∇e + xᵀHx turns it into
[[d(d−1)e, (d−1)gᵀ],[(d−1)g, H]]; Cauchy's bordered expansion finishes.
Verified fully symbolically for generic e of degrees 2 and 3.
Cert: `homogenization_identity.py`.

**THEOREM F [essentially KNOWN]** It is a three-line corollary of
de Bondt–van den Essen, *Singular Hessians*, J. Algebra 282 (2004) 195–204,
Thm 3.3 (the n = 3 classification of zero-Hessian polynomials), which this
project failed to cite; our homogenization proof was unnecessary. Retained
below only because it is self-contained. Statement: if
e ∈ C[x1,x2,x3] has B(e) ≡ 0 AND det Hess₃(e) ≡ 0, then e is affinely
2-variable. Proof: Identity H gives det Hess₄(E) ≡ 0; Gordan–Noether in
FOUR variables makes E a cone, D_vE ≡ 0 for some v = (v₀,v′) ≠ 0.
If v₀ = 0, restricting to x0 = 1 gives D_{v′}e ≡ 0. If v₀ ≠ 0, Euler
(E₀|_{x0=1} = de − x·∇e) turns D_vE ≡ 0 into de = (x−v′)·∇e, i.e.
ẽ(y) := e(y+v′) satisfies y·∇ẽ = dẽ, so ẽ is homogeneous (converse of
Euler); B is translation-invariant, so B(ẽ) ≡ 0, and Identity E gives
det Hess ẽ ≡ 0, whence Gordan–Noether (n = 3) makes ẽ a cone. ∎
**Corollary: Conjecture E holds whenever det Hess₃(e1) ≡ 0 — in
particular for every homogeneous pivot coefficient.**
Cert: `theorem_F.py` (translation-invariance of B fully symbolic;
converse-of-Euler step; end-to-end examples).
Adversarial: GN(n = 4) needs alg. closed char 0 ✓; "a form vanishing on
the chart x0 = 1 vanishes identically" ✓ (dehomogenization is injective on
forms of fixed degree); v′ is a constant vector so the translation is
affine ✓.

**THEOREM G [proved here] [machine-supported] (dimension-free rigidity)**
For e ∈ C[x_1..x_n], if B(e) := ∇eᵀadj(Hess e)∇e ≡ 0 then det Hess e ≡ 0.
Geometrically: **if every level hypersurface of a polynomial is developable
(identically zero Gaussian curvature), its Hessian determinant vanishes
identically.** The converse is false (e = x1²/2 + x2 has det Hess = 0,
B = 1), so the rigidity is one-directional.
Proof. Suppose det Hess e ≢ 0, so φ := ∇e is dominant. Where H := Hess e is
invertible set V := H⁻¹g, g := ∇e. Then (i) Dφ·V = g = φ, so V is φ-related
to the Euler field E on the target; (ii) D_V e = gᵀH⁻¹g = B/det H. So B ≡ 0
makes e a first integral of V. At a generic x₀ with det H(x₀) ≠ 0 and
p₀ := φ(x₀) ≠ 0, let ψ be the local inverse and L(p) := ⟨ψ(p),p⟩ − e(ψ(p))
the local Legendre transform, so ∇L = ψ and EL − L = e∘ψ. Then B ≡ 0 gives
E(e∘ψ) = 0, i.e. (E² − E)L = 0. Along rays p = tσ this is t²∂²_t L = 0, so
L is affine in t; hence L = A + B̃ with A, B̃ holomorphic homogeneous of
Euler-degree 0 and 1, and this extension IS the analytic continuation of L
over the cone. Therefore ψ(tp) = t⁻¹∇A(p) + ∇B̃(p) =: t⁻¹a + b, and the
continued identity ∇e(ψ(q)) = q gives ∇e(t⁻¹a + b) = tp for all t. The left
side is a polynomial in t⁻¹ and has no positive power of t; comparing
coefficients of t¹ forces p = 0, contradicting p₀ ≠ 0. ∎
Certs: `theorem_G.py` (structural identities fully symbolic; Legendre
identities on exact examples; the homogeneous eigenvalue d/(d−1) ∉ {0,1}
sanity check; a radical-membership test in 3 variables),
`theorem_G_n2_elementary.py` (the n = 2 case proved by an independent
elementary argument — B = Tᵀ(Hess e)T with T the level-curve tangent, so
B ≡ 0 means all level curves are lines; CORRECTED Aug 20: the "parallel"
step needs the generic fibre to be SMOOTH — it avoids the finitely many
critical values — not merely reduced, and the accompanying machine scan
DROPPED LINEAR TERMS, on which B depends (B(x1²/2+x2)=1 vs B(x1²/2)=0),
so it was NOT exhaustive as first claimed; the genuinely
full-coefficient-space decision, linear terms included, degrees ≤ 4, is
`pzG_n2_decision.py`, found and written by the independent pzG audit).
Adversarial (RESOLVED Aug 18): the load-bearing step above — "this
extension IS the analytic continuation of L over the cone" — was judged a
NON-SEQUITUR by the independent audit (`advG_flow_proof.py`): ψ is only a
local inverse, det H may vanish on the ray, and no continuation theorem
was actually invoked. The theorem STANDS, with two independent repaired
proofs, neither using analytic continuation:
(1) **Eigensplit proof** (`theorem_G_paper.tex` §2, audited by the
project lead step by step): near p₀ split L = α + β with α := −e∘ψ,
β := EL; B ≡ 0 gives Eα = 0, Eβ = β. Extend α, β to the punctured cone
Γ = C^×·W by their OWN homogeneity via explicit formulas (well-defined
because W∩ℓ is connected for every line ℓ through 0); Ψ := ∇(α̃+β̃)
satisfies Ψ(tp₀) = t⁻¹a + b, and ∇e∘Ψ = id propagates from W to Γ by the
identity theorem. Laurent comparison of t¹ kills p₀.
(2) **Flow proof** (advG, elementary, any characteristic-0 field):
B ≡ 0 ⟹ MV = 2g (from ∇(gᵀV) = 2g − MV) ⟹ Jac(V)V = V − 2H⁻¹g = −V, so
the local flow of V solves x″ + x′ = 0: x(s) = ae^{−s} + b with
∇e(x(s)) = e^s·p₀. Then Q(τ) := τ·∇e(aτ+b) − p₀ vanishes at the
infinitely many values τ = e^{−s} yet Q(0) = −p₀ ≠ 0. Formal-power-series
version avoids all analysis. Certs: J1/J2 identities fully generic
n = 2,3,4, random exact jets n = 5,6; non-vacuity on rational
Euler-degree-0 witnesses. Sharpness: x1 − sqrt(x1²−2x2) has B ≡ 0,
det Hess ≠ 0, and realises ψ(tp) = t⁻¹a + b exactly
(`tgc_geometry_and_sharpness.py`) — polynomiality is the only step that
fails off polynomials, so neither proof can be weakened.
(3) **Puiseux proof** (pzG, second independent auditor): embed C(x) into
the Puiseux field at t = ∞ via dominance of ∇e; the formal Legendre value
Λ = ⟨χ, tp⟩ − e(χ) satisfies Λ″ = t⁻²(B/D)(χ) = 0, so Λ is affine in t and
the same no-positive-powers comparison closes it. No analysis; any
algebraically closed characteristic-0 field. `pzG_audit_identities.py`.
(4) **Derivation-eigenvector proof** (atkG, third independent auditor):
with δ := D_V on K = k(x), (J1)/(J2) give δg = g and δV = −V; put
b := x + V (δb = 0), a := V (δa = −a). Finite Taylor expansion
g_i(b − a) = Σ_k c_{i,k} splits g_i into δ-eigencomponents with
eigenvalues −k; δg = g forces Σ(k+1)c_{i,k} = 0; eigenvectors for
distinct eigenvalues are independent over the constants, so all
c_{i,k} = 0, g = 0, H = 0, D = 0 — contradiction. Polynomiality enters
only as finiteness of the expansion. `atkG_algebraic_proof.py` (pillars
fully generic, ALL PASS). Also confirms the paper's eigensplit proof
line by line and classifies the original defect as expositional.
FOUR proofs, three of them by independent auditors instructed to refute.
n = 2 COMPLETELY DECIDED by machine for degree ≤ 5 over ALL polynomials
(linear terms included), via a certified normalisation e = x2^d + lower
terms — `atkG_n2_decision.py`; this supersedes both earlier n = 2 scans,
whose coverage claims were overstated (linear terms dropped in
`theorem_G_n2_elementary.py`; homogeneous-only slice in
`theorem_G_n2_deg45.py`).
Independent second proof of the n ≤ 3 consequence (Corollary E) via
de Bondt–van den Essen: `tgc_debondt_chain.py`.

**CONJECTURE E — NOW A THEOREM [proved here]** Every e ∈ C[x1,x2,x3] with
B(e) ≡ 0 is affinely 2-variable. Proof: Theorem G gives det Hess₃ e ≡ 0;
Theorem F then applies. Note the conclusion is strictly stronger than
det Hess e ≡ 0 alone: de Bondt–van den Essen's n = 3 classification
contains zero-Hessian polynomials that are not affinely 2-variable
(h = x1x2 + x1²x3), and indeed B(h) = −x1⁴ ≠ 0 as it must be.

**MAIN THEOREM (pivot trichotomy) [proved here] [machine]** Let
f ∈ C[x1..x4] with det Hess f ∈ C^× admit a pivot v; put γ = D_v²f,
e1 = D_vf. Then exactly one of (the cases are indexed by the invariants
(γ, deg e1) of the pair (f,v), hence genuinely exclusive):
(i) γ ≠ 0 (quadratic pivot) ⟹ ∇f is a polynomial automorphism
(unconditional, Theorem A);
(ii) γ = 0, deg e1 = 1 ⟹ unconditional (AP0, fiberwise Dillen); NOT
always a doubling — see the retraction under Theorem C;
(iii) γ = 0, deg e1 ≥ 2 ⟹ after an affine change,
f = a(x1,x2) + x3·b(x1,x2) + x4·e(x1,x2) — a Meng doubling of the planar
Keller map (b,e) plus an inert planar potential — with ∇f injective iff
(b,e) is.
**Hence HC_4 restricted to potentials admitting a pivot is EQUIVALENT to
JC_2.** The only case left open is potentials with no pivot at all.
Proof: (c2) + Corollary E make e1 = e(x1,x2); then
det Hess f is affine in x4 with
[x4]det Hess f = −(e0)_{x3x3}·B₂(e), B₂ the 2-variable bordered Hessian of
e. Constancy forces (e0)_{x3x3}·B₂(e) ≡ 0. If B₂(e) ≡ 0 then by Theorem G
in 2 variables e = φ(L); then det Hess f = −φ′(x1)²·det₂Hess_{(x2,x3)}e0,
so φ′² divides a nonzero constant, φ is affine, deg e1 = 1: case (ii). If
deg e1 ≥ 2 then B₂(e) ≢ 0, so (e0)_{x3x3} ≡ 0 ⟹ doubling: case (iii).
NOTE (earlier drafts): the first write-up asserted (e0)_{x3x3} ≡ 0
outright (missing case ii); a later one claimed the affine-pivot class
equals the doublings (false, witness w — see Theorem C entry). Both
corrections are unconditional-side, so the JC_2 equivalence is unchanged.
Cert: `pivot_dichotomy.py` (both identities fully symbolic + an instance of
each branch).
Adversarial note: branch (iii) was MISSED in the first write-up, which
asserted (e0)_{x3x3} ≡ 0 outright. It does not weaken the equivalence
(being unconditional) but its omission was a real gap in the argument as
first stated. Found by re-deriving the proof from scratch.
Proved for: deg e1 ≤ 2 (AP1–AP3); e1 homogeneous of any degree (Identity E
+ Gordan–Noether n = 3, since then (c2) ⟺ det Hess e1 ≡ 0 ⟺ cone).
Degree 3 tested by radical-membership certificates
(`conjE_degree3.py`). **If true, the affine-pivot class is exactly the
Meng-doubling class and HC_4 restricted to it is EQUIVALENT to JC_2** — an
exact localization of the planar Jacobian conjecture inside HC_4.
Route to a proof: the homogenization identity (derived and used here)
det Hess₄(E)|_{x0=1} = d(d−1)·e1·det Hess₃(e1) − (d−1)²·(bordered Hessian
of e1), which under (c2) makes det Hess₄(E) divisible by E, plus the
degree cascade forced by D′ (the top parts of det Hess₃ e1 must vanish
successively).

## Open / in progress

**E1 [RESOLVED — see Theorem B]** HC_4 for deg f ≤ 4: PROVED.
Next frontier: deg f = 5 (pivot coefficients of degree 3 appear; by
Theorem D′ their leading forms are cones; by Conjecture E in degree 3 they
would be affinely planar, landing in the doubling class where the relevant
planar Keller maps have degree ≤ 4 — invertible by Moh). So Conjecture E in
degree 3 plus the tower would give HC_4 in degree 5.
**C2 [conjectured]** No 5-variable potential Φ = tA(w) + B(w) with
det Hess Φ ∈ C^×, det(Hess_w B + s·Hess_w A) ≡ 0 (in s, w), possesses a
gradient collision over a common pivot value, unless JC_2 is false.
(Would show: the Schur-descent mechanism cannot settle HC_4 without first
refuting JC_2.) Structure notes: admissibility forces corank-1 pencils
M(s,w) with adj(M) = κ·kk^T and ⟨k, ∇_w A⟩ nowhere zero — being developed.

---

## R. Derivation reformulation of HC (new, 2026-08-27)

**Lemma R (the polynomial-in-partials criterion) [proved here] [machine
(structural identities)]**
Let f ∈ C[x_1..x_n] with det Hess f = c ∈ C^×, g = ∇f, H = Hess f,
A = adj H, and Λ := ⟨x, g⟩ − f (the pulled-back Legendre function). Define
the polynomial derivations W_i := A e_i (= c·(inverse-Hessian coordinate
field)) and the Euler pullback V := A g. The following are equivalent:
 (a) the formal Legendre transform of f is a polynomial (Meng's HC_n holds
     for f);
 (b) f ∈ C[g_1,…,g_n] — **f is a polynomial in its own first partials**;
 (c) Λ ∈ C[g_1,…,g_n];
 (d) each W_i is a locally nilpotent derivation of C[x].
Proof: (a)⇔(b) from the Legendre identity EL − L = f∘ψ and invertibility of
(E − 1) off degree 1; (b)⇔(c) from Λ = ⟨x,g⟩ − f with W_iΛ = x_i (identity
(I4)); (d)⇔(a) is the classical LND criterion for the Keller map g
(van den Essen, *Polynomial Automorphisms*, Ch. 1: an étale map is an
automorphism iff the derivations ∂/∂g_i are locally nilpotent, and here
∂/∂g_i = W_i). ∎
Cert (structural identities, fully generic ⇒ proof): `euler_pullback_reformulation.py`.

**Structural identities behind Lemma R [proved here] [machine].** For any
f (any n), with A = adj Hess f, D = det Hess f, B = g^T A g:
 (I1) div(A g) = n·D  — the Euler pullback V = Ag has divergence nc, so
      when D = c it is a polynomial vector field of constant divergence nc.
 (I2) W_i(g_j) = D·δ_ij and W_i(x_j) = A_{ji} is symmetric in (i,j).
 (I3) [W_i, W_j] = 0 and [V, W_i] = −W_i: the pullback φ* along φ = ∇f is a
      Lie-algebra homomorphism from {∂_i, Euler} on the target. (So
      {W_1,…,W_n} is a commuting family of derivations — an integrable
      polynomial framing of C^n off {D = 0}.)
 (I4) W_iΛ = x_i, W_i f = V_i, VΛ = Λ + f, Vf = B/D, Vg = g.
All verified with fully generic coefficients (n = 2,3,4). Fail-closed.
Adversarial: (I1)–(I4) are identities for ALL f; the equivalences in
Lemma R need det Hess f ∈ C^× (else W_i are not polynomial). NECESSARY
consequence of (b): fibre constancy — ∇f(p) = ∇f(q) ⟹ f(p) = f(q). The
Meng–Yang Ψ recovers this via non-injectivity: ∇Ψ agrees at the two
DISTINCT points (1,−3/2,0,0,0), (−1,3/2,0,0,0), so ∇Ψ is not an
automorphism and Ψ ∉ C[∇Ψ] by (a)⇔(b). CAUTION (self-corrected
2026-08-27, caught by a fail-closed assert): at this collision the
Ψ-values COINCIDE (both 0), so fibre-constancy alone does NOT witness the
failure — it is necessary but not sufficient, and a given collision need
not detect it. Cert: `euler_pullback_reformulation.py` (B2).

**Why Lemma R matters for HC_4.** It turns HC_4 into: *every
f ∈ C[x_1..x_4] with det Hess f ∈ C^× is a polynomial in its four first
partials.* This is a membership question in the subalgebra C[∇f] ⊆ C[x],
and — via (d) — a local-nilpotency question about four commuting
polynomial derivations. It is the gradient specialization of the LND form
of JC_2 under the bridge, and it is the natural language in which the
pivot-free residual class (no v with D_v²f constant) should be attacked:
a pivot is exactly a coordinate in which one W_i has an especially simple
(triangular) form. Open: use (I3)'s commuting framing to force local
nilpotency in the pivot-free case.

**R-D5 (degree-5 y1-slice reduction) [proved here] [machine]** For the
weighted leading form F = a5(y) + x4 b3(y) + (1/2)x4² y1 of a pivot-free
degree-5 rank-3 potential (branch 1 of E4 = 0, c5 ≠ 0):
det Hess₄ F ≡ 0 forces **a5|_{y1=0} to be a 5th power of a linear form**
in (y2,y3). Proof (certified, `d5_y1slice_reduction.py`): (R1) the
bordered-determinant identity det Hess₄F|_{y1=0} =
−(c5 y3²+x4)²·[adj(Hess₃a5 + x4 Hess₃b3)]₁₁|_{y1=0} (fully generic a5 ⇒
proof); (R2) its x4⁰ part equals det Hess₂(a5|_{y1=0}); so det Hess₂ of
the binary quintic ā = a5|_{y1=0} vanishes; (R3) by Hesse (binary
Gordan–Noether, certified via the apolar catalecticant radical test) ā is
a 5th power. (R4, correction: the x4¹ part is identically 0 for branch 1 —
Hess(b3)|_{y1=0} has (3,3) entry 2c5·y1 → 0 — so it gives nothing; an
earlier draft wrongly concluded ā = κy3⁵.) Consequence: Hess₃(a5)|_{y1=0}
has rank 1, so det₃ Hess₃ a5 vanishes on {y1=0}; the rank-3 hypothesis
must be carried entirely by the y1-dependent part of a5. This is a sharp
structural obstruction obtained WITHOUT a many-unknown Gröbner basis; the
residual decision lives in the y1 ≥ 1 graded pieces of det Hess₄ F. Open:
push the y1-graded chain (or the exported Singular decision) to close or
populate the branch.

**R-D5-GEN (degree-5 rank-3 pivot-free, generic-b3 emptiness) [proved here]
[machine, independently certified]** For the weighted leading form
F = a5(y) + x4·b3(y) + (1/2)x4²·y1 of the degree-5 rank-3
no-isotropic-direction pivot-free branch, E4 = 0 splits the generic cubic b3
into exactly two rational families (branch 0, branch 1; `d5_rank3_pivotfree_
decision.py`). Treat the free cubic coefficients c = (c0,…,c9) as PARAMETERS
in the field ℚ(c); the a5-system
    J = ⟨ coeff_y(E0), coeff_y(E1), coeff_y(E2), coeff_y(E3) ⟩ ⊆ ℚ(c)[a0..a20]
(21 unknowns) satisfies **std(J) = (1) for each branch** — the unit ideal
over ℚ(c), computed instantly (whereas the same system over ℚ with c as
unknowns is out of Gröbner reach). Hence for GENERIC b3 there is NO a5
solving the graded system at all: **both branches are empty off a proper
closed c-locus.** The emptiness was made independently checkable by extracting
a Nullstellensatz representation Σ gᵢ·Jᵢ = 1 via Singular `lift(J, ideal(1))`
and re-verifying Σ gᵢ·Jᵢ == 1 exactly in sympy (`verify_paramcert.py`). The
denominators of the lift are the exceptional loci off which emptiness holds:
    Z0 = V( 972·c8⁷·c9⁴·(3c0c9 − c1c8)⁶ )   [branch 0],
    Z1 = V( 3·c5·(2c0c5 − c1c4)⁵ )            [branch 1],
i.e. branch 0 is empty whenever c8·c9·(3c0c9−c1c8) ≠ 0, branch 1 whenever
c5·(2c0c5−c1c4) ≠ 0. Scripts: `_d5_param.py` (parametric std over 𝔽_p(c) and
ℚ(c)), `_d5_paramcert.py` (lift → `cert_b{0,1}_gens.txt`), `verify_paramcert.py`
(independent sympy check).
Adversarial: std(J) = (1) over ℚ(c) is a fully rigorous characteristic-0
statement (Gröbner over the transcendental field ℚ(c0,…,c9)); the extracted
certificate removes any need to trust Singular's std — a single sympy
polynomial identity Σ gᵢJᵢ = 1 witnesses 1 ∈ J. What remains is EMPTINESS OFF
Z; the loci Z0, Z1 themselves are a finite stratification tail (R-D5-TAIL),
and match exactly the y-degenerations b3 → (y1²)·(linear) reached there.

**R-D5-TAIL (degree-5 rank-3 pivot-free, the exceptional-locus tail)
[proved here; deductive core + multi-prime modular certificates]** The
recursive stratification of Z0, Z1 (`_d5_close.py`; each c-substratum decided
by std(J) = (1) or a cone/isotropy test, Noetherian) closes every branch
EMPTY, and every surviving leaf collapses to a SINGLE degenerate family
    b3 = y1²·(c2 y1 + c0 y2 + c1 y3).
Structure (DEDUCTIVE, fully generic; `d5_survivor_family.py` S1–S4):
B := Hess₃ b3 has det B ≡ 0 and constant kernel v* = (0,c1,−c0); adj B ≡
−4y1²·v*v*ᵀ; with q = ∇b3 one has v*·q ≡ 0, hence qᵀadjB q ≡ 0 and
tr(A·adjB) ≡ −4y1²·(v*ᵀA v*) = −4y1²·D²_{v*}a5 where A := Hess₃ a5; and the
y1⁰-piece det A|_{y1=0} equals the y1¹-part of E0 coefficient-by-coefficient,
so det Hess₃a5 vanishes on {y1=0} over V(J) (the R-D5 obstruction, recovered
for this family). Decision of the a-system J = ⟨coeff_y E0,E1,E2,E3⟩
(the E3 coefficients are ≡ 0 here):
  • det Hess₃a5 ∉ √J — the family is NOT a rank-collapse cone; the E-system
    genuinely HAS rank-3 solutions (unbiased full-variety single-combination
    Rabinowitsch 1 ∈ ⟨J, (Σλᵢ Tᵢ)w − 1⟩, 12/12 across primes 32003, 40009,
    15013 and independent λ; `_d5_surv_rabin.py`). A slice-only test gives a
    FALSE cone here — the det A ≠ 0 locus is a lower-dimensional component
    that generic slicing misses.
  • D²_{v*}a5 = v*ᵀA v* ∈ √J — EVERY solution's a5 has the isotropic/linear
    direction v* (full-variety radical test, all VANISH across the three
    primes, `_d5_surv_target.py iso`), with the explicit power certificate
    (v*ᵀA v*)² ∈ ⟨J⟩ (reduces to 0 mod std(J); v*ᵀA v* itself ∉ J;
    `_iso_power.py`).
CONCLUSION: every rank-3 solution of the E-system for b3 = y1²·(…) carries a
linear direction v*, so it FAILS the branch's "f5 has no linear direction"
hypothesis. **The degree-5 rank-3 NO-isotropic-direction pivot-free branch is
therefore EMPTY.** Adversarial / honesty: the two radical facts are certified
MODULARLY (multi-prime, unbiased full-variety Rabinowitsch) plus the char-0
power certificate (v*ᵀAv*)² ∈ ⟨J⟩ over F_p; the char-0 PARAMETRIC lift over
ℚ(c) is beyond this environment's compute ceiling (std/sat over ℚ(c) times
out; long jobs are killed near 10 min), so it is recorded as machine-certified
(modular) rather than a hand proof. What this does NOT do: it does not settle
HC_4 (JC_2 remains the blocker via the trichotomy), and the residual
"rank-3 WITH a linear direction" degree-5 case — where such a5 route to the
pivot theorem (Theorem A) — is a separate branch. Scripts: `_d5_close.py`,
`_d5_survivors.py`, `_d5_surv_rabin.py`, `_d5_surv_target.py`, `_iso_power.py`,
`d5_survivor_family.py`.
