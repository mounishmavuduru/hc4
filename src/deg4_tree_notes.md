# Scaffold: the degree-4 case tree for HC_4 (project-internal derivations)

These derivations were produced by the lead researcher and MUST be
re-verified symbolically by any agent that uses them (fail-closed asserts).
Notation: f = f4 + f3 + f2 in C[x1..x4] (WLOG f1 = f0 = 0; f0 irrelevant,
f1 only translates the gradient's image — drop), det Hess f == c != 0.
x' = (x1,x2,x3). Verified lemmas available (see lemma_ledger.md): T0-T4.

## 0. Global normalizations

- Affine changes act by f -> f o T (+affine); Hess(f o T) = T^T (Hess f o T) T;
  det scales by det(T)^2; gradient collisions correspond bijectively.
- T0: det Hess f2 = c, so S := Hess f2 is invertible.
- T1+T2 (GN_4): WLOG f4 in C[x1,x2,x3]. Let r = minimal number of variables
  of f4 after linear changes (essential rank), r in {1,2,3} (r=0 i.e. f4=0
  means deg f <= 3: HC holds by T4/Wang).
- T3: det3(Hess3 f4) * (f3)_{x4x4} == 0.

## 1. Case r = 3 (f4 essential ternary, det3 Hess3 f4 != 0)

T3 forces (f3)_{x4x4} = 0, i.e. f3 = A3(x') + x4*B2(x') (deg B2 <= 2).
Then the x4^2-coefficient of f is e2 = sigma/2 with sigma := S_44 CONSTANT,
and deg_{x4} f <= 2:
    f = e0(x') + x4*e1(x') + (sigma/2)*x4^2,
    e0 = f4 + A3 + (1/2) x'^T S' x',   e1 = B2 + s.x'   (S = [[S', s],[s^T, sigma]]).

### 1a. sigma != 0  —  PROVED (Theorem Q2, quadratic pivot; ANY degree)

Let g := grad' e1, K := Hess3 e1 (CONSTANT since deg e1 <= 2), and
etilde0 := e0 - e1^2/(2 sigma)  (polynomial).  In coordinates (x4, x'):
Hess f = [[sigma, g^T],[g, Hess3(e0 + x4 e1)]], and the Schur identity gives
    det Hess f = sigma * det3( Hess3 etilde0 + (x4 + e1(x')/sigma) * K ).
For fixed x', u := x4 + e1(x')/sigma sweeps all of C, hence the PENCIL IDENTITY
    det3( Hess3 etilde0 + s K ) == c/sigma   for ALL s in C.        (P)
Collision transfer: grad f(p) = grad f(q) gives (4th component) a common
lambda := p4 + e1(p')/sigma = q4 + e1(q')/sigma, and then (first 3 components)
grad'(etilde0 + lambda e1)(p') = grad'(etilde0 + lambda e1)(q').
   [CHECK: grad' f = grad' e0 + x4 grad' e1 = grad' etilde0 +
    (x4 + e1/sigma) grad' e1; with the common lambda this is
    grad'(etilde0 + lambda e1) evaluated at p', q'.]
By (P) with s = lambda, etilde0 + lambda e1 is a 3-variable polynomial with
constant Hessian determinant c/sigma != 0, so by de Bondt's HC_3 theorem its
gradient is injective: p' = q', then p4 = q4. Hence grad f injective; by
injectivity=>automorphism, f satisfies HC_4.  QED (modulo L-HC3).
AGENT TASK Q2V: verify the two bracketed identities symbolically for generic
e0 (deg <= 4), e1 (deg <= 2), symbolic sigma; verify (P) is forced; write
certificates/q2_verification.py, fail-closed.

### 1b. sigma = 0  —  the AFFINE-PIVOT class AP4 (the hard core)

f = e0(x') + x4 e1(x').  Bordered determinant (corner 0):
    det Hess f = - g^T adj3( Hess3 e0 + x4 K ) g   == c,             (AP)
identically in (x', x4), where g = grad' e1 (affine entries), K = Hess3 e1.
Expanding adj3(P + tK) = adj P + t L(P,K) + t^2 adj K
(L(P,K) := adj(P+K) - adj P - adj K), (AP) is equivalent to the three identities
    (c0)  g^T (adj3 Hess3 e0) g == -c
    (c1)  g^T L(Hess3 e0, K) g == 0
    (c2)  g^T (adj3 K) g == 0.
Injectivity endgame: grad f = (grad' e0 + x4 grad' e1,  e1).  Non-injectivity
<=> exist p' != q' with e1(p') = e1(q') and
    grad' e0(p') - grad' e0(q') = q4 grad' e1(q') - p4 grad' e1(p')
for some p4, q4  (i.e. the difference lies in the span of the two gradients;
NOTE p4, q4 are FREE — eliminate them: the condition is
    grad'e0(p') - grad'e0(q') IN span{ grad'e1(p'), grad'e1(q') } ).
AGENT TASK AP4: with deg e0 <= 4, deg e1 <= 2: use the residual group
( x' -> A x' + b;  x4 -> alpha x4 + ell(x');  actions:
  e1 -> alpha^{-1} * e1 o (A,b);  e0 -> e0 o (A,b) - (ell o ...) e1 ... —
  DERIVE the exact action and verify it preserves (AP) ) to normalize K by
rank: rank K in {0,1,2,3}; solve (c2),(c1),(c0) per rank case; enumerate the
resulting families; for each family either exhibit an explicit polynomial
inverse (triangular structure) or run the collision Groebner
(with e1-level and span conditions, saturating p' != q') to decide
injectivity. ANY consistent solution with p' != q' would be an HC_4
COUNTEREXAMPLE: report immediately with exact points and re-verify.
Warning: rank K = 0 means e1 affine, g constant: then (c0): g^T adj(Hess e0) g
== -c with g constant != 0 (g = 0 makes det Hess f = 0, impossible);
normalize g = e3-direction: adj(Hess3 e0)_{33} = -c: a 2x2-minor condition
det2(Hess_{(x1,x2)} e0) == -c: an MA_2-type constraint on e0 restricted to
the x3-fibers — classify via de Bondt HC_2 methods fiberwise.

## 2. Case r = 2 (f4 binary, det2 Hess2 f4 != 0 as binary quartic)

[det]_6 piece forces  det2( Hess_{(x3,x4)} f3 ) == 0.
LEMMA (verify symbolically): a symmetric 2x2 matrix M of LINEAR forms with
det M == 0 is  M = ell(x) * u u^T  for a linear form ell and a CONSTANT
vector u (possibly 0).  [Proof: l11 l22 = l12^2 in the UFD C[x]; linear
factors force proportionality; handle l12 = 0 cases separately.]
Consequence: there is v in span(e3,e4)\{0} with u.v = 0, and for that v:
  v^T Hess f4 v = 0 (f4 in C[x1,x2]),  v^T Hess f3 v = ell (u.v)^2 = 0,
  v^T Hess f2 v = v^T S v = CONSTANT.
So D_v^2 f == v^T S v =: tau, a constant.
  - tau != 0: quadratic pivot in direction v => Theorem Q2 applies => HC holds.
  - tau = 0: affine pivot: rotate v to e4-direction => f lands in AP4 form
    (with the residual (x3-mixing) freedom; the resulting (e0, e1) also
    satisfy (AP)).  Hand to the AP4 analysis.
AGENT TASK BR2: verify the 2x2 lemma and the consequence chain symbolically;
work out whether the r=2 tower imposes EXTRA constraints beyond AP4-form
(it does: f4 binary is a special leading structure inside AP4) — i.e., the
r=2/tau=0 families are a SUBFAMILY of AP4; parameterize them explicitly and
pass them through the same (c0)-(c2) + collision endgame. Also handle the
deg_{x4} f = 3 possibilities (x4^3 with constant coefficient beta from f3:
show they are absorbed or excluded — check [det]_m pieces).

## 3. Case r = 1 (f4 = x1^4 after scaling)

[det]_5 forces det3( Hess_{(x2,x3,x4)} f3 ) == 0  (3x3, linear entries).
This does NOT force the ell*uu^T shape (3x3 singular symmetric linear-form
matrices are richer — Perazzo-flavored).  AGENT TASK BR1: classify the
possibilities for the 3x3 kernel structure enough to either (i) find a
constant isotropic direction v in span(e2,e3,e4) with v^T M v == 0
(then D_v^2 f = v^T S v = const and the Q2/AP4 dichotomy repeats), or
(ii) exhibit subfamilies with NO constant pivot direction — these are the
potentially genuinely-new configurations: analyze their [det]_m constraints
directly and run collision Groebner. Be exhaustive about subcases; report
any subcase you cannot close as OPEN with its exact parameterization.

## 4. Bookkeeping demands (all agents)

- Verify EVERY scaffold identity before use; write fail-closed .py
  certificates under certificates/ named deg4_<branch>_*.py.
- Exact arithmetic only. Saturate v != 0 / p' != q' via a linear functional
  trick: add new symbols w with w . (p'-q') = 1 (this is equivalent to
  p' != q' over C — some coordinate differs — CHECK: it asserts existence of
  a linear functional valued 1 on the difference; solvable iff p' != q').
- If Groebner does not terminate in reasonable time, SPLIT the case further
  (rank conditions, vanishing/nonvanishing of specific coefficients) rather
  than reporting failure.
- Output: a summary JSON (via your structured return), files written, and
  for each closed family: the explicit inverse or the emptiness certificate;
  for each open family: exact parameterization + what was tried.
