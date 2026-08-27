# ============================================================================
# pzG_audit_identities.py -- ADVERSARIAL AUDIT of THEOREM G (agent key pzG)
#
# Independent second-party verification. This agent audited every step of the
# project's Legendre/analytic proof of
#
#   THEOREM G. e in C[x_1..x_n], B(e) := grad(e)^T adj(Hess e) grad(e) == 0
#              ==>  det Hess(e) == 0.
#
# AUDIT RESULT (details in the final report):
#   * Steps (1)-(4) of the write-up are correct (re-derived by hand; the
#     machine-checkable parts are certified below).
#   * The flagged continuation step (5)-(7) is REPAIRABLE: no continuation of
#     psi through the branch locus is needed. Correct mechanism: on a ball
#     U0 around p0 set Btld := EL, A := L - EL (holomorphic); (E^2-E)L = 0
#     makes L affine along every ray-segment inside the convex ball, whence
#     A, Btld are homogeneous of degree 0, 1 along such segments; the formula
#     L_ext(t p) := A(p) + t Btld(p) is then a WELL-DEFINED single-valued
#     holomorphic function on the punctured cone Gamma = C^* . U0 (consistency
#     of overlapping representations t p = t' p' follows from convexity of the
#     ball: the ray segment [p, (t/t')p] stays in it). phi(grad L_ext(q)) = q
#     then holds on all of Gamma by the IDENTITY THEOREM (both sides are
#     holomorphic on the connected open Gamma and agree on U0). det Hess e
#     vanishing somewhere on the ray is irrelevant: psi_ext is given by an
#     explicit formula, never by inverting phi along the ray.
#   * In addition this agent constructed a PURELY ALGEBRAIC proof (Puiseux
#     series at t = infinity; no analytic function theory at all), valid over
#     every algebraically closed field of characteristic 0:
#
# ---------------------------------------------------------------------------
# ALGEBRAIC PROOF (pzG, Puiseux at infinity). Suppose B(e) == 0 and
# D := det Hess e =/= 0 in C[x]. Then g := grad e is dominant (standard:
# if h of minimal degree vanished on the image, then h(g(x)) == 0; chain rule
# gives H . (grad h)(g(x)) == 0, so on the dense set {D != 0} each
# (d_i h)(g(x)) = 0, hence identically; minimality forces grad h = 0, h
# constant -- contradiction). Since both sides have dimension n, C(x) is a
# FINITE algebraic extension of C(g_1,...,g_n).
#
# Let F := C(p_1,...,p_n) (indeterminates), Fbar its algebraic closure, and
#   P := U_m  Fbar((t^{-1/m}))     (Puiseux field at t = infinity)
# which is algebraically closed (char 0). The n coordinates t*p_i are
# algebraically independent over C [if Q(t p) = 0 in C[p][t], put t = 1 to
# get Q == 0], so q_i |-> t p_i embeds C(q) into P; since C(x)/C(q) (via
# q_i = g_i(x)) is algebraic and P is algebraically closed, this extends to
# an embedding  C(x) -> P,  x_j |-> chi_j. Consequences:
#   (E1)  g(chi) = t p              (vector identity in P^n);
#   (E2)  h(chi) != 0 in P for EVERY nonzero h in C[x]   (field embedding!)
#         -- in particular D(chi) != 0, so H(chi) is invertible over P.
# Derivations: d/dt acts formally on P (Fbar-linearly, d/dt t^{1/m} =
# (1/m) t^{1/m - 1}) with kernel exactly Fbar; each d/dp_i extends uniquely
# from C(p) to Fbar (char 0) and acts coefficientwise on P with
# d/dp_i(t^{1/m}) = 0; all these derivations commute.
#
# Step 1. Apply d/dt to (E1):  H(chi) . chi' = p, so  chi' = H(chi)^{-1} p.
# Step 2. Legendre value  Lam := <chi, t p> - e(chi) in P.  Then
#           Lam' = <chi', t p> + <chi, p> - <g(chi), chi'> = <chi, p>
#         using (E1), and
#           Lam'' = <chi', p> = p^T H(chi)^{-1} p
#                 = t^{-2} (g^T H^{-1} g)(chi) = t^{-2} B(chi)/D(chi) = 0
#         because B is the zero polynomial. Since ker(d/dt) = Fbar,
#           Lam = alpha + beta t   with  alpha, beta in Fbar.
# Step 3. Apply D_i := d/dp_i to Lam = <chi, t p> - e(chi):
#           D_i Lam = <D_i chi, t p - g(chi)> + t chi_i = t chi_i.
#         Also D_i Lam = D_i alpha + t D_i beta with D_i alpha, D_i beta in
#         Fbar. Hence   chi_i = (D_i beta) + (D_i alpha) t^{-1}.
# Step 4. Substitute into (E1): with a := D alpha, b := D beta in Fbar^n
#         (t-free!),   g(b + a t^{-1}) = t p.  The left side is a polynomial
#         evaluated at b + a t^{-1}: a Laurent polynomial in t with exponents
#         <= 0. The right side has the term t^1 . p. Comparing Puiseux
#         coefficients of t^1:  p = 0  in Fbar^n -- absurd, the p_i are
#         indeterminates. Contradiction, so D == 0.                       []
#
# The proof uses no topology, no convergence, no continuation: only formal
# Puiseux series and derivations. It is the exact algebraic skeleton of the
# project's Legendre argument (Lam plays the role of L(t p), Step 3 is
# "grad L = psi", Step 4 is the "no positive powers of t" comparison).
#
# WHAT THIS SCRIPT CERTIFIES (fail-closed, exact arithmetic):
#  (A) the universal identity H adj(H) = det(H) I (generic symmetric matrix,
#      n = 2,3,4) -- underlies chi' = H^{-1}p and B/D = g^T H^{-1} g;
#  (B) the Legendre identities grad L = psi and E L - L = e o psi FULLY
#      SYMBOLICALLY for a generic nondegenerate quadratic in n = 3 (proof for
#      all quadratics), plus the ODE identity (E^2 - E)L = (B/det H) o psi on
#      an exact non-quadratic example;
#  (C) the whole Puiseux mechanism EXACTLY on a solvable branch (e = x1^3/3
#      + x2^2/2, chi = (s q, s^2 p2), t = s^2, p1 = q^2): Lam' = <chi, p>,
#      Lam'' = t^{-2} B(chi)/D(chi) != 0, and Lam is NOT affine in t --
#      confirming that the affineness of Lam is exactly the B == 0 input;
#  (D) the final comparison step as a THEOREM for generic shapes: for fully
#      symbolic t-free vectors a, b and fully generic e (n = 2 deg 5;
#      n = 3 deg 4), grad e(b + a/t) contains NO positive power of t;
#  (E) consistency witnesses: B(x1 x2 + x1^2 x3) = -x1^4 != 0 (de Bondt-
#      van den Essen), B(phi(L)) == 0 and det Hess == 0, and the n = 2
#      linear-term sensitivity B(x1^2/2 + x2) = 1 (see pzG report: the
#      exhaustive claim of theorem_G_n2_elementary.py dropped linear terms).
# ============================================================================

import itertools
import sympy as sp

print('pzG audit: machine-checkable pillars of the Theorem G proofs')
print('=' * 72)

# ---------------------------------------------------------------- (A)
print('(A) H*adj(H) == det(H)*I for a generic symmetric matrix (proof, n=2,3,4)')
for n in (2, 3, 4):
    ent = sp.symbols(f'h{n}_0:{n*(n+1)//2}')
    Hm = sp.zeros(n, n)
    k = 0
    for i in range(n):
        for j in range(i, n):
            Hm[i, j] = Hm[j, i] = ent[k]
            k += 1
    assert sp.expand(Hm*Hm.adjugate() - Hm.det(method='berkowitz')*sp.eye(n)) == sp.zeros(n, n)
    print(f'    n={n}: PASS (fully generic symmetric)')

# ---------------------------------------------------------------- (B)
print('(B) Legendre identities, fully symbolic for generic quadratics (n=3)')
xs = sp.Matrix(sp.symbols('x1 x2 x3'))
ps = sp.Matrix(sp.symbols('p1 p2 p3'))
s6 = sp.symbols('s11 s12 s13 s22 s23 s33')
S = sp.Matrix([[s6[0], s6[1], s6[2]], [s6[1], s6[3], s6[4]], [s6[2], s6[4], s6[5]]])
cv = sp.Matrix(sp.symbols('c1 c2 c3'))
e_q = (sp.Rational(1, 2)*(xs.T*S*xs)[0, 0] + (cv.T*xs)[0, 0])
# psi = S^{-1}(p - c) as exact rational expression
Sinv = S.adjugate()/S.det(method='berkowitz')
psi = Sinv*(ps - cv)
L = sp.together((psi.T*ps)[0, 0] - e_q.subs(dict(zip(list(xs), list(psi)))))
gradL = sp.Matrix([sp.diff(L, p) for p in ps])
for i in range(3):
    assert sp.simplify(sp.together(gradL[i] - psi[i])) == 0
print('    grad L == psi                     PASS (fully generic S, c)')
EL = sum(p*sp.diff(L, p) for p in ps)
e_psi = e_q.subs(dict(zip(list(xs), list(psi))))
assert sp.simplify(sp.together(EL - L - e_psi)) == 0
print('    E L - L == e o psi                PASS (fully generic S, c)')

# non-quadratic ODE identity: (E^2-E)L = (B/det H) o psi for e = x1^3/3+x2^2/2
p1, p2 = sp.symbols('p1 p2', positive=True)
psi2 = (sp.sqrt(p1), p2)
L2 = sp.Rational(2, 3)*p1**sp.Rational(3, 2) + p2**2/2
E = lambda F: p1*sp.diff(F, p1) + p2*sp.diff(F, p2)
x1, x2 = sp.symbols('x1 x2')
e_c = x1**3/3 + x2**2/2
Hc = sp.Matrix([[sp.diff(e_c, a, b) for b in (x1, x2)] for a in (x1, x2)])
gc = sp.Matrix([sp.diff(e_c, a) for a in (x1, x2)])
BoverD = sp.cancel((gc.T*Hc.adjugate()*gc)[0, 0]/Hc.det())
lhs = sp.simplify(E(E(L2)) - E(L2))
rhs = sp.simplify(BoverD.subs({x1: psi2[0], x2: psi2[1]}))
assert sp.simplify(lhs - rhs) == 0
print('    (E^2-E)L == (B/det H) o psi       PASS (exact non-quadratic example)')

# ---------------------------------------------------------------- (C)
print('(C) the Puiseux mechanism, exactly, on the branch of e = x1^3/3 + x2^2/2')
s, q, P2 = sp.symbols('s q P2')          # t = s^2, p1 = q^2, p2 = P2
chi = sp.Matrix([s*q, s**2*P2])
t = s**2
pvec = sp.Matrix([q**2, P2])
# (E1): g(chi) == t p exactly
gchi = sp.Matrix([chi[0]**2, chi[1]])
assert sp.expand(gchi - t*pvec) == sp.zeros(2, 1)
ddt = lambda F: sp.cancel(sp.diff(F, s)/(2*s))    # d/dt = (1/2s) d/ds
Lam = sp.expand((chi.T*(t*pvec))[0, 0] - e_c.subs({x1: chi[0], x2: chi[1]}))
assert sp.simplify(ddt(Lam) - (chi.T*pvec)[0, 0]) == 0
print("    Lam' == <chi, p>                  PASS (exact)")
Lam2 = sp.simplify(ddt(ddt(Lam)))
BD_chi = sp.cancel(BoverD.subs({x1: chi[0], x2: chi[1]}))
assert sp.simplify(Lam2 - BD_chi/t**2) == 0
print("    Lam'' == t^{-2} (B/det H)(chi)    PASS (exact)")
assert Lam2 != 0
print("    Lam'' != 0 here (B != 0): affineness of Lam is exactly the B==0 input")

# ---------------------------------------------------------------- (D)
print('(D) final step: grad e(b + a/t) has NO positive power of t (generic proof)')
tt = sp.Symbol('t')
for n, dmax in ((2, 5), (3, 4)):
    vs = sp.symbols(f'y1:{n+1}')
    mons = []
    for d in range(1, dmax+1):
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'e{n}_0:{len(mons)}')
    e_gen = sum(c*m for c, m in zip(cs, mons))
    av = sp.symbols(f'a1:{n+1}')
    bv = sp.symbols(f'b1:{n+1}')
    sub = {v: bv[i] + av[i]/tt for i, v in enumerate(vs)}
    maxdeg = 0
    for v in vs:
        comp = sp.expand(sp.diff(e_gen, v).subs(sub) * tt**dmax)  # clear denominators
        d_t = sp.Poly(comp, tt).degree()
        maxdeg = max(maxdeg, d_t - dmax)
    assert maxdeg <= 0, f'positive t-power appeared: {maxdeg}'
    print(f'    n={n}, deg<={dmax}, fully generic e, a, b: max t-exponent = {maxdeg}  PASS')
print('    => the identity g(b + a/t) = t p forces p = 0: the contradiction')
print('       closing both the analytic and the Puiseux proof.')

# ---------------------------------------------------------------- (E)
print('(E) consistency witnesses')
w1, w2, w3 = W = sp.symbols('w1 w2 w3')


def Bof(e, vs):
    H = sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T*H.adjugate()*g)[0, 0])


def dethess(e, vs):
    return sp.expand(sp.Matrix([[sp.diff(e, a, b) for b in vs]
                                for a in vs]).det(method='berkowitz'))


h = w1*w2 + w1**2*w3
assert dethess(h, W) == 0 and Bof(h, W) == sp.expand(-w1**4)
print('    B(x1x2 + x1^2x3) = -x1^4 != 0, det Hess == 0    PASS (dBvdE witness)')
for k in (3, 4, 5):
    e = (2*w1 - w2 + 3*w3)**k - 5*(2*w1 - w2 + 3*w3)
    assert Bof(e, W) == 0 and dethess(e, W) == 0
print('    e = phi(L): B == 0 and det Hess == 0, k = 3,4,5  PASS (positive control)')
e_lin = x1**2/2 + x2
assert Bof(e_lin, (x1, x2)) == 1
e_nolin = x1**2/2
assert Bof(e_nolin, (x1, x2)) == 0
print('    B(x1^2/2 + x2) = 1 but B(x1^2/2) = 0: B DEPENDS on linear terms.')
print('    (Hence dropping linear terms, as theorem_G_n2_elementary.py does,')
print('     is a real restriction of its "exhaustive" scan; repaired in')
print('     pzG_n2_decision.py.)')

print()
print('pzG AUDIT: all machine-checkable pillars PASS. Theorem G verdict: TRUE;')
print('proof status: original analytic proof repairable as stated above, and a')
print('purely algebraic Puiseux proof (header) now stands independently.')
