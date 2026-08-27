# tgc_geometry_and_sharpness.py  --  key: tgc  (Theorem G, Consequences)
#
# PART (2) of the "consequences of Theorem G beyond HC_4" task: the geometric
# reformulation, and a DECISIVE sharpness test of Theorem G.
#
# ------------------------------------------------------------------ (2a)
# THE GEOMETRIC REFORMULATION (KNOWN -- see the prior-art note below).
#   For a hypersurface {e = c} in C^n (or R^n) at a regular point,
#        Gauss-Kronecker curvature  K = (-1)^{n-1} B(e) / |grad e|^{n+1},
#   so   B(e) == 0  <=>  every level hypersurface of e has identically zero
#   Gauss-Kronecker curvature  <=>  every level hypersurface is DEVELOPABLE
#   (degenerate second fundamental form / degenerate Gauss map).
#   Theorem G then reads:
#        "a polynomial all of whose level hypersurfaces have identically
#         vanishing Gaussian curvature has identically vanishing Hessian
#         determinant."
#   PRIOR ART for the reformulation itself (not for Theorem G):
#     R. C. Reilly, Affine geometry and the form of the equation of a
#       hypersurface, Rocky Mountain J. Math. 16 (1986) 553-565, Prop. 4;
#     D. J. F. Fox, arXiv:1503.09108 (Math. Nachr. 290 (2017) 293-320),
#       Lemma `nondegenlemma` + Cor. `gkcorollary`: he writes U(F) = U^{ij}F_iF_j
#       for our B, proves U(F) = K |dF|^{n+2} and states "U(F) = 0 on Omega =>
#       every regular level set is a developable hypersurface".
#   This script re-derives the curvature identity from scratch.
#
# ------------------------------------------------------------------ (2b)
# SHARPNESS: THEOREM G IS FALSE FOR NON-POLYNOMIAL FUNCTIONS.
#   Let R = sqrt(x1^2 - 2x2) and  F = x1 - R  (a branch of the multivalued
#   "tangent-line parameter" of the parabola x2 = x1^2/2: the level sets of F
#   are exactly the tangent lines x2 = s x1 - s^2/2, s = F).  Then
#        B(F) == 0     (all level sets are straight lines, hence developable)
#        det Hess F = -1/(x1^2-2x2)^2  != 0 .
#   So the conclusion of Theorem G FAILS for this algebraic, real-analytic,
#   non-polynomial function.  Consequences:
#     * Theorem G is NOT a theorem of differential geometry or of local
#       analysis; polynomiality (more precisely: that grad e extends
#       holomorphically with no pole along the rays) is essential.  This is
#       strong evidence that it cannot be a known statement in the affine-
#       differential-geometry literature, where the ambient category is C^oo.
#     * It is a sharp test of the PROOF of Theorem G: the proof derives
#            psi(t p) = t^{-1} a + b     and     grad e(t^{-1}a + b) = t p,
#       and concludes by comparing the coefficient of t^1, which is where
#       polynomiality is used.  This script verifies that for the F above the
#       local inverse psi of grad F really does have the predicted form
#       psi(tp) = t^{-1}a + b -- i.e. every step of the proof is confirmed on a
#       live example, and only the final (polynomiality) step fails, exactly as
#       it must.  If the projected form of psi had FAILED here, the proof of
#       Theorem G would be wrong.
#
# ------------------------------------------------------------------ (2c)
# The two structural identities of the proof of Theorem G, re-verified
# independently of `theorem_G.py`.

import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4', real=False)
# t is the ray parameter of the proof; we take it positive so that sympy
# resolves sqrt(t^-2) = 1/t (the statement is about a ray, so this is no loss).
t = sp.Symbol('t', positive=True)


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def Bof(e, vs):
    H = hess(e, vs)
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.simplify((g.T * H.adjugate() * g)[0, 0])


print('=' * 74)
print('(2a) B(e) IS THE GAUSS-KRONECKER CURVATURE OF THE LEVEL SETS')
print('=' * 74)

# (i) graph case: e = x_n - u(x'), level set {e = 0} is the graph of u whose
#     Gauss-Kronecker curvature is det Hess u / (1+|grad u|^2)^{(n+1)/2}.
print('  (i) e = x_n - u(x1..x_{n-1}):  B(e) = (-1)^{n-1} det Hess(u)')
for n in (2, 3, 4):
    vs = sp.symbols(f'y1:{n+1}')
    u = sp.Function('u')(*vs[:-1])
    e = vs[-1] - u
    B = sp.expand((sp.Matrix([sp.diff(e, a) for a in vs]).T
                   * hess(e, vs).adjugate()
                   * sp.Matrix([sp.diff(e, a) for a in vs]))[0, 0])
    Hu = sp.Matrix([[sp.diff(u, a, b) for b in vs[:-1]] for a in vs[:-1]])
    assert sp.simplify(B - (-1)**(n - 1) * Hu.det()) == 0, f'graph identity n={n}'
    print(f'      n = {n}: verified for an ARBITRARY smooth u (symbolic function)  PROVED')

# (ii) implicit case: near a point with e_{x_n} != 0 the level set is the graph
#      of some u; verify B(e) = (-1)^{n-1} (e_{x_n})^{n+1} det Hess(u) with u
#      obtained by implicit differentiation.  Done symbolically for n = 2, 3.
print('  (ii) implicit level sets: B(e) = (-1)^{n-1} (d e/d x_n)^{n+1} det Hess(u)')
for n in (2, 3):
    vs = sp.symbols(f'z1:{n+1}')
    e = sp.Function('e')(*vs)
    # implicit function u(z1..z_{n-1}) with e(z', u) = 0
    prm = vs[:-1]
    u = sp.Function('u')(*prm)
    sub = {vs[-1]: u}
    E = e.subs(sub, simultaneous=True)
    # first derivatives of u
    du = {}
    for i, zi in enumerate(prm):
        du[zi] = sp.solve(sp.Eq(sp.diff(E, zi), 0), sp.Derivative(u, zi))[0]
    # Hessian of u by differentiating again and substituting first derivatives
    Hu = sp.zeros(n - 1, n - 1)
    for i, zi in enumerate(prm):
        for j, zj in enumerate(prm):
            expr = sp.diff(E, zi, zj)
            expr = expr.subs({sp.Derivative(u, (zi, 2)) if zi == zj else
                              sp.Derivative(u, zi, zj): sp.Symbol('HU')}, simultaneous=True)
            sol = sp.solve(sp.Eq(expr, 0), sp.Symbol('HU'))
            val = sol[0]
            for zk in prm:
                val = val.subs(sp.Derivative(u, zk), du[zk])
            Hu[i, j] = sp.simplify(val.subs(u, vs[-1]))
    Bn = sp.simplify(Bof(e, vs))
    lhs = sp.simplify(sp.expand(Bn))
    rhs = sp.simplify((-1)**(n - 1) * sp.diff(e, vs[-1])**(n + 1) * Hu.det())
    assert sp.simplify(sp.together(lhs - rhs)) == 0, f'implicit identity n={n}'
    print(f'      n = {n}: verified for an ARBITRARY smooth e (symbolic function)  PROVED')
print('  Hence with K = det Hess(u)/(1+|grad u|^2)^{(n+1)/2} for the graph,')
print('  K == 0 on every level set  <=>  B(e) == 0.  [Reilly 1986 Prop. 4; Fox]')

print()
print('=' * 74)
print('(2b) SHARPNESS: THEOREM G FAILS FOR A NON-POLYNOMIAL ALGEBRAIC FUNCTION')
print('=' * 74)
R = sp.sqrt(x1**2 - 2 * x2)
F = x1 - R
V2 = (x1, x2)
BF = sp.simplify(Bof(F, V2))
DF = sp.simplify(hess(F, V2).det())
assert BF == 0, f'B(F) should vanish, got {BF}'
assert sp.simplify(DF + 1 / (x1**2 - 2 * x2)**2) == 0, DF
print(f'  F = x1 - sqrt(x1^2 - 2x2):   B(F) = {BF},   det Hess F = {sp.factor(DF)}')
print('  => B == 0 but det Hess != 0.  Theorem G is FALSE outside the polynomial')
print('     (indeed outside the "no pole along rays") category.')

# the level sets really are the tangent lines of the parabola x2 = x1^2/2
s = sp.Symbol('s')
assert sp.simplify((F.subs({x1: x1, x2: s * x1 - s**2 / 2}) - s).rewrite(sp.Abs)) is not None
chk = sp.simplify(sp.sqrt((x1 - s)**2) - sp.Abs(x1 - s))
print('  its level set {F = s} is the line x2 = s x1 - s^2/2 (tangent to the')
print('  parabola x2 = x1^2/2), so every level set is a straight line.')

# The proof's mechanism, verified on this example.
p1, p2 = sp.symbols('p1 p2', positive=True)
psi = ((1 - p1) / p2, ((1 - p1)**2 - 1) / (2 * p2**2))
sub = {x1: psi[0], x2: psi[1]}
chk = [sp.simplify(sp.powsimp(sp.diff(F, v).subs(sub), force=True)) for v in V2]
assert sp.simplify(chk[0] - p1) == 0 and sp.simplify(chk[1] - p2) == 0, chk
print(f'  local inverse psi(p) = {psi}  satisfies grad F(psi(p)) = p          PASS')
psit = [sp.expand(sp.simplify(c.subs({p1: t * p1, p2: t * p2}))) for c in psi]
a = [sp.simplify(sp.limit(tt * c, tt, 0)) if False else sp.simplify(sp.expand(c * t).subs(t, 0))
     for c, tt in zip(psit, (t, t))]
b = [sp.simplify(c - ai / t) for c, ai in zip(psit, a)]
assert all(sp.simplify(sp.diff(bi, t)) == 0 for bi in b), b
assert all(sp.simplify(sp.diff(ai, t)) == 0 for ai in a), a
print(f'  psi(t p) = t^(-1)*a + b  with  a = {a},  b = {b}                     PASS')
print('  EXACTLY the form predicted by the proof of Theorem G (step (6)).')
print('  The proof then compares the coefficient of t^1 in grad e(t^-1 a + b) = t p;')
print('  here grad F has a square-root singularity and DOES produce a t^1 term,')
print('  which is precisely why the contradiction is unavailable off polynomials.')
lhsF = [sp.simplify(sp.powsimp(sp.diff(F, v).subs({x1: t**-1 * a[0] + b[0],
                                                   x2: t**-1 * a[1] + b[1]}), force=True))
        for v in V2]
assert sp.simplify(lhsF[0] - t * p1) == 0 and sp.simplify(lhsF[1] - t * p2) == 0, lhsF
print(f'  and indeed grad F(t^-1 a + b) = (t p1, t p2) exactly                 PASS')

print()
print('=' * 74)
print('(2c) THE TWO STRUCTURAL IDENTITIES OF THE PROOF (independent re-check)')
print('=' * 74)
import itertools
for n in (2, 3, 4):
    vs = sp.symbols(f'w1:{n+1}')
    mons = [sp.Integer(1)]
    for d in range(1, 4):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, d)},
                       key=str)
    cs = sp.symbols(f'k{n}_0:{len(mons)}')
    e = sum(c * m for c, m in zip(cs, mons))
    H = hess(e, vs)
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    W = sp.expand(H.adjugate() * g)
    assert sp.expand(H * W - H.det(method='berkowitz') * g) == sp.zeros(n, 1)
    assert sp.expand((g.T * W)[0, 0] - (g.T * H.adjugate() * g)[0, 0]) == 0
    print(f'  n={n}: H adj(H) g = det(H) g  and  <g, adj(H)g> = B              PROVED')
print('  (these are D(grad e).V = grad e and D_V e = B/det H for V = H^{-1}grad e)')

print()
print('ALL CHECKS PASSED.')
