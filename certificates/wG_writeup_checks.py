# wG_writeup_checks.py -- fail-closed certificate for theorem_G_paper.tex
#
# Verifies, in exact symbolic arithmetic, every algebraic identity quoted in
# the referee-ready write-up ../theorem_G_paper.tex that is not already covered
# verbatim by an existing project certificate. All checks with FULLY GENERIC
# symbolic coefficients are polynomial identities in those coefficients, hence
# complete proofs for all polynomials of the stated shape.
#
# Checks:
#  (1) bordered-determinant formula det[[M,g],[g^T,0]] = -g^T adj(M) g
#      (generic symmetric 3x3 M, generic g) -- proof.
#  (2) for f = e0(x1,x2,x3) + x4*e(x1,x2) with e0, e fully generic of degree 3:
#      det Hess f is AFFINE in x4 and
#      [x4] det Hess f = -(e0)_{x3x3} * B2(e),
#      B2(e) = e22 e1^2 - 2 e12 e1 e2 + e11 e2^2 -- proof (deg-3 generic).
#  (3) doubling: det Hess(a + x3 b + x4 e) = Jac(b,e)^2, a,b,e generic deg 3
#      in (x1,x2) -- proof (deg-3 generic).
#  (4) Schur identity for a quadratic pivot: with P = Hess3 e0, K = Hess3 e1,
#      g = grad3 e1, gamma a symbol:
#      gamma^2 * det Hess(e0 + x4 e1 + gamma/2 x4^2)
#            = det3( gamma*(P + x4 K) - g g^T ) -- proof (deg-3 generic).
#      (Equivalent, division-free form of
#       det Hess f = gamma * det3(Hess3(e0 - e1^2/(2 gamma)) + u K),
#       u = x4 + e1/gamma.)
#  (5) case-(iii) formulas: det Hess(e0 + x4*x3) = -det2 Hess_{(x1,x2)} e0 and
#      det Hess(e0 + x4*phi(x1)) = -phi'(x1)^2 * det2 Hess_{(x2,x3)} e0,
#      e0 generic deg 3, phi generic univariate deg 4 -- proof.
#  (6) the cautionary example for the dichotomy's branch (b): h = x1 + x2^2
#      has det Hess2 h == 0 identically, yet h is NOT a function of one affine
#      form (its two partials are linearly independent over C) and B2(h) = 2.
#      So "det Hess2 e == 0 => e = phi(L)" is FALSE as stated in
#      main_result.tex's dichotomy proof; the corrected write-up derives
#      e = phi(L) from B2(e) == 0 (the n = 2 case of Corollary E).
#  (7) the collision-transfer gradient identity for the quadratic pivot:
#      grad3(e0 - e1^2/(2 gamma) + lambda e1)
#            = grad3 e0 + (lambda - e1/gamma) grad3 e1 -- proof.
#
# Fail-closed: any assertion failure exits nonzero.

import itertools
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X2 = (x1, x2)
X3 = (x1, x2, x3)
X4 = (x1, x2, x3, x4)


def gen_poly(vs, deg, tag):
    mons = [sp.Integer(1)]
    for d in range(1, deg + 1):
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons))


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def grad(e, vs):
    return sp.Matrix([sp.diff(e, a) for a in vs])


def det(M):
    return M.det(method='berkowitz')


def border(M, g):
    n = M.shape[0]
    top = M.row_join(g)
    bot = g.T.row_join(sp.Matrix([[0]]))
    return top.col_join(bot)


# ------------------------------------------------------------------ (1)
print('(1) bordered-determinant formula, fully symbolic 3x3 ... ', end='')
m = sp.symbols('m11 m12 m13 m22 m23 m33')
gsym = sp.Matrix(sp.symbols('g1 g2 g3'))
M = sp.Matrix([[m[0], m[1], m[2]], [m[1], m[3], m[4]], [m[2], m[4], m[5]]])
lhs = det(border(M, gsym))
rhs = -sp.expand((gsym.T * M.adjugate() * gsym)[0, 0])
assert sp.expand(lhs - rhs) == 0
print('PASS (proof)')

# ------------------------------------------------------------------ (2)
print('(2) affine-pivot potential: det Hess f affine in x4 and')
print('    [x4] det Hess f = -(e0)_x3x3 * B2(e), fully generic deg 3 ... ',
      end='', flush=True)
e0 = gen_poly(X3, 3, 'a')
e = gen_poly(X2, 3, 'b')
f = e0 + x4 * e
D = sp.expand(det(hess(f, X4)))
Dp = sp.Poly(D, x4)
assert Dp.degree() <= 1, 'det Hess f not affine in x4'
c1 = sp.expand(Dp.coeff_monomial(x4))
e_1, e_2 = sp.diff(e, x1), sp.diff(e, x2)
e_11, e_12, e_22 = sp.diff(e, x1, x1), sp.diff(e, x1, x2), sp.diff(e, x2, x2)
B2 = sp.expand(e_22 * e_1**2 - 2 * e_12 * e_1 * e_2 + e_11 * e_2**2)
pred = sp.expand(-sp.diff(e0, x3, x3) * B2)
assert sp.expand(c1 - pred) == 0
print('PASS (proof)')

# ------------------------------------------------------------------ (3)
print('(3) doubling: det Hess(a + x3 b + x4 e) = Jac(b,e)^2, generic deg 3 ... ',
      end='', flush=True)
aa = gen_poly(X2, 3, 'p')
bb = gen_poly(X2, 3, 'q')
ee = gen_poly(X2, 3, 'r')
fd = aa + x3 * bb + x4 * ee
jac = sp.expand(sp.diff(bb, x1) * sp.diff(ee, x2)
                - sp.diff(bb, x2) * sp.diff(ee, x1))
assert sp.expand(det(hess(fd, X4)) - jac**2) == 0
print('PASS (proof)')

# ------------------------------------------------------------------ (4)
print('(4) quadratic-pivot Schur identity, generic deg 3, symbolic gamma ... ',
      end='', flush=True)
gam = sp.Symbol('gamma')
E0 = gen_poly(X3, 3, 's')
E1 = gen_poly(X3, 3, 't')
fq = E0 + x4 * E1 + gam / 2 * x4**2
P = hess(E0, X3)
K = hess(E1, X3)
gg = grad(E1, X3)
lhs4 = sp.expand(gam**2 * det(hess(fq, X4)))
rhs4 = sp.expand(det(gam * (P + x4 * K) - gg * gg.T))
assert sp.expand(lhs4 - rhs4) == 0
print('PASS (proof)')

# ------------------------------------------------------------------ (5)
print('(5) case-(iii) determinant formulas, generic deg 3 / deg 4 ... ',
      end='', flush=True)
E0b = gen_poly(X3, 3, 'u')
f5 = E0b + x4 * x3
lhs5 = sp.expand(det(hess(f5, X4)))
rhs5 = sp.expand(-det(hess(E0b, X2)))          # det2 Hess_{(x1,x2)} e0
assert sp.expand(lhs5 - rhs5) == 0
phi_c = sp.symbols('v0:5')
phi = sum(c * x1**k for k, c in enumerate(phi_c))
f5b = E0b + x4 * phi
lhs5b = sp.expand(det(hess(f5b, X4)))
rhs5b = sp.expand(-sp.diff(phi, x1)**2 * det(hess(E0b, (x2, x3))))
assert sp.expand(lhs5b - rhs5b) == 0
print('PASS (proof)')

# ------------------------------------------------------------------ (6)
print('(6) det Hess2 == 0 does NOT imply e = phi(affine form): h = x1 + x2^2')
h = x1 + x2**2
assert sp.expand(det(hess(h, X2))) == 0
h1, h2 = sp.diff(h, x1), sp.diff(h, x2)
B2h = sp.expand(sp.diff(h, x2, x2) * h1**2
                - 2 * sp.diff(h, x1, x2) * h1 * h2
                + sp.diff(h, x1, x1) * h2**2)
assert B2h == 2
# h = phi(l) for an affine form l would force c1*h_x1 + c2*h_x2 == 0 for some
# (c1,c2) != 0 (the partials proportional with CONSTANT ratio). Here
# h_x1 = 1, h_x2 = 2*x2 are linearly independent over C:
c1s, c2s = sp.symbols('c1s c2s')
comb = sp.expand(c1s * h1 + c2s * h2)
sols = sp.solve([sp.Poly(comb, x1, x2).coeff_monomial(1),
                 sp.Poly(comb, x1, x2).coeff_monomial(x2)], [c1s, c2s],
                dict=True)
assert sols == [{c1s: 0, c2s: 0}] or sols == [{c2s: 0, c1s: 0}]
print('    det Hess2 h == 0, B2(h) = 2 != 0, partials independent: h is not')
print('    a univariate composition -- the corrected dichotomy proof must use')
print('    B2(e) == 0 (n = 2 Corollary E), not det Hess2 e == 0 alone.  PASS')

# ------------------------------------------------------------------ (7)
print('(7) collision-transfer gradient identity, generic deg 3 ... ',
      end='', flush=True)
lam = sp.Symbol('lam')
h_lam = E0 - E1**2 / (2 * gam) + lam * E1
lhs7 = grad(h_lam, X3)
rhs7 = grad(E0, X3) + (lam - E1 / gam) * grad(E1, X3)
diff7 = sp.expand(gam * (lhs7 - rhs7))
assert diff7 == sp.zeros(3, 1)
print('PASS (proof)')

print()
print('ALL CHECKS PASSED -- every identity quoted in theorem_G_paper.tex that')
print('is not covered by an existing project certificate is verified above.')
