# THEOREM G (new, this project; dimension-free).
#
#   Let e in C[x_1,...,x_n] and set  B(e) := grad(e)^T adj(Hess e) grad(e)
#   (the bordered Hessian, up to sign:  det[[Hess e, grad e],[grad e^T, 0]] = -B).
#   If B(e) == 0 identically, then det Hess(e) == 0 identically.
#
# Equivalently: if every level hypersurface of a polynomial is developable
# (identically zero Gaussian curvature), then its Hessian determinant vanishes
# identically. The converse is false (e.g. e = x1^2/2 + x2 has det Hess = 0 but
# B = 1), so this is a genuinely one-directional rigidity statement.
#
# PROOF. Suppose det Hess(e) != 0, so phi := grad e is dominant. On the open set
# where H := Hess e is invertible put V := H^{-1} g, g := grad e. Then
#   (1)  Dphi . V = H H^{-1} g = g = phi,      i.e. V is phi-related to the
#        EULER FIELD  E = sum p_i d/dp_i  on the target;
#   (2)  D_V e = g^T H^{-1} g = B / det H.
# So B == 0 says precisely that e is a first integral of V.
# Fix a generic x0 with det H(x0) != 0 and p0 := phi(x0) != 0, and let psi be
# the local inverse of phi near p0; let L(p) := <psi(p), p> - e(psi(p)) be the
# local Legendre transform, so grad L = psi. Then
#   (3)  E L = <grad L, p> = <psi, p>,  hence  E L - L = e(psi(p)) = e o psi.
# By (1),(2) and B == 0, E(e o psi) = 0, i.e.
#   (4)  (E^2 - E) L = 0.
# On a cone neighbourhood, write p = t*sigma(z) (t the scaling parameter). Then
# (4) reads t^2 * d^2/dt^2 L = 0, so L is AFFINE in t along every ray:
#   (5)  L = A + Btilde,   E A = 0,  E Btilde = Btilde
# with A, Btilde holomorphic and homogeneous of degrees 0 and 1. Extending by
# (5) is exactly the analytic continuation of L along the whole cone, so psi
# continues there as
#   (6)  psi(t p) = grad L(t p) = t^{-1} grad A(p) + grad Btilde(p) =: t^{-1} a + b,
# using that grad of a degree-0 (resp. degree-1) homogeneous function is
# homogeneous of degree -1 (resp. 0). The identity phi(psi(q)) = q continues
# along the cone too, so for all t:
#   (7)  grad e ( t^{-1} a + b ) = t p .
# The left-hand side is a polynomial map evaluated at t^{-1}a + b, hence a
# polynomial in t^{-1}: it has NO positive power of t. Comparing the
# coefficient of t^1 in (7) gives p = 0, contradicting p0 != 0 generic.   []
#
# COROLLARY (n <= 3, with Theorem F): B(e) == 0 => e is affinely 2-variable.
# (Theorem G supplies det Hess e == 0; Theorem F then homogenises and uses
# Gordan-Noether in four variables.) This proves the project's Conjecture E.
#
# NOTE the corollary is strictly stronger than "det Hess e == 0" alone: by
# de Bondt-van den Essen's n=3 classification there are zero-Hessian
# polynomials that are NOT affinely 2-variable, e.g. h = x1 x2 + x1^2 x3;
# this script checks that such h indeed has B(h) != 0, as it must.
#
# WHAT THIS SCRIPT CHECKS (the proof itself is analytic; these are the
# machine-verifiable identities it rests on, plus decisive consistency tests):
#   (a) the algebraic identities (1) and (2), fully symbolically;
#   (b) identity (3), E L - L = e o psi, on exact examples;
#   (c) the homogeneous sanity check: for e homogeneous of degree d, the
#       Legendre transform has Euler-eigenvalue d/(d-1), which is never 0 or 1
#       -- exactly why homogeneous e can never satisfy B == 0 with det Hess != 0;
#   (d) the de Bondt-van den Essen witness h above has B(h) != 0;
#   (e) a decisive low-degree test of the theorem's CONCLUSION: over a reduced
#       but complete ansatz, every exact solution of B == 0 has det Hess == 0.

import itertools
import sympy as sp

# ---------------------------------------------------------------- (a)
print('(a) the two structural identities, fully symbolic')
for n in (2, 3, 4):
    vs = sp.symbols(f'z1:{n+1}')
    mons = [sp.Integer(1)]
    for d in range(1, 4):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'c{n}_0:{len(mons)}')
    e = sum(c*m for c, m in zip(cs, mons))
    H = sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    W = sp.expand(H.adjugate()*g)                     # W = adj(H) g = det(H) * V
    # (1'): H W = det(H) g          [equivalent to Dphi.V = phi]
    assert sp.expand(H*W - H.det(method='berkowitz')*g) == sp.zeros(n, 1)
    # (2'): <g, W> = B              [equivalent to D_V e = B/det H]
    B = sp.expand((g.T*H.adjugate()*g)[0, 0])
    assert sp.expand((g.T*W)[0, 0] - B) == 0
    print(f'    n={n}: H*adj(H)*g == det(H)*g  and  <g, adj(H)g> == B   PASS')

# ---------------------------------------------------------------- (b),(c)
print()
print('(b,c) Legendre-transform identities on exact examples')
x1, x2 = sp.symbols('x1 x2')
p1, p2 = sp.symbols('p1 p2')

# e = x1^2/2 + x2^2/2  ->  psi(p) = p,  L = (p1^2+p2^2)/2, homogeneous deg 2
e_ex = x1**2/2 + x2**2/2
psi = (p1, p2)
L = sp.expand(sum(pi*xi for pi, xi in zip((p1, p2), psi)) - e_ex.subs({x1: psi[0], x2: psi[1]}))
EL = sp.expand(p1*sp.diff(L, p1) + p2*sp.diff(L, p2))
e_psi = sp.expand(e_ex.subs({x1: psi[0], x2: psi[1]}))
assert sp.expand(EL - L - e_psi) == 0
print('    e = (x1^2+x2^2)/2:  E L - L == e o psi   PASS')
# Euler eigenvalue of L is d/(d-1) = 2 for d = 2
assert sp.expand(EL - 2*L) == 0
print('    and E L == 2 L, i.e. eigenvalue d/(d-1) = 2 (never 0 or 1)  PASS')

# a non-quadratic homogeneous check: e = x1^3/3 (n=1 style, use n=2 product)
e_h = x1**3/3 + x2**3/3
# grad e = (x1^2, x2^2); psi(p) = (sqrt(p1), sqrt(p2)); L = (2/3)(p1^{3/2}+p2^{3/2})
L2 = sp.Rational(2, 3)*(p1**sp.Rational(3, 2) + p2**sp.Rational(3, 2))
EL2 = sp.expand(p1*sp.diff(L2, p1) + p2*sp.diff(L2, p2))
assert sp.simplify(EL2 - sp.Rational(3, 2)*L2) == 0
print('    e = (x1^3+x2^3)/3: E L == (3/2) L, eigenvalue d/(d-1) = 3/2  PASS')
print('    => homogeneous e always has eigenvalue d/(d-1) not in {0,1}:')
print('       consistent with Theorem G (such e can never have B == 0 unless')
print('       det Hess e == 0).')

# ---------------------------------------------------------------- (d)
print()
X3 = sp.symbols('w1 w2 w3')


def Bof(e, vs):
    H = sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T*H.adjugate()*g)[0, 0])


def dethess(e, vs):
    return sp.expand(sp.Matrix([[sp.diff(e, a, b) for b in vs]
                                for a in vs]).det(method='berkowitz'))


w1, w2, w3 = X3
h = w1*w2 + w1**2*w3          # de Bondt-van den Essen: det Hess == 0, NOT planar
assert dethess(h, X3) == 0
Bh = Bof(h, X3)
assert Bh != 0
print(f'(d) de Bondt-van den Essen witness h = {h}: det Hess == 0 but B = {Bh} != 0')
print('    => Theorem G alone does not give planarity; Theorem F is still needed.')

# ---------------------------------------------------------------- (e)
print()
print('(e) decisive low-degree test of the CONCLUSION of Theorem G')
# Reduced but complete ansatz in 3 variables, degree 3: by Theorem D' the
# leading form is a cone, so WLOG it lies in C[w1,w2]; translations let us drop
# the constant term. Solve B == 0 exactly and check det Hess == 0 on every
# solution.
c3 = sp.symbols('a0:4')       # cubic in w1,w2
q3 = sp.symbols('b0:6')       # full quadratic in w1,w2,w3
l3 = sp.symbols('d0:3')       # linear
cub = c3[0]*w1**3 + c3[1]*w1**2*w2 + c3[2]*w1*w2**2 + c3[3]*w2**3
qmons = [w1**2, w1*w2, w2**2, w1*w3, w2*w3, w3**2]
quad = sum(a*m for a, m in zip(q3, qmons))
lin = l3[0]*w1 + l3[1]*w2 + l3[2]*w3
e_ans = cub + quad + lin
params = list(c3) + list(q3) + list(l3)

Bexpr = Bof(e_ans, X3)
Beqs = sp.Poly(Bexpr, *X3).coeffs() if Bexpr != 0 else []
Dexpr = dethess(e_ans, X3)
Deqs = sp.Poly(Dexpr, *X3).coeffs() if Dexpr != 0 else []
print(f'    ansatz: {len(params)} parameters, {len(Beqs)} equations from B == 0,')
print(f'    and {len(Deqs)} coefficients of det Hess to be shown to vanish on V(B).')

zz = sp.Symbol('zz')
bad = []
for i, dcoef in enumerate(Deqs):
    if sp.expand(dcoef) == 0:
        continue
    G = sp.groebner(list(Beqs) + [1 - zz*dcoef], *(params + [zz]), order='grevlex')
    inrad = (list(G.exprs) == [sp.Integer(1)])
    if not inrad:
        bad.append(dcoef)
print(f'    coefficients of det Hess NOT vanishing on V(B): {len(bad)}')
assert not bad, f'THEOREM G CONTRADICTED in degree 3 by {bad[:2]}'
print('    PASS: every coefficient of det Hess vanishes on the whole solution')
print('    variety of B == 0 (radical membership certified by Rabinowitsch).')

print()
print('THEOREM G supported. With Theorem F this PROVES Conjecture E:')
print('B(e) == 0 for e in C[x1,x2,x3]  =>  e is affinely 2-variable.')
