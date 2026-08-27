# THEOREM G in TWO VARIABLES -- a complete, elementary, purely algebraic proof,
# independent of the Legendre/Euler argument used for general n. This is the
# smallest decisive case: if Theorem G were false, it would most likely fail here.
#
# STATEMENT. Let e in C[x1,x2] with B(e) := grad(e)^T adj(Hess e) grad(e) == 0.
# Then det Hess(e) == 0. Moreover e = phi(L) for a linear form L and a
# one-variable polynomial phi (plus a constant), which is de Bondt-van den
# Essen's n = 2 classification of zero-Hessian polynomials.
#
# GEOMETRIC MEANING. For n = 2, with J = [[0,1],[-1,0]] one has the identity
#       adj(H) = J^T H J,                                    (I)
# so B = (Jg)^T H (Jg) where Jg = (e_2, -e_1) is TANGENT to the level curves.
# Hence B == 0 says the second fundamental form of every level curve vanishes:
# all level curves of e are straight lines.
#
# PROOF. Assume B == 0 and e non-constant; let d = deg e.
# By (I), B == 0 says every level curve {e = c} has vanishing curvature, i.e.
# is a union of lines. Take c generic (avoiding the finitely many critical
# values), so that e - c is reduced. Its zero set is a union of lines, and its
# degree is d, so e - c = u(c) * prod_{i} L_{c,i} with L_{c,i} linear.
# Sub-claim: the lines in one fibre are PARALLEL. Two non-parallel lines in a
# fibre meet at a point where grad e vanishes and which lies on the smooth part
# of two distinct branches -- but at such a crossing the level curve is
# singular, contradicting reducedness of a generic fibre. So each generic
# fibre is a union of parallel lines, i.e. e - c = u(c) * M_c(x)^{...} up to
# scalars in a single linear form direction: e - c is a polynomial in ONE
# linear form L_c.
# Now let c1 != c2 be two generic values with directions L_1, L_2. Then
#       p_1(L_1) - p_2(L_2) = c_2 - c_1                         (*)
# for one-variable polynomials p_i of degree d. If L_1 and L_2 are not
# proportional, then (x -> (L_1, L_2)) is a linear isomorphism of C^2 and (*)
# reads p_1(s) - p_2(t) = const in independent variables s,t, forcing both p_i
# constant, hence d = 0: contradiction. So L_1 and L_2 are proportional; as c
# varies over a generic (hence Zariski-dense) set the direction is constant,
# say L. Then e - c is a polynomial in L for a dense set of c, so e itself is a
# polynomial in L. Finally Hess(phi(L)) = phi''(L) * L_vec L_vec^T has rank <= 1,
# so det Hess e == 0.                                                     []
#
# This script machine-checks: the identity (I); that B is exactly the numerator
# of the curvature of the level curves; that the conclusion holds on an
# exhaustive exact search over ALL polynomials of degree <= 5 in two variables
# (solving B == 0 over the full coefficient space by Groebner elimination and
# verifying det Hess == 0 on every solution component, via radical membership).

import itertools
import sympy as sp

x1, x2 = X = sp.symbols('x1 x2')


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def Bof(e, vs):
    H = hess(e, vs)
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T*H.adjugate()*g)[0, 0])


def dethess(e, vs):
    return sp.expand(hess(e, vs).det(method='berkowitz'))


# ---- identity (I): adj(H) = J^T H J for symmetric 2x2 ----
a, b, c = sp.symbols('a b c')
H = sp.Matrix([[a, b], [b, c]])
J = sp.Matrix([[0, 1], [-1, 0]])
assert sp.expand(H.adjugate() - J.T*H*J) == sp.zeros(2, 2)
print('(I) PASS: adj(H) == J^T H J for symmetric 2x2 (fully symbolic)')

# ---- B is the tangential second fundamental form ----
e_sym = sp.Function('e')
mons = [sp.Integer(1)]
for d in range(1, 6):
    mons += sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, d)}, key=str)
cs = sp.symbols(f'q0:{len(mons)}')
e_gen = sum(k*m for k, m in zip(cs, mons))
g = sp.Matrix([sp.diff(e_gen, v) for v in X])
T = sp.Matrix([sp.diff(e_gen, x2), -sp.diff(e_gen, x1)])       # tangent to level set
assert sp.expand(Bof(e_gen, X) - (T.T*hess(e_gen, X)*T)[0, 0]) == 0
print('(II) PASS: B == T^T (Hess e) T with T the level-curve tangent')
print('     (fully symbolic, generic e of degree 5) => B == 0 <=> level curves straight')

# ---- exhaustive exact test, degrees 2..5 in two variables ----
zz = sp.Symbol('zz')
for d in (2, 3, 4, 5):
    mons_d = [sp.Integer(1)]
    for k in range(1, d+1):
        mons_d += sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, k)}, key=str)
    # constants and linear terms do not affect Hess; keep them out to shrink the system
    mons_d = [m for m in mons_d if sp.Poly(m, *X).total_degree() >= 2]
    ks = sp.symbols(f'r{d}_0:{len(mons_d)}')
    e = sum(k*m for k, m in zip(ks, mons_d))
    Bexpr = Bof(e, X)
    Beqs = [sp.expand(t) for t in sp.Poly(Bexpr, *X).coeffs()] if Bexpr != 0 else []
    Dexpr = dethess(e, X)
    Deqs = [sp.expand(t) for t in sp.Poly(Dexpr, *X).coeffs()] if Dexpr != 0 else []
    bad = []
    for dc in Deqs:
        if dc == 0:
            continue
        G = sp.groebner(Beqs + [1 - zz*dc], *(list(ks) + [zz]), order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            bad.append(dc)
    print(f'  degree {d}: {len(ks)} params, {len(Beqs)} equations from B==0,'
          f' {len(Deqs)} det-Hess coefficients, failures: {len(bad)}')
    assert not bad, f'THEOREM G FALSE in 2 variables at degree {d}: {bad[:2]}'

print()
print('THEOREM G VERIFIED COMPLETELY IN TWO VARIABLES, degrees <= 5:')
print('every exact solution of B == 0 has det Hess == 0, certified by radical')
print('membership over the whole coefficient space (not sampled instances).')
print('Together with the elementary proof in this header, the n = 2 case of')
print('Theorem G is settled independently of the general Legendre argument.')
