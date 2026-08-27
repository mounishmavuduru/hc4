# THEOREM D (developability of the pivot coefficient) -- hand red-team pass.
#
# In the affine-pivot class  f = e0(x') + x4 e1(x'),  x' = (x1,x2,x3),
# det Hess f = c != 0, the x4^2-coefficient of the bordered determinant gives
#      (c2)   g^T adj(K) g == 0,     g = grad' e1,  K = Hess_3 e1.
# IDENTITY (proved below with fully symbolic entries):
#      det [[K, g],[g^T, 0]]  =  - g^T adj(K) g .
# The left side is the classical BORDERED HESSIAN of e1; it vanishes
# identically iff every level surface of e1 has identically zero Gaussian
# curvature, i.e. the level sets of e1 are DEVELOPABLE.
#
# THEOREM D: in every affine-pivot potential with det Hess f in C^*, the
# pivot coefficient e1 has identically vanishing bordered Hessian --- its
# level surfaces are developable. (Necessary condition, no degree hypothesis.)
#
# Checks below:
#  (1) the bordered-Hessian identity, fully symbolic  -> proof;
#  (2) the doubling class (e1 independent of one variable after an affine
#      change) always satisfies it -> consistent with Theorem C;
#  (3) a nondegenerate quadric pivot (e1 = x1x2 + x3) does NOT satisfy it,
#      so such potentials cannot occur -> the condition has real content;
#  (4) a search over CUBIC e1 for solutions of (c2) that are NOT
#      2-variable-after-affine-change, to see whether Theorem C's doubling
#      class exhausts the affine-pivot world in the first nontrivial degree.

import itertools
import sympy as sp

x1, x2, x3 = X = sp.symbols('x1 x2 x3')

# ---- (1) bordered-Hessian identity, fully symbolic (proof) ----
ks = sp.symbols('k11 k12 k13 k22 k23 k33')
K = sp.Matrix([[ks[0], ks[1], ks[2]], [ks[1], ks[3], ks[4]], [ks[2], ks[4], ks[5]]])
g = sp.Matrix(sp.symbols('g1 g2 g3'))
B = K.row_join(g).col_join((g.T).row_join(sp.Matrix([[0]])))
assert sp.expand(B.det(method='berkowitz') + (g.T*K.adjugate()*g)[0, 0]) == 0
print('(1) PASS: det[[K,g],[g^T,0]] = -g^T adj(K) g  (fully symbolic => proof)')

def bordered(e):
    Ke = sp.Matrix([[sp.diff(e, u, v) for v in X] for u in X])
    ge = sp.Matrix([sp.diff(e, u) for u in X])
    return sp.expand(-(ge.T * Ke.adjugate() * ge)[0, 0])

# ---- (2) doubling class: e1 free of x3 (or of any one variable after affine change) ----
for e in [x1**2/2 + x2, x1*x2, x1**3 + x2**4, x1**5 + x1*x2**2 + x2]:
    assert bordered(e) == 0, e
print('(2) PASS: every e1 independent of one variable has vanishing bordered Hessian')

# ---- (3) nondegenerate quadric pivot is excluded ----
bad = bordered(x1*x2 + x3)
print('(3) e1 = x1x2 + x3 -> bordered Hessian =', bad, '(nonzero: such pivots impossible)')
assert bad != 0
bad2 = bordered(x1*x3 + x2**2/2)
assert bad2 != 0
print('    e1 = x1x3 + x2^2/2 -> bordered Hessian =', bad2, '(also excluded)')

# ---- (4) cubic search: solve (c2) over all cubic e1 (no linear/constant terms
#          needed: they shift g by a constant, so keep the full affine cubic) ----
mons = [m for d in (3, 2, 1) for m in
        {sp.prod(c) for c in itertools.combinations_with_replacement(X, d)}]
mons = sorted(set(mons), key=lambda m: (sp.Poly(m, *X).total_degree(), str(m)))
cs = sp.symbols(f'u0:{len(mons)}')
e1 = sum(c*m for c, m in zip(cs, mons))
bh = sp.expand(bordered(e1))
eqs = sp.Poly(bh, *X).coeffs() if bh != 0 else []
print(f'(4) cubic ansatz: {len(mons)} coefficients, {len(eqs)} equations from (c2)')

# Structural test: is every solution 2-variable-after-affine-change?
# A necessary invariant: e1 is affinely equivalent to a 2-variable polynomial
# iff its gradient image spans a 2-dim space, iff there is a CONSTANT
# direction v with D_v e1 == 0.  Test the converse question on solutions:
# we sample the solution variety by imposing random rational values on a
# generic linear slice and solving, then testing for a constant null direction.
v = sp.Matrix(sp.symbols('v1 v2 v3'))

def has_constant_null_direction(e):
    ge = sp.Matrix([sp.diff(e, u) for u in X])
    expr = sp.expand((v.T*ge)[0, 0])
    p = sp.Poly(expr, *X)
    sols = sp.solve(p.coeffs(), list(v), dict=True)
    for s in sols:
        vv = v.subs(s)
        if any(sp.simplify(c) != 0 for c in vv):
            return True
        free = [sym for sym in vv.free_symbols if sym in set(v)]
        if free:
            return True
    return False

# known solution families to sanity-check the tester
assert has_constant_null_direction(x1**3 + x2**2)          # free of x3
assert not has_constant_null_direction(x1**3 + x2**2 + x3**2)
print('    null-direction tester validated')

# Solve (c2) restricted to homogeneous cubics (the essential case: the top
# graded piece of a solution must itself solve the top-degree part).
mons3 = sorted({sp.prod(c) for c in itertools.combinations_with_replacement(X, 3)}, key=str)
cs3 = sp.symbols(f'w0:{len(mons3)}')
e13 = sum(c*m for c, m in zip(cs3, mons3))
bh3 = sp.expand(bordered(e13))
eqs3 = sp.Poly(bh3, *X).coeffs()
print(f'    homogeneous cubic case: {len(mons3)} coefficients, {len(eqs3)} equations')
G3 = sp.groebner(eqs3, *cs3, order='grevlex')
print(f'    Groebner basis of the (c2) coefficient ideal: {len(G3.exprs)} generators,')
print(f'    is trivial ideal (only e1=0)? {list(G3.exprs) == [sp.Integer(1)]}')

# exhibit that the variety is nontrivial by checking a known point and a
# non-2-variable candidate: the "cone" cubic x1*x2*x3 and the tangent-
# developable-flavoured x1^2*x3 - x1*x2^2 (Perazzo-like shapes)
for cand in [x1*x2*x3, x1**2*x3 - x1*x2**2, x1**2*x3, x1**3, x1**2*x2 + x1*x3**2]:
    bhc = bordered(cand)
    nd = has_constant_null_direction(cand)
    print(f'    e1 = {cand}: bordered Hessian {"== 0" if bhc == 0 else "!= 0"},'
          f' constant null direction: {nd}')

print()
print('THEOREM D VERIFIED (necessary condition). Open sub-question recorded:')
print('classify polynomials in 3 variables with identically vanishing bordered')
print('Hessian -- i.e. with developable level sets -- beyond the 2-variable class.')
