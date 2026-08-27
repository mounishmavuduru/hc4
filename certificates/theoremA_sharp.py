# THEOREM A, SHARP FORM.
#
# The hostile-referee pass (agent-authored rt_novelty_q2_anydegree.py) showed
# that the hypothesis "deg e1 <= 2" is NOT needed in the sigma != 0 branch.
# This script re-derives that independently and records the corrected theorem.
#
# CLAIM. Let f = e0(x') + x4 e1(x') + (sigma/2) x4^2 on C^4, x' = (x1,x2,x3),
# sigma != 0, det Hess f = c != 0.  Then grad f is a polynomial automorphism,
# with NO hypothesis on deg e1.
#
# PROOF. With etilde0 = e0 - e1^2/(2 sigma), g = grad' e1 and K(x') = Hess_3 e1
# (in general NOT constant), the product rule gives
#     Hess_3(etilde0) + u K   ==   Hess_3(e0 + x4 e1) - g g^T / sigma,
#     u := x4 + e1(x')/sigma,
# so the bordered/Schur expansion gives
#     det Hess f = sigma * det_3( Hess_3 etilde0(x') + u * K(x') ).
# FIX x'.  Then K(x') is a constant matrix and u sweeps all of C as x4 does,
# so det_3( Hess_3 etilde0(x') + s K(x') ) = c/sigma for EVERY s in C and every
# x'.  Hence for each SCALAR lambda the three-variable polynomial
# h_lambda := etilde0 + lambda e1 has det_3 Hess_3 h_lambda == c/sigma, a
# nonzero constant.
# A collision grad f(p) = grad f(q) forces (4th component, since
# d f/dx4 = e1 + sigma x4) a common lambda = p4 + e1(p')/sigma = q4 + e1(q')/sigma,
# and grad' f = grad' etilde0 + u grad' e1 gives grad' h_lambda(p') =
# grad' h_lambda(q').  By de Bondt's HC_3, h_lambda has injective gradient, so
# p' = q' and then p4 = q4.                                                []
#
# The point the first version of the proof missed: the pencil sweep is
# PER-POINT in x', so K need only be constant for fixed x' -- which is
# automatic -- not globally constant.
#
# Consequence (sharp corollary): any counterexample to HC_4 admits NO
# quadratic pivot whatsoever, i.e. for every direction v in C^4 \ {0}, either
# D_v^2 f is non-constant or D_v^2 f == 0 with pivot coefficient of degree >= 3.

import itertools
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X3 = (x1, x2, x3)
X4 = (x1, x2, x3, x4)
s, lam, sigma = sp.symbols('s lam sigma')


def hess(f, vs):
    return sp.Matrix([[sp.diff(f, a, b) for b in vs] for a in vs])


def monomials(vs, dmax):
    out = []
    for d in range(dmax + 1):
        out += [sp.prod(m) for m in itertools.combinations_with_replacement(vs, d)]
    return out


# ---- (1) the Schur/product-rule identity, fully symbolic, deg e1 = 4 ----
mons = monomials(X3, 4)
ac = sp.symbols(f'aa0:{len(mons)}')
bc = sp.symbols(f'bb0:{len(mons)}')
e0 = sum(c*m for c, m in zip(ac, mons))
e1 = sum(c*m for c, m in zip(bc, mons))
g = sp.Matrix([sp.diff(e1, w) for w in X3])
K = hess(e1, X3)
et0 = sp.expand(e0 - e1**2/(2*sigma))
u = x4 + e1/sigma
lhs = hess(et0, X3) + u*K
rhs = hess(sp.expand(e0 + x4*e1), X3) - (g*g.T)/sigma
assert sp.simplify(sp.expand(sp.Matrix(lhs - rhs))) == sp.zeros(3, 3)
print('(1) PASS: Hess3(etilde0) + u*K == Hess3(e0 + x4 e1) - gg^T/sigma')
print('    fully symbolic, generic e0,e1 of degree 4 => K is NOT constant.')
print('    The identity never uses constancy of K.')

# ---- (2) collision-transfer identity, fully symbolic, deg e1 = 4 ----
for w in X3:
    lhs_i = sp.diff(sp.expand(e0 + x4*e1), w)
    rhs_i = sp.diff(sp.expand(et0 + lam*e1), w).subs(lam, u)
    assert sp.simplify(sp.expand(lhs_i - rhs_i)) == 0
print('(2) PASS: grad_i f == grad_i(etilde0 + lam e1)|_{lam = x4 + e1/sigma}')
print('    fully symbolic for generic e1 of degree 4.')

# ---- (3) explicit high-degree-pivot instances that Theorem A now covers ----
#
# To satisfy the hypothesis det Hess f == const we must build pairs
# (etilde0, e1) whose PENCIL det_3(Hess3 etilde0 + s Hess3 e1) is constant in
# BOTH s and x'.  An anti-triangular base supplies them: with
#     etilde0 = x1 x3 + x2^2/2 + psi(x1)   =>  Hess3 etilde0 =
#         [[psi'', 0, 1], [0, 1, 0], [1, 0, 0]],  det == -1,
# adding s*Hess3(e1) for  e1 = x2 h(x1) + k(x1)  keeps the last row (1,0,0),
# so the determinant stays -1 for every s and x' -- while
# Hess3 e1 = [[x2 h'' + k'', h', 0], [h', 0, 0], [0, 0, 0]] is genuinely
# NON-CONSTANT.  These instances are therefore exactly what the old
# "deg e1 <= 2" hypothesis excluded and the sharp theorem now covers.
print()
print('(3) explicit sigma != 0 potentials with deg e1 > 2 (pencil verified):')
instances = [
    (x1**3, sp.Integer(3), sp.Integer(0)),                      # K constant-free, deg 3
    (x1**4, sp.Rational(-5, 2), sp.Integer(0)),                 # deg 4
    (x2*x1**2 + x1**3, sp.Integer(2), sp.Integer(0)),           # deg 3, K NON-constant
    (x2*x1**3, sp.Integer(1), x1**5),                           # deg 4, K non-constant
]
for e1v, sig, psi in instances:
    et0v = sp.expand(x1*x3 + x2**2/2 + psi)
    e0v = sp.expand(et0v + e1v**2/(2*sig))
    f = sp.expand(e0v + x4*e1v + sp.Rational(1, 2)*sig*x4**2)
    c = sp.expand(hess(f, X4).det(method='berkowitz'))
    d1 = sp.Poly(e1v, *X3).total_degree()
    Kv = hess(e1v, X3)
    kconst = all(sp.diff(Kv[i, j], w) == 0 for i in range(3) for j in range(3) for w in X3)
    pen = sp.expand(hess(sp.expand(et0v + s*e1v), X3).det(method='berkowitz'))
    ok_pen = sp.simplify(pen - c/sig) == 0
    print(f'    e1 = {e1v} (deg {d1}, K constant: {kconst}), sigma = {sig}:'
          f' det Hess f = {c}, pencil == c/sigma: {ok_pen}')
    assert c.is_number and c != 0 and ok_pen
    # direct injectivity check, independent of HC_3
    P = sp.symbols('pz1:5'); Q = sp.symbols('qz1:5'); W = sp.symbols('wz1:5')
    gr = [sp.diff(f, w) for w in X4]
    eqs = [sp.expand(gr[i].subs(dict(zip(X4, P))) - gr[i].subs(dict(zip(X4, Q))))
           for i in range(4)]
    sat = sp.expand(sum(W[i]*(P[i] - Q[i]) for i in range(4)) - 1)
    GB = sp.groebner(eqs + [sat], *(list(P) + list(Q) + list(W)), order='grevlex')
    assert list(GB.exprs) == [sp.Integer(1)], 'COLLISION FOUND - investigate!'
    print('        saturated collision system empty => gradient injective')

print()
print('THEOREM A (sharp) CERTIFIED: the sigma != 0 branch holds for pivot')
print('coefficients of ARBITRARY degree. Corollary: an HC_4 counterexample has')
print('no quadratic pivot in any affine coordinates.')
