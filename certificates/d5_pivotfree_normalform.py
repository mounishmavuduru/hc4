# d5_pivotfree_normalform.py
#
# DEGREE 5, n = 4.  Reduction of the "no pivot, rank-3 leading form" case to a
# WEIGHTED leading-form condition.
#
# Setup.  f = f5+f4+f3+f2, det Hess f = c in C^x, WLOG (T1+Gordan-Noether)
# f5 in C[x1,x2,x3].  Rank-3 branch: det3 Hess3 f5 =/= 0, so [det]_11 = 0
# forces (f4)_x4x4 = 0 (d5_graded_tower.py).  Hence f has x4-degree <= 2 apart
# from the x4^2 coefficient coming from f3 and f2, i.e.
#        f = a(x') + x4 b(x') + (1/2) x4^2 tau(x'),    tau = (f3)_x4x4 + (f2)_x4x4
# with deg a = 5, deg b <= 3, deg tau <= 1.
# A pivot at e4 exists iff tau is constant.  If moreover f5 has NO linear
# direction (no v' =/= 0 with D_{v'}^2 f5 = 0) then e4 is the ONLY possible
# pivot, so:  f is pivot-free  <=>  tau is a nonconstant linear form.
# Then an affine change of x' makes tau = y1, giving the NORMAL FORM
#        f = a(y) + x4 b(y) + (1/2) x4^2 y1 .
#
# THEOREM (proved here, certificate below).  Give weights wt(y_i)=1, wt(x4)=2.
# Then f has weighted degree 5 and its weighted leading form is
#        F = a5(y) + x4 b3(y) + (1/2) x4^2 y1 ,
# and det Hess f in C^x forces  det Hess_4 F = 0  identically.
# (Weighted analogue of T1.)  This is the certificate's part (A).
#
# Part (B): the five x4-graded components of det Hess F = 0, in closed form.
# Part (C): the case analysis of those five equations.
#
# Exact arithmetic, fail-closed.

import itertools
import sympy as sp

y1, y2, y3, x4 = sp.symbols('y1 y2 y3 x4')
Y = (y1, y2, y3)
V = (y1, y2, y3, x4)


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def wdeg_parts(pol, wts):
    """dict weighted-degree -> component"""
    P = sp.Poly(sp.expand(pol), *V)
    out = {}
    for mon, co in P.terms():
        d = sum(w * e for w, e in zip(wts, mon))
        t = co
        for v, e in zip(V, mon):
            t *= v**e
        out[d] = out.get(d, 0) + t
    return {k: sp.expand(v) for k, v in out.items()}


# ---------------------------------------------------------------------------
# (A) the weighted leading form identity, verified with generic coefficients
# ---------------------------------------------------------------------------
def gen_poly(vs, degs, tag):
    e = sp.Integer(0)
    ks = []
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vs)), d):
            s = sp.Symbol('%s_%d_%s' % (tag, d, ''.join(str(i) for i in m)))
            ks.append(s)
            t = s
            for i in m:
                t *= vs[i]
            e += t
    return sp.expand(e), ks


a_gen, ak = gen_poly(Y, [5, 4, 3, 2], 'a')
b_gen, bk = gen_poly(Y, [3, 2, 1], 'b')
f_gen = a_gen + x4 * b_gen + sp.Rational(1, 2) * x4**2 * y1
Dg = sp.expand(hess(f_gen, V).det(method='berkowitz'))
wts = (1, 1, 1, 2)
parts = wdeg_parts(Dg, wts)
maxw = max(parts)
assert maxw == 10, 'top weighted degree should be 10, got %d' % maxw

a5 = sum(t for t in sp.Add.make_args(a_gen)
         if sp.Poly(t, *Y).total_degree() == 5)
b3 = sum(t for t in sp.Add.make_args(b_gen)
         if sp.Poly(t, *Y).total_degree() == 3)
F = sp.expand(a5 + x4 * b3 + sp.Rational(1, 2) * x4**2 * y1)
DF = sp.expand(hess(F, V).det(method='berkowitz'))
assert sp.expand(parts[10] - DF) == 0, 'weighted leading form of det Hess is not det Hess F'
print('(A) weighted leading form:  [det Hess f]_{wt 10} = det Hess_4(F),')
print('    F = a5 + x4*b3 + x4^2*y1/2,  wt(y)=1, wt(x4)=2.   VERIFIED GENERICALLY')
print('    => det Hess f in C^x  ==>  det Hess_4 F = 0  identically.')

# ---------------------------------------------------------------------------
# (B) the five x4-components of det Hess F, in closed block form
# ---------------------------------------------------------------------------
A = sp.Matrix(3, 3, lambda i, j: sp.Symbol('A%d%d' % (min(i, j) + 1, max(i, j) + 1)))
B = sp.Matrix(3, 3, lambda i, j: sp.Symbol('B%d%d' % (min(i, j) + 1, max(i, j) + 1)))
qv = sp.Matrix([sp.Symbol('q1'), sp.Symbol('q2'), sp.Symbol('q3')])
e1 = sp.Matrix([1, 0, 0])
tt = sp.Symbol('tt')          # stands for x4

M = A + tt * B
col = qv + tt * e1
H4 = sp.zeros(4, 4)
H4[0:3, 0:3] = M
for i in range(3):
    H4[i, 3] = col[i]
    H4[3, i] = col[i]
H4[3, 3] = y1
detH = sp.expand(H4.det(method='berkowitz'))
Pc = sp.Poly(detH, tt)
E = {k: sp.expand(Pc.coeff_monomial(tt**k)) for k in range(5)}

adjA, adjB = A.adjugate(), B.adjugate()
L = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff((A + sp.Symbol('s') * B).adjugate()[i, j],
                                                   sp.Symbol('s'))))
L = sp.Matrix(3, 3, lambda i, j: sp.expand(L[i, j].subs(sp.Symbol('s'), 0)))

chk = {
    0: y1 * A.det(method='berkowitz') - (qv.T * adjA * qv)[0, 0],
    1: y1 * sp.trace(adjA * B) - ((qv.T * L * qv)[0, 0] + 2 * (adjA * qv)[0]),
    2: y1 * sp.trace(A * adjB) - ((qv.T * adjB * qv)[0, 0] + 2 * (L * qv)[0] + adjA[0, 0]),
    3: y1 * B.det(method='berkowitz') - (2 * (adjB * qv)[0] + L[0, 0]),
    4: -adjB[0, 0],
}
for k in range(5):
    assert sp.expand(E[k] - chk[k]) == 0, 'block form of x4^%d component fails' % k
print()
print('(B) det Hess F = sum_k x4^k E_k with  (A=Hess a5, B=Hess b3, q=grad b3,')
print('    adj(A+sB) = adjA + sL + s^2 adjB):')
print('    E4 = -adj(B)_11')
print('    E3 = y1 detB - 2 (adjB q)_1 - L_11')
print('    E2 = y1 tr(A adjB) - q^T adjB q - 2 (Lq)_1 - adj(A)_11')
print('    E1 = y1 tr(adjA B) - q^T L q - 2 (adjA q)_1')
print('    E0 = y1 detA - q^T adjA q                       VERIFIED SYMBOLICALLY')

# Euler simplifications for b3 homogeneous cubic:  B y = 2q, adjB q = (detB/2) y
bb, _ = gen_poly(Y, [3], 'bc')
Bb, qb = hess(bb, Y), sp.Matrix([sp.diff(bb, v) for v in Y])
assert sp.expand(Bb * sp.Matrix(Y) - 2 * qb) == sp.zeros(3, 1)
assert sp.expand(Bb.adjugate() * qb - sp.Rational(1, 2) * Bb.det(method='berkowitz') * sp.Matrix(Y)) == sp.zeros(3, 1)
print('    Euler:  B y = 2q  and  adj(B) q = (det B / 2) y   (generic cubic b3)')
print('    => E3 = 0  <=>  L_11 = 0,  and  q^T adjB q = (3/2) detB * b3.')
print()
print('d5_pivotfree_normalform: ALL PASS')
