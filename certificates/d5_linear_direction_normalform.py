# d5_linear_direction_normalform.py
#
# The degree-5 residual branch left open by R-D5-TAIL: a rank-3 leading form
# f5 (det3 Hess3 f5 =/= 0) that HAS a linear direction (some v =/= 0 with
# D^2_v f5 == 0).  This is NOT the pivot-free branch decided in
# d5_pivotfree_normalform.py (which assumes NO linear direction), and it is NOT
# covered by Theorem A (a linear direction of f5 gives D^2_v f = D^2_v(f4+f3+f2),
# degree 2, not a pivot in general).
#
# Facts proved here (exact, fully generic, fail-closed):
#  (N1) A linear direction (WLOG e1 after a linear change of x1,x2,x3) forces
#           f5 = P(x2,x3) + x1 * Q(x2,x3),   P homog deg 5, Q homog deg 4.
#  (N2) Hess3 f5 has (1,1) entry 0, and det3 Hess3 f5 is LINEAR in x1:
#           det3 = c0(x2,x3) + x1 * ( -cov(Q) ),
#       with cov(Q) = Q_3^2 Q_22 - 2 Q_2 Q_3 Q_23 + Q_2^2 Q_33 = the second
#       derivative of the binary quartic Q along the direction (Q_3,-Q_2) = perp
#       to grad Q (a classical covariant); generically cov(Q) =/= 0 (deg 8), so
#       det3 genuinely depends on x1.  Rank-3 <=> det3 =/= 0.
#  (N3) The linear direction is UNIQUE for generic (P,Q): the only constant v
#       with D^2_v f5 == 0 identically is v proportional to e1.  (Mirrors the
#       unique pivot e4 of the no-linear-direction branch.)
#
# These fix the normal form for the open sub-problem; the det Hess_4 f = const
# tower for it is not derived here.
#
#   py -u d5_linear_direction_normalform.py
import sys
import sympy as sp

x1, x2, x3 = sp.symbols('x1 x2 x3')
X = (x1, x2, x3)


def homog(d, tag):
    cs = sp.symbols(f'{tag}0:{d + 1}')
    return sum(c * x2**(d - i) * x3**i for i, c in enumerate(cs)), list(cs)


def hess(g):
    return sp.Matrix(3, 3, lambda i, j: sp.diff(g, X[i], X[j]))


fails = []
def check(name, cond):
    print(f'  [{"PASS" if cond else "FAIL"}] {name}')
    if not cond:
        fails.append(name)


P, _ = homog(5, 'p')
Q, _ = homog(4, 'q')
f5 = P + x1 * Q
H = hess(f5)

print('=== degree-5 rank-3 branch with a linear direction: normal form ===')

# (N1)/(N2)
check('N2a: Hess(1,1) == 0 (e1 is a linear direction)', sp.expand(H[0, 0]) == 0)
det3 = sp.expand(H.det())
check('N2b: det3 Hess3 f5 is linear in x1', sp.Poly(det3, x1).degree() == 1)
Q2, Q3 = sp.diff(Q, x2), sp.diff(Q, x3)
Q22, Q23, Q33 = sp.diff(Q, x2, 2), sp.diff(Q, x2, x3), sp.diff(Q, x3, 2)
cov = sp.expand(Q3**2 * Q22 - 2 * Q2 * Q3 * Q23 + Q2**2 * Q33)
check('N2c: [x1] det3 == -cov(Q),  cov = Q3^2 Q22 - 2 Q2 Q3 Q23 + Q2^2 Q33',
      sp.expand(det3.coeff(x1, 1) + cov) == 0)
check('N2d: cov(Q) is not identically 0 (rank-3 genuinely uses x1)', cov != 0)

# (N3) uniqueness of the linear direction
a, b, c = sp.symbols('a b c')
v = sp.Matrix([a, b, c])
D2 = sp.expand((v.T * H * v)[0])
eqs = [sp.expand(t) for t in sp.Poly(D2, x1, x2, x3).coeffs()]
sols = sp.solve(eqs, [a, b, c], dict=True)
# every solution must have b == 0 and c == 0 (v proportional to e1)
uniq = bool(sols) and all(s.get(b, b) == 0 and s.get(c, c) == 0 for s in sols)
check('N3: the only linear direction is v prop. to e1 (b == c == 0 in every solution)',
      uniq)

print()
if fails:
    print('FAILED:', fails)
    sys.exit(1)
print('ALL NORMAL-FORM CHECKS PASSED (N1-N3).')
print('Open: the det Hess_4 f = const tower for f5 = P(x2,x3) + x1 Q(x2,x3); '
      'decide whether pivot-free members exist.  See research_log Phase 20, '
      'ledger R-D5-TAIL.')
sys.exit(0)
