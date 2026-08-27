# Control test: de Bondt's dimension-4 example (arXiv:1203.6605, eq. dillen4)
#     f = (x1 + x2^2) x3 + (x2 + (x1 + x2^2)^2) x4
# = the Meng doubling of the planar Keller map F = (x1+x2^2, x2+(x1+x2^2)^2).
# Known: invertible gradient (F invertible), NOT anti-triangularizable.
# It lies in our AP4 class (affine pivot x4, with x' = (x1,x2,x3)):
#     f = e0(x') + x4 * e1(x'),   e0 = (x1+x2^2) x3,  e1 = x2 + (x1+x2^2)^2.
# This script verifies, exactly:
#   (1) det Hess f == 1  (doubling: (-1)^2 (Jac F)^2 = 1);
#   (2) the AP identities (c0),(c1),(c2) hold with c = 1;
#   (3) the collision system of the AP4 endgame has NO solution with p' != q'
#       (via the known explicit inverse of F -- direct injectivity check of
#        grad f by solving grad f(p) = grad f(q) with a Groebner basis).
# Passing (2) validates the engine's AP4 identity derivation on a known-good
# instance; (3) validates the endgame formulation.

import sys
import sympy as sp

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
from hc4lib import hessian

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Xp = [x1, x2, x3]

e0 = (x1 + x2**2)*x3
e1 = x2 + (x1 + x2**2)**2
f = sp.expand(e0 + x4*e1)

# (1) det Hess f == 1
H = hessian(f, [x1, x2, x3, x4])
d = sp.expand(H.det(method='berkowitz'))
print('det Hess f =', d)
assert d == 1

# (2) AP identities with c = 1
g = sp.Matrix([sp.diff(e1, v) for v in Xp])
P = hessian(e0, Xp)
K = hessian(e1, Xp)
t = sp.Symbol('t')
adjP = sp.Matrix(3, 3, lambda i, j: sp.expand((P.adjugate())[i, j]))
adjK = sp.Matrix(3, 3, lambda i, j: sp.expand((K.adjugate())[i, j]))
adjPK = (P + K).adjugate()
L = sp.Matrix(3, 3, lambda i, j: sp.expand(adjPK[i, j] - adjP[i, j] - adjK[i, j]))
c0 = sp.expand((g.T*adjP*g)[0, 0])
c1 = sp.expand((g.T*L*g)[0, 0])
c2 = sp.expand((g.T*adjK*g)[0, 0])
print('(c0) g^T adjP g =', c0, ' (expect -c = -1)')
print('(c1) g^T L g    =', c1, ' (expect 0)')
print('(c2) g^T adjK g =', c2, ' (expect 0)')
assert c0 == -1 and c1 == 0 and c2 == 0

# cross-check: adj(P + tK) decomposition identity on this instance
lhs = (P + t*K).adjugate()
for i in range(3):
    for j in range(3):
        assert sp.expand(lhs[i, j] - (adjP[i, j] + t*L[i, j] + t**2*adjK[i, j])) == 0

# (3) injectivity of grad f, directly: grad f(p) = grad f(q), saturate p != q
p = sp.symbols('p1:5')
q = sp.symbols('q1:5')
w = sp.symbols('w1:5')
gradf = [sp.diff(f, v) for v in (x1, x2, x3, x4)]
eqs = [sp.expand(gf.subs(dict(zip((x1, x2, x3, x4), p)))
                 - gf.subs(dict(zip((x1, x2, x3, x4), q)))) for gf in gradf]
sat = sp.expand(sum(w[i]*(p[i]-q[i]) for i in range(4)) - 1)
G = sp.groebner(eqs + [sat], *(list(p)+list(q)+list(w)), order='grevlex')
print('collision-system Groebner basis =', list(G.exprs))
assert list(G.exprs) == [sp.Integer(1)], 'unexpected: collision system not empty!'
print('PASS: dillen4 is a verified AP4 member with det Hess = 1, satisfies')
print('(c0),(c1),(c2), and its gradient is injective (empty collision system).')
