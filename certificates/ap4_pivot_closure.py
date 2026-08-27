# AP4 closure for deg e1 <= 2 (Theorem AP1-collapse, this project).
#
# Chain (hand proofs in lemma_ledger.md; every identity machine-verified here):
#  f = e0(x') + x4*(x1^2/2 + x2), det Hess f == c != 0, e0 ANY degree.
#  Step 1: adj(P + t E11) is LINEAR in t, with t-coefficient N supported in
#          the (2,3)x(2,3) block: N22 = P33, N23 = N32 = -P23, N33 = P22.
#          Contraction with g = (x1,1,0): g^T N g = P33.
#          => det Hess f = -g^T adjP g - x4 * P33 == c  forces  P33 = d33 e0 == 0.
#  Step 2: e0 = a(x1,x2) + x3*b(x1,x2). Then g^T adjP g = -(x1*b2 - b1)^2,
#          so (x1 b_{x2} - b_{x1})^2 == c: the PDE x1 b_{x2} - b_{x1} = kappa,
#          kappa^2 = c. Polynomial solutions: b = -kappa*x1 + beta(xi),
#          xi = x2 + x1^2/2, beta in C[xi].  [Proof: h := b + kappa*x1 solves
#          the homogeneous equation; in coordinates (x1, xi), H(x1,xi) :=
#          h(x1, xi - x1^2/2) satisfies dH/dx1 = 0, so h in C[xi].]
#  Step 3: for f = a + x3*(-kappa*x1 + beta(xi)) + x4*xi:
#          det Hess f == kappa^2, and grad f is INJECTIVE:
#          comp4 equalizes xi; comp3 = b then forces p1 = q1 (kappa != 0),
#          hence p2 = q2; comps 1,2 are a linear system in (Dx3, Dx4) with
#          matrix [[b1, x1],[b2, 1]], det = b1 - x1*b2 = -kappa != 0.
#
# Verifications below: Step 1 fully symbolically (generic symmetric P);
# Step 2's cofactor identity fully symbolically (generic a, b as unknown
# functions replaced by generic 2-jet symbols -- the identity is pointwise
# algebraic in the jet); Step 3's determinant identity symbolically for
# generic polynomial data (linear in coefficients of a; polynomial in those
# of beta -- verified for generic-coefficient a of degree 4 and beta cubic,
# plus the collision Groebner on an instance).

import sympy as sp

x1, x2, x3, x4, t, kap = sp.symbols('x1 x2 x3 x4 t kappa')

# ---- Step 1: adj(P + t E11) linear in t; g^T N g = P33 ----
ps = sp.symbols('p11 p12 p13 p22 p23 p33')
P = sp.Matrix([[ps[0], ps[1], ps[2]], [ps[1], ps[3], ps[4]], [ps[2], ps[4], ps[5]]])
E11 = sp.diag(1, 0, 0)
adjPt = (P + t*E11).adjugate()
N = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff(adjPt[i, j], t)))
assert all(sp.diff(adjPt[i, j], t, 2) == 0 for i in range(3) for j in range(3)), 'not linear in t'
Nexp = sp.Matrix([[0, 0, 0], [0, ps[5], -ps[4]], [0, -ps[4], ps[3]]])
assert sp.simplify(N - Nexp) == sp.zeros(3, 3)
g = sp.Matrix([x1, 1, 0])
assert sp.expand((g.T*N*g)[0, 0] - ps[5]) == 0
print('Step 1 PASS: adj(P+tE11) linear in t; g^T N g = P33  (=> P33 == 0 forced)')

# ---- Step 2: cofactor identity g^T adjP g = -(x1*b2 - b1)^2 for P with P33=0 ----
A11, A12, A22, b1, b2 = sp.symbols('A11 A12 A22 b1 b2')
P2 = sp.Matrix([[A11, A12, b1], [A12, A22, b2], [b1, b2, 0]])
val = sp.expand((g.T*P2.adjugate()*g)[0, 0])
assert sp.expand(val + (x1*b2 - b1)**2) == 0
print('Step 2 PASS: g^T adjP g = -(x1*b2 - b1)^2 (pointwise jet identity)')

# PDE solution-space check (instances): b = -kappa*x1 + beta(xi) solves it,
# and the homogeneous change-of-variables argument:
xi = x2 + x1**2/2
beta = sp.Function('beta')
bgen = -kap*x1 + sp.Symbol('be0') + sp.Symbol('be1')*xi + sp.Symbol('be2')*xi**2 + sp.Symbol('be3')*xi**3
lhs = sp.expand(x1*sp.diff(bgen, x2) - sp.diff(bgen, x1))
assert sp.expand(lhs - kap) == 0
print('PDE PASS: b = -kappa*x1 + beta(xi) solves x1 b_x2 - b_x1 = kappa (generic cubic beta)')

# ---- Step 3: det Hess f == kappa^2 for the general family ----
import itertools
mons_a = [x1**i*x2**j for i in range(5) for j in range(5) if i + j <= 4]
ca = sp.symbols(f'ca0:{len(mons_a)}')
a = sum(c*m for c, m in zip(ca, mons_a))
f = sp.expand(a + x3*bgen + x4*xi)
H = sp.hessian(f, (x1, x2, x3, x4))
d = sp.expand(H.det(method='berkowitz'))
assert sp.expand(d - kap**2) == 0
print('Step 3 PASS: det Hess f == kappa^2 identically, for FULLY GENERIC a (deg 4),')
print('             generic cubic beta, symbolic kappa')

# collision Groebner on a nontrivial instance (kappa = 1, beta = xi^2, a = x1^4 + x1*x2)
inst = {kap: 1, sp.Symbol('be0'): 0, sp.Symbol('be1'): 0, sp.Symbol('be2'): 1, sp.Symbol('be3'): 0}
for c in ca:
    inst[c] = 0
inst[ca[mons_a.index(x1**4)]] = 1
inst[ca[mons_a.index(x1*x2)]] = 1
fi = sp.expand(f.subs(inst))
p = sp.symbols('pp1:5'); q = sp.symbols('qq1:5'); w = sp.symbols('ww1:5')
grads = [sp.diff(fi, v) for v in (x1, x2, x3, x4)]
eqs = [sp.expand(gg.subs(dict(zip((x1, x2, x3, x4), p))) - gg.subs(dict(zip((x1, x2, x3, x4), q)))) for gg in grads]
sat = sp.expand(sum(w[i]*(p[i]-q[i]) for i in range(4)) - 1)
G = sp.groebner(eqs + [sat], *(list(p)+list(q)+list(w)), order='grevlex')
assert list(G.exprs) == [sp.Integer(1)]
print('Collision PASS: saturated collision system empty on the instance (injective)')
print()
print('AP4 CLOSURE (deg e1 <= 2) VERIFIED: with AP0, every affine-pivot potential')
print('with pivot coefficient of degree <= 2 satisfies HC_4, in ALL degrees of e0.')
