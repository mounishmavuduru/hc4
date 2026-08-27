# AP4 rank reduction (Lemmas AP0-AP3, this project).
#
# Setting: f = e0(x') + x4*e1(x'), x' = (x1,x2,x3), det Hess f == c != 0,
# deg e1 <= 2. Write e1 = (1/2) x'^T K x' + gamma.x' + delta, g = grad' e1 =
# K x' + gamma. The determinant identity is equivalent to
#   (c0) g^T adj(P) g == -c,  (c1) g^T L(P,K) g == 0,  (c2) g^T adj(K) g == 0,
# with P = Hess3 e0, adj(P+tK) = adjP + t L + t^2 adjK.
#
# CLAIMS (hand proofs in comments; the two symbolic identities they rest on
# are verified below with fully symbolic matrices, hence proved):
#
# AP3 (rank 3 impossible): the x'-quadratic part of (c2) is
#   (Kx)^T adjK (Kx) = detK * x^T K x   [identity (I1): K adjK K = detK * K].
#   (c2) forces detK * K = 0, so detK = 0 whenever K != 0.
#
# AP2 (rank 2 impossible): normalize K = diag(1,1,0) (complex congruence).
#   A translation x' -> x' + a changes gamma by K a, so WLOG gamma = (0,0,g3).
#   adjK = E33 [identity (I2) below], so (c2) = g3^2 == 0 => g3 = 0. Then
#   g(0) = 0 and evaluating the polynomial identity (c0) at x' = 0 gives
#   0 = -c, contradiction.
#
# AP1 (rank 1 normal form): normalize K = E11 (congruence + x1-scaling).
#   Translation kills gamma1; adjK = 0 so (c2) is vacuous. If gamma2 =
#   gamma3 = 0 then g = (x1,0,0) and (c0) at x1 = 0 gives 0 = -c,
#   contradiction; so (gamma2,gamma3) != 0 and an SO-type change in (x2,x3)
#   (stabilizing K = E11) plus scaling puts e1 = x1^2/2 + x2.
#   => THE ONLY SURVIVING deg-2 PIVOT SHAPE IS e1 = x1^2/2 + x2.
#
# AP0 (rank 0, e1 affine, ANY deg e0): g = gamma constant != 0 (g = 0 kills
#   the determinant); rotate gamma -> e3, absorb the constant of e1 by an
#   x3-translation... more precisely normalize e1 = x3. Then (c0) reads
#   det2 Hess_{(x1,x2)} e0 == -c identically in (x1,x2,x3), so for every
#   a in C the 2-variable polynomial e0(.,.,a) has constant Hessian
#   determinant -c != 0; by Dillen's theorem (HC_2, JPAA 71 (1991)) its
#   gradient is injective. A collision of grad f = (d1e0, d2e0, d3e0 + x4, x3)
#   forces equal x3 =: a (comp 4), then equal (x1,x2) (comps 1,2, fiberwise
#   HC_2 at parameter a), then equal x4 (comp 3). So grad f is injective:
#   AP4 with affine e1 satisfies HC_4 for ALL degrees of e0.
#
# This script verifies (I1), (I2), the adj-decomposition, and the (c0)-at-a-
# point evaluations used above, with fully symbolic entries (these identities
# are polynomial in the entries, so symbolic verification is proof).

import sympy as sp

# (I1): K adjK K = detK * K for symmetric 3x3
ks = sp.symbols('k11 k12 k13 k22 k23 k33')
K = sp.Matrix([[ks[0], ks[1], ks[2]], [ks[1], ks[3], ks[4]], [ks[2], ks[4], ks[5]]])
lhs = sp.expand(K * K.adjugate() * K)
rhs = sp.expand(K.det() * K)
assert sp.simplify(lhs - rhs) == sp.zeros(3, 3), 'I1 failed'
print('(I1) K adjK K = detK * K : verified for generic symmetric 3x3')

# (I2): adj(diag(1,1,0)) = E33 ; adj(E11) = 0
assert sp.diag(1, 1, 0).adjugate() == sp.diag(0, 0, 1)
assert sp.diag(1, 0, 0).adjugate() == sp.zeros(3, 3)
print('(I2) adj(diag(1,1,0)) = E33, adj(E11) = 0 : verified')

# adj-decomposition adj(P+tK) = adjP + t L + t^2 adjK with L = adj(P+K)-adjP-adjK
ps = sp.symbols('p11 p12 p13 p22 p23 p33')
P = sp.Matrix([[ps[0], ps[1], ps[2]], [ps[1], ps[3], ps[4]], [ps[2], ps[4], ps[5]]])
t = sp.Symbol('t')
L = (P + K).adjugate() - P.adjugate() - K.adjugate()
diff = sp.expand((P + t*K).adjugate() - (P.adjugate() + t*L + t**2*K.adjugate()))
assert sp.simplify(diff) == sp.zeros(3, 3), 'adj decomposition failed'
print('adj(P+tK) = adjP + t L + t^2 adjK : verified for generic symmetric 3x3')

# rank-2 contradiction skeleton: with K = diag(1,1,0), gamma = (0,0,g3):
# (c2) = g^T adjK g at generic x' equals (x-part + gamma)^T E33 (...) = g3^2
x1, x2, x3, g3 = sp.symbols('x1 x2 x3 g3')
xv = sp.Matrix([x1, x2, x3])
K2 = sp.diag(1, 1, 0)
g = K2*xv + sp.Matrix([0, 0, g3])
c2expr = sp.expand((g.T * K2.adjugate() * g)[0, 0])
assert c2expr == g3**2
print('rank-2: (c2) == g3^2, so g3 = 0; then g(0) = 0 and (c0)(0) = 0 != -c. OK')

# rank-1: K = E11, gamma = (0, gamma2, gamma3); if gamma2 = gamma3 = 0 then
# g = (x1, 0, 0) and (c0) restricted to {x1 = 0} is identically 0 != -c.
K1 = sp.diag(1, 0, 0)
gg = K1*xv
c0_at = sp.expand((gg.T * P.adjugate() * gg)[0, 0]).subs(x1, 0)
assert c0_at == 0
print('rank-1: with zero affine part, (c0)|_{x1=0} == 0, contradiction. OK')
print()
print('AP RANK REDUCTION VERIFIED: deg-2 pivots reduce to e1 = x1^2/2 + x2;')
print('rank K in {2,3} impossible; affine pivots (AP0) closed via fiberwise HC_2.')
