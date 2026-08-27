# Conjecture E in degree 3.
#
# CONJECTURE E: if e1 in C[x1,x2,x3] has identically vanishing bordered
# Hessian  g^T adj(Hess e1) g == 0  (constraint (c2) of the affine-pivot
# class), then after an affine change of coordinates e1 depends on at most
# two variables. [Proved already: deg e1 <= 2 (AP1-AP3); e1 homogeneous of
# any degree (Identity E + Gordan-Noether n=3). Here: deg e1 = 3.]
#
# If Conjecture E holds in all degrees then, with Theorem C, the ENTIRE
# affine-pivot class consists of Meng doublings of planar Keller maps plus an
# inert planar potential, so HC_4 restricted to that class is EQUIVALENT to
# JC_2 -- an exact localization of the planar Jacobian conjecture inside HC_4.
#
# Setup for d = 3. By Theorem D' the leading form of e1 is a cone, so WLOG
# it lies in C[x1,x2]. Constants and the choice of origin do not matter, so
#   e1 = c(x1,x2) + q0(x1,x2) + x3*lam(x1,x2) + (mu/2)*x3^2 + l0(x1,x2) + nu*x3
# with c cubic, q0 quadratic, lam linear, l0 linear, mu, nu constants.
# e1 is affinely 2-variable  <=>  its three partials are linearly dependent
# over C  <=>  all 3x3 minors of the coefficient matrix of (d1e1,d2e1,d3e1)
# vanish.

import itertools
import sys
import sympy as sp

x1, x2, x3 = X = sp.symbols('x1 x2 x3')

c_c = sp.symbols('c0:4')      # cubic in x1,x2
q_c = sp.symbols('q0:3')      # quadratic in x1,x2
lam_c = sp.symbols('t0:2')    # linear in x1,x2 (coefficient of x3)
mu, nu = sp.symbols('mu nu')
l_c = sp.symbols('l0:2')      # linear in x1,x2
params = list(c_c) + list(q_c) + list(lam_c) + [mu, nu] + list(l_c)

c = c_c[0]*x1**3 + c_c[1]*x1**2*x2 + c_c[2]*x1*x2**2 + c_c[3]*x2**3
q0 = q_c[0]*x1**2 + q_c[1]*x1*x2 + q_c[2]*x2**2
lam = lam_c[0]*x1 + lam_c[1]*x2
l0 = l_c[0]*x1 + l_c[1]*x2
e1 = c + q0 + x3*lam + mu*x3**2/2 + l0 + nu*x3

K = sp.Matrix([[sp.diff(e1, u, v) for v in X] for u in X])
g = sp.Matrix([sp.diff(e1, u) for u in X])
c2 = sp.expand((g.T * K.adjugate() * g)[0, 0])
eqs = [sp.expand(t) for t in sp.Poly(c2, *X).coeffs()]
print(f'(c2) yields {len(eqs)} polynomial equations in {len(params)} parameters')

# The "affinely 2-variable" locus: partials linearly dependent.
mons = [m for d in (0, 1, 2) for m in
        sorted({sp.prod(cc) for cc in itertools.combinations_with_replacement(X, d)}, key=str)]
mons = [sp.Integer(1)] + [m for m in mons if m != 1]
rows = []
for v in X:
    p = sp.Poly(sp.diff(e1, v), *X)
    rows.append([p.coeff_monomial(m) for m in mons])
Mat = sp.Matrix(rows)
minors = []
for cols in itertools.combinations(range(len(mons)), 3):
    m = Mat[:, list(cols)].det()
    if m != 0:
        minors.append(sp.expand(m))
minors = list({sp.factor(m) for m in minors})
print(f'{len(minors)} distinct nonzero 3x3 minors define the 2-variable locus')

# Radical membership test (Rabinowitsch): minor vanishes on V(c2)  <=>
# 1 in (c2-ideal) + (1 - y*minor).  Test each minor; report failures.
y = sp.Symbol('y')
fails = []
for i, m in enumerate(minors):
    G = sp.groebner(eqs + [1 - y*m], *(params + [y]), order='grevlex')
    inrad = (list(G.exprs) == [sp.Integer(1)])
    if not inrad:
        fails.append(m)
    print(f'  minor {i+1}/{len(minors)}: vanishes on V(c2)? {inrad}')

print()
if not fails:
    print('CONJECTURE E VERIFIED FOR d = 3: every cubic pivot coefficient with')
    print('vanishing bordered Hessian is affinely 2-variable.')
    print('=> in degree 3 the affine-pivot class is exactly the Meng-doubling class.')
else:
    print(f'CONJECTURE E FAILS (or is undecided) for d = 3: {len(fails)} minors do')
    print('not vanish on V(c2). A solution outside the 2-variable locus may exist;')
    print('extract and inspect it. Offending minors:')
    for m in fails[:5]:
        print('   ', m)
