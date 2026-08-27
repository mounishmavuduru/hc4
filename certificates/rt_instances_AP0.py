# REFEREE INSTANCE ATTACK 2 (hostile referee, rt_instances_*):
# Theorem A, AP0 branch (gamma = 0, deg e1 <= 1, normalized e1 = x3):
#   f = e0(x1,x2,x3) + x4*x3,  constraint (c0): det2 Hess_(x1,x2) e0 == -c.
# The paper claims every x3-fiber of e0 is then a planar polynomial with
# constant Hessian determinant, and grad f is injective (Dillen fiberwise).
#
# SHAPE CHECK first (which e0-shapes satisfy (c0)?):
#   (S1) e0 = k*(x1+u(x3))*(x2+v(x3)) + g(x1+s(x3)) + H(x1,x3):
#        Hess2 = [[g''+H_11, k],[k, 0]]  => det2 == -k^2 for ARBITRARY u,v,s,g,H.
#   (S2) e0 = x1*x2 + x2*r(x1): det2 = -(1+r'(x1))^2, NOT constant unless
#        deg r <= 1  (negative control -- shape does NOT work).
#   (S3) e0 = x1*x2 + x2^2*x3: det2 = -1 + 2*x3*(g'' terms)... checked
#        explicitly below (negative control: det Hess f nonconstant).
# Then 10 random S1 instances (degrees up to 6, so beyond Theorem B range):
#   verify det Hess f == k^2 exactly, and injectivity by saturated collision
#   Groebner. FAIL-CLOSED.

import random
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
random.seed(1234)

def rpoly(var, maxdeg, dens=0.6):
    P = sp.Integer(0)
    for i in range(1, maxdeg + 1):
        if random.random() < dens:
            P += random.choice([-3, -2, -1, 1, 2, 3]) * var**i
    return P

# ---------- shape checks ----------
u3, v3, s3 = sp.symbols('u3 v3 s3')  # placeholders not needed; do symbolic jets
# S1 symbolic: generic u,v,s of degree 3, g of degree 4, H bivariate degree 4
uu = rpoly(x3, 3); vv = rpoly(x3, 3); ss = rpoly(x3, 3)
gsym = sp.symbols('g0:5'); g = sum(c * sp.Symbol('T')**i for i, c in enumerate(gsym))
T = sp.Symbol('T')
kk = sp.Symbol('kk')
Hmons = [x1**i * x3**j for i in range(5) for j in range(5) if i + j <= 4]
Hc = sp.symbols(f'hh0:{len(Hmons)}')
Hpoly = sum(c * m for c, m in zip(Hc, Hmons))
e0_S1 = sp.expand(kk * (x1 + uu) * (x2 + vv) + g.subs(T, x1 + ss) + Hpoly)
det2 = sp.expand(sp.diff(e0_S1, x1, 2) * sp.diff(e0_S1, x2, 2) - sp.diff(e0_S1, x1, x2)**2)
assert det2 == sp.expand(-kk**2), f'S1 shape check failed: {det2}'
print('S1 shape PASS: det2 Hess_(x1,x2) e0 == -k^2 for generic u,v,s (deg 3), g (deg 4),')
print('               H(x1,x3) (deg 4), symbolic k  -- a genuinely x3-mixed family.')

# S2 negative control
r = x1**2
e0_S2 = x1 * x2 + x2 * r
det2_S2 = sp.expand(sp.diff(e0_S2, x1, 2) * sp.diff(e0_S2, x2, 2) - sp.diff(e0_S2, x1, x2)**2)
assert det2_S2 != sp.expand(-(sp.Integer(1))), 'S2 unexpectedly constant'
assert sp.Poly(det2_S2, x1, x2, x3).total_degree() > 0
print(f'S2 negative control PASS: e0 = x1*x2 + x2*x1^2 gives det2 = {det2_S2} (nonconstant),')

f_S3 = sp.expand(x1 * x2 + x2**2 * x3 + x4 * x3)
d_S3 = sp.expand(sp.hessian(f_S3, X).det(method='berkowitz'))
assert sp.Poly(d_S3, *X).total_degree() > 0, 'S3 unexpectedly constant'
print(f'S3 negative control PASS: e0 = x1*x2 + x2^2*x3 gives det Hess f = {d_S3} (nonconstant).')
print('   => the (c0) constraint is tight: x2-degree >= 2 fiber coupling breaks constancy.')
print()

# ---------- 10 random S1 instances ----------
n_ok_det = 0
n_ok_inj = 0
N = 10
for j in range(N):
    k = [1, 2, -1, 3, -2, 1, 2, 5, -3, sp.Rational(1, 2)][j]
    uu = rpoly(x3, random.choice([2, 3]))
    vv = rpoly(x3, random.choice([2, 3]))
    ss = rpoly(x3, 2)
    gdeg = random.choice([3, 4, 5])
    gp = sum(random.choice([-2, -1, 0, 1, 2]) * (x1 + ss)**i for i in range(2, gdeg + 1))
    Hp = sp.Integer(0)
    for m in [x1**a * x3**b for a in range(7) for b in range(7) if 2 < a + b <= 6]:
        if random.random() < 0.2:
            Hp += random.choice([-3, -2, -1, 1, 2, 3]) * m
    e0 = sp.expand(k * (x1 + uu) * (x2 + vv) + gp + Hp)
    f = sp.expand(e0 + x4 * x3)
    ftot = sp.Poly(f, *X).total_degree()

    d = sp.expand(sp.hessian(f, X).det(method='berkowitz'))
    assert d == sp.expand(k**2), f'AP0 INSTANCE {j}: det Hess f = {d} != k^2 = {k**2}; f={f}'
    n_ok_det += 1

    p = sp.symbols('pp1:5'); q = sp.symbols('qq1:5'); w = sp.symbols('ww1:5')
    grads = [sp.diff(f, v) for v in X]
    eqs = [sp.expand(gg.subs(dict(zip(X, p))) - gg.subs(dict(zip(X, q)))) for gg in grads]
    sat = sp.expand(sum(w[i] * (p[i] - q[i]) for i in range(4)) - 1)
    G = sp.groebner(eqs + [sat], *(list(p) + list(q) + list(w)), order='grevlex')
    assert list(G.exprs) == [sp.Integer(1)], f'AP0 INSTANCE {j}: NON-INJECTIVE! f={f}'
    n_ok_inj += 1
    print(f'AP0 instance {j}: deg f = {ftot}, c = k^2 = {k**2}: det PASS, injective PASS')

print()
print(f'AP0 ATTACK: {N} random x3-mixed instances, {n_ok_det}/{N} determinant checks,')
print(f'{n_ok_inj}/{N} injectivity Groebner checks passed (degrees up to 7 > 4).')
