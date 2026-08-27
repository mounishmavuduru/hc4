# REFEREE INSTANCE ATTACK 1 (hostile referee, rt_instances_*):
# Theorem A classified family (main_result.tex, end of proof of thm:pivot):
#   f = a(x1,x2) + x3*(-kappa*x1 + beta(xi)) + x4*xi,  xi = x2 + x1^2/2,
#   claim: det Hess f == kappa^2 and grad f injective for ARBITRARY a, beta.
# We generate 12 random members with deg a up to 6 and beta up to degree 2 in xi
# (so f has degree up to 6 -- deliberately BEYOND degree 4, testing the
# all-degree claim of Theorem A), and for each:
#   (1) verify det Hess f == kappa^2 EXACTLY (berkowitz, expand);
#   (2) verify injectivity by the saturated collision Groebner basis:
#       grad f(p) - grad f(q) = 0,  sum_i w_i (p_i - q_i) = 1  ==> basis {1}.
# FAIL-CLOSED: any assertion failure aborts with nonzero exit.

import random
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
xi = x2 + sp.Rational(1, 2) * x1**2

random.seed(20260810)

def random_a(maxdeg, force_top=True):
    mons = [x1**i * x2**j for i in range(maxdeg + 1) for j in range(maxdeg + 1)
            if 0 < i + j <= maxdeg]
    a = sp.Integer(0)
    chosen = 0
    for m in mons:
        if random.random() < 0.35:
            c = random.choice([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])
            a += c * m
            chosen += 1
    if force_top:
        # force a genuine top-degree monomial so deg a = maxdeg
        i = random.randint(0, maxdeg)
        a += random.choice([-3, -2, -1, 1, 2, 3]) * x1**i * x2**(maxdeg - i)
    return sp.expand(a)

def random_beta():
    # beta(t) = b0 + b1 t + b2 t^2, b2 forced nonzero in most draws
    b0 = random.randint(-4, 4)
    b1 = random.randint(-4, 4)
    b2 = random.choice([-3, -2, -1, 1, 2, 3])
    return (b0, b1, b2)

n_ok_det = 0
n_ok_inj = 0
N = 12
for k in range(N):
    kappa = [1, -1, 2, -2, 3, sp.Rational(1, 2), -3, 5, 1, 2, -1, 4][k % 12]
    dega = [6, 6, 6, 5, 5, 4, 6, 3, 6, 5, 6, 4][k]
    a = random_a(dega)
    b0, b1, b2 = random_beta()
    if k in (5, 7):
        b2 = 0  # a couple of low-degree-beta members too
    beta = b0 + b1 * xi + b2 * xi**2
    f = sp.expand(a + x3 * (-kappa * x1 + beta) + x4 * xi)
    ftot = sp.Poly(f, *X).total_degree()

    # (1) exact determinant check
    H = sp.hessian(f, X)
    d = sp.expand(H.det(method='berkowitz'))
    assert d == sp.expand(kappa**2), (
        f'INSTANCE {k}: det Hess f != kappa^2!  det={d}, kappa={kappa}, f={f}')
    n_ok_det += 1

    # (2) saturated collision Groebner
    p = sp.symbols('pp1:5')
    q = sp.symbols('qq1:5')
    w = sp.symbols('ww1:5')
    grads = [sp.diff(f, v) for v in X]
    eqs = [sp.expand(g.subs(dict(zip(X, p))) - g.subs(dict(zip(X, q))))
           for g in grads]
    sat = sp.expand(sum(w[i] * (p[i] - q[i]) for i in range(4)) - 1)
    G = sp.groebner(eqs + [sat], *(list(p) + list(q) + list(w)), order='grevlex')
    assert list(G.exprs) == [sp.Integer(1)], (
        f'INSTANCE {k}: NON-INJECTIVE member of classified family! f={f}')
    n_ok_inj += 1
    print(f'instance {k:2d}: deg f = {ftot}, kappa = {kappa}, '
          f'beta deg(in xi) = {2 if b2 else (1 if b1 else 0)}: '
          f'det==kappa^2 PASS, injective PASS')

print()
print(f'CLASSIFIED-FAMILY ATTACK: {N} random instances, '
      f'{n_ok_det}/{N} exact determinant checks passed, '
      f'{n_ok_inj}/{N} saturated collision Groebner injectivity checks passed.')
print('Degrees of f used include 6 (beyond Theorem B range): Theorem A all-degree')
print('claim survives on these instances.')
