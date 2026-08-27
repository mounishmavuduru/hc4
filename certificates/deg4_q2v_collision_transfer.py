# deg4_q2v_collision_transfer.py -- branch q2v, TASK Q2V part (2).
#
# Collision-transfer step of Theorem Q2 (scaffold sec. 1a, lines 36-44), for
#   f = e0(x') + x4*e1(x') + (sigma/2)*x4^2,   x' = (x1,x2,x3),  sigma != 0.
#
# Verified here, fail-closed, with FULLY SYMBOLIC generic coefficients
# (e0: all monomials deg <= 4 in x', 35 symbols; e1: deg <= 2, 10 symbols),
# symbolic sigma, and a fully symbolic point p = (p1,p2,p3,p4):
#
#  (2a) [PROOF, linear in coefficients]  grad' f = grad' e0 + x4 * grad' e1
#       and  df/dx4 = e1(x') + sigma*x4  = sigma*(x4 + e1(x')/sigma).
#       Both sides linear in the coefficients of e0, e1 => symbolic identity
#       check is a proof for all such e0, e1.
#       CONSEQUENCE (uses sigma != 0): grad f(p) = grad f(q) forces, from the
#       4th component, sigma*(p4 + e1(p')/sigma) = sigma*(q4 + e1(q')/sigma),
#       i.e. a COMMON pivot value lambda := p4 + e1(p')/sigma = q4 + e1(q')/sigma.
#
#  (2b) [hand argument + symbolic check + instances]
#           grad' f = grad' etilde0 + (x4 + e1/sigma) * grad' e1,
#       with etilde0 := e0 - e1^2/(2 sigma).
#       NOTE (discrepancy with the task wording): this identity is linear in
#       the coefficients of e0 but QUADRATIC in the coefficients of e1 --
#       both grad'(e1^2/(2 sigma)) and (e1/sigma)*grad' e1 are quadratic in
#       them; the two quadratic contributions cancel.  HAND ARGUMENT
#       (product rule): grad(e1^2) = 2*e1*grad(e1), so
#       grad' etilde0 = grad' e0 - (e1/sigma) grad' e1 and adding
#       (x4 + e1/sigma) grad' e1 leaves grad' e0 + x4 grad' e1 = grad' f
#       by (2a).  Per project rules the fully-symbolic expand==0 below is
#       recorded as a CONSISTENCY CHECK backing this hand argument, and it is
#       re-verified on 3 random dense rational instances.
#
#  (2c) [PROOF given (2b); the substitution step itself is linear in the
#       coefficients of etilde0+lambda*e1 once lambda is a free symbol]
#       With lambda a FREE symbol, define  h_lambda := etilde0 + lambda*e1.
#       Then   [grad' h_lambda](x') |_{lambda -> x4 + e1(x')/sigma}
#            = grad' etilde0 + (x4 + e1/sigma) grad' e1  = grad' f.
#       Evaluated at a symbolic point p with lam_p := p4 + e1(p')/sigma:
#            grad' f (p) == [grad' h_lambda](p') |_{lambda -> lam_p}.
#       COLLISION TRANSFER (logic, using (2a)+(2c)): if grad f(p) = grad f(q)
#       then lam_p = lam_q =: lambda (common, by (2a) and sigma != 0), and the
#       first three components give
#            grad'(etilde0 + lambda e1)(p') = grad'(etilde0 + lambda e1)(q'),
#       i.e. a gradient collision of the 3-variable polynomial
#       etilde0 + lambda*e1 for the SINGLE scalar lambda.  (The pencil
#       identity -- certificate deg4_q2v_pencil_forcing.py -- makes this
#       polynomial have constant Hessian determinant c/sigma != 0, so de
#       Bondt's HC_3 (L-HC3) forces p' = q', and then the 4th component gives
#       p4 = q4.)
#
# All arithmetic exact.  No numerics.

import sys
import itertools
import random

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
import sympy as sp
from hc4lib import hessian  # noqa: F401  (import kept for env sanity)

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Xp = (x1, x2, x3)
sigma = sp.Symbol('sigma')
lam = sp.Symbol('lambda')
p1, p2, p3, p4 = sp.symbols('p1 p2 p3 p4')


def generic_poly(vars_, degs, prefix):
    p = sp.Integer(0)
    idx = 0
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            c = sp.Symbol(f'{prefix}{idx}')
            idx += 1
            t = c
            for i in m:
                t *= vars_[i]
            p += t
    return p


def rand_poly(vars_, degs, rnd):
    p = sp.Integer(0)
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            c = sp.Rational(rnd.choice([q for q in range(-9, 10) if q != 0]),
                            rnd.randint(1, 4))
            t = c
            for i in m:
                t *= vars_[i]
            p += t
    return sp.expand(p)


def run_checks(e0, e1, tag, symbolic_proof_labels):
    f = e0 + x4*e1 + sp.Rational(1, 2)*sigma*x4**2
    g = [sp.diff(e1, w) for w in Xp]
    et0 = e0 - e1**2/(2*sigma)

    # ---- (2a) grad' f decomposition and 4th component: LINEAR in coeffs ----
    for i, w in enumerate(Xp):
        d = sp.expand(sp.diff(f, w) - (sp.diff(e0, w) + x4*g[i]))
        assert d == 0, f'[{tag}] (2a) grad\'-component {i+1} FAILED'
    d4 = sp.expand(sp.diff(f, x4) - (e1 + sigma*x4))
    assert d4 == 0, f'[{tag}] (2a) df/dx4 FAILED'
    d4b = sp.expand(sigma*(sp.diff(f, x4) - sigma*(x4 + e1/sigma)))
    assert d4b == 0, f'[{tag}] (2a) df/dx4 pivot form FAILED'

    # ---- (2b) grad' f == grad' etilde0 + (x4 + e1/sigma) grad' e1 ----
    for i, w in enumerate(Xp):
        d = sp.expand(sigma*(sp.diff(f, w)
                             - (sp.diff(et0, w) + (x4 + e1/sigma)*g[i])))
        assert d == 0, f'[{tag}] (2b) transfer component {i+1} FAILED'

    # ---- (2c) evaluation at a symbolic point p, lambda substitution ----
    subp = {x1: p1, x2: p2, x3: p3, x4: p4}
    lam_p = p4 + e1.subs({x1: p1, x2: p2, x3: p3})/sigma
    h = et0 + lam*e1                      # 3-variable polynomial, scalar lambda
    for i, w in enumerate(Xp):
        lhs = sp.diff(f, w).subs(subp)
        rhs = sp.diff(h, w).subs({x1: p1, x2: p2, x3: p3}).subs(lam, lam_p)
        d = sp.expand(sigma*(lhs - rhs))
        assert d == 0, f'[{tag}] (2c) pointwise transfer component {i+1} FAILED'
    print(f'[{tag}] (2a),(2b),(2c) all hold' + symbolic_proof_labels)


# ---------- fully symbolic generic coefficients (35 + 10 symbols) ----------
e0g = generic_poly(Xp, [0, 1, 2, 3, 4], 'a')
e1g = generic_poly(Xp, [0, 1, 2], 'b')
run_checks(e0g, e1g, 'FULLY SYMBOLIC',
           '  [(2a): PROOF (linear in coeffs); (2b),(2c): symbolic identity, '
           'quadratic in e1-coeffs -> recorded as consistency backing the '
           'product-rule hand argument]')

# ---------- 3 random dense rational instances (symbolic sigma, p) ----------
rnd = random.Random(20260810)
for k in range(3):
    e0r = rand_poly(Xp, [2, 3, 4], rnd)
    e1r = rand_poly(Xp, [0, 1, 2], rnd)
    run_checks(e0r, e1r, f'instance {k+1}/3', '')

print('COLLISION TRANSFER: any grad f(p) = grad f(q) forces the common pivot '
      'lambda = p4 + e1(p\')/sigma = q4 + e1(q\')/sigma (sigma != 0), and then '
      'grad\'(etilde0 + lambda e1)(p\') = grad\'(etilde0 + lambda e1)(q\').')
print('deg4_q2v_collision_transfer: ALL CHECKS PASS')
