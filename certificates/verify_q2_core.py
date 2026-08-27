# Theorem Q2 core verification (lead-researcher copy; agent runs its own).
#
# Q2: f = e0(x') + x4 e1(x') + (sigma/2) x4^2, sigma != 0, det Hess f == c
#     ==> det3(Hess3 etilde0 + s K) == c/sigma for ALL s  (pencil identity),
#         where etilde0 = e0 - e1^2/(2 sigma), K = Hess3 e1;
#     and any gradient collision of f descends to a collision of
#     grad'(etilde0 + lambda e1) for a single lambda in C, contradicting
#     de Bondt's HC_3.  Hence grad f is injective.
#
# The hand proof is in src/deg4_tree_notes.md sec 1a (Schur complement +
# "u sweeps C for fixed x'").  This script machine-checks the two load-bearing
# identities on instances (labelled consistency checks) and the one
# fully-symbolic step that is linear in coefficients.

import sys
import sympy as sp

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
from hc4lib import hessian, X4, X3

x1, x2, x3, x4 = X4

def check_instance(e0, e1, sigma, expect_const=None, tag=''):
    f = sp.expand(e0 + x4*e1 + sp.Rational(1, 2)*sigma*x4**2)
    H = hessian(f, list(X4))
    d = sp.expand(H.det(method='berkowitz'))
    et = sp.expand(e0 - e1**2/(2*sigma))
    K = hessian(e1, list(X3))
    s = sp.Symbol('s')
    pencil = hessian(et, list(X3)) + s*K
    dp = sp.expand(pencil.det(method='berkowitz'))
    # Schur identity: det Hess f == sigma * det3(Hess3 et + (x4 + e1/sigma) K)
    rhs = sp.expand(sigma*dp.subs(s, x4 + e1/sigma))
    assert sp.simplify(d - rhs) == 0, f'[{tag}] Schur identity FAILED'
    print(f'[{tag}] Schur identity OK; det Hess f = {sp.simplify(d) if d.is_number else "(nonconstant)"}')
    if expect_const is not None:
        assert sp.simplify(d - expect_const) == 0, f'[{tag}] det value mismatch'
        # pencil identity forced:
        assert sp.simplify(dp - expect_const/sigma) == 0, f'[{tag}] pencil identity FAILED'
        print(f'[{tag}] pencil det3(Hess etilde0 + sK) == {sp.simplify(dp)} for ALL s  OK')

# instance 1: K = 0 (e1 affine), etilde0 = x1 x2 + x1^4 + x3^2/2 (MA_3 solution)
sig = sp.Integer(3)
et0 = x1*x2 + x1**4 + x3**2/2
e1a = x3
check_instance(sp.expand(et0 + e1a**2/(2*sig)), e1a, sig, expect_const=-3, tag='K=0')

# instance 2: nontrivial pencil K = E11: e1 = x1^2/2 + x3
e1b = x1**2/2 + x3
check_instance(sp.expand(et0 + e1b**2/(2*sig)), e1b, sig, expect_const=-3, tag='K=E11')

# instance 3: random NON-solution (det not constant) - Schur identity must still hold
import random
random.seed(42)
e0r = sum(sp.Rational(random.randint(-5, 5), random.randint(1, 3)) *
          m for m in [x1**4, x1**3*x2, x1*x2*x3**2, x2**2*x3, x1*x3, x2**2, x1**2*x3**2])
e1r = sum(sp.Rational(random.randint(-5, 5), random.randint(1, 3)) *
          m for m in [x1**2, x1*x3, x2**2, x2, x3, sp.Integer(1)])
check_instance(sp.expand(e0r), sp.expand(e1r), sp.Rational(7, 2), tag='random')

# fully-symbolic (linear-in-coefficients) step: grad'(e1^2/(2 sigma)) == (e1/sigma) grad' e1
cs = sp.symbols('a0:10')
sigS = sp.Symbol('sigma')
mons2 = [x1**2, x1*x2, x1*x3, x2**2, x2*x3, x3**2, x1, x2, x3, sp.Integer(1)]
e1g = sum(c*m for c, m in zip(cs, mons2))
for w in X3:
    assert sp.expand(sp.diff(e1g**2/(2*sigS), w) - (e1g/sigS)*sp.diff(e1g, w)) == 0
print('symbolic collision-transfer identity OK (all quadratics e1, symbolic sigma)')

print('Q2 CORE CHECKS: ALL PASS')
