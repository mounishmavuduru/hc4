# HOSTILE REFEREE, instance attack #1  (key: x2)
# TARGET: Theorem A (main_result.tex thm:pivot), the *complete classification*
# claimed at the end of its proof, gamma = 0 / deg e1 = 2 branch:
#
#     f = a(x1,x2) + x3*( -kappa*x1 + beta(xi) ) + x4*xi,   xi = x2 + x1^2/2,
#     CLAIM (i)  det Hess f == kappa^2   for ARBITRARY a and beta;
#     CLAIM (ii) grad f is injective.
#
# The claim is made for ALL degrees, so we deliberately push deg a up to 6 and
# beta up to degree 2 in xi (=> deg f up to 6), far outside Theorem B's range.
#
# Checks per instance (exact rational arithmetic only):
#   (1) det Hess f == kappa^2, Berkowitz determinant, fully expanded;
#   (2) injectivity: saturated collision Groebner.  Collision variables
#       q = (q1..q4), difference d = p - q; p != q  <=>  some d_i != 0, so we
#       run 4 Rabinowitsch saturations (d_i*z - 1) and require basis {1} in all
#       four.  This is *equivalent* to p != q and is much cheaper than a
#       4-multiplier saturation.
#   (3) NEGATIVE CONTROLS (anti-vacuity): the same machinery must DETECT
#       failures -- we perturb b off the classified shape and assert that
#       det Hess is then NOT constant; and we feed a known non-injective
#       gradient map (a doubled planar non-Keller collision) and assert the
#       Groebner test returns a nontrivial basis.
#
# FAIL-CLOSED: every check is an assert; nonzero exit on any failure.

import random
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
XI = x2 + sp.Rational(1, 2) * x1**2

random.seed(20260810 + 77)


def injective_saturated(f, label, expect_injective=True):
    """Saturated collision Groebner. Returns True iff grad f is injective."""
    q = sp.symbols('q1:5')
    d = sp.symbols('d1:5')
    z = sp.Symbol('zz')
    p = [q[i] + d[i] for i in range(4)]
    grads = [sp.diff(f, v) for v in X]
    eqs = []
    for g in grads:
        gp = g.subs(dict(zip(X, p)), simultaneous=True)
        gq = g.subs(dict(zip(X, q)), simultaneous=True)
        eqs.append(sp.expand(gp - gq))
    gens = list(q) + list(d) + [z]
    for i in range(4):
        G = sp.groebner(eqs + [d[i] * z - 1], *gens, order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            print(f'  [{label}] saturation d{i+1} != 0 has SOLUTIONS: '
                  f'basis size {len(G.exprs)}')
            return False
    return True


def rand_a(maxdeg):
    """random a(x1,x2), genuine total degree maxdeg."""
    a = sp.Integer(0)
    for i in range(maxdeg + 1):
        for j in range(maxdeg + 1):
            if 0 < i + j <= maxdeg and random.random() < 0.4:
                a += random.choice([-4, -3, -2, -1, 1, 2, 3, 4]) * x1**i * x2**j
    i = random.randint(0, maxdeg)
    a += random.choice([-3, -2, -1, 1, 2, 3]) * x1**i * x2**(maxdeg - i)
    return sp.expand(a)


def rand_beta(deg):
    t = sp.Symbol('t')
    b = sum(random.choice([-3, -2, -1, 0, 1, 2, 3]) * t**k for k in range(deg))
    b += random.choice([-3, -2, -1, 1, 2, 3]) * t**deg
    return sp.Lambda(t, sp.expand(b))


N = 12
ok_det = ok_inj = 0
records = []
for k in range(N):
    kappa = sp.Integer(random.choice([-3, -2, -1, 1, 2, 3, 5]))
    adeg = random.choice([2, 3, 4, 5, 6, 6])
    bdeg = random.choice([0, 1, 1, 2, 2])          # degree of beta in xi
    a = rand_a(adeg)
    beta = rand_beta(bdeg)
    b = -kappa * x1 + beta(XI)
    f = sp.expand(a + x3 * b + x4 * XI)
    dtot = sp.Poly(f, *X).total_degree()

    H = sp.hessian(f, X)
    det = sp.expand(H.det(method='berkowitz'))
    assert det == sp.expand(kappa**2), (
        f'FATAL: instance {k}: det Hess f = {det} != kappa^2 = {kappa**2}; f = {f}')
    ok_det += 1

    inj = injective_saturated(f, f'A{k}')
    assert inj, f'FATAL: instance {k} of the CLASSIFIED FAMILY is NON-INJECTIVE: f = {f}'
    ok_inj += 1
    records.append((k, dtot, int(kappa), adeg, bdeg))
    print(f'  A{k:2d}: deg f = {dtot}, kappa = {kappa}, deg a = {adeg}, '
          f'deg beta = {bdeg}  ->  det == kappa^2 PASS, injective PASS')

print()
print(f'[1] classified family: {ok_det}/{N} determinant checks, '
      f'{ok_inj}/{N} injectivity checks PASS.')
print('    degrees of f attained:', sorted({r[1] for r in records}))

# ---------------------------------------------------------------------------
# NEGATIVE CONTROL 1: perturbing b off the PDE solution set must destroy
# constancy of det Hess f  (so the determinant check above is not vacuous).
bad = 0
for k in range(6):
    kappa = sp.Integer(random.choice([-2, -1, 1, 2, 3]))
    a = rand_a(random.choice([3, 4]))
    beta = rand_beta(random.choice([1, 2]))
    # perturb: add a monomial in x1,x2 that is NOT a function of xi
    pert = random.choice([x1**2, x1 * x2, x2**2, x1**3, x1 * x2**2])
    b = -kappa * x1 + beta(XI) + pert
    f = sp.expand(a + x3 * b + x4 * XI)
    det = sp.expand(sp.hessian(f, X).det(method='berkowitz'))
    assert det.free_symbols, (
        f'FATAL(anti-vacuity): perturbed b still gives CONSTANT det = {det}; '
        f'b = {b}  -- either the classification is incomplete or the test is blind')
    bad += 1
print(f'[2] negative control (perturbed b): {bad}/6 gave non-constant det Hess '
      f'as required -- the determinant test discriminates.')

# NEGATIVE CONTROL 2: the injectivity machinery must be able to say NO.
# Meng doubling of the *non-injective* planar map (u,v) |-> (u + v^2 ...)?  Use
# a genuine explicit gradient collision instead: f = (x1^2)/2 ... take
# g(x) = x1^3/3 - x1  (grad has g'(x1) = x1^2 - 1, so x1 = 1 and x1 = -1 collide).
f_bad = x1**3 / 3 - x1 + x2**2 / 2 + x3**2 / 2 + x4**2 / 2
assert not injective_saturated(f_bad, 'NEG', expect_injective=False), (
    'FATAL(anti-vacuity): the collision Groebner declared a KNOWN non-injective '
    'gradient map injective -- the injectivity test is broken.')
print('[3] negative control (known collision at x1 = +-1): correctly DETECTED. '
      'Injectivity test is not vacuous.')

print()
print('RESULT rt_instances_x2_A_classified: Theorem A classified family '
      'SURVIVED all instance attacks (degrees up to 6).')
