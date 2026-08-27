# HOSTILE REFEREE, instance attack #2  (key: x2)
# TARGET: Theorem A / branch AP0 (lemma_ledger AP0; main_result.tex, the
# "deg e1 <= 1 (K = 0)" paragraph):
#     f = e0(x1,x2,x3) + x4*x3     (after the claimed normalization e1 = x3),
#     determinant identity  <=>  (c0)  det2 Hess_(x1,x2) e0 == -c,
#     i.e. every x3-fiber of e0 is a planar polynomial with constant Hessian
#     determinant; the paper then invokes Dillen FIBERWISE and concludes
#     injectivity FOR EVERY DEGREE of e0.
#
# WHICH SHAPES ACTUALLY SATISFY (c0)?  We construct and test three:
#
#  FamA  (fixed isotropic direction):
#        e0 = x1*x2 + G(x1,x3) + x1*A(x3) + x2*B(x3) + C(x3)
#        Hess2 = [[G_11, 1],[1, 0]]  =>  det2 == -1  for ARBITRARY G,A,B,C.
#
#  FamB  (ROTATING isotropic direction -- the hostile case):
#        e0 = x2^2/2 + x2*(x1*s(x3) + t(x3)) + (s(x3)^2 - c)*x1^2/2
#             + x1*r(x3) + n(x3)
#        Hess2 = [[s^2 - c, s],[s, 1]],  det2 == -c  for ARBITRARY s,t,r,n.
#        Here BOTH diagonal entries are nonzero and the kernel direction of
#        Hess2 + c*(...) depends on x3, so this is NOT reducible to FamA by a
#        single (x3-independent) affine change -- exactly the configuration the
#        fiberwise argument must survive.  (Derived, not guessed: with
#        e0_{x2x2} = 1 the integrability of (c0) FORCES e0 = this shape.)
#
#  FamC  NEGATIVE CONTROL: e0 = x1*x2 + x1^2*x3 + x2^2*x3  has
#        det2 = 4*x3^2 - 1, NOT constant -- det Hess f must come out nonconstant.
#
# Per instance: det Hess f == c exactly (Berkowitz), plus injectivity by the
# saturated collision Groebner (4 Rabinowitsch saturations d_i != 0).
# FAIL-CLOSED.

import random
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
random.seed(4242)


def injective_saturated(f, label):
    q = sp.symbols('q1:5')
    d = sp.symbols('d1:5')
    z = sp.Symbol('zz')
    p = [q[i] + d[i] for i in range(4)]
    eqs = []
    for v in X:
        g = sp.diff(f, v)
        eqs.append(sp.expand(g.subs(dict(zip(X, p)), simultaneous=True)
                             - g.subs(dict(zip(X, q)), simultaneous=True)))
    gens = list(q) + list(d) + [z]
    for i in range(4):
        G = sp.groebner(eqs + [d[i] * z - 1], *gens, order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            print(f'  [{label}] d{i+1}-saturation NONEMPTY (basis size '
                  f'{len(G.exprs)})')
            return False
    return True


def rp(var, deg, dens=0.7):
    e = sp.Integer(0)
    for k in range(deg + 1):
        if random.random() < dens:
            e += random.choice([-3, -2, -1, 1, 2, 3]) * var**k
    return sp.expand(e)


def rp2(v, w, deg, dens=0.5):
    e = sp.Integer(0)
    for i in range(deg + 1):
        for j in range(deg + 1):
            if i + j <= deg and random.random() < dens:
                e += random.choice([-3, -2, -1, 1, 2, 3]) * v**i * w**j
    return sp.expand(e)


results = []

# ---------------------------- FamA ----------------------------------------
print('=== FamA: e0 = x1*x2 + G(x1,x3) + x1*A(x3) + x2*B(x3) + C(x3) ===')
nA = 6
for k in range(nA):
    deg = random.choice([3, 4, 5, 6])
    G = rp2(x1, x3, deg)
    A = rp(x3, random.choice([1, 2, 3]))
    B = rp(x3, random.choice([1, 2, 3]))
    C = rp(x3, random.choice([0, 2, 3]))
    e0 = sp.expand(x1 * x2 + G + x1 * A + x2 * B + C)
    f = sp.expand(e0 + x4 * x3)
    det = sp.expand(sp.hessian(f, X).det(method='berkowitz'))
    assert det == sp.Integer(1), f'FATAL FamA{k}: det Hess f = {det} != 1; e0={e0}'
    d2 = sp.expand(sp.Matrix([[sp.diff(e0, x1, 2), sp.diff(e0, x1, x2)],
                              [sp.diff(e0, x1, x2), sp.diff(e0, x2, 2)]]).det())
    assert d2 == sp.Integer(-1), f'FATAL FamA{k}: (c0) violated: det2 = {d2}'
    assert injective_saturated(f, f'FamA{k}'), (
        f'FATAL: AP0 FamA instance {k} is NON-INJECTIVE -- HC_4 COUNTEREXAMPLE?  '
        f'f = {f}')
    dt = sp.Poly(f, *X).total_degree()
    results.append(('FamA', k, dt))
    print(f'  FamA{k}: deg f = {dt}, det Hess == 1 PASS, (c0) PASS, injective PASS')

# ---------------------------- FamB ----------------------------------------
print()
print('=== FamB: rotating isotropic direction ===')
nB = 6
for k in range(nB):
    c = sp.Integer(random.choice([-3, -2, -1, 1, 2, 3]))
    s = rp(x3, random.choice([1, 2, 3]))
    if s == 0:
        s = x3
    t = rp(x3, random.choice([0, 1, 2, 3]))
    r = rp(x3, random.choice([0, 1, 2, 3]))
    n = rp(x3, random.choice([0, 2, 4]))
    e0 = sp.expand(x2**2 / 2 + x2 * (x1 * s + t)
                   + (s**2 - c) * x1**2 / 2 + x1 * r + n)
    assert all(co.q == 1 for co in sp.Poly(e0 * 2, x1, x2, x3).coeffs()), 'rational ok'
    f = sp.expand(e0 + x4 * x3)
    # (c0)
    d2 = sp.expand(sp.Matrix([[sp.diff(e0, x1, 2), sp.diff(e0, x1, x2)],
                              [sp.diff(e0, x1, x2), sp.diff(e0, x2, 2)]]).det())
    assert d2 == -c, f'FATAL FamB{k}: (c0) violated: det2 = {d2}, expected {-c}'
    det = sp.expand(sp.hessian(f, X).det(method='berkowitz'))
    assert det == c, f'FATAL FamB{k}: det Hess f = {det} != {c}; e0 = {e0}'
    # genuinely rotating: the Hess2 kernel-of-(Hess2 - ...) direction moves with x3
    assert sp.degree(sp.Poly(s, x3), x3) >= 1, 'want nonconstant s'
    assert injective_saturated(f, f'FamB{k}'), (
        f'FATAL: AP0 FamB instance {k} is NON-INJECTIVE -- HC_4 COUNTEREXAMPLE?  '
        f'f = {f}')
    dt = sp.Poly(f, *X).total_degree()
    results.append(('FamB', k, dt))
    print(f'  FamB{k}: deg f = {dt}, c = {c}, deg s = {sp.degree(sp.Poly(s,x3),x3)}, '
          f'det Hess == c PASS, (c0) PASS, injective PASS')

# ---------------------------- FamC (negative control) ----------------------
print()
e0 = x1 * x2 + x1**2 * x3 + x2**2 * x3
f = sp.expand(e0 + x4 * x3)
det = sp.expand(sp.hessian(f, X).det(method='berkowitz'))
assert det.free_symbols, f'FATAL(anti-vacuity): FamC det = {det} is constant?!'
print(f'=== FamC negative control: det Hess f = {det} (non-constant, as required) ===')

print()
print(f'RESULT rt_instances_x2_AP0: {len(results)} AP0 instances '
      f'({nA} FamA + {nB} FamB), degrees {sorted({r[2] for r in results})}, '
      f'ALL passed det-Hess and saturated-collision injectivity. '
      f'AP0 branch of Theorem A SURVIVED.')
