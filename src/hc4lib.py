# hc4lib -- exact tools for the HC_4 project.
# All arithmetic exact (QQ). No numerics anywhere.

import itertools
import sympy as sp

X4 = sp.symbols('x1 x2 x3 x4')
X3 = X4[:3]


def hessian(f, vars_):
    return sp.Matrix([[sp.diff(f, a, b) for b in vars_] for a in vars_])


def graded_parts(p, vars_):
    """dict: total degree -> homogeneous component (expression)."""
    d = {}
    P = sp.Poly(p, *vars_)
    for monom, coeff in P.terms():
        deg = sum(monom)
        term = coeff
        for v, e in zip(vars_, monom):
            term *= v**e
        d[deg] = d.get(deg, 0) + term
    return {k: sp.expand(v) for k, v in d.items()}


def graded_part(p, vars_, m):
    return graded_parts(p, vars_).get(m, sp.Integer(0))


def det_hess(f, vars_):
    return sp.expand(hessian(f, vars_).det(method='berkowitz'))


def random_form(vars_, deg, rng, dens=1):
    """dense random homogeneous form with rational coefficients from rng list."""
    mons = [m for m in itertools.combinations_with_replacement(range(len(vars_)), deg)]
    f = sp.Integer(0)
    for m in mons:
        c = rng.pop() if rng else 0
        term = sp.Integer(c)
        for i in m:
            term *= vars_[i]
        f += term
    return sp.expand(f)


def collision_system(f, vars_):
    """T4 system for deg <= 4: Hess f(m) v + grad f4(v) = 0 (components), in
    unknowns m1..m4, v1..v4.  Caller adds v != 0 saturation."""
    n = len(vars_)
    ms = sp.symbols(f'm1:{n+1}')
    vs = sp.symbols(f'v1:{n+1}')
    H = hessian(f, vars_)
    Hm = H.subs(dict(zip(vars_, ms)))
    f4 = graded_part(f, vars_, 4)
    g4 = [sp.diff(f4, w).subs(dict(zip(vars_, vs))) for w in vars_]
    eqs = [sp.expand(sum(Hm[i, j]*vs[j] for j in range(n)) + g4[i]) for i in range(n)]
    return eqs, ms, vs


# ----------------------------- self-tests -----------------------------------

def _selftest():
    import random
    random.seed(20260808)

    def coeffs(k):
        return [sp.Rational(random.randint(-9, 9), random.randint(1, 4)) for _ in range(k)]

    x1, x2, x3, x4 = X4

    # --- T0/T1 spot checks on a random quartic in 4 vars ---
    f4 = random_form(X3, 4, coeffs(15))          # leading form x4-free (GN-normalized shape)
    f3 = random_form(X4, 3, coeffs(20))
    f2 = random_form(X4, 2, coeffs(10))
    f = f4 + f3 + f2
    D = det_hess(f, X4)
    parts = graded_parts(D, X4)
    # T0: degree-0 piece is det Hess f2
    assert sp.expand(parts.get(0, 0) - det_hess(f2, X4)) == 0
    # T1: degree-8 piece is det Hess f4 (which vanishes identically here: f4 is x4-free)
    assert sp.expand(parts.get(8, 0) - det_hess(f4, X4)) == 0
    assert det_hess(f4, X4) == 0
    print('T0, T1 spot checks PASS')

    # --- T3 spot check: [det]_7 = det3(Hess3 f4) * (f3)_{x4x4} ---
    lhs = parts.get(7, 0)
    rhs = sp.expand(hessian(f4, X3).det(method='berkowitz') * sp.diff(f3, x4, x4))
    assert sp.expand(lhs - rhs) == 0
    print('T3 spot check PASS')

    # --- T4: FULLY SYMBOLIC proof for a generic quartic in 4 variables ---
    # the identity is linear in the coefficients, so symbolic verification is a proof.
    mons = list(itertools.combinations_with_replacement(range(4), 4)) \
         + list(itertools.combinations_with_replacement(range(4), 3)) \
         + list(itertools.combinations_with_replacement(range(4), 2))
    cs = sp.symbols(f'c0:{len(mons)}')
    fgen = sp.Integer(0)
    for c, m in zip(cs, mons):
        term = c
        for i in m:
            term *= X4[i]
        fgen += term
    ms = sp.symbols('mm1:5')
    vs = sp.symbols('vv1:5')
    grad = [sp.diff(fgen, w) for w in X4]
    H = hessian(fgen, X4)
    f4gen = sp.Integer(0)
    for c, m in zip(cs, mons):
        if len(m) == 4:
            term = c
            for i in m:
                term *= X4[i]
            f4gen += term
    sub_p = dict(zip(X4, [a + b for a, b in zip(ms, vs)]))
    sub_m = dict(zip(X4, [a - b for a, b in zip(ms, vs)]))
    sub_mm = dict(zip(X4, ms))
    sub_v = dict(zip(X4, vs))
    for i in range(4):
        lhs_i = sp.expand(grad[i].subs(sub_p) - grad[i].subs(sub_m))
        rhs_i = sp.expand(2*(sum(H[i, j].subs(sub_mm)*vs[j] for j in range(4))
                             + sp.diff(f4gen, X4[i]).subs(sub_v)))
        assert sp.expand(lhs_i - rhs_i) == 0
    print('T4 symbolic proof PASS (identity holds for ALL quartics)')

    print('hc4lib selftest: ALL PASS')


if __name__ == '__main__':
    _selftest()
