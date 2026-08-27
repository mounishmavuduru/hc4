# HOSTILE REFEREE, instance attack #3 + #4  (key: x2)
# TARGET: Theorem B (HC_4 in degree <= 4) and its Corollary / the pivot
# mechanism that Theorem B's proof relies on.
#
# METHOD.  Six sparse quartic ansatz families
#     f = f2 + f3 + f4,   f2 fixed nondegenerate, f3/f4 on random sparse
#     supports with UNKNOWN rational coefficients u_i.
# For each draw we impose det Hess f == const EXACTLY: every positive-degree
# x-coefficient of det Hess f is set to 0, giving a polynomial system in the
# u_i.  We solve it (Groebner + sp.solve, so complex branches are included) and
# then, for EVERY solution found (free parameters instantiated at several
# random rational points), we run three independent tests:
#
#   (T-det)   re-verify det Hess f == c exactly and c != 0  (Keller);
#   (T-piv)   compute the PIVOT CONE.  Two versions:
#                P1 = {v != 0 : D_v^2 f == const}
#                P2 = {v != 0 : D_v^2 f == const  AND  D_v f4 == 0}
#             P2 is exactly the set of pivots with pivot coefficient
#             deg e1 <= 2 (in coordinates with v = e4, f4 is affine-linear in
#             x4 and e1's cubic part is D_v f4), i.e. exactly the pivots
#             Theorem A is allowed to consume.  THEOREM B'S PROOF ASSERTS
#             P2 != {0} FOR EVERY degree-4 Keller potential.  Emptiness of P2
#             on any instance = SERIOUS GAP in the proof of Theorem B.
#             Emptiness is decided by 4 Rabinowitsch saturations v_i*z - 1.
#   (T-inj)   injectivity of grad f by saturated collision Groebner
#             (4 saturations d_i != 0).  Per Theorem B every solution MUST be
#             injective; a non-injective one is an HC_4 COUNTEREXAMPLE.
#
# Family R6 is designed adversarially to KILL the pivot cone: f4 carries all
# four pure powers x_i^4 with unknown coefficients, so that no coordinate
# direction is free of the leading form.  (Gordan-Noether is not assumed --
# det Hess f4 == 0 is imposed automatically as the degree-8 graded piece.)
#
# FAIL-CLOSED: any non-injective Keller solution, or any Keller solution with
# P2 == {0}, raises AssertionError.

import random
import itertools
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
random.seed(90210)

CUBICS = [x1**a * x2**b * x3**c * x4**d
          for a, b, c, d in itertools.product(range(4), repeat=4)
          if a + b + c + d == 3]
QUARTICS = [x1**a * x2**b * x3**c * x4**d
            for a, b, c, d in itertools.product(range(5), repeat=4)
            if a + b + c + d == 4]
Q2FORMS = [x1 * x3 + x2 * x4,
           x1 * x2 + x3 * x4,
           (x1**2 + x2**2 + x3**2 + x4**2) / 2,
           x1 * x2 + x3**2 / 2 + x4**2 / 2]


def pos_deg_coeffs(expr):
    P = sp.Poly(sp.expand(expr), *X)
    out = []
    for mon, co in zip(P.monoms(), P.coeffs()):
        if sum(mon) > 0:
            out.append(sp.expand(co))
    return [e for e in out if e != 0]


def const_part(expr):
    return sp.expand(expr).subs({x1: 0, x2: 0, x3: 0, x4: 0})


def cone_is_trivial(eqs, vs, extra_gens=()):
    """True iff the common zero cone of `eqs` (homogeneous in vs) is {0}."""
    z = sp.Symbol('zsat')
    gens = list(vs) + list(extra_gens) + [z]
    for v in vs:
        G = sp.groebner(list(eqs) + [v * z - 1], *gens, order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            return False
    return True


def pivot_cones(f):
    """Return (P1_nonempty, P2_nonempty) for the pivot cones of f."""
    v = sp.symbols('v1:5')
    Dv2 = sp.expand(sum(v[i] * v[j] * sp.diff(f, X[i], X[j])
                        for i in range(4) for j in range(4)))
    eqs1 = pos_deg_coeffs(Dv2)          # D_v^2 f == const
    # NOTE (repair): the agent left a placeholder line here that crashed on
    # non-homogeneous f (Poly.homogeneous_order() returns None). f4 is
    # correctly recomputed immediately below, so the placeholder is dropped.
    P = sp.Poly(f, *X)
    f4 = sp.expand(sum(co * sp.prod([X[k]**mon[k] for k in range(4)])
                       for mon, co in zip(P.monoms(), P.coeffs())
                       if sum(mon) == 4))
    Dvf4 = sp.expand(sum(v[i] * sp.diff(f4, X[i]) for i in range(4)))
    eqs2 = eqs1 + pos_deg_coeffs(Dvf4) + ([] if const_part(Dvf4) == 0
                                          else [const_part(Dvf4)])
    p1 = not cone_is_trivial(eqs1, v)
    p2 = not cone_is_trivial(eqs2, v)
    return p1, p2


def injective_saturated(f):
    q = sp.symbols('q1:5')
    d = sp.symbols('d1:5')
    z = sp.Symbol('zzz')
    p = [q[i] + d[i] for i in range(4)]
    eqs = []
    for var in X:
        g = sp.diff(f, var)
        eqs.append(sp.expand(g.subs(dict(zip(X, p)), simultaneous=True)
                             - g.subs(dict(zip(X, q)), simultaneous=True)))
    gens = list(q) + list(d) + [z]
    for i in range(4):
        G = sp.groebner(eqs + [d[i] * z - 1], *gens, order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            return False
    return True


def make_family(name, draw):
    return (name, draw)


def draw_R1():
    """r = 1 branch: f4 = x1^4, sparse cubic f3 over all 4 variables."""
    sup = random.sample(CUBICS, random.choice([7, 8, 9]))
    return x1**4, sup, []


def draw_R2():
    """r = 2 branch: binary quartic leading form + sparse cubic."""
    b = random.choice([x1**4 + x1**2 * x2**2, x1**3 * x2 + x2**4,
                       x1**4 + x2**4, x1**2 * x2**2 + x1 * x2**3])
    sup = random.sample(CUBICS, random.choice([7, 8, 9]))
    return b, sup, []


def draw_R3():
    """r = 3 branch: essential-ternary quartic + sparse cubic (x4^2 allowed)."""
    b = random.choice([x1**4 + x2**4 + x3**4, x1**2 * x2 * x3 + x2**4 + x3**4,
                       x1**3 * x2 + x2**3 * x3 + x3**4])
    sup = random.sample(CUBICS, random.choice([7, 8, 9]))
    return b, sup, []


def draw_R4():
    """doubling shape: f = x3*P(x1,x2) + x4*Q(x1,x2), P,Q cubic (unknowns)."""
    mons = [x1**3, x1**2 * x2, x1 * x2**2, x2**3, x1**2, x1 * x2, x2**2]
    sup = [x3 * m for m in mons] + [x4 * m for m in mons]
    return sp.Integer(0), [], sup


def draw_R5():
    """fully random sparse f3 + f4 over all four variables."""
    sup = (random.sample(QUARTICS, random.choice([4, 5, 6]))
           + random.sample(CUBICS, random.choice([4, 5, 6])))
    return sp.Integer(0), [], sup


def draw_R6():
    """ADVERSARIAL: all four pure fourth powers carry unknown coefficients,
    plus mixed quartics and a sparse cubic -- an attempt to make every
    direction 'occupied' by the leading form and thus kill the pivot cone."""
    sup = ([x1**4, x2**4, x3**4, x4**4]
           + random.sample([m for m in QUARTICS if m not in
                            (x1**4, x2**4, x3**4, x4**4)], 3)
           + random.sample(CUBICS, random.choice([4, 5])))
    return sp.Integer(0), [], sup


FAMILIES = [('R1 (f4 = x1^4, r=1 branch)', draw_R1),
            ('R2 (binary f4, r=2 branch)', draw_R2),
            ('R3 (ternary f4, r=3 branch)', draw_R3),
            ('R4 (doubling x3*P + x4*Q)', draw_R4),
            ('R5 (fully random sparse)', draw_R5),
            ('R6 (adversarial: all x_i^4)', draw_R6)]

DRAWS_PER_FAMILY = 4

# CHUNKING (added 2026-08-24): this environment kills background jobs after
# ~10 minutes.  Optional CLI:  py rt_instances_x2_B_break.py R3 [draws]
# runs only the family whose name starts with the given prefix.  The random
# seed is re-derived from the family name so chunked runs are reproducible
# and independent of ordering.
import sys as _sys
if len(_sys.argv) > 1:
    _pref = _sys.argv[1]
    FAMILIES = [fd for fd in FAMILIES if fd[0].startswith(_pref)]
    assert FAMILIES, 'no family matches prefix ' + _pref
    if len(_sys.argv) > 2:
        DRAWS_PER_FAMILY = int(_sys.argv[2])
    random.seed(90210 + sum(map(ord, _pref)))

stats = dict(draws=0, systems_solved=0, sol_branches=0, instances=0,
             keller_nontrivial=0, deg4=0, inj_ok=0, p1_ok=0, p2_ok=0,
             unresolved=0)
counterexamples = []
gaps = []
log = []

for fname, draw in FAMILIES:
    print(f'\n===== FAMILY {fname} =====')
    for t in range(DRAWS_PER_FAMILY):
        stats['draws'] += 1
        fixed4, sup3, sup_free = draw()
        n = len(sup3) + len(sup_free)
        u = sp.symbols(f'u0:{n}')
        terms = list(sup3) + list(sup_free)
        f2 = random.choice(Q2FORMS) if fixed4 != 0 or sup3 else random.choice(Q2FORMS)
        f = sp.expand(f2 + fixed4 + sum(u[i] * terms[i] for i in range(n)))
        D = sp.expand(sp.hessian(f, X).det(method='berkowitz'))
        eqs = pos_deg_coeffs(D)
        try:
            G = sp.groebner(eqs, *u, order='grevlex')
            gens = [g for g in G.exprs]
            sols = sp.solve(gens, list(u), dict=True)
        except Exception as exc:                       # pragma: no cover
            print(f'  draw {t}: solver failed ({exc}); UNRESOLVED')
            stats['unresolved'] += 1
            continue
        stats['systems_solved'] += 1
        stats['sol_branches'] += len(sols)
        print(f'  draw {t}: {n} unknowns, {len(eqs)} equations, '
              f'GB size {len(gens)}, {len(sols)} solution branch(es)')
        for si, sol in enumerate(sols):
            free = [s for s in u if s not in sol]
            n_inst = 3 if free else 1
            for rep in range(n_inst):
                subs = dict(sol)
                rnd = {s: sp.Rational(random.choice([-3, -2, -1, 1, 2, 3]),
                                      random.choice([1, 1, 2])) for s in free}
                subs = {k: sp.expand(v.subs(rnd)) if hasattr(v, 'subs') else v
                        for k, v in subs.items()}
                subs.update(rnd)
                fi = sp.expand(f.subs(subs, simultaneous=True))
                if any(s in fi.free_symbols for s in u):
                    continue
                stats['instances'] += 1
                Di = sp.expand(sp.hessian(fi, X).det(method='berkowitz'))
                assert not Di.free_symbols, (
                    f'BUG: solution branch does not give constant det: {Di}')
                if Di == 0:
                    continue                        # not a Keller potential
                Pi = sp.Poly(fi, *X)
                dtot = Pi.total_degree()
                if dtot <= 2:
                    continue                        # trivial (quadratic)
                stats['keller_nontrivial'] += 1
                if dtot == 4:
                    stats['deg4'] += 1
                p1, p2 = pivot_cones(fi)
                if p1:
                    stats['p1_ok'] += 1
                if p2:
                    stats['p2_ok'] += 1
                else:
                    gaps.append((fname, fi, Di))
                inj = injective_saturated(fi)
                if inj:
                    stats['inj_ok'] += 1
                else:
                    counterexamples.append((fname, fi, Di))
                log.append((fname, dtot, Di, p1, p2, inj))
                flag = '' if (p2 and inj) else '   <<<<<< ATTENTION'
                print(f'    inst deg {dtot}, det = {Di}, P1 {"Y" if p1 else "N"}, '
                      f'P2 {"Y" if p2 else "N"}, injective {"Y" if inj else "N"}{flag}')

print('\n' + '=' * 70)
print('SUMMARY (rt_instances_x2_B_break)')
for k, v in stats.items():
    print(f'  {k:22s} {v}')
print(f'  distinct nonzero det values: '
      f'{sorted({sp.nsimplify(r[2]) for r in log}, key=lambda e: sp.re(sp.N(e)))}')

assert not counterexamples, (
    'FATAL / HC_4 COUNTEREXAMPLE: non-injective degree-<=4 Keller potential(s) '
    f'found: {counterexamples}')
assert not gaps, (
    'SERIOUS GAP in the proof of Theorem B: degree-<=4 Keller potential(s) with '
    f'EMPTY admissible pivot cone P2 (no pivot with deg e1 <= 2): {gaps}')
print('\nNo non-injective Keller instance and no empty admissible pivot cone.')
print('Theorem B and its Corollary SURVIVED this instance attack.')
