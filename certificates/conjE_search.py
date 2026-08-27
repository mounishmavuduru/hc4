# Fast targeted search for a counterexample to Conjecture E.
#
# CONJECTURE E: e1 in C[x1,x2,x3] with identically vanishing bordered Hessian
#   g^T adj(Hess e1) g == 0   (g = grad e1)
# is affinely 2-variable (i.e. some constant direction v != 0 has D_v e1 == 0).
#
# Strategy: by Theorem D' the leading form is a cone, WLOG in C[x1,x2].
# Any genuine counterexample must therefore hide its third variable in the
# lower-order terms. We solve (c2) exactly on structured slices where the
# system is small enough to solve outright, and test each solution for a
# constant null direction. Exact arithmetic throughout; a solution without a
# null direction would REFUTE Conjecture E (and would be a genuinely new
# object: a polynomial with developable level surfaces that is not planar).

import itertools
import sympy as sp

x1, x2, x3 = X = sp.symbols('x1 x2 x3')


def bordered(e):
    K = sp.Matrix([[sp.diff(e, u, v) for v in X] for u in X])
    g = sp.Matrix([sp.diff(e, u) for u in X])
    return sp.expand((g.T * K.adjugate() * g)[0, 0])


def null_direction(e):
    """Return a nonzero constant v with D_v e == 0, or None."""
    v = sp.symbols('n1 n2 n3')
    expr = sp.expand(sum(vi*sp.diff(e, u) for vi, u in zip(v, X)))
    p = sp.Poly(expr, *X)
    sols = sp.linsolve(p.coeffs(), list(v))
    for s in sols:
        s = list(s)
        free = set()
        for comp in s:
            free |= (set(comp.free_symbols) & set(v))
        if free:
            sub = {f: 1 for f in free}
            vv = [sp.simplify(comp.subs(sub)) for comp in s]
            if any(c != 0 for c in vv):
                return vv
        elif any(sp.simplify(c) != 0 for c in s):
            return s
    return None


counterexamples = []
checked = 0

# ---------- Slice 1: cubic, leading form a cone, general lower terms ----------
# e1 = c(x1,x2) + q(x1,x2,x3) + linear, solved exactly on sub-slices.
c_c = sp.symbols('c0:4')
q_c = sp.symbols('q0:6')
l_c = sp.symbols('l0:3')
c = c_c[0]*x1**3 + c_c[1]*x1**2*x2 + c_c[2]*x1*x2**2 + c_c[3]*x2**3
qmons = [x1**2, x1*x2, x2**2, x1*x3, x2*x3, x3**2]
q = sum(a*m for a, m in zip(q_c, qmons))
lin = l_c[0]*x1 + l_c[1]*x2 + l_c[2]*x3
e_gen = c + q + lin

# Solve on slices: fix the cubic part to each of the GL2-orbit
# representatives of binary cubics, then solve (c2) for the rest.
binary_cubics = [x1**3, x1**2*x2, x1**3 + x2**3, x1*x2*(x1 + x2), sp.Integer(0)]
for cub in binary_cubics:
    e = cub + q + lin
    bh = sp.expand(bordered(e))
    eqs = sp.Poly(bh, *X).coeffs() if bh != 0 else []
    unknowns = list(q_c) + list(l_c)
    sols = sp.solve(eqs, unknowns, dict=True) if eqs else [{}]
    for s in sols or []:
        es = sp.expand(e.subs(s))
        # residual free parameters: specialize them over a small grid
        free = sorted({sym for sym in es.free_symbols if sym not in set(X)}, key=str)
        grids = itertools.product([0, 1, -1, 2], repeat=len(free)) if free else [()]
        for vals in grids:
            ei = sp.expand(es.subs(dict(zip(free, vals))))
            # sympy's solve can return solutions with denominators; a grid
            # value hitting a pole yields zoo/nan -- skip those specializations.
            if ei.has(sp.zoo) or ei.has(sp.nan) or ei.has(sp.oo):
                continue
            if ei == 0 or sp.Poly(ei, *X).total_degree() < 2:
                continue
            checked += 1
            assert bordered(ei) == 0, f'solver error on {ei}'
            if null_direction(ei) is None:
                counterexamples.append(ei)

print(f'slice 1 (cubic, cone leading form): {checked} exact solutions tested')

# ---------- Slice 2: quartic pivot coefficients, cone leading form ----------
quartic_leads = [x1**4, x1**3*x2, x1**2*x2**2, x1**4 + x2**4, x1**2*x2**2 + x1**4]
c2_c = sp.symbols('d0:6')
c3_c = sp.symbols('e0:10')
cub_mons = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, 3)}, key=str)
quad_mons = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, 2)}, key=str)
checked2 = 0
for lead in quartic_leads:
    cu = sum(a*m for a, m in zip(c3_c, cub_mons))
    qd = sum(a*m for a, m in zip(c2_c, quad_mons))
    e = lead + cu + qd
    bh = sp.expand(bordered(e))
    eqs = sp.Poly(bh, *X).coeffs() if bh != 0 else []
    unknowns = list(c3_c) + list(c2_c)
    try:
        sols = sp.solve(eqs, unknowns, dict=True) if eqs else [{}]
    except Exception as exc:                       # noqa: BLE001
        print(f'  lead {lead}: solver gave up ({type(exc).__name__})')
        continue
    for s in sols or []:
        es = sp.expand(e.subs(s))
        free = sorted({sym for sym in es.free_symbols if sym not in set(X)}, key=str)
        if len(free) > 6:
            free = free[:6]
        grids = itertools.product([0, 1, -1], repeat=len(free)) if free else [()]
        for vals in grids:
            ei = sp.expand(es.subs(dict(zip(free, vals))))
            if ei.has(sp.zoo) or ei.has(sp.nan) or ei.has(sp.oo) or ei == 0:
                continue
            if bordered(ei) != 0:
                continue
            checked2 += 1
            if null_direction(ei) is None:
                counterexamples.append(ei)
print(f'slice 2 (quartic, cone leading form): {checked2} exact solutions tested')

# ---------- Slice 3: structured "developable" ansaetze that are NOT planar ----
# natural candidates: tangent-developable / ruled shapes in which x3 enters
# multiplicatively; each is checked exactly.
cands = [
    x1*x3 + x2, x1*x3 + x2**2, x1**2*x3 + x2, x1*x2*x3, x3*(x1**2 + x2**2),
    x1**3 + x1*x3 + x2, x1**3 + x2*x3, (x1 + x2*x3)**2, (x1 + x2*x3)**3,
    x1 + x2*x3 + x2**2, x3**2 + x1*x2, x1**2 + x2*x3 + x3**2,
    x1*x2 + x3*(x1 + x2), (x1*x2 + x3)**2, x1**2*x2 + x2**2*x3,
]
zero_bh_nonplanar = []
for e in cands:
    if bordered(e) == 0 and null_direction(e) is None:
        zero_bh_nonplanar.append(e)
print(f'slice 3 (structured candidates): {len(cands)} tested,'
      f' {len(zero_bh_nonplanar)} non-planar with vanishing bordered Hessian')
counterexamples.extend(zero_bh_nonplanar)

print()
if counterexamples:
    print('!!! CONJECTURE E REFUTED by the following exact examples:')
    for e in counterexamples[:10]:
        print('   e1 =', e)
        print('      bordered Hessian:', bordered(e), '  null direction:', null_direction(e))
else:
    print('NO COUNTEREXAMPLE FOUND. Every exact solution of the vanishing-bordered-')
    print('Hessian equation found in these slices is affinely 2-variable, consistent')
    print('with Conjecture E. (Evidence, not proof.)')
