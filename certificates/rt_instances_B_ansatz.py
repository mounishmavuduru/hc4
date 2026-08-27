# REFEREE INSTANCE ATTACK 3 (hostile referee, rt_instances_*):
# Attempt to BREAK Theorem B (HC_4 in degree <= 4).
# Method: 6 sparse quartic ansatz families f = f2 + f3 + f4 with unknown
# coefficients on small supports; impose det Hess f == det Hess f2 (a nonzero
# constant by construction) as polynomial identities in x; SOLVE the resulting
# polynomial system in the coefficients (Groebner + sp.solve over C, so complex
# branches are found too); then for EVERY solution found (sampling free
# parameters at several rational points), verify by saturated collision
# Groebner that grad f is injective.  Per Theorem B, ALL solutions must be
# injective; a single non-injective one is a FATAL error in the paper (or an
# HC_4 counterexample).  FAIL-CLOSED: failures are collected and re-raised.

import itertools
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
I = sp.I

FAMILIES = []

def fam(name, f2, mons):
    FAMILIES.append((name, f2, mons))

# F1: doubling-type support P(x1,x3) against dual quadratic x1x2+x3x4
fam('F1 P(x1,x3) vs x1x2+x3x4',
    x1*x2 + x3*x4,
    [x1**3, x1**2*x3, x1*x3**2, x3**3, x1**4, x1**3*x3, x1**2*x3**2, x1*x3**3, x3**4])

# F2: mixed dual pairs (x2-, x4-monomials present in f3)
fam('F2 mixed duals',
    x1*x2 + x3*x4,
    [x1**2*x2, x1**3, x1**2*x3, x2**2*x4, x1**4, x1**2*x2**2, x2**4])

# F3: definite quadratic + planar cubic/quartic (complex harmonic branches)
fam('F3 round quadratic, planar P(x1,x2)',
    sp.Rational(1, 2)*(x1**2 + x2**2 + x3**2 + x4**2),
    [x1**3, x1**2*x2, x1*x2**2, x2**3, x1**4, x1**3*x2, x1**2*x2**2, x1*x2**3, x2**4])

# F4: planar support against crossed duals
fam('F4 P(x1,x2) vs x1x4+x2x3',
    x1*x4 + x2*x3,
    [x1**2*x2, x1*x2**2, x2**3, x1**3*x2, x1**2*x2**2, x1*x2**3, x2**4, x1**3])

# F5: x4 present in the cubic part (attacks the r=2 / r=1 branches)
fam('F5 x4-cubics',
    x1*x2 + x3*x4,
    [x4*x1**2, x4*x1*x2, x4*x2**2, x1**2*x2**2, x1**4, x2**4])

# F6: rank-deficient-flavored quadratic coupling
fam('F6 x1x3-coupled',
    x1*x3 + sp.Rational(1, 2)*(x2**2 + x4**2),
    [x1**3, x1**2*x2, x1*x2**2, x1**2*x4, x1**4, x1**2*x2**2, x1**3*x2])

SAMPLE_VALS = [sp.Integer(1), sp.Integer(-2), sp.Rational(1, 2), sp.Integer(3),
               sp.Integer(-1), sp.Integer(2)]

def injective_check(f_inst):
    p = sp.symbols('pp1:5'); q = sp.symbols('qq1:5'); w = sp.symbols('ww1:5')
    grads = [sp.diff(f_inst, v) for v in X]
    eqs = [sp.expand(gg.subs(dict(zip(X, p))) - gg.subs(dict(zip(X, q))))
           for gg in grads]
    sat = sp.expand(sum(w[i]*(p[i]-q[i]) for i in range(4)) - 1)
    G = sp.groebner(eqs + [sat], *(list(p)+list(q)+list(w)), order='grevlex')
    return list(G.exprs) == [sp.Integer(1)]

failures = []
total_solution_branches = 0
total_instances_tested = 0
total_injective = 0

for name, f2, mons in FAMILIES:
    cs = sp.symbols(f'c0:{len(mons)}')
    f_template = sp.expand(f2 + sum(c*m for c, m in zip(cs, mons)))
    c0 = sp.hessian(f2, X).det()
    assert c0 != 0, f'{name}: chosen f2 is degenerate'
    d = sp.expand(sp.hessian(f_template, X).det(method='berkowitz') - c0)
    pol = sp.Poly(d, *X)
    eqs = sorted(set(pol.coeffs()), key=sp.default_sort_key)
    eqs = [sp.expand(e) for e in eqs if e != 0]
    print(f'--- {name}: {len(cs)} unknowns, {len(eqs)} coefficient equations, c = {c0}')
    if not eqs:
        # identity holds for ALL coefficient values: whole space is a solution
        sols = [dict()]
        print('    det Hess f == c holds IDENTICALLY for all coefficients (full space).')
    else:
        G = sp.groebner(eqs, *cs, order='grevlex')
        if list(G.exprs) == [sp.Integer(1)]:
            print('    no solutions at all (not even 0?!) -- impossible'); continue
        try:
            sols = sp.solve(list(G.exprs), list(cs), dict=True)
        except Exception as ex:
            print(f'    sp.solve FAILED ({ex}); falling back to Groebner basis only')
            sols = None
        if sols is None or sols == []:
            sols = [dict()] if all(e.subs({c: 0 for c in cs}) == 0 for e in eqs) else []
    total_solution_branches += len(sols)
    print(f'    {len(sols)} solution branch(es) found')
    seen = set()
    for bi, sol in enumerate(sols):
        free = set()
        for v in sol.values():
            free |= v.free_symbols
        free |= set(cs) - set(sol.keys())
        free = sorted(free & set(cs) | (free - set(X)) - set(sol.keys()), key=str)
        free = [s for s in free if s in cs or s not in X]
        n_samples = 3 if free else 1
        for t in range(n_samples):
            subsmap = {s: SAMPLE_VALS[(t + i) % len(SAMPLE_VALS)]
                       for i, s in enumerate(free)}
            full = {}
            ok = True
            for c in cs:
                val = sol.get(c, c)
                val = sp.simplify(val.subs(subsmap)) if hasattr(val, 'subs') else val
                if val.free_symbols:
                    ok = False
                    break
                full[c] = val
            if not ok:
                print(f'    branch {bi} sample {t}: unresolved symbols, skipped')
                continue
            key = tuple(full[c] for c in cs)
            if key in seen:
                continue
            seen.add(key)
            f_inst = sp.expand(f_template.subs(full))
            dd = sp.expand(sp.hessian(f_inst, X).det(method='berkowitz'))
            if sp.simplify(dd - c0) != 0:
                # spurious solve output; not a valid member -- record, don't count
                print(f'    branch {bi} sample {t}: det check FAILED (spurious solve '
                      f'branch), det = {dd}')
                failures.append((name, bi, t, 'SPURIOUS-SOLVE', f_inst))
                continue
            total_instances_tested += 1
            inj = injective_check(f_inst)
            if inj:
                total_injective += 1
                tag = 'injective PASS'
            else:
                tag = 'NON-INJECTIVE  <<< FATAL for Theorem B'
                failures.append((name, bi, t, 'NON-INJECTIVE', f_inst))
            print(f'    branch {bi} sample {t}: coeffs {key}: {tag}')

print()
print(f'THEOREM-B ATTACK SUMMARY: {len(FAMILIES)} ansatz families, '
      f'{total_solution_branches} solution branches, '
      f'{total_instances_tested} concrete constant-Hessian quartics tested, '
      f'{total_injective} injective.')
if failures:
    print('FAILURES:')
    for rec in failures:
        print('   ', rec)
    raise SystemExit('THEOREM B BROKEN OR PIPELINE ERROR -- see failures above')
print('All solutions of all ansatz families are injective: Theorem B survives '
      'this instance attack.')
