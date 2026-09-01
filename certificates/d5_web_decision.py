# d5_web_decision.py -- build CORRECT, ready-to-paste web-CAS decisions for the
# degree-5 rank-3 pivot-free branches.  Supersedes the malformed target in the
# earlier .sing export (that wrote det_3 Hess_3 a5 verbatim, which still
# contains y1,y2,y3 -- variables NOT in the coefficient ring).
#
# CORRECT DECISION.  A rank-3 survivor in a branch is a choice of coefficients
# with (i) E2 = E1 = E0 == 0 identically in y  (i.e. in V(J), J = residual
# ideal in the coefficient ring), and (ii) det_3 Hess_3 a5 NOT == 0 in y.
# Condition (ii) means: not every y-coefficient of det_3 Hess_3 a5 vanishes.
# Hence:  branch EMPTY  <=>  EVERY y-coefficient d_alpha of det_3 Hess_3 a5
# lies in radical( J : sat^inf ).  A single random rational combination
#   dcomb = sum_alpha r_alpha d_alpha
# is in the radical iff all d_alpha are (generic r_alpha); we test dcomb, a
# single radical-membership query cheap enough for a free web CAS, and print
# the individual-coefficient test as a Singular loop for a fully rigorous run.
#
# Fast (no local Groebner).  Writes, per branch:
#   d5_web_branch{K}.singular.txt   (sagecell.sagemath.org, language Singular)
#   d5_web_branch{K}.magma.txt      (magma.maths.usyd.edu.au/calc)

import itertools
import time
import sympy as sp

T0 = time.time()


def stamp(m):
    print(f'[{time.time()-T0:6.1f}s] {m}', flush=True)


y1, y2, y3, x4, s = sp.symbols('y1 y2 y3 x4 s')
Y = (y1, y2, y3)


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def grad(g, vs):
    return sp.Matrix([sp.diff(g, v) for v in vs])


def form(vs, deg, tag):
    mons = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(vs, deg)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sp.expand(sum(c * m for c, m in zip(cs, mons))), list(cs)


def coeffs_in_y(expr):
    e = sp.expand(expr)
    return [] if e == 0 else [sp.expand(c) for c in sp.Poly(e, *Y).coeffs()]


def E_components(a5, b3):
    A, B, q = hess(a5, Y), hess(b3, Y), grad(b3, Y)
    adjA, adjB = A.adjugate(), B.adjugate()
    L = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff((A + s * B).adjugate()[i, j], s).subs(s, 0)))
    E2 = y1 * sp.trace(A * adjB) - (q.T * adjB * q)[0, 0] - 2 * (L * q)[0] - adjA[0, 0]
    E1 = y1 * sp.trace(adjA * B) - (q.T * L * q)[0, 0] - 2 * (adjA * q)[0]
    E0 = y1 * A.det(method='berkowitz') - (q.T * adjA * q)[0, 0]
    E3 = y1 * B.det(method='berkowitz') - 2 * (adjB * q)[0] - L[0, 0]
    return [sp.expand(e) for e in (E3, E2, E1, E0)], A


def clear(exprs):
    nums, dens = [], set()
    for e in exprs:
        n, d = sp.fraction(sp.together(e))
        nums.append(sp.expand(n))
        for fac in sp.Mul.make_args(d):
            base = fac.base if fac.is_Pow else fac
            if base.free_symbols:
                dens.add(sp.expand(base))
    return nums, dens


def s_poly(p):
    return str(sp.expand(p)).replace('**', '^')


def emit(k, unknowns, Jpolys, sat, dcomb, dcoeffs):
    ring = ','.join(str(u) for u in unknowns)
    # Singular (SageMathCell)
    with open(f'd5_web_branch{k}.singular.txt', 'w') as f:
        f.write('// PASTE INTO https://sagecell.sagemath.org (Language: Singular)\n')
        f.write(f'// degree-5 rank-3 pivot-free branch {k}: is it EMPTY?\n')
        f.write('//   final print 0  => EMPTY  => no counterexample in this branch\n')
        f.write('//   final print !=0 => a pivot-free 4-variable potential EXISTS\n')
        f.write('LIB "primdec.lib";\n')
        f.write(f'ring R = 0,({ring}),dp;\n')
        f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jpolys) + ';\n')
        if sat != 1:
            # NB: 'sat' is a reserved proc name in elim.lib (pulled in by
            # primdec.lib); name the saturation poly 'satp' to avoid shadowing.
            f.write(f'poly satp = {s_poly(sat)};\n  J = sat(J, satp)[1];\n')
        f.write('ideal Jr = radical(J);\n')
        # RIGOROUS: reduce EVERY y-coefficient of det_3 Hess_3 a5; branch is
        # EMPTY iff all reduce to 0.  (A single random combination could give a
        # false 0; testing all coefficients is exact.)
        f.write('ideal T = ' + ',\n  '.join(s_poly(d) for d in dcoeffs) + ';\n')
        f.write('ideal red = reduce(T, std(Jr));\n')
        f.write('// EMPTY  <=>  red is the zero ideal (every entry 0):\n')
        f.write('size(simplify(red,2));   // prints 0  <=> branch EMPTY\n')
    # Magma
    with open(f'd5_web_branch{k}.magma.txt', 'w') as f:
        f.write('// PASTE INTO http://magma.maths.usyd.edu.au/calc/\n')
        f.write(f'// degree-5 rank-3 pivot-free branch {k}: is it EMPTY?\n')
        f.write('//   prints true  => EMPTY => no counterexample\n')
        f.write('//   prints false => a pivot-free 4-variable potential EXISTS\n')
        f.write(f'R<{ring}> := PolynomialRing(RationalField(), {len(unknowns)});\n')
        f.write('J := ideal<R | ' + ',\n  '.join(s_poly(p) for p in Jpolys) + '>;\n')
        if sat != 1:
            f.write(f'J := Saturation(J, {s_poly(sat)});\n')
        f.write('J := Radical(J);\n')
        # RIGOROUS: every y-coefficient of det_3 Hess_3 a5 must lie in J.
        f.write('T := [\n  ' + ',\n  '.join(s_poly(d) for d in dcoeffs) + '\n];\n')
        f.write('print &and[ t in J : t in T ];   // true <=> branch EMPTY\n')


# ---- branches of E4 == 0
b3g, cs = form(Y, 3, 'c')
eqs4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
sols_b = sp.solve(eqs4, cs, dict=True)
a5, as_ = form(Y, 5, 'a')
# fixed small integer weights for the coefficient combination (deterministic,
# avoids Math.random; a generic-enough combination)
WT = [1, -2, 3, -5, 7, -11, 13, -17, 19, -23, 29, -31, 37, -41, 43, -47, 53, -59, 61, -67, 71]

for k, sb in enumerate(sols_b):
    b3k = sp.expand(b3g.subs(sb))
    if b3k == 0:
        print(f'branch {k}: b3 == 0 -> EMPTY by the elementary domain argument (E0 = y1 det A); no web run needed.')
        continue
    Es, A = E_components(a5, b3k)
    sol3 = sp.solve(coeffs_in_y(Es[0]), as_, dict=True)   # E3 linear in a5
    if not sol3:
        sol3 = sp.solve(coeffs_in_y(Es[0]), as_ + [v for v in b3k.free_symbols if v not in set(Y)], dict=True)
    s3 = sol3[0]
    a5m, b3m = sp.expand(a5.subs(s3)), sp.expand(b3k.subs(s3))
    Es2, A2 = E_components(a5m, b3m)
    resid = [e for e in (coeffs_in_y(Es2[1]) + coeffs_in_y(Es2[2]) + coeffs_in_y(Es2[3])) if e != 0]
    unknowns = sorted({v for v in (a5m.free_symbols | b3m.free_symbols) if v not in set(Y)}, key=str)
    detA = sp.expand(A2.det(method='berkowitz'))
    # clear the c5 denominators in the det-A y-coefficients (safe: we saturate
    # by those denominators, under which they are units, so radical membership
    # is unchanged).  Drop duplicates / zeros.
    dcoeffs_raw = coeffs_in_y(detA)
    dcoeffs = []
    seen = set()
    for d in dcoeffs_raw:
        num = sp.expand(sp.fraction(sp.together(d))[0])
        if num != 0 and num not in seen:
            seen.add(num)
            dcoeffs.append(num)
    dcomb = sp.expand(sum(WT[i % len(WT)] * dcoeffs[i] for i in range(len(dcoeffs))))
    Jn, dens = clear(resid + [dcomb] + dcoeffs)
    sat = sp.prod(dens) if dens else sp.Integer(1)
    Jpolys = Jn[:len(resid)]
    dcomb_num = Jn[len(resid)]
    emit(k, unknowns, Jpolys, sat, dcomb_num, dcoeffs)
    stamp(f'branch {k}: {len(unknowns)} vars, {len(Jpolys)} generators, '
          f'{len(dcoeffs)} y-coeffs of det A -> single combined target')
    print(f'   wrote d5_web_branch{k}.singular.txt and d5_web_branch{k}.magma.txt')

print()
print('HOW TO RUN (no install; bypasses the blocked package manager):')
print(' * sagecell.sagemath.org  -> language dropdown "Singular" -> paste a')
print('   d5_web_branch*.singular.txt file -> Evaluate.  0 => branch EMPTY.')
print(' * magma.maths.usyd.edu.au/calc  -> paste d5_web_branch*.magma.txt ->')
print('   Submit.  true => branch EMPTY.  false => a counterexample exists.')
print(' A "branch EMPTY" verdict on every branch closes degree-5 HC_4 in the')
print(' rank-3 case; a "false"/nonzero is the first pivot-free 4-var potential.')
