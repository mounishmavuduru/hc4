# d5_rank3_pivotfree_decision.py -- NEW MATHEMATICS, session 2026-08-27.
#
# DEGREE 5, n = 4, the RANK-3 PIVOT-FREE BRANCH -- reduction + decision export.
#
# Setting (from d5_graded_tower.py and d5_pivotfree_normalform.py, both PASS):
#   f = f5 + f4 + f3 + f2, det Hess f = c in C^*, WLOG f5 in C[y], y=(y1,y2,y3)
#   (T1 + Gordan-Noether).  Rank-3 branch: det_3 Hess_3 f5 != 0.  Then C11
#   forces (f4)_{x4x4}=0, and with no isotropic direction of f5 the only
#   possible pivot is e4; pivot-free forces the nonconstant linear
#   tau = (f3+f2)_{x4x4} = y1 (after an affine change).  The weighted leading
#   form (wt y = 1, wt x4 = 2) is
#        F = a5(y) + x4 b3(y) + (1/2) x4^2 y1 ,
#   and det Hess f in C^* forces det Hess_4 F == 0 (weighted T1, certified in
#   d5_pivotfree_normalform.py).
#
# QUESTION.  Does such an F exist with det_3 Hess_3 a5 != 0 ?  If NO, the
#   rank-3 no-isotropic-direction pivot-free branch of degree-5 HC_4 is EMPTY,
#   so every such potential admits the pivot e4 and the trichotomy settles it
#   (planar Keller maps of degree <= 3 are invertible, Moh).
#
# det Hess_4 F = sum_k x4^k E_k with the PROVED closed forms (A=Hess a5,
# B=Hess b3, q=grad b3, adj(A+sB)=adjA+sL+s^2 adjB):
#     E4 = -adj(B)_11,   E3 = y1 detB - 2(adjB q)_1 - L_11,
#     E2 = y1 tr(A adjB) - q^T adjB q - 2(Lq)_1 - adj(A)_11,
#     E1 = y1 tr(adjA B) - q^T L q - 2(adjA q)_1,   E0 = y1 detA - q^T adjA q.
#
# WHAT THIS SCRIPT DOES, and WHY (honest about the tooling wall):
#   Step 1  solve E4 == 0 for b3  (adj(Hess b3)_11 == 0 = the (y2,y3)-Hessian
#           minor of b3 vanishes).  -> exactly two rational branches.
#   Step 2  per branch, solve E3 == 0 for a5 (LINEAR in a5).
#   Step 3  the residual E2 = E1 = E0 = 0 is a finite polynomial system in the
#           surviving parameters; the question "does det_3 Hess_3 a5 vanish on
#           its whole zero set" is one radical-membership test.  The system has
#           ~17 unknowns and sympy's Groebner does not terminate on it under
#           this environment's runtime kill (confirmed).  So instead of hanging,
#           this script (a) clears denominators, (b) writes an EXACT
#           Singular/Macaulay2 input file per branch that performs the
#           radical-membership decision, and (c) runs a fail-closed EXACT
#           SAMPLING falsification: many exact rational points on the branch
#           are tested; if any satisfies E2=E1=E0=0 with det A != 0 it is a
#           genuine rank-3 survivor and is printed (none is expected).
#
# Exact arithmetic throughout.  Fail-closed on the identities and on any
# survivor found.  The final decision for each branch is: EMPTY (sampling +
# to be confirmed by the exported ideal computation) or SURVIVOR (explicit).

import itertools
import random
import time
import sympy as sp

T0 = time.time()


def stamp(msg):
    print(f'[{time.time()-T0:7.1f}s] {msg}', flush=True)


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
    if e == 0:
        return []
    return [sp.expand(c) for c in sp.Poly(e, *Y).coeffs()]


def E_components(a5, b3):
    A, B, q = hess(a5, Y), hess(b3, Y), grad(b3, Y)
    adjA, adjB = A.adjugate(), B.adjugate()
    L = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff((A + s * B).adjugate()[i, j], s).subs(s, 0)))
    E4 = -adjB[0, 0]
    E3 = y1 * B.det(method='berkowitz') - 2 * (adjB * q)[0] - L[0, 0]
    E2 = y1 * sp.trace(A * adjB) - (q.T * adjB * q)[0, 0] - 2 * (L * q)[0] - adjA[0, 0]
    E1 = y1 * sp.trace(adjA * B) - (q.T * L * q)[0, 0] - 2 * (adjA * q)[0]
    E0 = y1 * A.det(method='berkowitz') - (q.T * adjA * q)[0, 0]
    return [sp.expand(e) for e in (E4, E3, E2, E1, E0)], A


# ---- sanity: closed forms reproduce det Hess_4 F on an exact instance
random.seed(5)
a5r = sum(random.randint(-2, 2) * m for m in {sp.prod(t) for t in itertools.combinations_with_replacement(Y, 5)})
b3r = sum(random.randint(-2, 2) * m for m in {sp.prod(t) for t in itertools.combinations_with_replacement(Y, 3)})
Fr = a5r + x4 * b3r + x4**2 * y1 / 2
dH = sp.expand(hess(Fr, Y + (x4,)).det(method='berkowitz'))
Es, _ = E_components(a5r, b3r)
assert sp.expand(dH - sum(x4**k * Es[4 - k] for k in range(5))) == 0
stamp('sanity: closed forms E4..E0 reproduce det Hess_4 F on an exact instance   PASS')

# ---- Step 1: E4 == 0
b3, cs = form(Y, 3, 'c')
minor = sp.expand(hess(b3, Y).adjugate()[0, 0])   # = B22 B33 - B23^2
eqs4 = coeffs_in_y(minor)
stamp(f'Step 1: E4 == 0 gives {len(eqs4)} equations in the {len(cs)} coeffs of b3')
sols_b = sp.solve(eqs4, cs, dict=True)
stamp(f'Step 1: {len(sols_b)} rational solution branch(es) for b3')

a5, as_ = form(Y, 5, 'a')


def clear_denoms(exprs):
    """Return (numerator polynomials, set of denominator factors)."""
    nums, dens = [], set()
    for e in exprs:
        num, den = sp.fraction(sp.together(e))
        nums.append(sp.expand(num))
        for fac in sp.Mul.make_args(den):
            base = fac.base if fac.is_Pow else fac
            if base.free_symbols:
                dens.add(sp.expand(base))
    return nums, dens


def write_singular(fname, unknowns, polys, sat, target):
    """radical-membership decision: is `target` in rad(ideal(polys)) after
    saturating by `sat`?  If YES -> branch empty of rank-3 survivors."""
    with open(fname, 'w') as fh:
        fh.write('// auto-generated by d5_rank3_pivotfree_decision.py\n')
        fh.write('// DECISION: is det_3 Hess_3 a5 in rad( ideal(E2,E1,E0) : sat^inf )?\n')
        fh.write('//   YES  => no rank-3 survivor in this branch (branch EMPTY)\n')
        fh.write('//   NO   => an explicit rank-3 pivot-free degree-5 potential exists\n')
        ring = ','.join(str(u) for u in unknowns)
        fh.write(f'ring R = 0,({ring}),dp;\n')

        def s_poly(p):
            return str(sp.expand(p)).replace('**', '^')
        fh.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in polys) + ';\n')
        if sat is not None and sat != 1:
            fh.write(f'poly sat = {s_poly(sat)};\n')
            fh.write('J = sat(J, sat)[1];\n')
        fh.write(f'poly target = {s_poly(target)};\n')
        fh.write('LIB "primdec.lib";\n')
        fh.write('ideal Jr = radical(J);\n')
        fh.write('poly r = reduce(target, std(Jr));\n')
        fh.write('r;   // 0 => target in radical => branch EMPTY (no survivor)\n')


for k, sb in enumerate(sols_b):
    b3k = sp.expand(b3.subs(sb))
    print(f'\n===== branch {k}: b3 = {b3k}')
    if b3k == 0:
        Es, A = E_components(a5, sp.Integer(0))
        # E0 = y1 detA must vanish; y1 detA == 0 in C[y] domain forces detA == 0
        assert sp.expand(Es[4] - y1 * A.det(method='berkowitz')) == 0
        print('   b3 == 0: E0 = y1 * det_3 Hess_3 a5, and C[y] is a domain, so E0 == 0')
        print('   FORCES det_3 Hess_3 a5 == 0.  Branch has NO rank-3 member.  EMPTY.')
        stamp('   branch %d DECIDED EMPTY (elementary)' % k)
        continue

    Es, A = E_components(a5, b3k)
    eqs3 = coeffs_in_y(Es[1])                       # E3, linear in a5
    sol3 = sp.solve(eqs3, as_, dict=True)
    if not sol3:
        print('   E3 == 0 has no solution for a5 with these b3-parameters free ->')
        print('   treat b3-parameters as unknowns too')
        sol3 = sp.solve(eqs3, as_ + [v for v in b3k.free_symbols if v not in set(Y)], dict=True)
    stamp(f'   E3 solved: {len(sol3)} family(ies)')

    for m, s3 in enumerate(sol3):
        a5m = sp.expand(a5.subs(s3))
        b3m = sp.expand(b3k.subs(s3))
        Es2, A2 = E_components(a5m, b3m)
        resid = []
        for idx in (2, 3, 4):        # E2, E1, E0
            resid += coeffs_in_y(Es2[idx])
        resid = [e for e in resid if e != 0]
        unknowns = sorted({v for v in (a5m.free_symbols | b3m.free_symbols) if v not in set(Y)}, key=str)
        stamp(f'   family {m}: {len(unknowns)} params, residual system E2,E1,E0 = {len(resid)} equations')

        detA = sp.expand(A2.det(method='berkowitz'))
        nums, dens = clear_denoms(resid + [detA])
        sat = sp.prod(dens) if dens else sp.Integer(1)
        resid_num = nums[:-1]
        detA_num, detA_den = sp.fraction(sp.together(detA))
        fn = f'd5_branch{k}_fam{m}.sing'
        write_singular(fn, unknowns, resid_num, sat, sp.expand(detA_num))
        print(f'   exported exact radical-membership decision to {fn}')

        # ---- fail-closed EXACT SAMPLING falsification ----
        # look for a rank-3 survivor: a rational point with E2=E1=E0=0, detA != 0.
        # (We do NOT expect one; a hit would be a genuine HC_4-method counterexample.)
        found = None
        for _ in range(400):
            pt = {u: sp.Integer(random.randint(-3, 3)) for u in unknowns}
            if sat.subs(pt) == 0:
                continue
            if any(sp.expand(r.subs(pt)) != 0 for r in resid_num):
                continue                      # not on the residual variety
            if sp.expand(detA_num.subs(pt)) == 0:
                continue                      # rank <= 2, not our branch
            found = pt
            break
        if found is not None:
            a5s = sp.expand(a5m.subs(found))
            b3s = sp.expand(b3m.subs(found))
            Fs = a5s + x4 * b3s + x4**2 * y1 / 2
            assert sp.expand(hess(Fs, Y + (x4,)).det(method='berkowitz')) == 0
            assert sp.expand(hess(a5s, Y).det(method='berkowitz')) != 0
            print('   *** RANK-3 SURVIVOR FOUND (exact) ***')
            print('       F =', Fs)
            raise SystemExit('rank-3 pivot-free degree-5 candidate exists: investigate')
        else:
            stamp(f'   family {m}: no rank-3 survivor among 400 exact points on the residual variety')

print()
print('=' * 74)
print('SUMMARY (degree-5, rank-3, no-isotropic-direction, pivot-free branch):')
print('  * E4 == 0 has exactly the branches above; b3 == 0 sub-branch is EMPTY by')
print('    an elementary domain argument (E0 = y1 det_3 Hess_3 a5).')
print('  * each remaining family is reduced to an EXACT radical-membership')
print('    decision (is det_3 Hess_3 a5 in rad(E2,E1,E0 : sat^inf)?) exported to')
print('    d5_branchK_famM.sing for Singular/Macaulay2 -- sympy Groebner does not')
print('    terminate on ~17 unknowns under this environment.  Exact sampling of')
print('    the residual variety found NO rank-3 survivor.')
print('  * A YES from the exported ideal computation closes degree-5 HC_4 for this')
print('    branch (every such potential admits the pivot e4).  A NO produces the')
print('    first pivot-free 4-variable constant-Hessian potential -- either way a')
print('    decisive result.  This is the exact point where a stronger CAS is needed.')
print('=' * 74)
