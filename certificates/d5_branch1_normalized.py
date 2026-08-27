# d5_branch1_normalized.py -- close the E4-branch-1 of the degree-5 rank-3
# pivot-free question by NORMALISING first, then deciding.  Session 2026-08-27.
#
# From d5_rank3_pivotfree_decision.py: E4 == 0 (the (y2,y3)-Hessian minor of
# b3 vanishes) has two rational branches.  Branch 1 is
#     b3 = y1 * ( c2 y1^2 + c0 y1 y2 + c1 y1 y3 + (c4^2/(4 c5)) y2^2
#                 + c4 y2 y3 + c5 y3^2 )
# whose (y2,y3)-quadratic part (c4^2/4c5) y2^2 + c4 y2 y3 + c5 y3^2
#     = c5 ( (c4/(2 c5)) y2 + y3 )^2   is a PERFECT SQUARE of a linear form.
#
# RESIDUAL SYMMETRY that preserves F = a5 + x4 b3 + x4^2 y1/2 and tau = y1:
#   (S1) GL on (y2,y3) fixing y1  -- move the square's linear form to y3;
#   (S2) shear y2 -> y2 + s y1, y3 -> y3 + t y1 -- weight-preserving;
#   (S3) scale y1 -> u y1 (with compensating x4-scale) and overall scale of F;
#   (S4) x4 -> x4 + Q(y), Q a quadratic form -- weight-2 "translation" of x4,
#        which changes b3 -> b3 + Q_{y1}?  NO: it shifts a5 and b3 by
#        derivatives; we DO NOT use (S4) here to keep the reduction transparent.
# Using (S1) we set the perfect-square linear form to y3, i.e. c4 = 0, and
# (S3) scales c5 = 1.  Then
#     b3 = y1 ( c2 y1^2 + c0 y1 y2 + c1 y1 y3 + y3^2 ) .
# Finally (S2) with y3 -> y3 + t y1 kills c1 (choose t = -c1/2), and the
# y1 y2-shear y2 -> y2 + r y1 leaves the y3^2 term intact while shifting the
# y1 y2 and y1^2 coefficients; we keep c0, c2 general (the shear that removes
# c0 also changes a5, tracked below).  NORMAL FORM used here:
#     b3 = y1 ( c2 y1^2 + c0 y1 y2 + y3^2 ) .            (2 parameters)
# This is a genuine sub-slice; if it contains NO rank-3 survivor we have NOT
# closed branch 1 (the full branch has the extra shears), but if it DOES we
# would have found one.  To keep the decision COMPLETE we instead run the
# 4-parameter form b3 = y1(c2 y1^2 + c0 y1 y2 + c1 y1 y3 + c5 y3^2) with c5
# free (only c4 = 0 imposed, which is a genuine GL_2 normalisation that loses
# nothing), and decide THAT.
#
# Exact, fail-closed.  Prints a definite verdict for the c4 = 0 normal form.

import itertools
import time
import sympy as sp

T0 = time.time()


def stamp(m):
    print(f'[{time.time()-T0:7.1f}s] {m}', flush=True)


y1, y2, y3, x4 = sp.symbols('y1 y2 y3 x4')
Y = (y1, y2, y3)


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def coeffs_in_y(expr):
    e = sp.expand(expr)
    return [] if e == 0 else [sp.expand(c) for c in sp.Poly(e, *Y).coeffs()]


# normalised branch-1 b3 (c4 = 0 by GL_2; c5 kept free)
c0, c1, c2, c5 = sp.symbols('c0 c1 c2 c5')
b3 = y1 * (c2 * y1**2 + c0 * y1 * y2 + c1 * y1 * y3 + c5 * y3**2)
# check E4 == 0 on the nose
assert sp.expand(hess(b3, Y).adjugate()[0, 0]) == 0
stamp('normalised branch-1 b3 satisfies E4 == 0 exactly (c4 = 0 GL_2 normal form)')

# generic quintic a5
mons5 = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(Y, 5)}, key=str)
acoef = sp.symbols(f'a0:{len(mons5)}')
a5 = sum(a * m for a, m in zip(acoef, mons5))

s = sp.Symbol('s')
A = hess(a5, Y)
B = hess(b3, Y)
q = sp.Matrix([sp.diff(b3, v) for v in Y])
adjA, adjB = A.adjugate(), B.adjugate()
L = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff((A + s * B).adjugate()[i, j], s).subs(s, 0)))
E3 = sp.expand(y1 * B.det(method='berkowitz') - 2 * (adjB * q)[0] - L[0, 0])
E2 = sp.expand(y1 * sp.trace(A * adjB) - (q.T * adjB * q)[0, 0] - 2 * (L * q)[0] - adjA[0, 0])
E1 = sp.expand(y1 * sp.trace(adjA * B) - (q.T * L * q)[0, 0] - 2 * (adjA * q)[0])
E0 = sp.expand(y1 * A.det(method='berkowitz') - (q.T * adjA * q)[0, 0])
stamp('E0..E3 built')

# solve E3 (linear in a5)
sol3 = sp.solve(coeffs_in_y(E3), list(acoef), dict=True)
assert len(sol3) == 1
a5 = sp.expand(a5.subs(sol3[0]))
A = hess(a5, Y)
stamp(f'E3 solved; remaining a-parameters: {len(sorted(a5.free_symbols - set(Y) - {c0,c1,c2,c5}, key=str))}')

# residual system
E2 = sp.expand(E2.subs(sol3[0]))
E1 = sp.expand(E1.subs(sol3[0]))
E0 = sp.expand(E0.subs(sol3[0]))
resid = [e for e in (coeffs_in_y(E2) + coeffs_in_y(E1) + coeffs_in_y(E0)) if e != 0]
unknowns = sorted({v for v in set().union(*[r.free_symbols for r in resid]) if v not in set(Y)}, key=str)
stamp(f'residual E2=E1=E0: {len(resid)} equations in {len(unknowns)} unknowns')

detA = sp.expand(A.det(method='berkowitz'))

# clear denominators, saturate against c5 (the only denominator)
zz = sp.Symbol('zz')
nums, dens = [], set()
for e in resid + [detA]:
    n, d = sp.fraction(sp.together(e))
    nums.append(sp.expand(n))
    for fac in sp.Mul.make_args(d):
        base = fac.base if fac.is_Pow else fac
        if base.free_symbols:
            dens.add(sp.expand(base))
sat = sp.prod(dens) if dens else sp.Integer(1)
resid_num, detA_num = nums[:-1], nums[-1]
stamp(f'denominators cleared; saturation poly = {sat}')

# DECISION via radical membership: is detA_num in rad( <resid_num> : sat^inf )?
# Rabinowitsch for saturation + power test.  Bounded: try Groebner; if it does
# not return quickly the environment kill ends it and we still have the export.
try:
    G = sp.groebner(resid_num + [sat * zz - 1], *(unknowns + [zz]), order='grevlex')
    stamp(f'saturated Groebner: {len(G.exprs)} elements')
    if list(G.exprs) == [1]:
        print('RESULT: residual ideal is EMPTY after saturation (no b3,a5 with det A != 0).')
        print('  => normalised branch 1 has NO rank-3 survivor: EMPTY.')
    else:
        red = G.reduce(detA_num)[1]
        # power test up to 4
        inrad = (red == 0)
        p = detA_num
        for _ in range(4):
            if G.reduce(sp.expand(p))[1] == 0:
                inrad = True
                break
            p = sp.expand(p * detA_num)
        if inrad:
            print('RESULT: det_3 Hess_3 a5 lies in the radical of the saturated residual ideal.')
            print('  => every solution has det_3 Hess_3 a5 == 0: normalised branch 1 has NO')
            print('     rank-3 survivor.  EMPTY.  (Degree-5 rank-3 pivot-free, c4=0 slice, closed.)')
        else:
            print('RESULT: det_3 Hess_3 a5 NOT detected in the radical by power test <= 4.')
            print('  => a rank-3 survivor MAY exist; run the exported Singular decision.')
except Exception as ex:
    print('Groebner did not terminate here:', repr(ex)[:120])
    print('  exact system exported by d5_rank3_pivotfree_decision.py stands.')

stamp('done')
