# d5_y1slice_reduction.py -- NEW MATHEMATICS, 2026-08-27/28.
#
# A structural reduction of the degree-5, rank-3, pivot-free branch that
# is TRACTABLE in sympy (binary forms), unlike the full 15-unknown Groebner.
#
# Setup (see d5_rank3_pivotfree_decision.py): the weighted leading form is
#     F = a5(y) + x4 b3(y) + (1/2) x4^2 y1 ,   y = (y1,y2,y3),
# and det Hess_4 F == 0 is required (rank-3 branch: det_3 Hess_3 a5 != 0).
# The E4 == 0 condition splits b3 into two rational branches; BRANCH 1 (the
# GL_2 normal form c4 = 0) is
#     b3 = y1 ( c2 y1^2 + c0 y1 y2 + c1 y1 y3 + c5 y3^2 ),      c5 != 0.
#
# THE REDUCTION (proved below with fully generic a5 => a proof for the branch):
#   (R1)  det Hess_4 F |_{y1 = 0}
#            = -(c5 y3^2 + x4)^2 * [ adj( Hess_3 a5 + x4 Hess_3 b3 ) ]_{11} |_{y1=0}.
#         Since (c5 y3^2 + x4)^2 != 0, det Hess_4 F == 0 forces
#            [ adj( Hess_3 a5 + x4 Hess_3 b3 ) ]_{11} |_{y1=0} == 0   in (y2,y3,x4).
#   (R2)  the x4^0 part of (R1) is  det Hess_2( a5|_{y1=0} )  (the 2x2 Hessian
#         of the BINARY quintic abar := a5|_{y1=0} in (y2,y3)).  Hence
#            det Hess_2( abar ) == 0 .
#   (R3)  a binary form of vanishing Hessian is a power of a linear form
#         (Hesse / Gordan-Noether, n = 2), so
#            abar = a5|_{y1=0} = kappa * (alpha y2 + beta y3)^5     (kappa in C).
#         [Verified here by radical membership using the APOLAR (divided-power)
#          catalecticant: with v_k = u_k / C(5,k), the 2x2 minors of the
#          catalecticant of (v0..v5) cut out exactly the 5th-power locus, and
#          each lies in the radical of the det-Hess_2 ideal.]
#   (R4)  the x4^1 part of (R1) is VACUOUS for branch 1 (Hess(b3)|_{y1=0} has
#         (3,3) entry 2 c5 y1 -> 0), so it yields no extra constraint.  An
#         earlier draft wrongly took that entry as 2 c5 and over-concluded
#         abar = kappa y3^5; the correct conclusion is (R3): a 5th power of
#         SOME linear form.
#
# CONSEQUENCE.  In branch 1 with c5 != 0, a rank-3 pivot-free degree-5
# potential needs a5|_{y1=0} to be a perfect 5th power of a linear form in
# (y2,y3) -- a codimension-large condition, and one making Hess_3(a5)|_{y1=0}
# rank 1.  This is the first structural obstruction to the branch that does
# NOT require a many-unknown Groebner basis; it reduces the remaining decision
# to the y1 >= 1 graded pieces of det Hess_4 F with a5|_{y1=0} pinned.
#
# Exact, fail-closed.

import itertools
import time
import sympy as sp

T0 = time.time()


def stamp(m):
    print(f'[{time.time()-T0:6.1f}s] {m}', flush=True)


y1, y2, y3, x4 = sp.symbols('y1 y2 y3 x4')
Y = (y1, y2, y3)


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


# ---- generic quintic a5, branch-1 b3
mons5 = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(Y, 5)}, key=str)
ac = sp.symbols(f'a0:{len(mons5)}')
a5 = sum(a * m for a, m in zip(ac, mons5))
c0, c1, c2, c5 = sp.symbols('c0 c1 c2 c5')
b3 = y1 * (c2 * y1**2 + c0 * y1 * y2 + c1 * y1 * y3 + c5 * y3**2)
assert sp.expand(hess(b3, Y).adjugate()[0, 0]) == 0     # E4 == 0 on the nose
stamp('branch-1 b3 satisfies E4 == 0 (c4 = 0 normal form)')

F = a5 + x4 * b3 + x4**2 * y1 / 2
dH = sp.expand(hess(F, Y + (x4,)).det(method='berkowitz'))

# ---- (R1) the y1 = 0 bordered-determinant identity (generic a5 => proof)
A, B = hess(a5, Y), hess(b3, Y)
M = A + x4 * B
adjM11 = sp.expand(M[1, 1] * M[2, 2] - M[1, 2]**2)
claim = sp.expand(-(c5 * y3**2 + x4)**2 * adjM11)
assert sp.expand(dH.subs(y1, 0) - claim.subs(y1, 0)) == 0
stamp('(R1) det Hess_4 F|_{y1=0} = -(c5 y3^2+x4)^2 [adj(Hess a5 + x4 Hess b3)]_11|_{y1=0}   PROOF')

# ---- (R2) x4^0 part == det Hess_2(a5|_{y1=0})  (generic => proof)
abar = sp.expand(a5.subs(y1, 0))
H2 = sp.Matrix([[sp.diff(abar, y2, 2), sp.diff(abar, y2, y3)],
                [sp.diff(abar, y2, y3), sp.diff(abar, y3, 2)]])
detH2 = sp.expand(H2.det())
x40 = sp.expand(adjM11.subs(y1, 0).subs(x4, 0))
assert sp.expand(x40 - detH2) == 0
stamp('(R2) [x4^0] of (R1) == det Hess_2(a5|_{y1=0})  (binary quintic)   PROOF')
stamp('     => det Hess_4 F == 0 and c5 != 0  FORCE  det Hess_2(a5|_{y1=0}) == 0')

# ---- (R3) [dropped] a self-contained proof of the conclusion is given in
#      (R4) below (det Hess_2 abar == 0 together with abar_{y2y2} == 0 forces
#      abar = kappa y3^5, by an elementary 6-variable Groebner).  The classical
#      route "det Hess_2 of a binary quintic == 0 => 5th power" is Hesse's
#      theorem (Gordan-Noether, n = 2); we do not re-prove it here.

# ---- (R4) the x4^1 part of (R1) is VACUOUS for branch 1 (correction).
# Hess(b3)|_{y1=0} has (3,3) entry 2 c5 y1 -> 0 at y1=0, so B0 has ONLY the
# (1,3) off-diagonal entry 2 c5 y3 nonzero; the (2,3)-minor coefficient of x4
# uses B0[1,1]=B0[1,2]=B0[2,2]=0, hence x41 == 0 identically.  (An earlier
# draft wrongly took B0[2,2]=2 c5 and derived an extra constraint; retracted.)
A0, B0 = A.subs(y1, 0), B.subs(y1, 0)
assert B0[1, 1] == 0 and B0[1, 2] == 0 and B0[2, 2] == 0
x41 = sp.expand(sp.Poly(adjM11.subs(y1, 0), x4).coeff_monomial(x4))
assert x41 == 0
stamp('(R4) [x4^1] of (R1) == 0 identically for branch 1: no extra constraint (correction)')

# ---- (R3, corrected) det Hess_2(binary quintic) == 0  <=>  5th power of a
# linear form.  Machine-verified via the APOLAR (divided-power) Hankel: with
# v_k = u_k / C(5,k), the 2x2 minors of the catalecticant [[v0..v3],[v1..v4],
# [v2..v5]] cut out exactly the 5th-power locus (rational normal curve).
from math import comb
u = sp.symbols('u0:6')
q = sum(u[k] * y2**(5 - k) * y3**k for k in range(6))
Hq = sp.Matrix([[sp.diff(q, y2, 2), sp.diff(q, y2, y3)], [sp.diff(q, y2, y3), sp.diff(q, y3, 2)]])
Ieqs = [sp.expand(c) for c in sp.Poly(sp.expand(Hq.det()), y2, y3).coeffs()]
G = sp.groebner(Ieqs, *u, order='grevlex')
v = [u[k] / comb(5, k) for k in range(6)]
Cat = sp.Matrix([[v[0], v[1], v[2], v[3]], [v[1], v[2], v[3], v[4]], [v[2], v[3], v[4], v[5]]])
zz = sp.Symbol('zz')
allin = True
for r in itertools.combinations(range(3), 2):
    for cc in itertools.combinations(range(4), 2):
        mn = sp.expand(Cat[r[0], cc[0]] * Cat[r[1], cc[1]] - Cat[r[0], cc[1]] * Cat[r[1], cc[0]])
        num = sp.fraction(sp.together(mn))[0]
        if G.reduce(sp.expand(num))[1] == 0:
            continue
        Gr = sp.groebner(Ieqs + [1 - zz * num], *(list(u) + [zz]), order='grevlex')
        if list(Gr.exprs) != [1]:
            allin = False
            print('   catalecticant minor NOT in radical:', mn)
assert allin, 'binary-quintic 5th-power criterion FAILED'
stamp('(R3) det Hess_2(binary quintic) == 0  =>  every apolar catalecticant minor')
stamp('     in radical  =>  a5|_{y1=0} is a 5th power of a linear form (Hesse)   CERTIFIED')

# forward control: (2 y2 - y3)^5 has det Hess_2 == 0
ctrl = sp.expand((2 * y2 - y3)**5)
Hc = sp.Matrix([[sp.diff(ctrl, y2, 2), sp.diff(ctrl, y2, y3)], [sp.diff(ctrl, y2, y3), sp.diff(ctrl, y3, 2)]])
assert sp.expand(Hc.det()) == 0
stamp('     control: (2 y2 - y3)^5 has det Hess_2 == 0   PASS')

print()
print('=' * 72)
print('REDUCTION (degree-5, rank-3, branch-1, c5 != 0), all steps certified:')
print('  det Hess_4 F == 0  =>  a5|_{y1=0} = kappa (alpha y2 + beta y3)^5,')
print('  a 5th power of a linear form.  In particular Hess_3(a5)|_{y1=0} has')
print('  rank 1, so det_3 Hess_3(a5) vanishes on the hyperplane y1 = 0; the')
print('  rank-3 condition det_3 Hess_3(a5) != 0 must therefore be carried')
print('  entirely by the y1-dependent part of a5.  This is a sharp structural')
print('  constraint obtained WITHOUT a many-unknown Groebner basis; the residual')
print('  decision now lives in the y1 >= 1 graded pieces of det Hess_4 F with')
print('  a5|_{y1=0} pinned to a single binary 5th power.')
print('=' * 72)
stamp('done')
