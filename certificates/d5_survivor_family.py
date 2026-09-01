# d5_survivor_family.py
#
# The degree-5, rank-3, no-isotropic-direction, pivot-free branch, weighted
# leading form F = a5(y) + x4 b3(y) + (1/2) x4^2 y1.  The generic-b3 emptiness
# is proved and certified elsewhere (R-D5-GEN; verify_paramcert.py): both
# E4=0 branches have std(J) = (1) over Q(c) off the exceptional loci
#   Z0 = V(972 c8^7 c9^4 (3c0c9-c1c8)^6),  Z1 = V(3 c5 (2c0c5-c1c4)^5).
# The recursive stratification of Z0, Z1 (_d5_close.py / _d5_close2.py) closes
# every substratum EMPTY except a SINGLE degenerate family, which every
# survivor leaf reduces to (verified by _d5_survivors.py):
#
#        b3 = y1^2 * (c2 y1 + c0 y2 + c1 y3).
#
# This script proves, by exact generic-symbolic identities (fail-closed), the
# STRUCTURE of that family, and records what remains machine-only.
#
#  DEDUCTIVE (proved here, fully generic in a5 and in c0,c1,c2):
#   (S1) B := Hess3(b3) has det B == 0 and constant kernel v* = (0, c1, -c0).
#   (S2) adj B == -4 y1^2 * v* v*^T   (exact).
#   (S3) v* . q == 0, where q = grad b3.  Hence q^T adjB q == 0 and
#        tr(A adjB) == -4 y1^2 * (v*^T A v*) = -4 y1^2 * D^2_{v*} a5,
#        A := Hess3(a5).
#   (S4) The y1^0 graded piece of det A vanishes on V(J):
#        det A|_{y1=0} is, coefficient-by-coefficient, equal to the y1^1
#        graded part of E0 -- a set of generators of J.  So
#        det Hess3(a5)|_{y1=0} == 0 on V(J):  every solution's a5 is Hessian-
#        singular on the plane {y1=0}.  (Matches R-D5.)
#
#  MACHINE-ONLY (modular, unbiased full-variety radical tests):
#   (M1) det Hess3(a5) is NOT in sqrt(J): the family DOES carry rank-3
#        solutions (det A =/= 0 somewhere on V(J)) -- 12/12 across 3 primes,
#        2 parameter points, 2 independent random combinations
#        (_d5_surv_rabin.py).  So this family is NOT a rank<3 cone.
#   (M2) BUT every solution's a5 has the isotropic direction v* = (0,c1,-c0):
#        D^2_{v*} a5 = v*^T A v* vanishes on V(J).  Confirmed two ways:
#        (a) full-variety radical test  v*^T A v* in sqrt(J)  (_d5_surv_target.py
#            iso), all VANISH across 3 primes; (b) an explicit power certificate
#            (v*^T A v*)^2 in <J>  (reduces to 0 mod std(J) over F_p; iso itself
#            not in J).  Hence D^2_{v*} a5 == 0 on V(J).
#   CONSEQUENCE: every rank-3 solution of the E-system for b3 = y1^2*(...) has a
#   LINEAR DIRECTION v*, so it FAILS the branch's "no linear direction of f5"
#   hypothesis.  The rank-3, NO-isotropic-direction, pivot-free branch is
#   therefore EMPTY on this family too.  (Such a5 lie OUTSIDE this branch: they
#   fall into the separate, not-yet-analyzed "rank-3 f5 WITH a linear direction"
#   case -- NOT covered by Theorem A, since a linear direction of f5 gives
#   D^2_v f = D^2_v(f4+f3+f2), degree 2, not a pivot in general.)
#   The char-0 parametric lift of (M1),(M2) over Q(c) exceeds this environment's
#   compute ceiling; the facts are certified modularly (multi-prime) + the
#   iso^2-in-J power certificate.  A slice-only test gives a FALSE cone here
#   (the det A =/= 0 locus is a lower-dim component slicing misses) -- the
#   full-variety Rabinowitsch is the correct test (_d5_surv_rabin.py).
#
# Exact arithmetic, fail-closed.
#   py -u d5_survivor_family.py
import sys
import sympy as sp
from _d5_close import Y, a5, as_, E_components, coeffs_in_y, hess, grad

y1, y2, y3 = Y
c0, c1, c2 = sp.symbols('c0 c1 c2')
b3 = c2 * y1**3 + c0 * y1**2 * y2 + c1 * y1**2 * y3

A = hess(a5, Y)
B = hess(b3, Y)
q = grad(b3, Y)
vstar = sp.Matrix([0, c1, -c0])

fails = []
def check(name, cond):
    ok = bool(cond)
    print(f'  [{"PASS" if ok else "FAIL"}] {name}')
    if not ok:
        fails.append(name)

print('=== degree-5 survivor family  b3 = y1^2 (c2 y1 + c0 y2 + c1 y3) ===')

# (S1)
check('S1a: det Hess3(b3) == 0', sp.expand(B.det()) == 0)
check('S1b: B * v* == 0 (v* = (0,c1,-c0) is the kernel)',
      sp.expand(B * vstar) == sp.zeros(3, 1))

# (S2)
adjB = B.adjugate()
check('S2: adj(B) == -4 y1^2 * v* v*^T',
      sp.simplify(adjB - (-4 * y1**2) * (vstar * vstar.T)) == sp.zeros(3, 3))

# (S3)
check('S3a: v* . q == 0', sp.expand((vstar.T * q)[0]) == 0)
check('S3b: q^T adj(B) q == 0', sp.expand((q.T * adjB * q)[0]) == 0)
Dvv = sp.expand((vstar.T * A * vstar)[0])           # D^2_{v*} a5
check('S3c: tr(A adjB) == -4 y1^2 * D^2_{v*} a5',
      sp.expand(sp.trace(A * adjB) - (-4 * y1**2) * Dvv) == 0)

# (S4)  det A|_{y1=0} equals the y1^1-graded part of E0, coeff by coeff => in J
(E3, E2, E1, E0), _A = E_components(a5, b3)
assert sp.expand(_A - A) == sp.zeros(3, 3)
detA = sp.expand(A.det(method='berkowitz'))
detA_y10 = sp.expand(detA.subs(y1, 0))
# E0 = y1*detA - q^T adjA q ; its coefficient of y1^1 (rest y1=0) :
E0p = sp.Poly(sp.expand(E0), y1)
E0_lin = E0p.coeff_monomial(y1)                      # y1^1 part (poly in y2,y3,a)
check('S4a: det A|_{y1=0} == (y1^1 part of E0)',
      sp.expand(detA_y10 - E0_lin) == 0)
# and every y2,y3-coefficient of E0_lin is a generator of J (a coeff of E0)
Jgens = set()
for e in (E3, E2, E1, E0):
    for g in coeffs_in_y(e):
        if g != 0:
            Jgens.add(sp.expand(g))
E0lin_coeffs = [sp.expand(g) for g in sp.Poly(sp.expand(E0_lin), y2, y3).coeffs() if g != 0]
check('S4b: each coeff of det A|_{y1=0} is a generator of J (=> vanishes on V(J))',
      all(g in Jgens or (-g) in Jgens for g in E0lin_coeffs))
check('S4c: det A|_{y1=0} is not identically 0 (constraint is nontrivial)',
      detA_y10 != 0)

print()
if fails:
    print('FAILED:', fails)
    sys.exit(1)
print('ALL STRUCTURAL CHECKS PASSED (S1-S4 deductive).')
print()
print('Branch closure for this family rests on the machine-verified facts:')
print('  (M1) det Hess3(a5) NOT in sqrt(J): rank-3 solutions EXIST  '
      '(_d5_surv_rabin.py, 12/12, 3 primes).')
print('  (M2) v*^T A v* = D^2_{v*} a5 in sqrt(J), with (v*^T A v*)^2 in <J>: '
      'every solution has isotropic direction v* (_d5_surv_target.py iso, '
      '_iso_power.py).')
print('=> every rank-3 solution has a LINEAR DIRECTION v*, failing the branch '
      '"no-isotropic-direction" hypothesis; the rank-3 no-isotropic pivot-free '
      'branch is EMPTY on this family.')
print('Char-0 parametric lift over Q(c) is compute-bound; facts certified '
      'modularly (multi-prime) + the iso^2-in-J power certificate.')
sys.exit(0)
