# r5_rank21_topforms.py  --  key: r5  (HC_4 in DEGREE 5: rank-2 and rank-1
#                                       leading-form branches, top graded pieces)
#
# SETTING.  f = f5 + f4 + f3 + f2 in C[x1..x4], det Hess f = c in C^x, f5 in
# C[x'] (T1 + Gordan-Noether).  Let r = essential rank of f5.  The rank-3
# branch is closed in r5_rank3_topform.py + r5_rank3_cascade.py (e4 is always a
# pivot there).  Here we record the exact top graded pieces of det Hess f in
# the r = 2 and r = 1 branches, proved by the UNIVERSAL mixed-determinant
# expansion (independent symbols for every Hessian block entry; a t-graded
# determinant; the coefficient extraction is a proof), and read off the
# forced singular-matrix conditions.
#
# RANK 2 (WLOG f5 in C[x1,x2], det Hess_{(x1,x2)} f5 != 0):
#     [det Hess f]_10 = det Hess_{(x1,x2)}(f5) * det Hess_{(x3,x4)}(f4).
#   Proof: a weight-10 permutation product of Hess entries (deg f5 = 3,
#   f4 = 2, f3 = 1, f2 = 0 per entry) uses exactly two entries from Hess f5
#   and two from Hess f4; Hess f5 is supported on the (x1,x2) 2x2 block, so the
#   two f5-entries occupy rows/cols {1,2} and the two f4-entries rows/cols
#   {3,4}.  Hence det Hess_{(x3,x4)} f4 == 0.
#   BUT the entries of Hess_{(x3,x4)} f4 are QUADRICS, so (unlike the deg-4
#   case, S1) this does NOT force a constant kernel: no pivot follows from the
#   top piece alone.  This is the essential new difficulty in degree 5.
#
# RANK 1 (WLOG f5 = x1^5):
#     [det Hess f]_9 = 20 x1^3 * det Hess_{(x2,x3,x4)}(f4).
#   Proof: the only weight-9 product using three entries of Hess f5 needs all
#   three in position (1,1) (the sole nonzero entry, = 20 x1^3) -- impossible
#   for distinct rows/cols -- so weight 9 = 3 + 2 + 2 + 2: ONE f5-entry, forced
#   to (1,1), and three f4-entries in rows/cols {2,3,4}.  Hence
#   det Hess_{(x2,x3,x4)} f4 == 0.  Again the entries are QUADRICS.
#
# CONSEQUENCE.  In both low-rank branches the forced condition is the singular
# locus of a symmetric matrix of QUADRATIC forms, not linear forms; the deg-4
# closers (2x2 lemma S1, isotropic-vector lemma S3) do not apply verbatim, and
# a pivot is not produced by the top piece.  Whether a pivot is nonetheless
# forced (by the whole graded tower) or a pivot-free witness exists is the
# residual degree-5 question; see r5_rank2_probe.py / r5_summary.
#
# Exact arithmetic, fail-closed.

import itertools
import sys
import time
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
V = (x1, x2, x3, x4)
t = sp.Symbol('t')
T0 = time.time()


def say(m):
    print('[%5.1fs] %s' % (time.time() - T0, m))
    sys.stdout.flush()


def sym4(name, support):
    """generic 4x4 symmetric matrix with nonzero entries only on `support`
    (a set of frozensets/pairs of indices), else 0."""
    M = sp.zeros(4, 4)
    for i in range(4):
        for j in range(i, 4):
            if (i, j) in support or (j, i) in support:
                s = sp.Symbol('%s_%d%d' % (name, i + 1, j + 1))
                M[i, j] = s
                M[j, i] = s
    return M


def full_sym4(name):
    return sym4(name, {(i, j) for i in range(4) for j in range(4)})


def tgraded_det(H, k):
    P = sp.Poly(sp.expand(H.det(method='berkowitz')), t)
    return sp.expand(P.coeff_monomial(t**k))


# ============================================================ RANK 2
say('RANK 2: f5 in C[x1,x2]')
# Hess f5 supported on (x1,x2) block; Hess f4, f3, f2 full symmetric.
H5 = sym4('P5', {(0, 0), (0, 1), (1, 1)})          # deg-3 entries
H4 = full_sym4('P4')                               # deg-2 entries
H3 = full_sym4('P3')                               # deg-1 entries
H2 = full_sym4('P2')                               # deg-0 entries
H = t**3 * H5 + t**2 * H4 + t * H3 + H2
top = tgraded_det(H, 10)
d2f5 = sp.expand(H5[0, 0] * H5[1, 1] - H5[0, 1]**2)
d2f4_34 = sp.expand(H4[2, 2] * H4[3, 3] - H4[2, 3]**2)
assert sp.expand(top - d2f5 * d2f4_34) == 0
say('  [det]_10 = det Hess_{(x1,x2)} f5 * det Hess_{(x3,x4)} f4        PROVED')
say('  => det Hess_{(x3,x4)} f4 == 0 (a 2x2 symmetric matrix of QUADRICS)')

# demonstrate the obstruction: a singular 2x2 matrix of quadrics need NOT have
# a constant kernel vector (contrast S1 for linear entries).
a2, b2 = sp.symbols('a b')                          # candidate v = (0,0,a,b)
# take f4_33 = x1^2, f4_34 = x1 x2, f4_44 = x2^2  => det = 0 but Q(a,b) =
# (a x1 + b x2)^2 has NO constant (a,b) making it identically 0.
Q = sp.expand(a2**2 * x1**2 + 2 * a2 * b2 * x1 * x2 + b2**2 * x2**2)
assert sp.expand(Q - (a2 * x1 + b2 * x2)**2) == 0
sols = sp.solve([sp.Poly(Q, x1, x2).coeff_monomial(m) for m in (x1**2, x1 * x2, x2**2)],
                [a2, b2], dict=True)
assert sols == [{a2: 0, b2: 0}]
say('  witness: Hess_{(x3,x4)}f4 = [[x1^2,x1x2],[x1x2,x2^2]] is singular but')
say('           has only the trivial constant kernel: NO pivot from top piece')

# ============================================================ RANK 1
say('RANK 1: f5 = x1^5')
H5b = sym4('R5', {(0, 0)})                          # only (1,1) entry, = 20 x1^3
H4b = full_sym4('R4')
H3b = full_sym4('R3')
H2b = full_sym4('R2')
Hb = t**3 * H5b + t**2 * H4b + t * H3b + H2b
top9 = tgraded_det(Hb, 9)
M3 = sp.Matrix(3, 3, lambda i, j: H4b[i + 1, j + 1])
assert sp.expand(top9 - H5b[0, 0] * M3.det(method='berkowitz')) == 0
say('  [det]_9 = (f5)_x1x1 * det Hess_{(x2,x3,x4)} f4                  PROVED')
# concretely (f5)_x1x1 = 20 x1^3:
assert sp.expand(sp.diff(x1**5, x1, 2) - 20 * x1**3) == 0
say('  = 20 x1^3 * det Hess_{(x2,x3,x4)} f4  => det Hess_{(x2,x3,x4)} f4 == 0')
say('    (a 3x3 symmetric matrix of QUADRICS -- Perazzo-rich, not linear)')

print()
print('r5_rank21_topforms: ALL PASS')
print('Top-piece structure of the low-rank degree-5 branches established.')
print('Both reduce to a symmetric matrix of QUADRATIC forms being singular;')
print('unlike degree 4, the top piece alone does NOT force a pivot.')
