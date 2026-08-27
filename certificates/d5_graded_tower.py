# d5_graded_tower.py -- graded pieces of det Hess f for deg f = 5, n = 4.
#
# Setting: f = f5 + f4 + f3 + f2, det Hess f = c in C^x.
# By T1 + Gordan-Noether (n=4) WLOG f5 in C[x1,x2,x3], so Hess f5 has zero
# 4th row/column.  Write
#   H5 = [[A5,0],[0,0]],  H4 = [[A4,q],[q^T,d4]],
#   H3 = [[A3,p],[p^T,r]], H2 = [[S3,u],[u^T,w]]
# where the entries of H_k are homogeneous of degree k-2.  Hence introducing a
# formal grading variable t and computing det(t^3 H5 + t^2 H4 + t H3 + H2) with
# INDEPENDENT symbols for every entry yields, coefficient by coefficient in t,
# universal polynomial expressions C_k in the block entries, and
#      [det Hess f]_k = C_k(actual blocks).
# This is a PROOF (the determinant is a universal polynomial in the entries and
# the substitution is degree-preserving), not a spot check.
#
# All arithmetic exact.  Fail-closed.

import itertools
import random
import sympy as sp

t = sp.Symbol('t')


def sym3(name):
    M = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            s = sp.Symbol('%s%d%d' % (name, i + 1, j + 1))
            M[i, j] = s
            M[j, i] = s
    return M


def vec3(name):
    return sp.Matrix([sp.Symbol('%s%d' % (name, i + 1)) for i in range(3)])


def block(M3, v, sc):
    M = sp.zeros(4, 4)
    M[0:3, 0:3] = M3
    for i in range(3):
        M[i, 3] = v[i]
        M[3, i] = v[i]
    M[3, 3] = sc
    return M


A5 = sym3('A')
A4 = sym3('B')
q = vec3('q')
d4 = sp.Symbol('d4')
A3 = sym3('C')
p = vec3('p')
r = sp.Symbol('r')
S3 = sym3('S')
u = vec3('u')
w = sp.Symbol('w')

H5 = block(A5, sp.zeros(3, 1), sp.Integer(0))
H4 = block(A4, q, d4)
H3 = block(A3, p, r)
H2 = block(S3, u, w)

H = t**3 * H5 + t**2 * H4 + t * H3 + H2
D = sp.expand(H.det(method='berkowitz'))
Dp = sp.Poly(D, t)
C = {k: sp.expand(Dp.coeff_monomial(t**k)) for k in range(0, 13)}

detA5 = sp.expand(A5.det(method='berkowitz'))
adjA5 = A5.adjugate()

print('--- universal graded coefficients (n=4, deg f = 5, f5 x4-free) ---')

# [det]_12  =  det H5 = 0
assert C[12] == 0, 'C12 must vanish (T1 for f5 x4-free)'
print('C12 = 0                                        [T1]')

# [det]_11 = d4 * det A5
assert sp.expand(C[11] - d4 * detA5) == 0
print('C11 = (f4)_x4x4 * det3(Hess3 f5)               [T3, d = 5]')

# [det]_10 : general form, then the d4 = 0 specialization
gen10 = sp.expand(C[10])
claim10 = sp.expand(r * detA5 - (q.T * adjA5 * q)[0, 0])
rem10 = sp.expand(gen10 - claim10)
assert sp.expand(rem10.subs(d4, 0)) == 0, 'C10 closed form fails at d4=0'
# and the d4-part is d4 * (sum over i of det(A5 with row i replaced by A4 row i))
tr_term = sp.expand(sum(sum(A4[i, j] * adjA5[i, j] for j in range(3)) for i in range(3)))
assert sp.expand(rem10 - d4 * tr_term) == 0, 'C10 d4-part closed form fails'
print('C10 = (f3)_x4x4*det(A5) - q^T adj(A5) q + d4*tr(A4 adj A5)')
print('      with q = grad_x\' (f4)_x4 ;  at d4 = 0:')
print('C10 = (f3)_x4x4 * det(A5)  -  q^T adj(A5) q')

# [det]_9 at d4 = 0
C9 = sp.expand(C[9].subs(d4, 0))
claim9 = sp.expand(
    w * detA5                                        # (f2)_x4x4 * det A5
    - 2 * (p.T * adjA5 * q)[0, 0]                    # cross term f3/f4
    + r * tr_term                                    # (f3)_x4x4 * tr(A4 adj A5)
)
rest9 = sp.expand(C9 - claim9)
print()
print('C9 - [ w*detA5 - 2 p^T adj(A5) q + r*tr(A4 adjA5) ] has %d terms'
      % len(sp.Poly(rest9, *sorted(rest9.free_symbols, key=str)).terms()))

# identify the remaining piece: it is the "one H5 row + three H4 rows" family.
# express it as  - sum_{i} q_i * det( A5 row i replaced ... )  -- verify by an
# independent construction: the mixed determinant with rows (H5,H4,H4,H4).
def mixed_det(rows):
    """rows: list of 4 matrices; take row i from rows[i]."""
    M = sp.Matrix(4, 4, lambda i, j: rows[i][i, j])
    return M.det(method='berkowitz')


H5g, H4g = H5, H4.subs(d4, 0)
caseC = sp.expand(sum(mixed_det([H4g if k != i else H5g for k in range(4)])
                      for i in range(3)))
assert sp.expand(rest9 - caseC) == 0, 'C9 residual is not the (H5,H4,H4,H4) family'
print('C9 = (f2)_x4x4*det(A5) - 2 p^T adj(A5) q + (f3)_x4x4*tr(A4 adj A5)')
print('     + sum_{i<=3} det(row i from H5, other rows from H4)   [d4 = 0]')

# ---------------------------------------------------------------------------
# Validation of the grading argument on a random exact instance.
# ---------------------------------------------------------------------------
X = sp.symbols('x1 x2 x3 x4')
x1, x2, x3, x4 = X
Xp = X[:3]
random.seed(20260810)


def rnd():
    return sp.Rational(random.randint(-5, 5), random.randint(1, 3))


def rand_form(vars_, deg):
    e = sp.Integer(0)
    for m in itertools.combinations_with_replacement(range(len(vars_)), deg):
        term = rnd()
        for i in m:
            term *= vars_[i]
        e += term
    return sp.expand(e)


f5 = rand_form(Xp, 5)
f4 = rand_form(X, 4)
f3 = rand_form(X, 3)
f2 = rand_form(X, 2)
f = f5 + f4 + f3 + f2
Hf = sp.Matrix(4, 4, lambda i, j: sp.diff(f, X[i], X[j]))
Df = sp.expand(Hf.det(method='berkowitz'))


def graded(pol, m):
    P = sp.Poly(pol, *X)
    out = sp.Integer(0)
    for mon, co in P.terms():
        if sum(mon) == m:
            trm = co
            for v, e in zip(X, mon):
                trm *= v**e
            out += trm
    return sp.expand(out)


def hess_part(g):
    return sp.Matrix(4, 4, lambda i, j: sp.diff(g, X[i], X[j]))


sub = {}
Hk = {5: hess_part(f5), 4: hess_part(f4), 3: hess_part(f3), 2: hess_part(f2)}
for i in range(3):
    for j in range(i, 3):
        sub[A5[i, j]] = Hk[5][i, j]
        sub[A4[i, j]] = Hk[4][i, j]
        sub[A3[i, j]] = Hk[3][i, j]
        sub[S3[i, j]] = Hk[2][i, j]
for i in range(3):
    sub[q[i]] = Hk[4][i, 3]
    sub[p[i]] = Hk[3][i, 3]
    sub[u[i]] = Hk[2][i, 3]
sub[d4] = Hk[4][3, 3]
sub[r] = Hk[3][3, 3]
sub[w] = Hk[2][3, 3]
assert Hk[5][3, 3] == 0 and all(Hk[5][i, 3] == 0 for i in range(3))

for k in range(0, 13):
    assert sp.expand(graded(Df, k) - C[k].subs(sub)) == 0, 'grading mismatch at k=%d' % k
print()
print('grading argument validated on a random exact instance, degrees 0..12  OK')
print('d5_graded_tower: ALL PASS')
