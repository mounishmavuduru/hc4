# r5_rank3_topform.py  --  key: r5  (HC_4 in DEGREE 5, rank-3 branch,
#                                      top weighted form)
#
# SETTING.  f = f5 + f4 + f3 + f2 in C[x1..x4], det Hess f = c in C^x.  By
# T1 + Gordan-Noether, WLOG f5 in C[x1,x2,x3] =: C[x'].  RANK-3 BRANCH:
# det_3 Hess_3 f5 =/= 0.  Then (d5_graded_tower.py, re-verified in part 0)
#     [det]_11 = (f4)_{x4x4} det_3 Hess_3 f5      =>  (f4)_{x4x4} = 0,
#     [det]_10 = (f3)_{x4x4} det_3 Hess_3 f5 - q^T adj(Hess_3 f5) q,
#                q = grad' (f4)_{x4},   everything except (f3)_{x4x4} x4-free
#                                              =>  (f3)_{x4x4} is x4-FREE.
# Hence  f = a(x') + x4 b(x') + x4^2 tau(x')/2,  deg a <= 5, deg b <= 3,
# tau = (f3)_{x4x4} + (f2)_{x4x4} AFFINE in x'.
#   * tau constant  =>  e4 is a pivot  =>  Theorem A / pivot trichotomy.
#   * tau non-constant: affine change (x1 -> tau, translation) makes tau = x1.
#     THIS FILE analyses that case.
#
# WEIGHTED LEADING FORM.  Weights wt(x1,x2,x3) = 1, wt(x4) = 2.  f has
# weighted degree 5 with leading form  G = a5 + x4 b3 + x1 x4^2/2,  and
# [det Hess f]_{wt 10} = det Hess_4 G  (d5_pivotfree_normalform.py (A);
# re-verified in part 0 on a generic instance of the shape).  So
#          det Hess f in C^x   ==>   det Hess_4 G == 0.
#
# THEOREM R3-TOP (proved here).  Let a5 be a ternary quintic with
# det Hess_3 a5 =/= 0 and b3 a ternary cubic with
#          det Hess_4( a5 + x4 b3 + x1 x4^2/2 ) == 0.
# Then after a linear change of (x2,x3) (fixing x1 and x4),
#          a5 = x1^3 x2^2/2 + x1^4 x3,     b3 = x1^2 x2,
# i.e.  G = G_I := x1 * ( (x4 + x1 x2)^2/2 + x1^3 x3 ).
# Conversely det Hess G_I == 0 (G_I is x2-free in the non-affine coordinates
# (x1,x2,x3,u = x4 + x1 x2)) and det Hess_3(x1^3 x2^2/2 + x1^4 x3) = -16 x1^9.
#
# PROOF (each machine step below is an identity in generic coefficients,
# hence a proof for all forms of the stated shape).
#   Write E_k := [x4^k] det Hess_4 G  (k = 0..4; k = 5 vanishes).
#   E4 = -det Hess_{(x2,x3)} b3.  [Z2] (de Bondt, JPAA 219 (2015) Thm
#   "zerohess"(i) = de Bondt-van den Essen 2004: a 2-variable polynomial with
#   zero Hessian is phi(l) + (affine), l a linear form) applied to
#   b3(1,x2,x3) and re-homogenised gives  b3 = Phi(x1,l) + x1^2 M  with l, M
#   linear in (x2,x3); a linear change of (x2,x3) makes
#          b3 = B(x1,x2) + beta x1^2 x3,   B a binary cubic.           (N1)
#   E3 = -(a5)_{x3x3} B_{x2x2}.                                        (N2)
#   CASE I: B_{x2x2} = 0, i.e. b3 = x1^2 (c0 x1 + al x2 + be x3).
#     I.0 (al,be) = 0:  E2 = -adj(A)_11 and E0 = x1 det A - 9c0^2 x1^4 adj(A)_11
#         (A = Hess_3 a5)  =>  det A = 0: contradiction.
#     I.1 (al,be) =/= 0: WLOG b3 = x1^2 x2.  E2 = -adj(A)_11, so
#         det Hess_{(x2,x3)} a5 = 0 and [Z2] gives a5 = Psi(x1,L) + x1^4 M,
#         L = l2 x2 + l3 x3, M = m2 x2 + m3 x3, Psi a binary quintic.  Then
#           det A = -16 x1^6 (l2 m3 - l3 m2)^2 Psi_LL,
#           E1 = 2 l3 x1^4 (4(l2 m3 - l3 m2) x1 + 3 l3 x2) Psi_LL,
#         so rank 3 forces l3 = 0, and then
#           E0 = -16 x1^7 m3^2 ( d^2/dx2^2 [Psi(x1,l2 x2) + m2 x1^4 x2] - x1^3 ),
#         so A5 := Psi(x1,l2x2) + m2 x1^4 x2 has A5_{x2x2} = x1^3, i.e.
#         a5 = x1^3 x2^2/2 + p1 x1^4 x2 + p0 x1^5 + m3 x1^4 x3, and
#         x3 -> (x3 - p1 x2 - p0 x1)/m3 gives the normal form.
#   CASE II: B_{x2x2} =/= 0, so a5 = A5(x1,x2) + x3 C4(x1,x2).  Then
#           E2 = (C4_{x2})^2 + beta B_{x2x2} x1 (6 beta x1 b3 - 8 C4);
#         its x3-coefficient is 6 beta^3 x1^4 B_{x2x2}, so beta = 0, then
#         C4_{x2} = 0, C4 = gamma x1^4, and
#           det Hess_4 G = -16 gamma^2 x1^6 ( x1 A5_{x2x2} + x1 x4 B_{x2x2} - B_{x2}^2 ),
#         whose x4-coefficient forces B_{x2x2} = 0: contradiction.
#
# All arithmetic exact (sympy 1.14, rationals).  Fail-closed.

import itertools
import random
import sys
import time
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Y = (x1, x2, x3)
V = (x1, x2, x3, x4)
W = (1, 1, 1, 2)
T0 = time.time()


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def gen(vs, d, tag):
    e = sp.Integer(0)
    ks = []
    for m in itertools.combinations_with_replacement(range(len(vs)), d):
        s = sp.Symbol(tag + ''.join(str(i + 1) for i in m))
        ks.append(s)
        t = s
        for i in m:
            t *= vs[i]
        e += t
    return sp.expand(e), ks


def x4coeffs(G):
    D = sp.expand(hess(G, V).det(method='berkowitz'))
    P = sp.Poly(D, x4)
    return {k: sp.expand(P.coeff_monomial(x4**k)) for k in range(6)}


def wparts(pol):
    P = sp.Poly(sp.expand(pol), *V)
    out = {}
    for mon, co in P.terms():
        d = sum(w * e for w, e in zip(W, mon))
        t = co
        for v, e in zip(V, mon):
            t *= v**e
        out[d] = out.get(d, 0) + t
    return {k: sp.expand(v) for k, v in out.items()}


def graded(pol, m):
    P = sp.Poly(sp.expand(pol), *V)
    out = sp.Integer(0)
    for mon, co in P.terms():
        if sum(mon) == m:
            t = co
            for v, e in zip(V, mon):
                t *= v**e
            out += t
    return sp.expand(out)


def say(msg):
    print('[%6.1fs] %s' % (time.time() - T0, msg))
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# PART 0.  Tower facts used (rank-3 branch), re-verified generically.
# ---------------------------------------------------------------------------
say('PART 0: tower facts for deg f = 5, f5 in C[x\']  (universal block expansion)')
# Universal method (as in d5_graded_tower.py): independent symbols for every
# block entry of Hess f_k, grading variable t; the t^k coefficient of
# det(t^3 H5 + t^2 H4 + t H3 + H2) is a universal polynomial in the block
# entries, and [det Hess f]_k is its evaluation.  This is a proof.
tt = sp.Symbol('t')


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


uA5, uA4, uA3, uS3 = sym3('uA'), sym3('uB'), sym3('uC'), sym3('uS')
uq, up, uu = vec3('uq'), vec3('up'), vec3('uu')
ud4, ur, uw = sp.symbols('ud4 ur uw')
uH = (tt**3 * block(uA5, sp.zeros(3, 1), 0) + tt**2 * block(uA4, uq, ud4)
      + tt * block(uA3, up, ur) + block(uS3, uu, uw))
uD = sp.Poly(sp.expand(uH.det(method='berkowitz')), tt)
uC = {k: sp.expand(uD.coeff_monomial(tt**k)) for k in range(13)}
udet = sp.expand(uA5.det(method='berkowitz'))
uadj = uA5.adjugate()
assert uC[12] == 0
assert sp.expand(uC[11] - ud4 * udet) == 0
say('  [det]_11 = (f4)_x4x4 * det3(Hess3 f5)                          PROVED')
utr = sp.expand(sum(uA4[i, j] * uadj[i, j] for i in range(3) for j in range(3)))
assert sp.expand(uC[10] - (ur * udet - (uq.T * uadj * uq)[0, 0] + ud4 * utr)) == 0
say('  [det]_10 = (f3)_x4x4 det3 - q^T adj3 q + (f4)_x4x4 tr(A4 adj3),')
say('             q = grad\'(f4)_x4, A4 = Hess3 f4                     PROVED')
# In the rank-3 branch det3 != 0 forces (f4)_x4x4 = 0 (so d4 = 0 and every
# entry of A5, A4, q is x4-free), and then the x4-dependence of [det]_10 is
# only through (f3)_x4x4 = r'(x') + rho x4:  [x4][det]_10 = rho det3.
# Instance check of exactly this on a random exact f with (f4)_x4x4 = 0:
random.seed(20260825)


def rnd():
    return sp.Rational(random.randint(-7, 7), random.randint(1, 3))


def rand_poly(vs, degs):
    e = sp.Integer(0)
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vs)), d):
            t = rnd()
            for i in m:
                t *= vs[i]
            e += t
    return sp.expand(e)


f5i = rand_poly(Y, [5])
f4i = rand_poly(Y, [4]) + x4 * rand_poly(Y, [3])
rho = sp.Symbol('rho')
f3i = rand_poly(Y, [3]) + x4 * rand_poly(Y, [2]) + x4**2 * rand_poly(Y, [1]) / 2 + rho * x4**3 / 6
Di = sp.expand(hess(f5i + f4i + f3i, V).det(method='berkowitz'))
d10i = graded(Di, 10)
det3i = sp.expand(hess(f5i, Y).det())
assert det3i != 0
assert sp.expand(sp.diff(d10i, x4) - rho * det3i) == 0
say('  d/dx4 [det]_10 = rho * det3(Hess3 f5) on an exact rank-3 instance  OK')
say('  => rank 3 forces (f4)_x4x4 = 0 and rho = 0:')
say('     f = a(x\') + x4 b(x\') + x4^2 tau(x\')/2 with tau AFFINE.           PROVED')

# weighted leading form identity on the shape  a + x4 b + x1 x4^2/2.
# GRADING ARGUMENT (proof).  Under the scaling D_s = diag(s,s,s,s^2),
#   Hess_x( f o D_s ) = D_s (Hess f)(D_s x) D_s,  so
#   det Hess_x( f o D_s ) = s^{2(1+1+1+2)} det Hess f (D_s x) = s^10 * det Hess f(D_s x).
# Writing f = sum_j f_[j] in weighted-homogeneous parts (f_[5] = G), the left
# side is det Hess_x( sum_j s^j f_[j] ), a polynomial in s whose top power is
# s^{10} det Hess G (only j = 5 pieces can reach weighted degree 10 in each of
# the four rows).  Matching the coefficient of s^{10} on both sides gives
#   det Hess G = [ det Hess f (D_s x) ]_{s^0-in-the-rescaled-sense}
#             = [det Hess f]_{wt 10} (evaluated at x).
# This is d5_pivotfree_normalform.py part (A), already certified generically
# (symbolic a5, b3).  Here we re-validate it on an exact instance AND verify
# the scaling identity det Hess_x(f o D_s) = s^10 det Hess f(D_s x) itself
# fully generically (which IS the grading argument).
say('  weighted leading form: [det Hess f]_{wt10} = det Hess_4(a5 + x4 b3 + x1 x4^2/2)')
s = sp.Symbol('s')
Ds = {x1: s * x1, x2: s * x2, x3: s * x3, x4: s**2 * x4}
# generic quartic-and-below f in the shape, symbolic top, symbolic-ish lower
ftest = (rand_poly(Y, [5]) + x4 * rand_poly(Y, [3]) + x1 * x4**2 / 2
         + rand_poly(Y, [4, 3, 2]) + x4 * rand_poly(Y, [2, 1]))
lhs = sp.expand(hess(ftest.subs(Ds), V).det(method='berkowitz'))
rhs = sp.expand(s**10 * hess(ftest, V).det(method='berkowitz').subs(Ds))
assert sp.expand(lhs - rhs) == 0
say('  scaling identity det Hess(f o D_s) = s^10 (det Hess f)(D_s x)   PROVED (generic)')
# the weight-10 part of det Hess f equals det Hess G (the actual claim):
Gpart = wparts(ftest)[5]                              # = a5 + x4 b3 + x1 x4^2/2
detf10 = wparts(sp.expand(hess(ftest, V).det(method='berkowitz')))[10]
assert sp.expand(detf10 - hess(Gpart, V).det(method='berkowitz')) == 0
say('  [wt 10] det Hess f = det Hess G on an exact instance; generic proof in')
say('  d5_pivotfree_normalform.py (A).                                         OK')

# ---------------------------------------------------------------------------
# PART 1.  E4 and E3 in closed form (generic a5, b3).
# ---------------------------------------------------------------------------
say('PART 1: E4..E0 closed forms (block form, generic symmetric matrices)')
a5s, _ = gen(Y, 5, 'A')          # generic ternary quintic
b3s, _ = gen(Y, 3, 'B')          # generic ternary cubic
Am = hess(a5s, Y)                # A = Hess_3 a5  (used symbolically below)
Bm = hess(b3s, Y)                # B = Hess_3 b3
# BLOCK FORM of Hess_4 G, G = a5 + x4 b3 + x1 x4^2/2:
#   Hess_4 G = [[A + x4 B, q + x4 e1],[(q+x4 e1)^T, x1]],  q = grad b3, e1 = (1,0,0).
# The E_k are UNIVERSAL polynomials in (A, B, q, y1); we prove them with
# GENERIC constant symmetric A, B and generic q, scalar y1 (fast: single-symbol
# entries), then specialise.  (This is d5_pivotfree_normalform.py (B).)
tt = sp.Symbol('t')


def _sym3(name):
    M = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            s = sp.Symbol('%s%d%d' % (name, i + 1, j + 1))
            M[i, j] = s
            M[j, i] = s
    return M


Ag, Bg = _sym3('AA'), _sym3('BB')
qg = sp.Matrix([sp.Symbol('qq%d' % i) for i in range(3)])
y1s = sp.Symbol('y1s')
e1v = sp.Matrix([1, 0, 0])
Mblk = sp.zeros(4, 4)
Mblk[0:3, 0:3] = Ag + tt * Bg
colb = qg + tt * e1v
for i in range(3):
    Mblk[i, 3] = colb[i]
    Mblk[3, i] = colb[i]
Mblk[3, 3] = y1s
detblk = sp.Poly(sp.expand(Mblk.det(method='berkowitz')), tt)
Eg = {k: sp.expand(detblk.coeff_monomial(tt**k)) for k in range(5)}
adjAg, adjBg = Ag.adjugate(), Bg.adjugate()
Lg = sp.Matrix(3, 3, lambda i, j: sp.expand(
    sp.diff((Ag + sp.Symbol('sv') * Bg).adjugate()[i, j], sp.Symbol('sv')).subs(sp.Symbol('sv'), 0)))
assert sp.expand(Eg[4] + adjBg[0, 0]) == 0
assert sp.expand(Eg[3] - (y1s * Bg.det() - 2 * (adjBg * qg)[0] - Lg[0, 0])) == 0
assert sp.expand(Eg[2] - (y1s * sp.trace(Ag * adjBg) - (qg.T * adjBg * qg)[0]
                          - 2 * (Lg * qg)[0] - adjAg[0, 0])) == 0
assert sp.expand(Eg[1] - (y1s * sp.trace(adjAg * Bg) - (qg.T * Lg * qg)[0] - 2 * (adjAg * qg)[0])) == 0
assert sp.expand(Eg[0] - (y1s * Ag.det() - (qg.T * adjAg * qg)[0])) == 0
say('  E4=-adjB_11, E3=y1 detB-2(adjB q)_1-L_11, E2=y1 tr(A adjB)-q^T adjB q')
say('   -2(Lq)_1-adjA_11, E1=y1 tr(adjA B)-q^T L q-2(adjA q)_1, E0=y1 detA-q^T adjA q  PROVED')
# Euler relations for a REAL cubic b3 (B = Hess b3, q = grad b3):
qb = sp.Matrix([sp.diff(b3s, v) for v in Y])
assert sp.expand(Bm * sp.Matrix(Y) - 2 * qb) == sp.zeros(3, 1)          # B y = 2q
assert sp.expand(Bm.adjugate() * qb - sp.Rational(1, 2) * Bm.det() * sp.Matrix(Y)) == sp.zeros(3, 1)
say('  Euler (real cubic b3): B y = 2q, adj(B) q = (detB/2) y            PROVED')
# specialise E4, E3 with y1 = x1, q = grad b3, and Euler => E4 = -det Hess_{(x2,x3)}b3,
# E3 = -L_11 = -(A22 B33 + A33 B22 - 2 A23 B23):
E4s = -(Bm[1, 1] * Bm[2, 2] - Bm[1, 2]**2)
L11 = Am[1, 1] * Bm[2, 2] + Am[2, 2] * Bm[1, 1] - 2 * Am[1, 2] * Bm[1, 2]
# check the specialisation is consistent with the block formula (adjB q)_1 = detB/2 * x1:
assert sp.expand((Bm.adjugate() * qb)[0] - sp.Rational(1, 2) * Bm.det() * x1) == 0
E3s = sp.expand(x1 * Bm.det() - 2 * sp.Rational(1, 2) * Bm.det() * x1 - L11)
assert sp.expand(E3s + L11) == 0
say('  => E4 = -det Hess_{(x2,x3)} b3,  E3 = -(A22 B33 + A33 B22 - 2 A23 B23)   PROVED')

# (N1) shape consequences: with b3 = B(x1,x2) + beta x1^2 x3, generic a5 (21 sym)
c0, al, mu, nu, be = sp.symbols('c0 alpha mu nu beta')
B3 = c0 * x1**3 + al * x1**2 * x2 + mu * x1 * x2**2 + nu * x2**3
b3N = B3 + be * x1**2 * x3
EN = x4coeffs(a5s + x4 * b3N + x1 * x4**2 / 2)
assert EN[4] == 0
assert sp.expand(EN[3] + sp.diff(a5s, x3, 2) * sp.diff(B3, x2, 2)) == 0
say('  (N1) => E4 = 0 and E3 = -(a5)_x3x3 * B_x2x2                      PROVED')
# converse sanity for [Z2]: Phi(x1,l) + x1^2 M has zero (x2,x3)-Hessian
l2, l3, m2, m3 = sp.symbols('l2 l3 m2 m3')
Lf = l2 * x2 + l3 * x3
Mf = m2 * x2 + m3 * x3
ph = sp.symbols('ph0:4')
PhiL = sum(ph[j] * x1**(3 - j) * Lf**j for j in range(4))
bz = sp.expand(PhiL + x1**2 * Mf)
Hz = hess(bz, (x2, x3))
assert sp.expand(Hz.det()) == 0
say('  [Z2]-shape b3 = Phi(x1,l) + x1^2 M indeed has det Hess_{(x2,x3)} = 0 (converse)')

# ---------------------------------------------------------------------------
# PART 2.  CASE I:  B_x2x2 = 0,  b3 = x1^2 (c0 x1 + al x2 + be x3)
# ---------------------------------------------------------------------------
say('PART 2: CASE I')
adjA = Am.adjugate()
detA = sp.expand(Am.det(method='berkowitz'))
# I.0
E0c = x4coeffs(a5s + x4 * c0 * x1**3 + x1 * x4**2 / 2)
assert E0c[4] == 0 and E0c[3] == 0 and E0c[1] == 0
assert sp.expand(E0c[2] + adjA[0, 0]) == 0
assert sp.expand(E0c[0] - x1 * detA + 9 * c0**2 * x1**4 * adjA[0, 0]) == 0
say('  I.0: E2 = -adj(A)_11,  E0 = x1 det A - 9 c0^2 x1^4 adj(A)_11  PROVED')
say('       => adj(A)_11 = 0 => det A = 0: CONTRADICTION with rank 3')
# I.1: b3 = x1^2 x2
E1c = x4coeffs(a5s + x4 * x1**2 * x2 + x1 * x4**2 / 2)
assert E1c[4] == 0 and E1c[3] == 0
assert sp.expand(E1c[2] + adjA[0, 0]) == 0
say('  I.1: E2 = -adj(A)_11 = -det Hess_{(x2,x3)} a5                   PROVED')
# [Z2] => a5 = Psi(x1,L) + x1^4 M
ps = sp.symbols('p0:6')
Psi = sum(ps[j] * x1**(5 - j) * Lf**j for j in range(6))
a5z = sp.expand(Psi + x1**4 * Mf)
Az = hess(a5z, Y)
detAz = sp.expand(Az.det(method='berkowitz'))
PsiLL = sp.expand(sum(ps[j] * j * (j - 1) * x1**(5 - j) * Lf**(j - 2) for j in range(2, 6)))
assert sp.expand(detAz + 16 * x1**6 * (l2 * m3 - l3 * m2)**2 * PsiLL) == 0
say('  I.1: det A = -16 x1^6 (l2 m3 - l3 m2)^2 Psi_LL                  PROVED')
Ez = x4coeffs(a5z + x4 * x1**2 * x2 + x1 * x4**2 / 2)
assert Ez[4] == 0 and Ez[3] == 0 and Ez[2] == 0
assert sp.expand(Ez[1] - 2 * l3 * x1**4 * (4 * (l2 * m3 - l3 * m2) * x1 + 3 * l3 * x2) * PsiLL) == 0
say('  I.1: E1 = 2 l3 x1^4 (4(l2 m3 - l3 m2) x1 + 3 l3 x2) Psi_LL      PROVED')
say('       rank 3 (Psi_LL != 0, l2 m3 - l3 m2 != 0) => l3 = 0')
Ez0 = sp.expand(Ez[0].subs(l3, 0))
A5z = sp.expand((Psi + x1**4 * m2 * x2).subs(l3, 0))
assert sp.expand(Ez0 + 16 * x1**7 * m3**2 * (sp.diff(A5z, x2, 2) - x1**3)) == 0
say('  I.1: E0|_{l3=0} = -16 x1^7 m3^2 (A5_x2x2 - x1^3),  A5 = a5 - m3 x1^4 x3   PROVED')
say('       => A5_x2x2 = x1^3 => a5 = x1^3 x2^2/2 + p1 x1^4 x2 + p0 x1^5 + m3 x1^4 x3')
# the normal form and its properties
GI = x1 * ((x4 + x1 * x2)**2 / 2 + x1**3 * x3)
assert sp.expand(GI - (x1**3 * x2**2 / 2 + x1**4 * x3 + x4 * x1**2 * x2 + x1 * x4**2 / 2)) == 0
assert sp.expand(hess(GI, V).det(method='berkowitz')) == 0
assert sp.expand(hess(x1**3 * x2**2 / 2 + x1**4 * x3, Y).det() + 16 * x1**9) == 0
# full normal-form family (before the x3-shear) also has zero Hessian:
p0, p1 = sp.symbols('p0 p1')
Gfam = x1**3 * x2**2 / 2 + p1 * x1**4 * x2 + p0 * x1**5 + m3 * x1**4 * x3 + x4 * x1**2 * x2 + x1 * x4**2 / 2
assert sp.expand(hess(Gfam, V).det(method='berkowitz')) == 0
say('  G_I = x1((x4+x1x2)^2/2 + x1^3 x3): det Hess_4 = 0, det Hess_3 a5 = -16 x1^9   OK')

# ---------------------------------------------------------------------------
# PART 3.  CASE II:  B_x2x2 != 0  =>  (a5)_x3x3 = 0,  a5 = A5(x1,x2) + x3 C4(x1,x2)
# ---------------------------------------------------------------------------
say('PART 3: CASE II')
A5s, _ = gen((x1, x2), 5, 'S')
C4s, _ = gen((x1, x2), 4, 'C')
a5II = sp.expand(A5s + x3 * C4s)
EII = x4coeffs(a5II + x4 * b3N + x1 * x4**2 / 2)
assert EII[4] == 0 and EII[3] == 0
B22 = sp.diff(B3, x2, 2)
claimE2 = sp.expand(sp.diff(C4s, x2)**2 + be * B22 * x1 * (6 * be * x1 * b3N - 8 * C4s))
assert sp.expand(EII[2] - claimE2) == 0
say('  II: E2 = (C4_x2)^2 + beta B_x2x2 x1 (6 beta x1 b3 - 8 C4)        PROVED')
assert sp.expand(sp.diff(EII[2], x3) - 6 * be**3 * x1**4 * B22) == 0
say('  II: d/dx3 E2 = 6 beta^3 x1^4 B_x2x2  =>  beta = 0')
E2b0 = sp.expand(EII[2].subs(be, 0))
assert sp.expand(E2b0 - sp.diff(C4s, x2)**2) == 0
say('  II: beta = 0 => E2 = (C4_x2)^2 => C4 = gamma x1^4')
gam = sp.Symbol('gamma')
a5IIb = sp.expand(A5s + gam * x1**4 * x3)
DII = sp.expand(hess(a5IIb + x4 * B3 + x1 * x4**2 / 2, V).det(method='berkowitz'))
claimD = sp.expand(-16 * gam**2 * x1**6 * (x1 * sp.diff(A5s, x2, 2) + x1 * x4 * B22 - sp.diff(B3, x2)**2))
assert sp.expand(DII - claimD) == 0
say('  II: det Hess_4 G = -16 gamma^2 x1^6 (x1 A5_x2x2 + x1 x4 B_x2x2 - B_x2^2)   PROVED')
say('      rank 3 needs gamma != 0; the x4-coefficient forces B_x2x2 = 0: CONTRADICTION')

print()
print('THEOREM R3-TOP: in the rank-3 branch with non-constant tau, the weighted')
print('leading form of f is, after a linear change of (x2,x3),')
print('      G_I = x1^3 x2^2/2 + x1^4 x3 + x1^2 x2 x4 + x1 x4^2/2')
print('          = x1 * ( (x4 + x1 x2)^2/2 + x1^3 x3 ).')
print('r5_rank3_topform: ALL PASS')
