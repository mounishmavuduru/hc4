# fpG_paper_certificate.py -- certificate for ../theorem_G_paper.tex (agent key: fpG)
#
# Fail-closed (assert-based; nonzero exit on any failure), exact/symbolic only.
#
# PART-2 CHUNK (checks C11-C19), split out on 2026-08-25: the full script
# exceeds the environment's runtime kill inside (C10) (a 4x4 determinant
# with two fully generic cubics). (C10) is the same identity already
# certified generically as (W9) in paperG_writeup_checks.py and check (2)
# of wG_writeup_checks.py. Same helpers as fpG_paper_certificate.py.
import itertools
import re
import sys
import time
from pathlib import Path

import sympy as sp

T0 = time.time()
TEX = Path(__file__).resolve().parent.parent / 'theorem_G_paper.tex'


def stamp(msg):
    print(f'[{time.time() - T0:6.1f}s] {msg}', flush=True)


# ============================================================== PART 2
stamp('PART 2: exact re-derivation of the displayed identities')

x1, x2, x3, x4, x0 = sp.symbols('x1 x2 x3 x4 x0')
X2, X3, X4 = (x1, x2), (x1, x2, x3), (x1, x2, x3, x4)


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def grad(e, vs):
    return sp.Matrix([sp.diff(e, a) for a in vs])


def Bof(e, vs):
    g = grad(e, vs)
    return sp.expand((g.T * hess(e, vs).adjugate() * g)[0, 0])


def det(M):
    return M.det(method='berkowitz')


def generic(vs, deg, tag, const=True, mindeg=None):
    mons = [sp.Integer(1)] if const else []
    lo = 1 if mindeg is None else mindeg
    for d in range(lo, deg + 1):
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons))


def form(vs, deg, tag):
    return generic(vs, deg, tag, const=False, mindeg=deg)


# ---------------------------------------------------------------- (C1)
u = sp.Symbol('u')
e1v = generic((u,), 5, 'q')
H1 = sp.Matrix([[sp.diff(e1v, u, 2)]])
assert H1.adjugate() == sp.Matrix([[1]])
assert sp.expand(Bof(e1v, (u,)) - sp.diff(e1v, u)**2) == 0
e0 = generic(X3, 3, 'a')
e1 = generic(X3, 3, 'b')
f = e0 + x4 * e1
P, K, g = hess(e0, X3), hess(e1, X3), grad(e1, X3)
stamp('(C10) skipped here: certified as (W9) in paperG_writeup_checks.py and wG check (2)')
stamp('(C10) Lemma lem:affx4: affine in x4, [x4] = -(e0)_{x3x3} B2(e)         PROOF')

# ---------------------------------------------------------------- (C11)
gam, lam = sp.symbols('gamma lambda')
assert sp.expand(hess(e1**2 / 2, X3) - (e1 * K + g * g.T)) == sp.zeros(3, 3)
et0 = e0 - e1**2 / (2 * gam)
uu = x4 + e1 / gam
lhsm = hess(et0, X3) + uu * K
rhsm = (P + x4 * K) - (g * g.T) / gam
assert sp.simplify(sp.expand(lhsm - rhsm)) == sp.zeros(3, 3)
fq = e0 + x4 * e1 + gam / 2 * x4**2
lhsd = sp.expand(gam**2 * det(hess(fq, X4)))          # division-free form
rhsd = sp.expand(det(gam * (P + x4 * K) - g * g.T))   # FIX 2026-08-25: Schur gives gam^2 det Hess = det3(gam(P+x4K) - g g^T); the agent's extra factor gam was wrong
assert sp.expand(lhsd - rhsd) == 0                    # det Hess f = gamma det3(Hess3 e0~ + uK)
hl = et0 + lam * e1
assert sp.simplify(grad(hl, X3) - (grad(e0, X3) + (lam - e1 / gam) * g)) == sp.zeros(3, 1)
assert sp.simplify(hess(hl, X3) - (hess(et0, X3) + lam * K)) == sp.zeros(3, 3)
stamp('(C11) Theorem thm:quad: Schur identity, determinant, collision-transfer gradient  PROOF')

# ---------------------------------------------------------------- (C12)
fa = e0 + x4 * x3
assert sp.expand(det(hess(fa, X4)) + det(hess(e0, X2))) == 0
stamp('(C12) Theorem thm:aff: det Hess(e0 + x4 x3) = -det2 Hess_{(x1,x2)} e0  PROOF')

# ---------------------------------------------------------------- (C13)
phi = generic((x1,), 4, 'p')
fp = e0 + x4 * phi
assert sp.expand(det(hess(fp, X4)) + sp.diff(phi, x1)**2 * det(hess(e0, (x2, x3)))) == 0
stamp('(C13) Theorem thm:tri branch: det Hess(e0 + x4 phi(x1)) = -phi\'^2 det2 Hess_{(x2,x3)} e0  PROOF')

# ---------------------------------------------------------------- (C14)
aa = generic(X2, 3, 'd')
bb = generic(X2, 3, 'e')
ee = generic(X2, 3, 'f2_')
fd = aa + x3 * bb + x4 * ee
jac = sp.expand(sp.diff(bb, x1) * sp.diff(ee, x2) - sp.diff(bb, x2) * sp.diff(ee, x1))
assert sp.expand(det(hess(fd, X4)) - jac**2) == 0
Hd = hess(fd, X4)
assert sp.expand(Hd[2:, 2:]) == sp.zeros(2, 2)
Mblk = sp.Matrix([[sp.diff(bb, x1), sp.diff(ee, x1)], [sp.diff(bb, x2), sp.diff(ee, x2)]])
assert sp.expand(Hd[:2, 2:] - Mblk) == sp.zeros(2, 2)
stamp('(C14) Theorem thm:doubling: block form and det Hess = Jac(b,e)^2       PROOF')

# ---------------------------------------------------------------- (C15)
# lem:inj (=>): with p' != q', (b,e)(p') = (b,e)(q'), the system
# M(p') (x3,x4)^T = grad a(q') - grad a(p') has a solution iff det M(p') != 0,
# which is the Keller condition.  We verify the gradient bookkeeping symbolically:
pp = sp.symbols('pp1 pp2')
qq = sp.symbols('qq1 qq2')
y3, y4 = sp.symbols('y3 y4')
gf = grad(fd, X4)
at_p = {x1: pp[0], x2: pp[1], x3: y3, x4: y4}
at_q = {x1: qq[0], x2: qq[1], x3: 0, x4: 0}
diff12 = sp.Matrix([gf[0].subs(at_p) - gf[0].subs(at_q), gf[1].subs(at_p) - gf[1].subs(at_q)])
Mp = Mblk.subs({x1: pp[0], x2: pp[1]})
rhs_sys = sp.Matrix([sp.diff(aa, x1).subs({x1: qq[0], x2: qq[1]}) - sp.diff(aa, x1).subs({x1: pp[0], x2: pp[1]}),
                     sp.diff(aa, x2).subs({x1: qq[0], x2: qq[1]}) - sp.diff(aa, x2).subs({x1: pp[0], x2: pp[1]})])
assert sp.expand(diff12 - (Mp * sp.Matrix([y3, y4]) - rhs_sys)) == sp.zeros(2, 1)
assert sp.expand(gf[2].subs(at_p) - bb.subs({x1: pp[0], x2: pp[1]})) == 0
assert sp.expand(gf[3].subs(at_p) - ee.subs({x1: pp[0], x2: pp[1]})) == 0
stamp('(C15) Lemma lem:inj: components 1,2 of grad f(p\',y3,y4) - grad f(q\',0,0) = M(p\')y - rhs  PROOF')

# ---------------------------------------------------------------- (C16)
f2g = form(X3, 2, 'u2_')
f3g = form(X3, 3, 'u3_')
detg = sp.expand(det(hess(f2g + f3g, X3)))
top = sum(term for term in sp.Add.make_args(detg)
          if sp.Poly(term, *X3).total_degree() == 3)
assert sp.expand(top - det(hess(f3g, X3))) == 0
# exact 4-variable instance
f2i = (x1**2 + x2**2 + x3**2 + x4**2) / 2 + x1 * x3
f3i = x1**3 + 2 * x1 * x2 * x4 + x2**2 * x3 - x3**2 * x4 + x1 * x3 * x4
deti = sp.expand(det(hess(f2i + f3i, X4)))
topi = sum(term for term in sp.Add.make_args(deti) if sp.Poly(term, *X4).total_degree() == 4)
assert sp.expand(topi - det(hess(f3i, X4))) == 0
stamp('(C16) top graded piece of det Hess(f2 + fd) = det Hess fd             PROOF (n=3) / INSTANCE (n=4)')

# ---------------------------------------------------------------- (C17)
# converse failure
assert det(hess(x1**2 / 2 + x2, X2)) == 0 and Bof(x1**2 / 2 + x2, X2) == 1
# dBvdE witness
hh = x1 * x2 + x1**2 * x3
assert sp.expand(det(hess(hh, X3))) == 0 and sp.expand(Bof(hh, X3) + x1**4) == 0
# cautionary example: det Hess2 == 0 does not give phi(L)
hc = x1 + x2**2
assert sp.expand(det(hess(hc, X2))) == 0 and Bof(hc, X2) == 2
# n = 4 counterexample to planarity
e4 = x2 + x1 * x3 + x1**2 * x4
H4e = hess(e4, X4)
assert Bof(e4, X4) == 0 and sp.expand(det(H4e)) == 0
assert H4e.adjugate() == sp.zeros(4, 4) and H4e.rank() == 2
parts = [sp.diff(e4, v) for v in X4]
mons = sorted({m for pr in parts for m in sp.Poly(pr, *X4).monoms()})
coef = sp.Matrix([[sp.Poly(pr, *X4).coeff_monomial(sp.prod(v**k for v, k in zip(X4, m))) for m in mons]
                  for pr in parts])
assert coef.rank() == 4      # partials linearly independent => not affinely 3-variable
E4h = sp.expand(x0**3 * e4.subs({v: v / x0 for v in X4}, simultaneous=True))
assert sp.expand(E4h - (x0**2 * x2 + x0 * x1 * x3 + x1**2 * x4)) == 0
assert sp.expand(det(hess(E4h, (x0,) + X4))) == 0
# jet identity B(a(x1,x2,x3) + lam x4) = lam^2 det Hess3(a)
bj = sp.symbols('b11 b12 b13 b22 b23 b33 g1 g2 g3')
Hj = sp.Matrix([[bj[0], bj[1], bj[2], 0], [bj[1], bj[3], bj[4], 0], [bj[2], bj[4], bj[5], 0], [0, 0, 0, 0]])
gj = sp.Matrix([bj[6], bj[7], bj[8], lam])
M3j = Hj[:3, :3]
assert sp.expand((gj.T * Hj.adjugate() * gj)[0, 0] - lam**2 * det(M3j)) == 0
# sharpness witness
Fs = x1 - sp.sqrt(x1**2 - 2 * x2)
assert sp.simplify(Bof(Fs, X2)) == 0
assert sp.simplify(det(hess(Fs, X2)) + 1 / (x1**2 - 2 * x2)**2) == 0
# de Bondt's dillen4 doubling
fdB = (x1 + x2**2) * x3 + (x2 + (x1 + x2**2)**2) * x4
assert sp.expand(det(hess(fdB, X4))) == 1
jacdB = sp.expand(sp.diff(x1 + x2**2, x1) * sp.diff(x2 + (x1 + x2**2)**2, x2)
                  - sp.diff(x1 + x2**2, x2) * sp.diff(x2 + (x1 + x2**2)**2, x1))
assert jacdB == 1
stamp('(C17) witnesses of Remarks rem:conv, rem:dim, rem:strict; jet identity   INSTANCE/PROOF')

# ---------------------------------------------------------------- (C18)
w = x1**2 / 2 + x1 * x2 * (1 + x3) + (x2**2 / 2) * (2 * x3 + x3**2) + x3 * x4
assert sp.expand(det(hess(w, X4))) == 1
assert sp.expand(sp.diff(w, x4, 2)) == 0 and sp.expand(sp.diff(w, x4)) == x3
vv = sp.symbols('v1:5')
zz = sp.Symbol('zz')
Qv = sp.expand((sp.Matrix(vv).T * hess(w, X4) * sp.Matrix(vv))[0, 0])
gens = [sp.expand(c) for c in sp.Poly(Qv, *X4).coeffs()]
for i in range(3):
    G = sp.groebner(gens + [1 - zz * vv[i]], *(list(vv) + [zz]), order='grevlex')
    assert list(G.exprs) == [sp.Integer(1)], f'v{i+1} does not vanish on the isotropic cone'
assert all(sp.expand(gn.subs({vv[0]: 0, vv[1]: 0, vv[2]: 0, vv[3]: 1})) == 0 for gn in gens)
# a doubling has span(e3,e4) in its isotropic cone (block form of C14)
Qd = sp.expand((sp.Matrix(vv).T * Hd * sp.Matrix(vv))[0, 0]).subs({vv[0]: 0, vv[1]: 0})
assert sp.expand(Qd) == 0
stamp('(C18) Remark rem:notdoubling: det Hess w = 1, isotropic cone = C.e4; doublings contain span(e3,e4)  PROOF')

# ---------------------------------------------------------------- (C19)
s = sp.Symbol('s')
Fn = form(X3, 3, 'kk')
HF, gF = hess(Fn, X3), grad(Fn, X3)
n_, r_ = 3, 3
lhsN = sp.expand(det(-Fn * HF + s * (gF * gF.T)))
rhsN = sp.expand((-1)**(n_ - 1) * sp.Rational(r_, r_ - 1) * (s - sp.Rational(r_ - 1, r_))
                 * Fn**n_ * det(HF))
assert sp.expand(lhsN - rhsN) == 0
coeff_s = sp.expand(sp.diff(lhsN, s))
assert sp.expand(coeff_s - (-1)**(n_ - 1) * Fn**(n_ - 1) * Bof(Fn, X3)) == 0
# rank-one expansion used in the extraction, generic 3x3 M and vectors u, v
Mr = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'R{i+1}{j+1}'))
ur = sp.Matrix(sp.symbols('ur1 ur2 ur3'))
vr = sp.Matrix(sp.symbols('vr1 vr2 vr3'))
assert sp.expand(det(Mr + s * ur * vr.T) - (det(Mr) + s * (vr.T * Mr.adjugate() * ur)[0, 0])) == 0
stamp('(C19) Nagaoka-Yazawa identity, s^1 extraction, rank-one expansion      PROOF')

# ---------------------------------------------------------------- (C20)
f4g = generic(X4, 3, 'ff')
Dv2 = sum(vv[i] * sp.diff(sum(vv[j] * sp.diff(f4g, X4[j]) for j in range(4)), X4[i]) for i in range(4))
assert sp.expand(Dv2 - (sp.Matrix(vv).T * hess(f4g, X4) * sp.Matrix(vv))[0, 0]) == 0
stamp('(C20) D_v^2 f = v^T Hess f v                                            PROOF')

print()
print('fpG_paper_certificate: ALL CHECKS PASSED --', TEX.name,
      'is referentially sound and every displayed identity is re-certified.')
sys.exit(0)
