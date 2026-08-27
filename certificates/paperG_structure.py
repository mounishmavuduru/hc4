# paperG_structure.py  --  certificates for theorem_G_paper.tex (agent key: paperG)
#
# Fail-closed, exact/symbolic only.  These certify the NEW material of
# ../theorem_G_paper.tex that is not already covered by another agent's script:
#
#   (P1) bordered identity  det[[M,g],[g^T,0]] = -g^T adj(M) g   (3x3 M, generic)  PROOF
#   (P2) affine-pivot potential with PLANAR pivot coefficient: the x4-coefficient of
#        det Hess f equals  -(d^2 e0/dx3^2) * B_2(e1)            (generic)          PROOF
#   (P3) affine-pivot potential with LINEAR-VARIABLE pivot coefficient e1 = h(x1):
#        det Hess f = -(h')^2 (P22 P33 - P23^2), free of x4      (generic)          PROOF
#   (P4) COUNTEREXAMPLE to the claim "the affine-pivot class consists exactly of Meng
#        doublings plus an inert planar potential": an explicit w with det Hess w = 1,
#        an affine pivot, and NO 2-dimensional totally isotropic subspace of Hess w
#        (hence not affinely equivalent to any doubling).                           PROOF
#   (P5) doubling determinant  det[[A,B],[B^T,0]] = det(B)^2     (generic 2x2 blocks) PROOF
#   (P6) Theorem C converse for ARBITRARY inert potential a: a collision of the planar
#        Keller map lifts to a collision of grad f for every a.  (generic)          PROOF
#   (P7) homogenisation identity (Identity H) in n = 2, degrees 2,3,4 (generic)     PROOF
#   (P8) translation invariance of B in n = 2 (generic cubic)                       PROOF
#   (P9) Identity E (Nagaoka-Yazawa) for binary forms of degrees 2,3,4 (generic)    PROOF
#  (P10) decisive reduced-ansatz test of Theorem G in n = 2 (consistency check).
#
# "PROOF" = the checked expression is a polynomial identity in independent symbolic
# entries of the jets involved, so the check settles all instances of that shape.

import itertools
import sympy as sp

FAILS = []


def ok(tag, msg):
    print(f'  [OK] {tag}: {msg}')


def sym_sym_matrix(n, tag):
    """generic symmetric n x n matrix with independent entries"""
    e = {}
    for i in range(n):
        for j in range(i, n):
            e[(i, j)] = sp.Symbol(f'{tag}{i+1}{j+1}')
            e[(j, i)] = e[(i, j)]
    return sp.Matrix(n, n, lambda i, j: e[(i, j)])


def hess(f, vs):
    return sp.Matrix([[sp.diff(f, a, b) for b in vs] for a in vs])


def grad(f, vs):
    return sp.Matrix([sp.diff(f, a) for a in vs])


def Bform(f, vs):
    H = hess(f, vs)
    g = grad(f, vs)
    return sp.expand((g.T * H.adjugate() * g)[0, 0])


def dethess(f, vs):
    return sp.expand(hess(f, vs).det(method='berkowitz'))


def gen_poly(vs, deg, tag, with_const=True):
    mons = [sp.Integer(1)] if with_const else []
    for d in range(1, deg + 1):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, d)},
                       key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons)), list(cs)


def gen_form(vs, deg, tag):
    mons = sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, deg)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons)), list(cs)


x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X4 = (x1, x2, x3, x4)
X3 = (x1, x2, x3)
X2 = (x1, x2)

print('=' * 78)
print('(P1) bordered identity for the 4-variable affine-pivot Hessian  [PROOF]')
M = sym_sym_matrix(3, 'm')
gg = sp.Matrix(sp.symbols('gg1 gg2 gg3'))
Big = sp.Matrix(sp.BlockMatrix([[M, gg], [gg.T, sp.zeros(1, 1)]]))
lhs = sp.expand(Big.det(method='berkowitz'))
rhs = sp.expand(-(gg.T * M.adjugate() * gg)[0, 0])
assert sp.expand(lhs - rhs) == 0
ok('P1', 'det[[M,g],[g^T,0]] == -g^T adj(M) g for generic symmetric 3x3 M')

print()
print('=' * 78)
print('(P2) planar pivot coefficient: x4-coefficient of det Hess f  [PROOF]')
# f = e0(x1,x2,x3) + x4*e1(x1,x2), so
#   Hess f = [[P + x4 K, g],[g^T, 0]],  K = Hess3(e1) has zero 3rd row/col, g3 = 0.
P = sym_sym_matrix(3, 'p')
k11, k12, k22 = sp.symbols('k11 k12 k22')
K = sp.Matrix([[k11, k12, 0], [k12, k22, 0], [0, 0, 0]])
g1, g2 = sp.symbols('g1 g2')
g = sp.Matrix([g1, g2, 0])
Mx = P + x4 * K
det4 = sp.expand(sp.Matrix(sp.BlockMatrix([[Mx, g], [g.T, sp.zeros(1, 1)]])).det(method='berkowitz'))
assert sp.degree(sp.Poly(det4, x4), x4) <= 1, 'det should be affine in x4'
coef_x4 = sp.expand(sp.Poly(det4, x4).coeff_monomial(x4))
q = sp.expand(g1**2 * k22 - 2 * g1 * g2 * k12 + g2**2 * k11)      # = B_2(e1)
assert sp.expand(coef_x4 + P[2, 2] * q) == 0
ok('P2', 'coeff of x4 in det Hess f is  -P33 * (g1^2 K22 - 2 g1 g2 K12 + g2^2 K11)')
# and that q really is the 2-variable bordered Hessian of a planar e1
e1gen, _ = gen_poly(X2, 3, 'e1c')
H2 = hess(e1gen, X2)
g2v = grad(e1gen, X2)
B2 = sp.expand((g2v.T * H2.adjugate() * g2v)[0, 0])
q_sub = sp.expand(q.subs({g1: sp.diff(e1gen, x1), g2: sp.diff(e1gen, x2),
                         k11: sp.diff(e1gen, x1, 2), k12: sp.diff(e1gen, x1, x2),
                         k22: sp.diff(e1gen, x2, 2)}))
assert sp.expand(q_sub - B2) == 0
ok('P2', 'that coefficient q equals B_2(e1), the 2-variable bordered Hessian (generic cubic e1)')

print()
print('=' * 78)
print('(P3) pivot coefficient in one variable, e1 = h(x1)  [PROOF]')
hp, hpp = sp.symbols('hp hpp')                # h'(x1), h''(x1)
Kh = sp.Matrix([[hpp, 0, 0], [0, 0, 0], [0, 0, 0]])
gh = sp.Matrix([hp, 0, 0])
deth = sp.expand(sp.Matrix(sp.BlockMatrix([[P + x4 * Kh, gh], [gh.T, sp.zeros(1, 1)]])).det(
    method='berkowitz'))
assert sp.expand(sp.diff(deth, x4)) == 0
assert sp.expand(deth + hp**2 * (P[1, 1] * P[2, 2] - P[1, 2]**2)) == 0
ok('P3', "det Hess f = -(h')^2 (P22 P33 - P23^2), independent of x4  =>  (h')^2 | c")

print()
print('=' * 78)
print('(P4) an affine-pivot potential that is NOT a Meng doubling  [PROOF]')
# w = x1^2/2 + x1 x2 (1+x3) + (x2^2/2)(2 x3 + x3^2) + x3 x4
w = x1**2 / 2 + x1 * x2 * (1 + x3) + (x2**2 / 2) * (2 * x3 + x3**2) + x3 * x4
dw = dethess(w, X4)
assert dw == 1, f'det Hess w = {dw}'
ok('P4', 'det Hess w == 1 exactly (w is a Hessian-Keller potential)')
# affine pivot: D_{e4}^2 w == 0, and the pivot coefficient e1 = x3 has degree 1
assert sp.expand(sp.diff(w, x4, 2)) == 0
e1_of_w = sp.expand(sp.diff(w, x4))
assert e1_of_w == x3
ok('P4', 'w = e0(x1,x2,x3) + x4*e1 with e1 = x3: an affine pivot with deg e1 = 1')
# the isotropic cone {v : D_v^2 w == 0 identically} is the LINE span(e4)
v = sp.symbols('v1:5')
Hw = hess(w, X4)
Qv = sp.expand((sp.Matrix(v).T * Hw * sp.Matrix(v))[0, 0])
gens = [sp.expand(c) for c in sp.Poly(Qv, *X4).coeffs()]
G = sp.groebner(gens, *v, order='grevlex')
for i in range(3):
    assert G.reduce(v[i]**2)[1] == 0, f'v{i+1}^2 not in the ideal'
ok('P4', 'every v with D_v^2 w == 0 satisfies v1^2 = v2^2 = v3^2 = 0, hence v in span(e4)')
assert all(sp.expand(gn.subs({v[0]: 0, v[1]: 0, v[2]: 0, v[3]: 1})) == 0 for gn in gens)
ok('P4', 'and e4 itself is isotropic: the cone is exactly the line C.e4 (dim 1)')
print('       => Hess w admits NO 2-dimensional totally isotropic subspace, so w is')
print('          NOT affinely equivalent to  a(y1,y2) + y3 b(y1,y2) + y4 e(y1,y2).')
print('          The "affine-pivot class = Meng doublings" claim is therefore FALSE;')
print('          the correct dichotomy is deg e1 = 1 (class AP0) vs deg e1 >= 2.')

print()
print('=' * 78)
print('(P5) doubling determinant  [PROOF]')
A2 = sym_sym_matrix(2, 'a')
Bm = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'b{i+1}{j+1}'))
Dbl = sp.Matrix(sp.BlockMatrix([[A2, Bm.T], [Bm, sp.zeros(2, 2)]]))
assert sp.expand(Dbl.det(method='berkowitz') - Bm.det()**2) == 0
ok('P5', 'det[[A, B^T],[B, 0]] == (det B)^2 for generic symmetric A and generic B')

print()
print('=' * 78)
print('(P6) converse of Theorem C for an ARBITRARY inert potential a  [PROOF]')
# f = a(y) + y3 b(y) + y4 e(y),  y = (y1,y2).  Suppose (b,e)(p) = (b,e)(qq), p != qq.
# Claim: for every u in C^2 there is a w in C^2 with grad f (p,u) = grad f (qq,w),
# namely w = (Bq^T)^{-1} ( grad a(p) - grad a(qq) + Bp^T u ),  Bx = jacobian of (b,e).
ap = sp.Matrix(sp.symbols('ap1 ap2'))          # grad a at p
aq = sp.Matrix(sp.symbols('aq1 aq2'))          # grad a at qq
Bp = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'Bp{i}{j}'))
Bq = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'Bq{i}{j}'))
u = sp.Matrix(sp.symbols('u1 u2'))
wv = Bq.T.inv() * (ap - aq + Bp.T * u)
lhs6 = ap + Bp.T * u
rhs6 = aq + Bq.T * wv
assert sp.simplify(sp.expand(lhs6 - rhs6)) == sp.zeros(2, 1)
ok('P6', 'the lifted collision exists for every a (det Bq != 0 is all that is used)')

print()
print('=' * 78)
print('(P7) homogenisation identity (Identity H) in n = 2  [PROOF, degrees 2,3,4]')
x0 = sp.Symbol('x0')
for d in (2, 3, 4):
    e, _ = gen_poly(X2, d, f'h{d}_')
    E = sp.expand(x0**d * e.subs({x1: x1 / x0, x2: x2 / x0}))
    E = sp.simplify(sp.expand(E))
    lhs7 = sp.expand(dethess(E, (x0, x1, x2)).subs({x0: 1}))
    rhs7 = sp.expand(d * (d - 1) * e * dethess(e, X2) - (d - 1)**2 * Bform(e, X2))
    assert sp.expand(lhs7 - rhs7) == 0, f'Identity H failed in n=2, d={d}'
    ok('P7', f'd = {d}: det Hess_3(E)|_(x0=1) == d(d-1) e det Hess_2 e - (d-1)^2 B(e)')

print()
print('=' * 78)
print('(P8) translation invariance of B in n = 2  [PROOF, generic cubic]')
t1, t2 = sp.symbols('t1 t2')
e, _ = gen_poly(X2, 3, 'tr_')
et = sp.expand(e.subs({x1: x1 + t1, x2: x2 + t2}, simultaneous=True))
assert sp.expand(Bform(et, X2) - Bform(e, X2).subs({x1: x1 + t1, x2: x2 + t2},
                                                   simultaneous=True)) == 0
ok('P8', 'B(e o tau) == B(e) o tau for every translation tau')

print()
print('=' * 78)
print('(P9) Identity E of Nagaoka-Yazawa (cited, not claimed) for binary forms  [PROOF]')
for d in (2, 3, 4):
    F, _ = gen_form(X2, d, f'ny{d}_')
    assert sp.expand(Bform(F, X2) - sp.Rational(d, d - 1) * F * dethess(F, X2)) == 0
    ok('P9', f'd = {d}: B(F) == d/(d-1) * F * det Hess F')

print()
print('=' * 78)
print('(P10) reduced-ansatz test of THEOREM G in n = 2 (consistency check)')
# If B(e) == 0 then the top-degree part of B is B(e_d) = d/(d-1) e_d det Hess e_d
# (Identity E), so det Hess(e_d) == 0, so by Hesse/Gordan-Noether for binary forms
# e_d = (linear)^d; after a linear change and scaling, e_d = x1^d.  Translations kill
# the constant term.  So this ansatz is complete up to affine equivalence.
zz = sp.Symbol('zz')
for d in (3, 4):
    low, pars = gen_poly(X2, d - 1, f'lw{d}_', with_const=False)
    e = x1**d + low
    Bexpr = Bform(e, X2)
    Beqs = [sp.expand(c) for c in sp.Poly(Bexpr, *X2).coeffs()] if sp.expand(Bexpr) != 0 else []
    Dexpr = dethess(e, X2)
    Deqs = [sp.expand(c) for c in sp.Poly(Dexpr, *X2).coeffs()] if sp.expand(Dexpr) != 0 else []
    bad = []
    for dc in Deqs:
        if dc == 0:
            continue
        Gr = sp.groebner(list(Beqs) + [1 - zz * dc], *(pars + [zz]), order='grevlex')
        if list(Gr.exprs) != [sp.Integer(1)]:
            bad.append(dc)
    assert not bad, f'THEOREM G CONTRADICTED in n=2, degree {d}: {bad[:2]}'
    ok('P10', f'n = 2, deg {d}: every coefficient of det Hess e vanishes on V(B(e)) '
              f'({len(pars)} params, {len(Beqs)} equations)')

print()
print('=' * 78)
print('ALL paperG CHECKS PASSED.')
print('Headline consequence certified here: the affine-pivot class of HC_4 is NOT')
print('exhausted by Meng doublings (P4); the correct statement is the degree dichotomy')
print('deg e1 = 1  (AP0, closed unconditionally)  vs  deg e1 >= 2  (doubling class),')
print('and the identities (P2),(P3) are exactly what forces it.')
