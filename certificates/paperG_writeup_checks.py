# paperG_writeup_checks.py -- certificates for ../theorem_G_paper.tex   (agent key: paperG)
#
# Fail-closed, exact/symbolic only.  Every claim of the write-up that is an
# ALGEBRAIC identity is certified here or in the project certificate named in
# the comment; the analytic core of Theorem G (local holomorphic inverse + the
# identity theorem) is a hand proof and is NOT claimed to be machine-checked.
#
# "PROOF" marks a check whose expression is a polynomial identity in independent
# symbolic entries of the jets involved, so the check settles every instance of
# that shape.  Heavier fully-generic versions of (W3), (W7) already exist in
# ../certificates/identityE_is_known.py, cone_leading_form.py and
# homogenization_identity.py; they are not duplicated here.
#
#  (W1) H adj(H) g = det(H) g  and  g^T adj(H) g = B                    PROOF (n=2,3,4)
#  (W2) affine covariance  B(e o T) = (det T)^2 (B e) o T               PROOF
#  (W3) THEOREM G for HOMOGENEOUS e, unconditionally, from Identity E   PROOF (n=3, d=2,3,4)
#  (W4) Nagaoka-Yazawa Prop. 2.3, and Identity E as its s-linear
#       coefficient (attribution certificate)                           PROOF (n=3, d=2,3)
#  (W5) the Legendre chain of the proof of Theorem G on an EXACT example exact example
#  (W6) the "no positive power of t" step                               PROOF (n=2, deg 4)
#  (W7) Identity H (homogenisation) in n = 3                            PROOF (d=2)
#  (W8) affine-pivot expansion; [x4^2] det Hess f = -B(e1)              PROOF
#  (W9) planar pivot coefficient: [x4] det Hess f = -(e0)_{x3x3} B_2(e1) PROOF
# (W10) one-variable pivot coefficient e1 = h(x1)                       PROOF
# (W11) AP0: f = e0(x') + x4 x3 has det Hess f = -det_2 Hess_{(x1,x2)}e0 PROOF
# (W12) doubling determinant + collision transfer both ways, any inert a PROOF
# (W13) the CORRECTION: an explicit affine-pivot potential w with
#       det Hess w = 1 that is NOT affinely equivalent to a doubling    PROOF
# (W14) two instances quoted in the paper                               exact instances

import itertools
import time
import sympy as sp

x0, x1, x2, x3, x4 = sp.symbols('x0 x1 x2 x3 x4')
X2 = (x1, x2)
X3 = (x1, x2, x3)
X4 = (x1, x2, x3, x4)
T0 = time.time()


def hess(f, vs):
    return sp.Matrix([[sp.diff(f, a, b) for b in vs] for a in vs])


def grad(f, vs):
    return sp.Matrix([sp.diff(f, a) for a in vs])


def Bform(f, vs):
    return sp.expand((grad(f, vs).T * hess(f, vs).adjugate() * grad(f, vs))[0, 0])


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


def sym_matrix(n, tag):
    e = {}
    for i in range(n):
        for j in range(i, n):
            e[(i, j)] = sp.Symbol(f'{tag}{i+1}{j+1}')
            e[(j, i)] = e[(i, j)]
    return sp.Matrix(n, n, lambda i, j: e[(i, j)])


def ok(tag, msg):
    print(f'  [OK] {tag}: {msg}   ({time.time()-T0:.1f}s)', flush=True)


BAR = '=' * 78

# ---------------------------------------------------------------- (W1)
print(BAR, flush=True)
print('(W1) the two structural identities behind V = H^{-1} grad e   [PROOF]', flush=True)
for n in (2, 3, 4):
    M = sym_matrix(n, f'M{n}_')
    g = sp.Matrix(sp.symbols(f'gg{n}_1:{n+1}'))
    assert sp.expand(M * M.adjugate() * g - M.det(method='berkowitz') * g) == sp.zeros(n, 1)
    assert sp.expand((g.T * (M.adjugate() * g))[0, 0] - (g.T * M.adjugate() * g)[0, 0]) == 0
    ok('W1', f'n={n}: M adj(M) g == det(M) g for a generic symmetric M and generic g')
print('       With M = Hess e, g = grad e this says: W := adj(H)g satisfies H W = det(H) g,')
print('       i.e. V := H^{-1}g solves (D grad e) V = grad e on {det H != 0}, and')
print('       <g, V> = B/det H.  Both are linear-algebra identities, valid for every e.',
      flush=True)

# ---------------------------------------------------------------- (W2)
print()
print(BAR)
print('(W2) affine covariance of B: B(e o T) == (det T)^2 (B e) o T   [PROOF]', flush=True)
T2 = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'T{i+1}{j+1}'))
t2 = sp.Matrix(sp.symbols('tt1 tt2'))
e2, _ = gen_poly(X2, 3, 'cov2_')
yv = T2 * sp.Matrix([x1, x2]) + t2
sub2 = {x1: yv[0], x2: yv[1]}
assert sp.expand(Bform(sp.expand(e2.subs(sub2, simultaneous=True)), X2)
                 - T2.det()**2 * Bform(e2, X2).subs(sub2, simultaneous=True)) == 0
ok('W2', 'n=2: fully symbolic affine T (matrix + translation), generic cubic e')
T3 = sp.Matrix([[1, 2, 0], [0, 1, 3], [1, 0, 1]])
assert T3.det() == 7
e3, _ = gen_poly(X3, 2, 'cov3_')
y3 = T3 * sp.Matrix(list(X3)) + sp.Matrix([2, -1, 5])
sub3 = {x1: y3[0], x2: y3[1], x3: y3[2]}
assert sp.expand(Bform(sp.expand(e3.subs(sub3, simultaneous=True)), X3)
                 - 49 * Bform(e3, X3).subs(sub3, simultaneous=True)) == 0
ok('W2', 'n=3: affine T with det T = 7, generic quadratic e')
print('       => "B(e) == 0" and "e is affinely k-variable" are both affine invariants.',
      flush=True)

# ---------------------------------------------------------------- (W3)
print()
print(BAR)
print('(W3) THEOREM G for HOMOGENEOUS e, unconditionally (Identity E)   [PROOF]', flush=True)
for d in (2, 3, 4):
    F, _ = gen_form(X3, d, f'hg{d}_')
    assert sp.expand(Bform(F, X3) - sp.Rational(d, d - 1) * F * dethess(F, X3)) == 0
    ok('W3', f'n=3, d={d}: B(F) == d/(d-1) F det Hess F')
print('       C[x] is a domain and F != 0, so B(F) == 0 <=> det Hess F == 0:')
print('       the homogeneous case of Theorem G is a corollary of the KNOWN')
print('       identity of Nagaoka-Yazawa.  (n = 4 and higher degrees: see')
print('       identityE_is_known.py and cone_leading_form.py.)', flush=True)

# ---------------------------------------------------------------- (W4)
print()
print(BAR)
print('(W4) attribution: Identity E is the s-linear coefficient of Nagaoka-Yazawa')
print('     arXiv:1904.01800, Prop. 2.3 (source label identity1)          [PROOF]', flush=True)
s = sp.Symbol('s')
n = 3
for d in (2, 3):
    F, _ = gen_form(X3, d, f'ny{d}_')
    H = hess(F, X3)
    grow = sp.Matrix([[sp.diff(F, a) for a in X3]])
    lhs = sp.expand((-F * H + s * (grow.T * grow)).det(method='berkowitz'))
    rhs = sp.expand((-1)**(n - 1) * sp.Rational(d, d - 1) * (s - sp.Rational(d - 1, d))
                    * F**n * H.det(method='berkowitz'))
    assert sp.expand(lhs - rhs) == 0
    assert sp.expand(sp.diff(lhs, s) - (-1)**(n - 1) * F**(n - 1) * Bform(F, X3)) == 0
    ok('W4', f'n=3, d={d}: NY identity holds; its s-coefficient is (-1)^(n-1) F^(n-1) B(F)')

# ---------------------------------------------------------------- (W5)
print()
print(BAR)
print('(W5) the Legendre chain of the proof of Theorem G, on an exact example', flush=True)
p1, p2 = sp.symbols('p1 p2')
e_ex = x1 * x2 + x1**3 / 3
H_ex = hess(e_ex, X2)
assert sp.expand(H_ex.det()) == -1
psi = sp.Matrix([p2, p1 - p2**2])
subpsi = {x1: psi[0], x2: psi[1]}
assert sp.expand(grad(e_ex, X2).subs(subpsi, simultaneous=True)
                 - sp.Matrix([p1, p2])) == sp.zeros(2, 1)
ok('W5', 'psi is a global polynomial inverse of grad e for e = x1 x2 + x1^3/3')
e_psi = sp.expand(e_ex.subs(subpsi, simultaneous=True))
L = sp.expand(p1 * psi[0] + p2 * psi[1] - e_psi)
assert sp.expand(sp.Matrix([sp.diff(L, p1), sp.diff(L, p2)]) - psi) == sp.zeros(2, 1)
ok('W5', 'grad L == psi  (L the Legendre transform)')
EL = sp.expand(p1 * sp.diff(L, p1) + p2 * sp.diff(L, p2))
assert sp.expand(EL - L - e_psi) == 0
ok('W5', 'E L - L == e o psi')
Bex = Bform(e_ex, X2)
assert sp.expand(p1 * sp.diff(e_psi, p1) + p2 * sp.diff(e_psi, p2)
                 - (Bex / H_ex.det()).subs(subpsi, simultaneous=True)) == 0
ok('W5', 'E(e o psi) == (B / det Hess e) o psi')
assert sp.expand(Bex) != 0
print(f'       B(e) = {sp.expand(Bex)} is NOT zero here, so this is a live test of the')
print('       identity that converts B == 0 into (E^2 - E) L = 0.', flush=True)

# ---------------------------------------------------------------- (W6)
print()
print(BAR)
print('(W6) "no positive power of t": grad e(a/t + b) lies in C[a,b][1/t]  [PROOF]', flush=True)
tt = sp.Symbol('tt')
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2')
e6, _ = gen_poly(X2, 4, 'w6_')
for comp in grad(e6, X2):
    expr = sp.expand(comp.subs({x1: a1 / tt + b1, x2: a2 / tt + b2}, simultaneous=True))
    assert sp.degree(sp.Poly(sp.expand(expr * tt**4), tt), tt) <= 4
ok('W6', 'generic quartic e in 2 variables: no positive power of t occurs')
print('       => comparing coefficients of t^1 in grad e(a/t + b) = t p forces p = 0.',
      flush=True)

# ---------------------------------------------------------------- (W7)
print()
print(BAR)
print('(W7) Identity H (homogenisation) in n = 3   [PROOF, d = 2; d = 3 in',
      'homogenization_identity.py]', flush=True)
for d in (2,):
    e, _ = gen_poly(X3, d, f'ih{d}_')
    E = sp.expand(sp.simplify(x0**d * e.subs({x1: x1 / x0, x2: x2 / x0, x3: x3 / x0},
                                             simultaneous=True)))
    lhs7 = sp.expand(dethess(E, (x0, x1, x2, x3)).subs({x0: 1}))
    rhs7 = sp.expand(d * (d - 1) * e * dethess(e, X3) - (d - 1)**2 * Bform(e, X3))
    assert sp.expand(lhs7 - rhs7) == 0
    ok('W7', f'd = {d}: det Hess_4(E)|_(x0=1) == d(d-1) e det Hess_3 e - (d-1)^2 B(e)')

# ---------------------------------------------------------------- (W8)
print()
print(BAR)
print('(W8) affine-pivot expansion and constraint (c2)   [PROOF]', flush=True)
Pm = sym_matrix(3, 'P')
Km = sym_matrix(3, 'K')
gv = sp.Matrix(sp.symbols('gv1 gv2 gv3'))
Big = sp.Matrix(sp.BlockMatrix([[Pm + x4 * Km, gv], [gv.T, sp.zeros(1, 1)]]))
detBig = sp.expand(Big.det(method='berkowitz'))
assert sp.expand(detBig + (gv.T * (Pm + x4 * Km).adjugate() * gv)[0, 0]) == 0
ok('W8', 'det[[P + x4 K, g],[g^T, 0]] == -g^T adj(P + x4 K) g (generic symmetric P, K)')
pol = sp.Poly(detBig, x4)
assert pol.degree() <= 2
assert sp.expand(pol.coeff_monomial(x4**2) + (gv.T * Km.adjugate() * gv)[0, 0]) == 0
ok('W8', '[x4^2] det Hess f == -g^T adj(K) g == -B(e1): constraint (c2) IS B(e1) == 0')
e1g, _ = gen_poly(X3, 3, 'e1g_')
assert sp.expand((grad(e1g, X3).T * hess(e1g, X3).adjugate() * grad(e1g, X3))[0, 0]
                 - Bform(e1g, X3)) == 0
ok('W8', 'and g^T adj(K) g really is B(e1) for a generic cubic e1 in 3 variables')

# ---------------------------------------------------------------- (W9)
print()
print(BAR)
print('(W9) planar pivot coefficient: [x4] det Hess f == -(e0)_{x3x3} B_2(e1)  [PROOF]',
      flush=True)
k11, k12, k22 = sp.symbols('k11 k12 k22')
Kp = sp.Matrix([[k11, k12, 0], [k12, k22, 0], [0, 0, 0]])
gq1, gq2 = sp.symbols('gq1 gq2')
gp = sp.Matrix([gq1, gq2, 0])
det4 = sp.expand(sp.Matrix(sp.BlockMatrix([[Pm + x4 * Kp, gp],
                                           [gp.T, sp.zeros(1, 1)]])).det(method='berkowitz'))
assert sp.Poly(det4, x4).degree() <= 1
qq = sp.expand(gq1**2 * k22 - 2 * gq1 * gq2 * k12 + gq2**2 * k11)
assert sp.expand(sp.Poly(det4, x4).coeff_monomial(x4) + Pm[2, 2] * qq) == 0
ok('W9', 'generic symmetric P, planar K and planar g')
e1p, _ = gen_poly(X2, 4, 'e1p_')
assert sp.expand(qq.subs({gq1: sp.diff(e1p, x1), gq2: sp.diff(e1p, x2),
                          k11: sp.diff(e1p, x1, 2), k12: sp.diff(e1p, x1, x2),
                          k22: sp.diff(e1p, x2, 2)}) - Bform(e1p, X2)) == 0
ok('W9', 'the coefficient q equals B_2(e1), the planar bordered Hessian (generic quartic e1)')
print('       => det Hess f constant forces (e0)_{x3x3} * B_2(e1) == 0 in the domain C[x].',
      flush=True)

# ---------------------------------------------------------------- (W10)
print()
print(BAR)
print("(W10) one-variable pivot coefficient e1 = h(x1)   [PROOF]", flush=True)
hp, hpp = sp.symbols('hp hpp')
Kh = sp.Matrix([[hpp, 0, 0], [0, 0, 0], [0, 0, 0]])
gh = sp.Matrix([hp, 0, 0])
deth = sp.expand(sp.Matrix(sp.BlockMatrix([[Pm + x4 * Kh, gh],
                                           [gh.T, sp.zeros(1, 1)]])).det(method='berkowitz'))
assert sp.expand(sp.diff(deth, x4)) == 0
assert sp.expand(deth + hp**2 * (Pm[1, 1] * Pm[2, 2] - Pm[1, 2]**2)) == 0
ok('W10', "det Hess f == -(h')^2 (P22 P33 - P23^2), independent of x4")
print("       => if det Hess f is a nonzero constant then (h')^2 is a unit of C[x], so h")
print('          is affine: the branch B_2(e1) == 0 collapses to deg e1 = 1.', flush=True)

# ---------------------------------------------------------------- (W11)
print()
print(BAR)
print('(W11) class AP0: f = e0(x1,x2,x3) + x4 x3   [PROOF]', flush=True)
g3 = sp.Matrix([0, 0, 1])
det11 = sp.expand(sp.Matrix(sp.BlockMatrix([[Pm, g3], [g3.T, sp.zeros(1, 1)]]))
                  .det(method='berkowitz'))
assert sp.expand(det11 + (Pm[0, 0] * Pm[1, 1] - Pm[0, 1]**2)) == 0
ok('W11', 'det Hess f == -(P11 P22 - P12^2) == -det_2 Hess_{(x1,x2)} e0, generic P')
e0a, _ = gen_poly(X3, 3, 'e0a_')
fa = sp.expand(e0a + x4 * x3)
assert sp.expand(dethess(fa, X4)
                 + (sp.diff(e0a, x1, 2) * sp.diff(e0a, x2, 2) - sp.diff(e0a, x1, x2)**2)) == 0
ok('W11', 'same on a generic cubic e0 (polynomial level)')
print('       => every x3-fibre of e0 is a planar potential with constant nonzero Hessian')
print('          determinant, and Dillen applies fibrewise.', flush=True)

# ---------------------------------------------------------------- (W12)
print()
print(BAR)
print('(W12) doubling: determinant and collision transfer both ways   [PROOF]', flush=True)
A2 = sym_matrix(2, 'A')
Bm = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'B{i+1}{j+1}'))
assert sp.expand(sp.Matrix(sp.BlockMatrix([[A2, Bm.T], [Bm, sp.zeros(2, 2)]]))
                 .det(method='berkowitz') - Bm.det()**2) == 0
ok('W12', 'det[[A, B^T],[B, 0]] == (det B)^2 (generic symmetric A, generic B)')
ag, _ = gen_poly(X2, 3, 'dba_')
bg, _ = gen_poly(X2, 2, 'dbb_')
eg, _ = gen_poly(X2, 2, 'dbe_')
fdb = sp.expand(ag + x3 * bg + x4 * eg)
Jbe = sp.Matrix([[sp.diff(bg, x1), sp.diff(bg, x2)], [sp.diff(eg, x1), sp.diff(eg, x2)]])
assert sp.expand(dethess(fdb, X4) - sp.expand(Jbe.det()**2)) == 0
ok('W12', 'polynomial level: det Hess (a + x3 b + x4 e) == Jac(b,e)^2, generic a, b, e')
ap = sp.Matrix(sp.symbols('ap1 ap2'))
aq = sp.Matrix(sp.symbols('aq1 aq2'))
Bp = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'Bp{i}{j}'))
Bq = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'Bq{i}{j}'))
uu = sp.Matrix(sp.symbols('uu1 uu2'))
wv = Bq.T.inv() * (ap - aq + Bp.T * uu)
assert sp.simplify(sp.expand((ap + Bp.T * uu) - (aq + Bq.T * wv))) == sp.zeros(2, 1)
ok('W12', 'a collision of (b,e) lifts to a collision of grad f for EVERY inert a')
print('       (only det B_q != 0 is used, which the Keller condition supplies); conversely')
print('       components 3 and 4 of grad f are b and e, so a collision of grad f gives one')
print('       of (b,e).  Hence grad f is injective IFF (b,e) is, for every a.', flush=True)

# ---------------------------------------------------------------- (W13)
print()
print(BAR)
print('(W13) CORRECTION: an affine-pivot potential that is NOT a Meng doubling  [PROOF]',
      flush=True)
w = x1**2 / 2 + x1 * x2 * (1 + x3) + (x2**2 / 2) * (2 * x3 + x3**2) + x3 * x4
assert dethess(w, X4) == 1
ok('W13', 'det Hess w == 1 for w = x1^2/2 + x1 x2 (1+x3) + (x2^2/2)(2x3+x3^2) + x3 x4')
assert sp.expand(sp.diff(w, x4, 2)) == 0 and sp.expand(sp.diff(w, x4)) == x3
ok('W13', 'w has an affine pivot with pivot coefficient e1 = x3 (degree 1): class AP0')
vv = sp.symbols('v1:5')
zz = sp.Symbol('zz')
Qv = sp.expand((sp.Matrix(vv).T * hess(w, X4) * sp.Matrix(vv))[0, 0])
gens = [sp.expand(c) for c in sp.Poly(Qv, *X4).coeffs()]
# NOTE: v1^2 lies in the RADICAL of this ideal but NOT in the ideal itself (a
# Groebner basis is [v1^2 + 2 v3 v4, v1 v2, v2^2, v1 v3, v2 v3, v3^2]), so the
# plain ideal-membership test used in paperG_structure.py (P4) FAILS here.  The
# correct test is radical membership, certified by Rabinowitsch below.
for i in range(3):
    R = sp.groebner(gens + [1 - zz * vv[i]], *(list(vv) + [zz]), order='grevlex')
    assert list(R.exprs) == [sp.Integer(1)], f'v{i+1} does not vanish on the isotropic cone'
assert all(sp.expand(gn.subs({vv[0]: 0, vv[1]: 0, vv[2]: 0, vv[3]: 1})) == 0 for gn in gens)
ok('W13', 'v1, v2, v3 all vanish on {v : D_v^2 w == 0} (Rabinowitsch radical membership),')
print('       and e4 is isotropic: the cone is exactly the line C.e4, of dimension 1.')
print('       A Meng doubling a(y1,y2) + y3 b + y4 e has the 2-PLANE span(e3,e4) inside its')
print('       isotropic cone, and that cone is an affine invariant; so w is NOT affinely')
print('       equivalent to any doubling.  Hence "the affine-pivot class consists exactly of')
print('       Meng doublings" is FALSE: the correct split is deg e1 = 1 vs deg e1 >= 2.',
      flush=True)

# ---------------------------------------------------------------- (W14)
print()
print(BAR)
print('(W14) two instances quoted in the paper', flush=True)
dillen4 = sp.expand((x1 + x2**2) * x3 + (x2 + (x1 + x2**2)**2) * x4)
assert dethess(dillen4, X4) == 1
bb = x1 + x2**2
ee = x2 + (x1 + x2**2)**2
Jd = sp.Matrix([[sp.diff(bb, x1), sp.diff(bb, x2)], [sp.diff(ee, x1), sp.diff(ee, x2)]])
assert sp.expand(Jd.det()) == 1
ok('W14', "de Bondt's dillen4 potential is the doubling of the planar Keller map "
          "(x1 + x2^2, x2 + (x1 + x2^2)^2); det Hess == 1")
hwit = x1 * x2 + x1**2 * x3
assert dethess(hwit, X3) == 0 and Bform(hwit, X3) == -x1**4
ok('W14', 'de Bondt-van den Essen witness x1 x2 + x1^2 x3: det Hess == 0 but B == -x1^4 != 0')
print('       => Corollary E is strictly stronger than "det Hess == 0".', flush=True)

# ---------------------------------------------------------------- (W15)
print()
print(BAR)
print('(W15) SHARPNESS of Theorem G: polynomiality is essential  [PROOF + witnesses]',
      flush=True)
# Euler for a function homogeneous of degree 0: Hess(e) x = -grad e and x^T Hess x = 0,
# hence B(e) = (Hess x)^T adj(H) (Hess x) = det(H) * x^T H x = 0 identically.
for n in (2, 3):
    vs = sp.symbols(f'y1:{n+1}')
    Hs = sym_matrix(n, f'S{n}_')
    xs = sp.Matrix(list(vs))
    gg = -Hs * xs                                   # = grad e for a degree-0 homogeneous e
    lhs15 = sp.expand((gg.T * Hs.adjugate() * gg)[0, 0])
    rhs15 = sp.expand(Hs.det(method='berkowitz') * (xs.T * Hs * xs)[0, 0])
    assert sp.expand(lhs15 - rhs15) == 0
    ok('W15', f'n={n}: with grad e = -(Hess e)x one has B(e) == det(Hess e) * x^T Hess e x')
print('       and Euler in degree 0 gives x^T Hess e x = -x . grad e = 0, so EVERY')
print('       function homogeneous of degree 0 satisfies B(e) == 0.', flush=True)
for expr, vs in ((x1 / x2, X2), (x1 * x2 / x3**2, X3)):
    H = hess(expr, vs)
    g = grad(expr, vs)
    Bv = sp.simplify((g.T * H.adjugate() * g)[0, 0])
    Dv = sp.simplify(H.det())
    assert Bv == 0 and sp.simplify(Dv) != 0
    ok('W15', f'witness e = {expr}: B(e) == 0 but det Hess e = {Dv} != 0')
print('       These are RATIONAL, homogeneous of degree 0.  A polynomial homogeneous of')
print('       degree 0 is constant, so there is no contradiction with Theorem G: the')
print('       theorem genuinely uses polynomiality, exactly at Step 4 of its proof.',
      flush=True)

print()
print(BAR)
print(f'ALL paperG_writeup CHECKS PASSED.  ({time.time()-T0:.1f}s)')
