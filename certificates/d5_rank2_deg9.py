# d5_rank2_deg9.py
#
# Degree-5 leading forms of essential rank <= 2: the graded pieces [det Hess f]_10,
# [det Hess f]_9, [det Hess f]_8 as FORMAL identities, and their consequences.
#
# Setting.  f = f5 + f4 + f3 + f2 + (affine) in C[x1..x4], f5 in C[x1,x2] (binary
# quintic; rank 2 means D := det2 Hess2 f5 =/= 0, rank 1 means f5 = x1^5 after a
# linear change).  Then Hess f (t x) = t^3 H5 + t^2 H4 + t H3 + H2 with
# H5 = Hess f5 supported in the (x1,x2) block, H4 = Hess f4, H3 = Hess f3,
# H2 = Hess f2, and [det Hess f]_d = coefficient of t^d in det(t^3 H5 + ...).
#
# PART 1 (formal, fully symbolic in the 33 matrix entries -- a PROOF, not an
# instance check):  with N4 = (x3,x4)-block of H4, R3 = (x3,x4)-block of H3,
# R2 = (x3,x4)-block of H2, D = det H5, cof_ij = cofactor of the 4x4 matrix,
#   [det]_10 = D * det N4
#   [det]_9  = D * tr(adj(N4) R3) + sum_{i,j<=2} (H5)_ij cof_ij(H4)
#   [det]_8  = D * (det R3 + tr(adj(N4) R2))
#              + sum_{i,j<=2} (H5)_ij d/de cof_ij(H4 + e H3)|_{e=0} + det H4
# and [det]_d = 0 for d > 10.
#
# PART 2 (formal): if N4 = l * u u^T (rank <= 1), with w := (u4, -u3) (so
# adj N4 = l w w^T) and g_i := u4 (H4)_i3 - u3 (H4)_i4 (i = 1,2):
#   [det]_9 = l * ( D * w^T R3 w  -  g^T adj(H5) g ).
# If N4 = 0:  [det]_10 = [det]_9 = 0 and, with Q := [(H4)_ij]_{i<=2, j>=3},
#   [det]_8 = D det R3 - tr( adj(H5) Q adj(R3) Q^T ) + (det Q)^2.
# If f5 = x1^5 (rank 1, H5 = 20 x1^3 e1 e1^T):
#   [det]_10 = 0,  [det]_9 = 20 x1^3 det3 Hess_{x2,x3,x4} f4,
#   [det]_8  = 20 x1^3 tr( adj(Hess_{234} f4) Hess_{234} f3 ) + det Hess f4.
#
# PART 3 (rank 2, the two rank-1 shapes of N4 = Hess_{x3,x4} f4, entries quadratic).
# [det]_10 = 0 and D =/= 0 force det N4 = 0, so (UFD) N4 = l u u^T with either
#   (a) u constant, l quadratic:  WLOG u = e3, f4 = A(x1,x2,x3) + x4 B(x1,x2),
#       l = A_x3x3;  then  [det]_9 = l * ( D f3_x4x4 - grad(B)^T adj(H5) grad(B) ).
#       Consequence (l =/= 0): f3_x4x4 =: lambda lies in C[x1,x2] and
#       D * lambda = grad(B)^T adj(H5) grad(B).  NOT a contradiction: witnesses
#       below have lambda = 0 (f5 = x1^4 x2, B = x1^3) AND lambda =/= 0
#       (f5 = x1^5 + x2^5, B = x1^3 + x2^3, lambda = (9/20)(x1 + x2)).  So [det]_9
#       does NOT force e4 to be a pivot in degree 5 (unlike degree 4, Theorem B).
#   (b) u = (L1, L2) linear, l constant =: alpha.  Integrability of N4 = Hess f4
#       forces L1, L2 in C[x1,x2] whenever L1, L2 are independent (the dependent
#       case is (a) with l = alpha L1^2).  Then f4 = (alpha/2) m^2 + x3 C + x4 D'
#       + E with m = L1 x3 + L2 x4, C, D' binary cubics, E a binary quartic, and
#       the part of [det]_9 of degree 2 in (x3, x4) equals
#             -20 alpha^3 delta^2 f5 m^2,    delta = det[L1; L2] =/= 0,
#       which is NONZERO.  Hence sub-case (b) is EMPTY:  no such f has
#       det Hess f in C^x.   [proved here, fully symbolic]
#   (c) N4 = 0: [det]_10 = [det]_9 = 0 automatically; the (x3,x4)-bidegree-2 part
#       of [det]_8 is D * det(R3'') with R3'' = Hess_{x3,x4} of the pure-(x3,x4)
#       part f3'' of f3; so det Hess2(f3'') == 0, i.e. f3'' is the cube of a linear
#       form (Hesse, binary cubics).  The rest of the branch is OPEN.
#
# Every identity below is checked with all coefficients symbolic (exact, char 0);
# the two instance checks only guard the transcription of the specializations.
# Fail-closed: any mismatch asserts.
#
#   py -u d5_rank2_deg9.py
import itertools, random, sys
import sympy as sp

x1, x2, x3, x4 = X = sp.symbols('x1 x2 x3 x4')
t, eps = sp.symbols('t eps')


def symm(name, n):
    M = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = sp.Symbol(f'{name}{i+1}{j+1}')
    return M


def check(label, expr):
    ok = sp.expand(expr) == 0
    print(f'  [{"ok" if ok else "FAIL"}] {label}')
    assert ok, label


def hess(g, vs=X):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def graded(poly, d):
    P = sp.Poly(sp.expand(poly), *X)
    return sp.expand(sum(c * sp.prod(v**e for v, e in zip(X, m)) for m, c in P.terms() if sum(m) == d))


def gen_form(vs, d, name):
    mons = [sp.prod(c) for c in itertools.combinations_with_replacement(vs, d)]
    cs = sp.symbols(f'{name}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons)), cs


# ---------------------------------------------------------------- PART 1
print('PART 1: formal graded pieces of det(t^3 H5 + t^2 H4 + t H3 + H2)')
H5 = symm('p', 2)
H5e = sp.zeros(4, 4); H5e[:2, :2] = H5
H4, H3, H2 = symm('q', 4), symm('r', 4), symm('s', 4)
M = t**3 * H5e + t**2 * H4 + t * H3 + H2
detM = sp.Poly(sp.expand(M.det(method='berkowitz')), t)
c = {d: detM.coeff_monomial(t**d) for d in range(0, 13)}
assert detM.degree() <= 10, 'degree > 10 appears'
print('  [ok] [det]_d = 0 for d > 10')
N4, R3, R2 = H4[2:, 2:], H3[2:, 2:], H2[2:, 2:]
D = H5.det()
check('[det]_10 = D det N4', c[10] - D * N4.det())
S9 = sum(H5[i, j] * H4.cofactor(i, j) for i in range(2) for j in range(2))
check('[det]_9 = D tr(adj N4 R3) + sum H5_ij cof_ij(H4)', c[9] - (D * (N4.adjugate() * R3).trace() + S9))
H4e = H4 + eps * H3
S8 = sum(H5[i, j] * sp.diff(H4e.cofactor(i, j), eps).subs(eps, 0) for i in range(2) for j in range(2))
check('[det]_8 = D (det R3 + tr(adj N4 R2)) + sum H5_ij d_e cof_ij(H4+eH3) + det H4',
      c[8] - (D * (R3.det() + (N4.adjugate() * R2).trace()) + S8 + H4.det()))

# ---------------------------------------------------------------- PART 2
print('PART 2: rank-1 / zero / rank-1-f5 specializations (formal)')
l, u3, u4 = sp.symbols('l u3 u4')
q33, q34, q44 = H4[2, 2], H4[2, 3], H4[3, 3]
sub1 = {q33: l * u3**2, q34: l * u3 * u4, q44: l * u4**2}
w = sp.Matrix([u4, -u3])
g = sp.Matrix([u4 * H4[i, 2] - u3 * H4[i, 3] for i in range(2)])
check('N4 = l u u^T:  [det]_9 = l (D w^T R3 w - g^T adj(H5) g)',
      c[9].subs(sub1) - l * (D * (w.T * R3 * w)[0] - (g.T * H5.adjugate() * g)[0]))
sub0 = {q33: 0, q34: 0, q44: 0}
check('N4 = 0:  [det]_10 = 0', c[10].subs(sub0))
check('N4 = 0:  [det]_9 = 0', c[9].subs(sub0))
Q = H4[:2, 2:]
check('N4 = 0:  [det]_8 = D det R3 - tr(adj(H5) Q adj(R3) Q^T) + (det Q)^2',
      c[8].subs(sub0) - (D * R3.det() - (H5.adjugate() * Q * R3.adjugate() * Q.T).trace() + Q.det()**2))
subr1 = {H5[0, 1]: 0, H5[1, 1]: 0}
check('f5 = x1^5 (H5 = p11 e1e1^T):  [det]_10 = 0', c[10].subs(subr1))
check('f5 = x1^5:  [det]_9 = p11 det3(H4[2:4,2:4])', c[9].subs(subr1) - H5[0, 0] * H4[1:, 1:].det())
check('f5 = x1^5:  [det]_8 = p11 tr(adj(H4_234) H3_234) + det H4',
      c[8].subs(subr1) - (H5[0, 0] * (H4[1:, 1:].adjugate() * H3[1:, 1:]).trace() + H4.det()))

# ---------------------------------------------------------------- PART 3
print('PART 3: rank-2 branch, sub-cases of N4 = Hess_{x3,x4} f4 = l u u^T')
f5, cf5 = gen_form((x1, x2), 5, 'f')
Hf5 = hess(f5, (x1, x2)); Df5 = Hf5.det()
f3, cf3 = gen_form(X, 3, 'k')
Rf3 = hess(f3, (x3, x4))

# --- (b) integrability: u = (L1, L2) linear in x1..x4, N4 = alpha u u^T = Hess_{34} f4
a1, a2, a3, a4, b1, b2, b3, b4 = ab = sp.symbols('a1:5 b1:5')
L1 = a1 * x1 + a2 * x2 + a3 * x3 + a4 * x4
L2 = b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4
N = sp.Matrix([[L1**2, L1 * L2], [L1 * L2, L2**2]])           # alpha = 1 WLOG for the ideal
integ = [sp.diff(N[0, 0], x4) - sp.diff(N[0, 1], x3), sp.diff(N[1, 1], x3) - sp.diff(N[0, 1], x4)]
I_int = [co for e in integ for co in sp.Poly(sp.expand(e), *X).coeffs()]
minors = [a1 * b2 - a2 * b1, a1 * b3 - a3 * b1, a1 * b4 - a4 * b1,
          a2 * b3 - a3 * b2, a2 * b4 - a4 * b2, a3 * b4 - a4 * b3]
zz = sp.Symbol('zz')
nchk = 0
for v in (a3, a4, b3, b4):
    for dl in minors:
        G = sp.groebner(I_int + [1 - zz * v * dl], *ab, zz, order='grevlex')
        assert list(G.exprs) == [1], f'{v}*({dl}) not in rad(I_int)'
        nchk += 1
print(f'  [ok] (b) integrability: v*minor in rad(I_int) for v in (a3,a4,b3,b4), all 6 minors ({nchk} Rabinowitsch certs)')
print('       => L1, L2 independent  =>  L1, L2 in C[x1,x2]')

# --- (b) emptiness: L1, L2 in C[x1,x2] independent, f4 = (alpha/2) m^2 + x3 C + x4 D' + E
alpha = sp.Symbol('alpha')
L1b = a1 * x1 + a2 * x2; L2b = b1 * x1 + b2 * x2
delta = a1 * b2 - a2 * b1
m = L1b * x3 + L2b * x4
Cc, _ = gen_form((x1, x2), 3, 'c'); Dd, _ = gen_form((x1, x2), 3, 'd'); Ee, _ = gen_form((x1, x2), 4, 'e')
f4b = alpha / 2 * m**2 + x3 * Cc + x4 * Dd + Ee
H4b = hess(f4b)
N4b = H4b[2:, 2:]
check('(b) N4 = alpha (L1,L2)(L1,L2)^T', sp.expand(N4b - alpha * sp.Matrix([L1b, L2b]) * sp.Matrix([L1b, L2b]).T).norm(1))
wb = sp.Matrix([L2b, -L1b])
gb = sp.Matrix([L2b * H4b[i, 2] - L1b * H4b[i, 3] for i in range(2)])
det9b = sp.expand(alpha * (Df5 * (wb.T * Rf3 * wb)[0] - (gb.T * Hf5.adjugate() * gb)[0]))
# part of (x3,x4)-degree exactly 2
P9 = sp.Poly(det9b, x3, x4)
top2 = sp.expand(sum(co * x3**e3 * x4**e4 for (e3, e4), co in P9.terms() if e3 + e4 == 2))
check('(b) (x3,x4)-degree-2 part of [det]_9 = -20 alpha^3 delta^2 f5 m^2', top2 + 20 * alpha**3 * delta**2 * f5 * m**2)
print('       f5 =/= 0, alpha =/= 0, delta =/= 0, m =/= 0  =>  [det]_9 =/= 0  =>  sub-case (b) EMPTY')
# instance guard: formula vs direct determinant (exact integers)
rng = random.Random(7)
inst = {s: rng.randint(-3, 3) for s in list(cf5) + list(cf3) + list(ab[:2]) + list(ab[4:6]) + [alpha]}
inst.update({s: rng.randint(-3, 3) for s in f4b.free_symbols if s not in X and s not in inst})
inst[alpha] = 2; inst[a1] = 1; inst[a2] = -2; inst[b1] = 3; inst[b2] = 1
fb_inst = (f5 + f4b + f3).subs(inst)
check('(b) instance: formula [det]_9 == graded piece of det Hess f',
      graded(hess(fb_inst).det(), 9) - det9b.subs(inst))

# --- (a) u = e3: f4 = A(x1,x2,x3) + x4 B(x1,x2), l = A_x3x3
A, cA = gen_form((x1, x2, x3), 4, 'A'); B, cB = gen_form((x1, x2), 3, 'B')
f4a = A + x4 * B
H4a = hess(f4a)
ell = sp.diff(A, x3, 2)
check('(a) N4 = l e3 e3^T', sp.expand(H4a[2:, 2:] - sp.Matrix([[ell, 0], [0, 0]])).norm(1))
gradB = sp.Matrix([sp.diff(B, x1), sp.diff(B, x2)])
det9a = sp.expand(ell * (Df5 * sp.diff(f3, x4, 2) - (gradB.T * Hf5.adjugate() * gradB)[0]))
wa = sp.Matrix([0, -1]); ga = sp.Matrix([-H4a[i, 3] for i in range(2)])
check('(a) rank-1 formula gives [det]_9 = l (D f3_x4x4 - grad(B)^T adj(H5) grad(B))',
      ell * (Df5 * (wa.T * Rf3 * wa)[0] - (ga.T * Hf5.adjugate() * ga)[0]) - det9a)
inst_a = {s: rng.randint(-3, 3) for s in list(cf5) + list(cf3) + list(cA) + list(cB)}
fa_inst = (f5 + f4a + f3).subs(inst_a)
check('(a) instance: formula [det]_9 == graded piece of det Hess f',
      graded(hess(fa_inst).det(), 9) - det9a.subs(inst_a))
# witnesses: [det]_9 = 0 is satisfiable with lambda = 0 and with lambda =/= 0
for (f5w, Bw, lamw) in [(x1**4 * x2, x1**3, 0),
                        (x1**5 + x2**5, x1**3 + x2**3, sp.Rational(9, 20) * (x1 + x2))]:
    Hw = hess(f5w, (x1, x2)); gw = sp.Matrix([sp.diff(Bw, x1), sp.diff(Bw, x2)])
    check(f'(a) witness f5={f5w}, B={Bw}: D*lambda == grad(B)^T adj(H5) grad(B) with lambda={lamw}',
          Hw.det() * lamw - (gw.T * Hw.adjugate() * gw)[0])
print('       => [det]_9 = 0 admits lambda = f3_x4x4 =/= 0: e4 NOT forced to be a pivot at this stage')

# --- (c) N4 = 0: f4 = x3 C + x4 D' + E; bidegree-(.,2) part of [det]_8 is D det(R3'')
f4c = x3 * Cc + x4 * Dd + Ee
H4c = hess(f4c)
assert H4c[2:, 2:] == sp.zeros(2, 2)
Qc = H4c[:2, 2:]
f2, cf2 = gen_form(X, 2, 's')
det8c = sp.expand(Df5 * Rf3.det() - (Hf5.adjugate() * Qc * Rf3.adjugate() * Qc.T).trace() + Qc.det()**2)
f3pp = sum(co * x3**e3 * x4**e4 for (e1, e2, e3, e4), co in sp.Poly(f3, *X).terms() if e1 == e2 == 0)
R3pp = hess(f3pp, (x3, x4))
P8 = sp.Poly(det8c, x3, x4)
top2c = sp.expand(sum(co * x3**e3 * x4**e4 for (e3, e4), co in P8.terms() if e3 + e4 == 2))
check('(c) (x3,x4)-degree-2 part of [det]_8 = D det Hess2(f3\'\')', top2c - Df5 * R3pp.det())
inst_c = {s: rng.randint(-3, 3) for s in list(cf5) + list(cf3) + list(cf2) + list(f4c.free_symbols - set(X))}
fc_inst = (f5 + f4c + f3 + f2).subs(inst_c)
check('(c) instance: formula [det]_8 == graded piece of det Hess f',
      graded(hess(fc_inst).det(), 8) - det8c.subs(inst_c))
print('       => D =/= 0 forces det Hess2(f3\'\') == 0: the pure-(x3,x4) part of f3 is a cube of a linear form')

# --- context: "quadratic in x4 with a LINEAR x4^2-coefficient" is not empty in general
check('context: det Hess(x1 x4^2/2 + x2 x4 + x1 x3) = 1 (lambda = x1 =/= 0, pivot e3 instead)',
      hess(x1 * x4**2 / 2 + x2 * x4 + x1 * x3).det() - 1)

print()
print('ALL CHECKS PASSED.')
print('PROVED HERE: formal [det]_10/9/8 identities; rank-2 sub-case (b) (u linear) EMPTY;')
print('  sub-case (a) [det]_9 constraint (lambda in C[x1,x2], D*lambda = grad(B)^T adj(H5) grad(B));')
print('  N4 = 0: f3\'\' is a cube; r = 1: [det]_9 = 20 x1^3 det3 Hess_{234} f4.')
print('OPEN: rank-2 sub-case (a) beyond [det]_9 (lambda =/= 0 survives), N4 = 0 beyond the')
print('  bidegree-2 part of [det]_8, and the r = 1 residual det3 Hess_{234} f4 == 0.')
sys.exit(0)
