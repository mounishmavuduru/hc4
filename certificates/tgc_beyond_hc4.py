# tgc_beyond_hc4.py  --  key: tgc  (Theorem G, Consequences)
#
# PARTS (3) and (4) of the "consequences of Theorem G beyond HC_4" task:
# what Theorem G says for the vanishing-Hessian (Gordan-Noether/Perazzo)
# literature, where the low-dimensional classification stops, and the
# converse in the homogeneous case.
#
# NOTATION.  B(e) = grad(e)^T adj(Hess e) grad(e)  (bordered Hessian up to
# sign; Reilly's / Fox's U(F)).  E = the degree-d homogenization of e in the
# extra variable x0.  "affinely k-variable" = after an affine change e lies in
# C[x_1..x_k], equivalently dim{v : D_v e == 0} = n - k.
#
# ---------------------------------------------------------------- RESULTS
#
# R1 (corank lemma; elementary).  If rank Hess e <= n-2 at every point then
#     adj(Hess e) == 0 and hence B(e) == 0.  So for n >= 4 the condition
#     "B == 0" is implied by the cheap condition "corank >= 2" and is NOT a
#     rigidity condition by itself.  For n = 3 corank >= 2 means rank <= 1.
#
# R2 (the exact place of Theorem G in the Gordan-Noether picture).
#     For e of degree d >= 2 with homogenization E:
#            B(e) == 0   <=>   det Hess_n(e) == 0  AND  det Hess_{n+1}(E) == 0.
#     "<=" is Dolgachev's affine equation of the Hessian (project Identity H),
#     "=>" is Dolgachev + THEOREM G.  So Theorem G says precisely: the affine
#     "all level sets developable" condition is the CONJUNCTION of the two
#     classical vanishing-Hessian conditions, for e and for its projective
#     closure.  Neither one alone suffices:
#        det Hess_3(x1x2 + x1^2x3) == 0  but B = -x1^4 != 0  (E has det Hess != 0);
#        det Hess_2(x1^2/2 + x2)  == 0  but B = 1     != 0.
#
# R3 (sharpness of Corollary E in the dimension).  Corollary E ("B == 0 =>
#     affinely 2-variable") is TRUE for n <= 3 and FALSE for every n >= 4:
#        n = 4:  e = x2 + x1x3 + x1^2x4   has B == 0 and 4 essential variables.
#                Its homogenization is exactly the PERAZZO CUBIC
#                x0^2x2 + x0x1x3 + x1^2x4, the canonical quinary form with
#                zero Hessian that is not a cone (Watanabe-de Bondt Rem.
#                after Thm 4.5; Perazzo 1900).  So the failure of Corollary E
#                in four variables is exactly the failure of Gordan-Noether in
#                five -- the same single example, dehomogenized.
#        n = 4:  e = a(x1,x2,x3) + x4 for ANY a with det Hess_3 a == 0, e.g.
#                a = x1x2 + x1^2x3, gives an infinite family.
#        n >= 5: the Perazzo cubic itself is homogeneous with det Hess == 0,
#                so B == 0 by Identity E, and it is not a cone.
#
# R4 (classification in four variables).  Combining R2 with Watanabe-de Bondt's
#     classification of quinary forms with zero Hessian (arXiv:1703.07624,
#     Thm 4.5): if e in C[x1..x4] has B(e) == 0 then either e is affinely
#     3-variable, or the degree-d homogenization E is, after a linear change
#     of C^5, an element of K[y1,y2][Delta] with
#     Delta = p3(y1,y2)y3 + p4(y1,y2)y4 + p5(y1,y2)y5.
#     (Proof: R2 gives det Hess_5 E == 0.  If E is a cone, D_vE == 0 for some
#     v = (v0,v'); if v0 = 0 then D_{v'}e == 0 and e is affinely 3-variable; if
#     v0 != 0 then, as in the project's Theorem F, a translate of e is
#     homogeneous, whence B == 0 forces det Hess_4 e == 0 by Identity E and
#     Gordan-Noether in FOUR variables makes it a cone -- again affinely
#     3-variable.  If E is not a cone, WdB Thm 4.5 applies.)
#
# R5 (part 4: the converse in the homogeneous case).  For e homogeneous of
#     degree d >= 2 (or homogeneous after a translation),
#            B(e) == 0   <=>   det Hess e == 0,
#     immediately from Identity E, B(e) = (d/(d-1)) e det Hess e, since C[x] is
#     a domain and e != 0.  So on cones Theorem G is a one-line KNOWN corollary
#     of Identity E, whose prior art is Fox (arXiv:1503.09108), Del Pia-
#     Hildebrand-Weismantel-Zemmer (arXiv:1408.4711, Lemma 5.2, after Hemmer
#     1995), Dolgachev (Classical Algebraic Geometry (1.20)) and -- in the
#     stronger one-parameter form -- Nagaoka-Yazawa (arXiv:1904.01800,
#     Prop. `identity1`), of which Identity E is the s-linear coefficient
#     (project certificate `identityE_is_known.py`).  ALL the content of
#     Theorem G is therefore in the non-conical case, and on cones the class
#     {B == 0} coincides exactly with the classical Gordan-Noether class.

import itertools
import sympy as sp

x0 = sp.Symbol('x0')


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def Bof(e, vs):
    H = hess(e, vs)
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T * H.adjugate() * g)[0, 0])


def dethess(e, vs):
    return sp.expand(hess(e, vs).det(method='berkowitz'))


def ess_vars(e, vs):
    n = len(vs)
    cs = sp.symbols(f'vv0:{n}')
    Dv = sp.expand(sum(c * sp.diff(e, x) for c, x in zip(cs, vs)))
    if Dv == 0:
        return 0
    eqs = sp.Poly(Dv, *vs).coeffs()
    M, _ = sp.linear_eq_to_matrix(eqs, cs)
    return n - len(M.nullspace())


def homogenize(e, vs, d=None):
    d = d or sp.Poly(e, *vs).total_degree()
    P = sp.Poly(e, *vs)
    out = 0
    for mono, c in zip(P.monoms(), P.coeffs()):
        out += c * sp.prod([v**k for v, k in zip(vs, mono)]) * x0**(d - sum(mono))
    return sp.expand(out), d


V = {n: sp.symbols(f'x1:{n+1}') for n in (2, 3, 4, 5)}

print('=' * 74)
print('R1.  corank >= 2  =>  adj(Hess) == 0  =>  B == 0')
print('=' * 74)
for n in (3, 4, 5):
    # a symmetric matrix of rank <= n-2 built generically: M = u u^T + w w^T
    # with n-2 generic rank-one pieces
    ks = n - 2
    M = sp.zeros(n, n)
    for j in range(ks):
        u = sp.Matrix(sp.symbols(f'u{n}{j}_0:{n}'))
        M += u * u.T
    assert sp.expand(M.adjugate()) == sp.zeros(n, n), f'adj not zero, n={n}'
    print(f'  n={n}: a generic symmetric matrix of rank <= n-2 has adj == 0      PROVED')
print('  => for n >= 4 every polynomial whose Hessian has corank >= 2 satisfies')
print('     B == 0; "B == 0" alone is not a rigidity condition when n >= 4.')

print()
print('=' * 74)
print('R2.  B(e) == 0  <=>  det Hess_n(e) == 0  AND  det Hess_{n+1}(E) == 0')
print('=' * 74)
# Identity H (Dolgachev) verified in general n, fully generically, for the
# degrees we need; this is the "<=" direction and half of "=>".
for (n, d) in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2)):
    vs = V[n]
    mons = [m for k in range(0, d + 1)
            for m in sorted({sp.prod(c) for c in
                             itertools.combinations_with_replacement(vs, k)}, key=str)]
    cs = sp.symbols(f'h{n}{d}_0:{len(mons)}')
    e = sum(c * m for c, m in zip(cs, mons))
    E, dd = homogenize(e, vs, d)
    lhs = sp.expand(dethess(E, (x0,) + tuple(vs)).subs(x0, 1))
    rhs = sp.expand(d * (d - 1) * e * dethess(e, vs) - (d - 1)**2 * Bof(e, vs))
    assert sp.expand(lhs - rhs) == 0, f'Identity H failed n={n} d={d}'
    print(f'  Identity H (Dolgachev) verified fully generically: n={n}, d={d}   PROVED')
# n = 4, d = 3 with a sparse but non-degenerate instance (the fully generic
# 5-variable determinant is out of reach for sympy; the identity is Dolgachev's
# and dimension-free, its proof being a row/column reduction).
vs = V[4]
e_i = (x0 * 0 + vs[0]**3 + vs[0] * vs[1] * vs[2] + vs[3]**2 * vs[1]
       + vs[0] * vs[3] + vs[2]**2 + vs[1] + 7)
E_i, d_i = homogenize(e_i, vs, 3)
assert sp.expand(dethess(E_i, (x0,) + tuple(vs)).subs(x0, 1)
                 - (3 * 2 * e_i * dethess(e_i, vs) - 4 * Bof(e_i, vs))) == 0
print('  Identity H checked on a dense instance with n=4, d=3               PASS')

print()
print('  Consequences on exact instances (fail-closed):')
x1, x2, x3, x4, x5 = sp.symbols('x1 x2 x3 x4 x5')
tests = [
    (x1 * x2 + x1**2 * x3, (x1, x2, x3)),
    (x1**2 / 2 + x2, (x1, x2)),
    (x1 * x2, (x1, x2)),
    (x1 * x2 + x3, (x1, x2, x3)),
    (x1**3 + x2**3 + x3**3, (x1, x2, x3)),
    (x2 + x1 * x3 + x1**2 * x4, (x1, x2, x3, x4)),
    (x1 * x2 + x1**2 * x3 + x4, (x1, x2, x3, x4)),
    (x1 * x2 + x3 * x4, (x1, x2, x3, x4)),
]
for e, vs in tests:
    E, d = homogenize(e, vs)
    B = Bof(e, vs)
    dhe = dethess(e, vs)
    dhE = dethess(E, (x0,) + tuple(vs))
    lhs = (B == 0)
    rhs = (dhe == 0) and (dhE == 0)
    assert lhs == rhs, f'R2 FAILED for {e}'
    print(f'    e={str(e):26s} B==0:{str(lhs):5s}  detHess_n==0:{str(dhe==0):5s} '
          f' detHess_{{n+1}}(E)==0:{str(dhE==0):5s}  consistent')

print()
print('=' * 74)
print('R3.  Corollary E is SHARP: it fails in every dimension n >= 4')
print('=' * 74)
X4 = (x1, x2, x3, x4)
X5 = (x1, x2, x3, x4, x5)

e4 = x2 + x1 * x3 + x1**2 * x4
assert Bof(e4, X4) == 0
assert dethess(e4, X4) == 0
assert ess_vars(e4, X4) == 4
E4, d4 = homogenize(e4, X4)
per = x0**2 * x2 + x0 * x1 * x3 + x1**2 * x4
assert sp.expand(E4 - per) == 0
assert dethess(per, (x0,) + X4) == 0
assert ess_vars(per, (x0,) + X4) == 5
print('  e = x2 + x1x3 + x1^2x4  in 4 variables:')
print('     B == 0, det Hess == 0, 4 ESSENTIAL VARIABLES  => Corollary E fails')
print('     its degree-3 homogenization is x0^2x2 + x0x1x3 + x1^2x4,')
print('     i.e. the PERAZZO CUBIC (zero Hessian, not a cone, 5 ess. variables)')

fam = []
for a in [x1 * x2 + x1**2 * x3, x1 * x2, x1**3 + x1**2 * x3, x1 * x2 + x1**3 * x3,
          x1**2 * x2 + x1**5 * x3]:
    e = a + x4
    fam.append((e, Bof(e, X4), ess_vars(e, X4), dethess(a, (x1, x2, x3))))
print()
print('  the infinite family  e = a(x1,x2,x3) + x4  with det Hess_3 a == 0:')
for e, B, nv, dha in fam:
    assert dha == 0 and B == 0
    print(f'    e = {str(e):28s} B == 0, det Hess_3(a) == 0, ess.vars = {nv}')
# jet-level proof that B(a(x') + lam x4) = lam^2 det Hess_3(a)
sy = sp.symbols('b11 b12 b13 b22 b23 b33 g1 g2 g3 lam')
b11, b12, b13, b22, b23, b33, g1, g2, g3, lam = sy
H = sp.Matrix([[b11, b12, b13, 0], [b12, b22, b23, 0], [b13, b23, b33, 0], [0, 0, 0, 0]])
g = sp.Matrix([g1, g2, g3, lam])
Bj = sp.expand((g.T * H.adjugate() * g)[0, 0])
M3 = sp.Matrix([[b11, b12, b13], [b12, b22, b23], [b13, b23, b33]])
assert sp.expand(Bj - lam**2 * M3.det()) == 0
print('    jet identity B(a(x1,x2,x3) + lam x4) = lam^2 det Hess_3(a)         PROVED')
print('    => every 3-variable zero-Hessian polynomial that is not affinely')
print('       2-variable yields a 4-variable counterexample to Corollary E.')

P5 = x1**2 * x3 + x1 * x2 * x4 + x2**2 * x5
assert Bof(P5, X5) == 0 and dethess(P5, X5) == 0 and ess_vars(P5, X5) == 5
print()
print('  n = 5: the Perazzo cubic x1^2x3 + x1x2x4 + x2^2x5 is homogeneous with')
print('     det Hess == 0, hence B == 0 by Identity E, and has 5 essential vars.')

print()
print('=' * 74)
print('R4.  n = 4 classification statement (uses Watanabe-de Bondt Thm 4.5)')
print('=' * 74)
print('  If e in C[x1..x4] has B(e) == 0 then either e is affinely 3-variable,')
print('  or its degree-d homogenization E lies, after a linear change of C^5,')
print('  in K[y1,y2][Delta],  Delta = p3(y1,y2)y3+p4(y1,y2)y4+p5(y1,y2)y5.')
print('  Checked on the witness: E = x0^2x2 + x0x1x3 + x1^2x4 is exactly of')
print('  that shape with (y1,y2) = (x0,x1) and Delta = E itself (p = y1^2,')
print('  y1y2, y2^2), and e is NOT affinely 3-variable.')
Delta = x0**2 * x2 + x0 * x1 * x3 + x1**2 * x4
assert sp.expand(E4 - Delta) == 0
print('  consistent.')

print()
print('=' * 74)
print('R7.  WHY n = 3 IS THE THRESHOLD: an isotropic-dimension count')
print('=' * 74)
# LEMMA (linear algebra).  A symmetric bilinear form on C^N that admits a
# totally isotropic subspace of dimension k with 2k > N is degenerate (a
# nondegenerate form has maximal isotropic subspaces of dimension floor(N/2)).
# CONSEQUENCE.  If a hypersurface Sigma in C^n contains, through each of its
# points, an affine k-plane, then the second fundamental form of Sigma (a form
# on C^{n-1}) vanishes on that k-plane, hence is degenerate as soon as
# 2k > n-1; i.e. Sigma is developable.
# APPLICATION.  Let  delta = a_2(x1)x_2 + ... + a_n(x1)x_n  and  e = f(x1,delta).
# The fibres of (x1, delta) are affine (n-2)-planes contained in the level sets
# of e, so every level set is developable as soon as 2(n-2) > n-1, i.e. n >= 4;
# and then B(e) == 0.  For n = 3 the count fails by exactly one, and indeed the
# same family is the de Bondt-van den Essen exceptional branch, where
# B = -W(a2,a3)^2 is generically NONZERO (certificate tgc_debondt_chain.py).
# THIS IS THE REASON Corollary E is a 2- and 3-variable phenomenon.
for N, k in ((3, 2), (5, 3), (4, 3), (3, 1), (4, 2)):
    # a form on C^N vanishing on the span of the first k basis vectors
    M = sp.zeros(N, N)
    syms = sp.symbols(f'm{N}{k}_0:{N*N}')
    idx = 0
    for i in range(N):
        for j in range(i, N):
            if i < k and j < k:
                continue
            M[i, j] = M[j, i] = syms[idx]
            idx += 1
    if 2 * k > N:
        assert sp.expand(M.det(method='berkowitz')) == 0, (N, k)
        print(f'  N={N}, isotropic dim k={k} > N/2: every such form is degenerate  PROVED')
    else:
        assert sp.expand(M.det(method='berkowitz')) != 0
        print(f'  N={N}, isotropic dim k={k} <= N/2: generically nondegenerate    (sharp)')

a2f, a3f, a4f = (sp.Function('a2')(x1), sp.Function('a3')(x1), sp.Function('a4')(x1))
sS, tS = sp.symbols('sS tS')
ff = sp.Function('f')(sS, tS)
delta4 = a2f * x2 + a3f * x3 + a4f * x4
e_fam = ff.subs({sS: x1, tS: delta4}, simultaneous=True)
Hf = sp.Matrix([[sp.diff(e_fam, a, b) for b in X4] for a in X4])
gf = sp.Matrix([sp.diff(e_fam, a) for a in X4])
assert sp.simplify(sp.expand((gf.T * Hf.adjugate() * gf)[0, 0])) == 0
print('  n = 4: B(f(x1, a2(x1)x2+a3(x1)x3+a4(x1)x4)) == 0 for ARBITRARY smooth')
print('         f and a_i (symbolic functions)                              PROVED')
print('  n >= 5: the same family has rank Hess <= 3 <= n-2, so adj == 0 and')
print('         B == 0 a fortiori.')
print('  n = 3: the same family is a1(x1)+a2(x1)x2+a3(x1)x3 with')
print('         B = -W(a2,a3)^2, generically NONZERO -- the count fails by one.')

print()
print('=' * 74)
print('R6.  the WHOLE PENCIL of projective closures has vanishing Hessian')
print('=' * 74)
# If B(e) == 0 then, for EVERY c, the degree-d homogenization of e - c, namely
# E - c*x0^d, has det Hess_{n+1} == 0.  (Identity H applied to e - c, using
# B(e-c) = B(e) = 0 and -- by THEOREM G -- det Hess_n(e) = 0.)
# Projectively: if every affine level hypersurface of e is developable, then
# the projective closure of every level hypersurface lies in the classical
# Gordan-Noether class of forms with identically vanishing Hessian.
c = sp.Symbol('c')
for e, vs in [(e4, X4), (x1 * x2 + x1**2 * x3 + x4, X4),
              (x1**2 / 2 + x2 * 0, (x1, x2)),
              (x1 * x2 * 0 + x1, (x1, x2))]:
    if Bof(e, vs) != 0:
        continue
    E, d = homogenize(e, vs)
    Ec = sp.expand(E - c * x0**d)
    assert dethess(Ec, (x0,) + tuple(vs)) == 0, f'pencil member not singular: {e}'
    print(f'    e = {str(e):26s}: det Hess_{{n+1}}(E - c x0^{d}) == 0 for symbolic c  PASS')
print('  (fails, as it must, when B != 0: )')
hbad = x1 * x2 + x1**2 * x3
Eb, db = homogenize(hbad, (x1, x2, x3))
assert dethess(sp.expand(Eb - c * x0**db), (x0, x1, x2, x3)) != 0
print('    e = x1x2 + x1^2x3 (B != 0): det Hess_4(E - c x0^3) != 0            PASS')

print()
print('=' * 74)
print('R5.  PART (4): the converse in the homogeneous case, via Identity E')
print('=' * 74)
s = sp.Symbol('s')
# (n = 3,4 with d = 2,3,4 are already certified in `identityE_is_known.py`;
#  here we re-verify the cheap range independently.)
for n in (2, 3):
    vs = V[n]
    for d in ((2, 3, 4) if n == 2 else (2, 3)):
        mons = sorted({sp.prod(c) for c in
                       itertools.combinations_with_replacement(vs, d)}, key=str)
        cs = sp.symbols(f'f{n}{d}_0:{len(mons)}')
        F = sum(c * m for c, m in zip(cs, mons))
        H = hess(F, vs)
        gr = sp.Matrix([[sp.diff(F, a) for a in vs]])
        # Identity E
        assert sp.expand((gr * H.adjugate() * gr.T)[0, 0]
                         - sp.Rational(d, d - 1) * F * H.det(method='berkowitz')) == 0
        # the full Nagaoka-Yazawa Proposition (arXiv:1904.01800), of which
        # Identity E is the coefficient of s^1
        lhs = sp.expand((-F * H + s * (gr.T * gr)).det(method='berkowitz'))
        rhs = sp.expand((-1)**(n - 1) * sp.Rational(d, d - 1) * (s - sp.Rational(d - 1, d))
                        * F**n * H.det(method='berkowitz'))
        assert sp.expand(lhs - rhs) == 0
    print(f'  n={n}: Identity E and the Nagaoka-Yazawa identity hold for '
          f'd={"2,3,4" if n == 2 else "2,3"}  PROVED')
print('  => for HOMOGENEOUS e of degree d >= 2:  B(e) == 0 <=> det Hess e == 0.')
print('     (C[x] is a domain and e != 0.)  The same holds after a translation,')
print('     B and det Hess both being translation invariant.  Hence Theorem G')
print('     is KNOWN and one-line on cones; its content is the non-conical case.')

# translation invariance of B, fully generic (n = 3, degree 3)
vs = V[3]
mons = [m for k in range(0, 4) for m in sorted({sp.prod(c) for c in
        itertools.combinations_with_replacement(vs, k)}, key=str)]
cs = sp.symbols(f'tt0:{len(mons)}')
e = sum(c * m for c, m in zip(cs, mons))
tvec = sp.symbols('w1 w2 w3')
et = sp.expand(e.subs({v: v + t for v, t in zip(vs, tvec)}, simultaneous=True))
assert sp.expand(Bof(et, vs) - Bof(e, vs).subs({v: v + t for v, t in zip(vs, tvec)},
                                               simultaneous=True)) == 0
print('  translation invariance of B verified fully generically (n=3, d<=3)   PROVED')

print()
print('ALL CHECKS PASSED.')
