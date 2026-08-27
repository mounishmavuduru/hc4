# gcv_ny_homogeneous_converse.py -- key: gcv (G Consequences, independent Verification)
#
# TASK PART (4): the converse-ish statement of Theorem G in the homogeneous
# case, via the Nagaoka-Yazawa identity, and the attribution audit.
#
# ATTRIBUTION FACTS (source-verified 2026-08-18).
#   * Nagaoka-Yazawa, "Strict log-concavity of the Kirchhoff polynomial and
#     its applications to the strong Lefschetz property", arXiv:1904.01800,
#     IS PUBLISHED:  J. Algebra 577 (2021) 175-202.  The identity is their
#     Proposition 2.3 (TeX label `identity1`, Section 2 "Homogeneous
#     polynomials", after Lemmas 2.1-2.2; their Euler lemma 2.2 is credited
#     to Anari-Gharan-Vinzant):
#        det( -F H_F + s (grad F)^T grad F )
#            = (-1)^{n-1} (r/(r-1)) (s - (r-1)/r) F^n det H_F,
#     for F homogeneous of degree r >= 2 (stated over R; both sides are
#     polynomial identities in the coefficients, so char 0 suffices).
#   * main_result.tex \bibitem{NY} lacked the title and journal data; the
#     statement "the identity of [NY] is an equivalent published form of
#     Lemma lem:euler (Identity E)" is CORRECT, and this script proves the
#     equivalence in BOTH directions (identityE_is_known.py proved one
#     direction, NY => Identity E, on heavier generic cases).
#
# WHAT IS PROVED HERE (fully generic symbolic checks = proofs for the shapes
# tested; the derivations are shape-independent formal algebra):
#   (1) Rank-one expansion  det(M + s u v^T) = det M + s * v^T adj(M) u
#       for ARBITRARY (nonsymmetric) 3x3 and 4x4 M -- the only linear-algebra
#       input needed to pass between NY and Identity E.
#   (2) Identity E (Euler contraction)  g^T adj(H) g = (d/(d-1)) F det H
#       for fully generic ternary forms of degree d = 2, 3.
#   (3) NY  =>  Identity E:  extract the s^1 coefficient of NY and cancel
#       F^{n-1} in the domain C[coefficients][x].
#   (4) Identity E  =>  NY:  by (1) with M = -F H, u = v = g:
#       det(-FH + s g^T g) = (-F)^n det H + s (-F)^{n-1} g^T adj(H) g,
#       then substitute (2).  Hence NY and Identity E are formally EQUIVALENT
#       given the rank-one expansion: "equivalent published form" is the
#       right description.
#   (5) THE HOMOGENEOUS CONVERSE.  For e homogeneous of degree d >= 2,
#           B(e) = (d/(d-1)) e det Hess e ,
#       so in the domain C[x] (e != 0):
#           B(e) == 0  <=>  det Hess e == 0 .
#       I.e. ON CONES Theorem G's implication is an equivalence -- and both
#       directions are one-line corollaries of the KNOWN identity.  The
#       equivalence persists for "homogeneous after a translation" (B and
#       det Hess are translation invariant -- re-proved generically here).
#       For general inhomogeneous e the converse FAILS:  e = x1^2/2 + x2 has
#       det Hess == 0 but B == 1.  Witnesses checked: Perazzo cubic
#       (both vanish), Fermat cubic (neither vanishes).
#
# Exact arithmetic only; fail-closed asserts.

import itertools
import sympy as sp

s = sp.Symbol('s')


def gen_form(vs, deg, tag):
    mons = sorted({sp.prod(c) for c in
                   itertools.combinations_with_replacement(vs, deg)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons))


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def Bof(e, vs):
    H = hess(e, vs)
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T * H.adjugate() * g)[0, 0])


def dethess(e, vs):
    return sp.expand(hess(e, vs).det(method='berkowitz'))


print('=' * 74)
print('(1) RANK-ONE EXPANSION  det(M + s u v^T) = det M + s v^T adj(M) u')
print('=' * 74)
for n in (3, 4):
    M = sp.Matrix(n, n, sp.symbols(f'm{n}_0:{n*n}'))
    u = sp.Matrix(sp.symbols(f'u{n}_0:{n}'))
    v = sp.Matrix(sp.symbols(f'v{n}_0:{n}'))
    lhs = sp.expand((M + s * u * v.T).det(method='berkowitz'))
    rhs = sp.expand(M.det(method='berkowitz') + s * (v.T * M.adjugate() * u)[0, 0])
    assert sp.expand(lhs - rhs) == 0, f'rank-one expansion failed, n={n}'
    print(f'  n={n}: proved for ARBITRARY M, u, v (fully symbolic).'
          '  In particular the')
    print('       determinant is AFFINE in s: no higher powers of s occur.')

print()
print('=' * 74)
print('(2)-(4) IDENTITY E  <=>  NAGAOKA-YAZAWA Prop. 2.3   (n = 3, d = 2, 3)')
print('=' * 74)
n = 3
vs = sp.symbols('x1:4')
for d in (2, 3):
    F = gen_form(vs, d, f'a{d}_')
    H = hess(F, vs)
    g = sp.Matrix([sp.diff(F, a) for a in vs])          # column
    detH = H.det(method='berkowitz')

    # (2) Identity E, fully generic
    IdE = sp.expand((g.T * H.adjugate() * g)[0, 0]
                    - sp.Rational(d, d - 1) * F * detH)
    assert IdE == 0, f'Identity E failed, d={d}'
    print(f'  d={d}:  Identity E  g^T adj(H) g = (d/(d-1)) F det H'
          '            PROVED (generic)')

    # NY identity, fully generic
    NY_l = sp.expand((-F * H + s * (g * g.T)).det(method='berkowitz'))
    NY_r = sp.expand((-1)**(n - 1) * sp.Rational(d, d - 1)
                     * (s - sp.Rational(d - 1, d)) * F**n * detH)
    assert sp.expand(NY_l - NY_r) == 0, f'NY identity failed, d={d}'
    print(f'        NY Prop. 2.3 identity                        '
          '            PROVED (generic)')

    # (3) NY => Identity E: s^1 coefficient, then cancel F^{n-1}
    c1 = sp.expand(sp.diff(NY_l, s))                    # affine in s by (1)
    # c1 = (-F)^{n-1} g^T adj(H) g must equal the s-coeff of NY_r:
    c1r = sp.expand(sp.diff(NY_r, s))
    assert sp.expand(c1 - c1r) == 0
    quot = sp.cancel(sp.together(c1 / (-F)**(n - 1)))
    assert sp.expand(quot - (g.T * H.adjugate() * g)[0, 0]) == 0
    assert sp.expand(quot - sp.Rational(d, d - 1) * F * detH) == 0
    print('        NY => Identity E  (s^1 coefficient / F^{n-1} in the domain)  OK')

    # (4) Identity E => NY: rank-one expansion with M = -F H, u = v = g
    lhs4 = sp.expand(((-F) * H + s * g * g.T).det(method='berkowitz'))
    rhs4 = sp.expand((-F)**n * detH
                     + s * (-F)**(n - 1) * sp.Rational(d, d - 1) * F * detH)
    assert sp.expand(lhs4 - rhs4) == 0
    assert sp.expand(rhs4 - NY_r) == 0
    print('        Identity E => NY  (rank-one expansion, adj(-FH) = (-F)^{n-1}adj H)  OK')
print('  => the two identities are FORMALLY EQUIVALENT: "equivalent published')
print('     form" in main_result.tex is accurate.')

print()
print('=' * 74)
print('(5) THE HOMOGENEOUS CONVERSE OF THEOREM G')
print('=' * 74)
x1, x2, x3, x4, x5 = sp.symbols('y1:6')
# generic binary forms d = 2..5: B = (d/(d-1)) e det Hess, so B==0 <=> det==0
for d in (2, 3, 4, 5):
    vs2 = (x1, x2)
    F = gen_form(vs2, d, f'b{d}_')
    assert sp.expand(Bof(F, vs2)
                     - sp.Rational(d, d - 1) * F * dethess(F, vs2)) == 0
print('  generic binary forms d = 2..5:  B = (d/(d-1)) e det Hess       PROVED')
print('  => for HOMOGENEOUS e (deg >= 2, any n where Identity E is certified):')
print('       B(e) == 0  <=>  det Hess e == 0        (C[x] is a domain, e != 0)')
print('     so on cones the one-directional Theorem G becomes an EQUIVALENCE,')
print('     and BOTH directions are corollaries of the KNOWN NY identity.')

# witnesses
X5 = (x1, x2, x3, x4, x5)
per = x1**2 * x3 + x1 * x2 * x4 + x2**2 * x5
assert dethess(per, X5) == 0 and Bof(per, X5) == 0
print('  Perazzo cubic (n=5):        det Hess == 0  and  B == 0          PASS')
X3v = (x1, x2, x3)
fer = x1**3 + x2**3 + x3**3
assert dethess(fer, X3v) != 0 and Bof(fer, X3v) != 0
assert sp.expand(Bof(fer, X3v) - sp.Rational(3, 2) * fer * dethess(fer, X3v)) == 0
print('  Fermat cubic (n=3):         det Hess != 0  and  B != 0          PASS')

# translation invariance (so "homogeneous after a translation" inherits it)
vs2 = (x1, x2)
e_g = gen_form(vs2, 3, 'tg3_') + gen_form(vs2, 2, 'tg2_') \
    + gen_form(vs2, 1, 'tg1_') + sp.Symbol('tg0')
t1, t2 = sp.symbols('t1 t2')
sub = {x1: x1 + t1, x2: x2 + t2}
assert sp.expand(Bof(sp.expand(e_g.subs(sub, simultaneous=True)), vs2)
                 - Bof(e_g, vs2).subs(sub, simultaneous=True)) == 0
assert sp.expand(dethess(sp.expand(e_g.subs(sub, simultaneous=True)), vs2)
                 - dethess(e_g, vs2).subs(sub, simultaneous=True)) == 0
print('  translation invariance of B and det Hess (generic n=2, deg<=3)  PROVED')

# the converse FAILS off cones
epar = x1**2 / 2 + x2
assert dethess(epar, vs2) == 0 and Bof(epar, vs2) == 1
print('  e = x1^2/2 + x2:  det Hess == 0 but B == 1: converse FAILS      PASS')
print('     for inhomogeneous e -- Theorem G is genuinely one-directional.')

print()
print('CONCLUSION (part 4).')
print('  * On homogeneous polynomials (and translates of such):')
print('        B(e) == 0  <=>  det Hess e == 0,')
print('    both directions being corollaries of Nagaoka-Yazawa Prop. 2.3')
print('    [J. Algebra 577 (2021) 175-202; arXiv:1904.01800].')
print('  * Identity E and the NY identity are formally equivalent (rank-one')
print('    expansion), so main_result.tex\'s description is accurate; its')
print('    \\bibitem{NY} needed the title + journal data (now supplied).')
print('  * Off cones the converse fails (x1^2/2 + x2), so the content of')
print('    Theorem G is exactly the non-conical case, where no prior art')
print('    was located.')
