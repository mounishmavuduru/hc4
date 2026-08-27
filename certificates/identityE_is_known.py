# ATTRIBUTION CERTIFICATE: the project's "Identity E" is NOT new.
#
# Project statement (lemma_ledger.md, Identity E; main_result.tex Lemma 4.3):
#     for e homogeneous of degree d >= 2 in n variables,
#         grad(e)^T adj(Hess e) grad(e) = (d/(d-1)) * e * det Hess e.
#
# PRIOR ART. Nagaoka & Yazawa, "Strict log-concavity of the Kirchhoff
# polynomial and its applications to the strong Lefschetz property",
# arXiv:1904.01800, Section 2 ("Homogeneous polynomials"), Proposition
# (label `identity1`), equation (label `The Hessian of a homogeneous
# polynomial`), verbatim from the arXiv TeX source:
#
#     det( -F H_F + s (grad F)^T grad F )
#         = (-1)^{n-1} * (r/(r-1)) * (s - (r-1)/r) * F^n * det H_F
#
# for F homogeneous of degree r >= 2 in n variables. Their proof uses the
# same two ingredients as ours: Euler's identities
# (r(r-1)F = x^T H_F x and (r-1)(grad F)^T = H_F x -- their Lemma `hom2`,
# itself credited to Anari-Gharan-Vinzant) plus rank-one determinant algebra.
#
# THIS SCRIPT PROVES that our Identity E is exactly the coefficient of s^1
# in their identity, hence a corollary of it. Therefore Identity E must be
# CITED, not claimed. (Their statement is over R; both sides are polynomial
# identities in the coefficients, so it holds over any field of char 0.)
#
# Derivation checked below symbolically for n = 3, 4 and degrees d = 2,3,4:
#   det(-e H + s g^T g) = det(-e H) + s * g adj(-e H) g            [rank-one]
#                       = (-e)^n det H + s (-e)^{n-1} g adj(H) g,
# and matching the s^1 coefficient with (-1)^{n-1}(r/(r-1)) F^n det H_F gives
#   g adj(H) g = (r/(r-1)) e det H.

import itertools
import sympy as sp

s = sp.Symbol('s')


def gen_form(vs, deg, tag):
    mons = sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, deg)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c*m for c, m in zip(cs, mons))


for n in (3, 4):
    vs = sp.symbols(f'y1:{n+1}')
    for d in (2, 3, 4):
        F = gen_form(vs, d, f'k{n}{d}_')
        H = sp.Matrix([[sp.diff(F, a, b) for b in vs] for a in vs])
        g = sp.Matrix([[sp.diff(F, a) for a in vs]])          # ROW vector, as in the paper
        # left-hand side of Nagaoka-Yazawa
        lhs = sp.expand((-F*H + s*(g.T*g)).det(method='berkowitz'))
        # right-hand side of Nagaoka-Yazawa
        rhs = sp.expand((-1)**(n-1) * sp.Rational(d, d-1) * (s - sp.Rational(d-1, d))
                        * F**n * H.det(method='berkowitz'))
        assert sp.expand(lhs - rhs) == 0, f'NY identity failed at n={n}, d={d}'

        # our Identity E = the s^1 coefficient of the same identity
        ours_l = sp.expand((g*H.adjugate()*g.T)[0, 0])
        ours_r = sp.expand(sp.Rational(d, d-1) * F * H.det(method='berkowitz'))
        assert sp.expand(ours_l - ours_r) == 0, f'Identity E failed at n={n}, d={d}'

        # and that ours really is the s-coefficient of theirs
        coeff_s = sp.expand(sp.diff(lhs, s))
        pred = sp.expand((-1)**(n-1) * F**(n-1) * ours_l)
        assert sp.expand(coeff_s - pred) == 0, f'extraction failed at n={n}, d={d}'
        print(f'n={n}, d={d}: NY identity OK; Identity E OK; and Identity E IS'
              f' the s-coefficient of NY (up to the factor (-1)^(n-1) F^(n-1)).')

print()
print('CONCLUSION: Identity E is a corollary of Nagaoka-Yazawa arXiv:1904.01800,')
print('Prop. `identity1`. It is KNOWN and must be cited, not claimed as new.')
print('Consequence for the project: Theorem D\' (cone leading form) rests on a')
print('KNOWN identity plus Gordan-Noether; only its APPLICATION to pivot')
print('coefficients of affine-pivot HC_4 potentials is ours.')
