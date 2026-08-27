# Certificate for the homogenization identity quoted in main_result.tex.
#
# Let e in C[x1,x2,x3] have degree d, and let E(x0,x1,x2,x3) = x0^d e(x/x0)
# be its degree-d homogenization. Then, with g = grad e and H = Hess_3 e,
#
#   det Hess_4(E) |_{x0 = 1}  =  d(d-1) * e * det H  -  (d-1)^2 * g^T adj(H) g.
#
# Derivation (row/column reduction of the bordered form of Hess E, using
# E_ij = e_ij, E_0i = (d-1)e_i - (H x)_i, E_00 = d(d-1)e - 2(d-1)x.g + x^T H x
# at x0 = 1): adding x^T times the lower block rows to row 0 and likewise for
# the columns turns Hess E into [[d(d-1)e, (d-1)g^T],[(d-1)g, H]]; Cauchy's
# bordered expansion then gives the stated formula.
#
# Consequence used in the paper: under the affine-pivot constraint (c2)
# (g^T adj(H) g == 0), det Hess_4(E) is divisible by E -- the Hessian
# hypersurface of E contains E itself. Combined with Theorem D' (the leading
# form's Hessian determinant vanishes) this is the proposed route to
# Conjecture E.

import itertools
import sympy as sp

x0, x1, x2, x3 = sp.symbols('x0 x1 x2 x3')
X = [x1, x2, x3]


def gen_poly(maxdeg, tag):
    mons = [sp.Integer(1)]
    for d in range(1, maxdeg + 1):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(X, d)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c*m for c, m in zip(cs, mons)), cs


for d in (2, 3):
    e, _ = gen_poly(d, f'a{d}_')
    # homogenize to degree d in (x0, x1, x2, x3)
    E = sp.expand(x0**d * e.subs({v: v/x0 for v in X}))
    E = sp.simplify(sp.expand(E))
    assert sp.Poly(E, x0, *X).is_homogeneous, 'homogenization failed'
    HE = sp.Matrix([[sp.diff(E, u, v) for v in (x0, *X)] for u in (x0, *X)])
    lhs = sp.expand(HE.det(method='berkowitz').subs(x0, 1))

    H = sp.Matrix([[sp.diff(e, u, v) for v in X] for u in X])
    g = sp.Matrix([sp.diff(e, u) for u in X])
    rhs = sp.expand(d*(d-1)*e*H.det(method='berkowitz')
                    - (d-1)**2*(g.T*H.adjugate()*g)[0, 0])
    assert sp.expand(lhs - rhs) == 0, f'identity failed at d = {d}'
    print(f'PASS d = {d}: det Hess_4(E)|_(x0=1) = {d}*{d-1}*e*det H - {(d-1)**2}*g^T adj(H) g')
    print(f'         (fully generic e of degree {d} in three variables => proof for that degree)')

print()
print('HOMOGENIZATION IDENTITY VERIFIED. Under (c2) the first term survives')
print('alone, so det Hess_4(E) is divisible by E.')
