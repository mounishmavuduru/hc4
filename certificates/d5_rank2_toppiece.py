# d5_rank2_toppiece.py
#
# Degree-5, rank-2 leading-form branch (first graded reduction).  In degree 5,
# Gordan-Noether puts f5 in C[x1,x2,x3]; the essential-rank-2 branch has
# f5 in C[x1,x2] a binary quintic with det2 Hess2 f5 =/= 0.  Write
# f = f5(x1,x2) + f4 + f3 + f2 with det Hess_4 f in C^x.
#
# CLAIM (this file):  the top graded piece of det Hess_4 f is
#     [det Hess_4 f]_10 = det2(Hess_{x1,x2} f5) * det2(Hess_{x3,x4} f4),
# so  det Hess f in C^x  forces  det2(Hess_{x3,x4} f4) == 0.
#
# PROOF (generalized Laplace expansion + a degree bound).  Block the 4x4 Hessian
# as [[P, Q],[Q^T, N]] with P = Hess_{x1,x2} f (top part Hess2 f5, entries deg 3),
# N = Hess_{x3,x4} f (top part Hess_{x3,x4} f4, entries deg 2), Q the (x1,x2)x
# (x3,x4) block (from f4,f3; entries deg <= 2).  Laplace-expand det along the
# first two rows: det = sum over column-pairs {i<j} of
#     (+/-) det M[{1,2},{i,j}] * det M[{3,4}, complement].
# Only the pair {i,j}={1,2} keeps both chosen columns inside P, giving
# det(P)*det(N), of total degree <= 6+4 = 10.  Every other pair uses at least one
# column from Q (deg <= 2) in place of a deg-3 P column, and its complementary
# minor then uses a column of Q^T (deg <= 2) in place of a deg-2 N column, so its
# total degree is <= 9 (checked term by term below).  Hence the degree-10 part
# comes ONLY from det(P)*det(N), and its top part is
# det2(Hess2 f5) * det2(Hess_{x3,x4} f4).  QED.
#
# The identity is polynomial (not linear) in the coefficients, so the machine
# check here is multi-instance (random exact-integer f5, f4; det2 Hess2 f5 =/= 0)
# rather than a single generic run (the fully symbolic 4x4 det is out of reach in
# this environment) -- the Laplace-degree argument above is the proof; the
# instances are corroboration.  Fail-closed: any mismatch asserts.
#
# SCOPE (honest): this is only the FIRST reduction of the rank-2 branch.  The
# completion -- det2(Hess_{x3,x4} f4) == 0 (QUADRATIC entries) => rank <= 1
# => Hess_{x3,x4} f4 = l * u u^T -- has u of degree 0 OR 1.  For u constant it
# yields a pivot exactly as in degree 4 (Theorem B, r=2); the u-LINEAR sub-case
# is the genuinely new degree-5 structure and is NOT settled here.  The rank-1
# branch (f5 = x1^5) is not addressed here either.  See research_log Phase 22.
#
#   py -u d5_rank2_toppiece.py
import itertools, random, sys
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X = (x1, x2, x3, x4)
NTRIALS = 5


def hess(g):
    return sp.Matrix(4, 4, lambda i, j: sp.diff(g, X[i], X[j]))


def topdeg(poly, d):
    P = sp.Poly(sp.expand(poly), *X)
    return sum(c * sp.prod(v**e for v, e in zip(X, m)) for m, c in P.terms()
               if sum(m) == d)


def det2(g, vs):
    H = sp.Matrix(2, 2, lambda i, j: sp.diff(g, vs[i], vs[j]))
    return sp.expand(H.det())


def main():
    rng = random.Random(2)
    mons5 = [x1**(5 - i) * x2**i for i in range(6)]
    mons4 = [sp.prod(t) for t in itertools.combinations_with_replacement(X, 4)]
    ok = True
    checked = 0
    for trial in range(NTRIALS):
        f5 = sum(rng.randint(-5, 5) * m for m in mons5)
        if det2(f5, (x1, x2)) == 0:      # keep only rank-2 f5
            continue
        f4 = sum(rng.randint(-4, 4) * m for m in mons4)
        f = f5 + f4
        top10 = sp.expand(topdeg(sp.expand(hess(f).det()), 10))
        claim = sp.expand(det2(f5, (x1, x2)) * det2(f4, (x3, x4)))
        match = sp.expand(top10 - claim) == 0
        ok = ok and match
        checked += 1
        print(f'  instance {trial}: [det Hess]_10 == det2(Hess2 f5)*det2(Hess_(x3,x4) f4)?  {match}')
        assert match, f'instance {trial}: factorization FAILED'
    assert checked >= 3, 'too few rank-2 instances sampled'
    print()
    print(f'ALL {checked} INSTANCES MATCH -- rank-2 top-piece reduction holds:')
    print('  det Hess f in C^x  =>  det2(Hess_(x3,x4) f4) == 0.')
    print('(Proof: Laplace-degree argument in the header; instances corroborate.)')
    print('OPEN: pivot-forcing completion (u-linear sub-case) and the r=1 branch '
          '(f5=x1^5) -- see research_log Phase 22.')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
