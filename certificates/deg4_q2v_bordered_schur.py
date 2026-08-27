# deg4_q2v_bordered_schur.py -- branch q2v, TASK Q2V part (1).
#
# Verifies, fail-closed, for  f = e0(x') + x4*e1(x') + (sigma/2)*x4^2,
# x' = (x1,x2,x3), sigma a symbol (generic, != 0 where it is inverted):
#
#   (1a) [PROOF, linear in coefficients] In coordinates (x4, x'),
#        Hess f = [[sigma, g^T], [g, Hess3(e0 + x4*e1)]],  g := grad' e1,
#        verified FULLY SYMBOLICALLY for generic e0 (all monomials of degree
#        <= 4 in x'), generic e1 (degree <= 2), symbolic sigma.  Every entry
#        of both sides is LINEAR in the coefficients of e0, e1, so the
#        symbolic identity check is a proof for ALL such e0, e1.
#
#   (1b) [PROOF, polynomial identity in matrix entries] The bordered/Schur
#        determinant identity for a scalar corner:
#            det [[sigma, g^T],[g, M]] = sigma * det3(M - g g^T / sigma)
#        for a FULLY SYMBOLIC symmetric 3x3 M, symbolic g, symbolic sigma
#        (identity in Q[sigma, g, m][1/sigma]; both sides multiplied by
#        sigma^2 to clear denominators -> polynomial identity, expand == 0).
#        Also the rank-1 form  det = sigma*det M - g^T adj(M) g.
#
#   (1c) [hand argument + symbolic check + instances] Block identification:
#            Hess3(etilde0) + (x4 + e1/sigma)*K
#                == Hess3(e0 + x4*e1) - g g^T / sigma,
#        where etilde0 := e0 - e1^2/(2 sigma), K := Hess3(e1).
#        HAND ARGUMENT (product rule): Hess(e1^2) = 2*e1*K + 2*g g^T, since
#        d_i d_j (e1^2) = d_i(2 e1 d_j e1) = 2 d_i e1 d_j e1 + 2 e1 d_i d_j e1.
#        Hence Hess3 etilde0 = Hess3 e0 - (e1*K + g g^T)/sigma, and adding
#        (x4 + e1/sigma)*K gives Hess3 e0 + x4*K - g g^T/sigma
#        = Hess3(e0 + x4 e1) - g g^T/sigma (K = Hess3 e1).
#        This identity is QUADRATIC in the coefficients of e1, so per project
#        rules the fully-symbolic expand==0 below is recorded as a
#        CONSISTENCY CHECK backing the stated hand argument, and it is
#        additionally verified on random dense rational instances.
#
#   (1d) [consistency check, >= 4 random dense rational instances, symbolic
#        sigma] The full Schur identity of scaffold sec. 1a / ledger Q2(i):
#            det Hess f == sigma * det3( Hess3(e0 - e1^2/(2 sigma))
#                                        + (x4 + e1/sigma) * Hess3(e1) ).
#        (det is cubic in matrix entries, hence NOT linear in coefficients;
#        instances + the proof chain (1a)+(1b)+(1c) give the general case:
#        (1a) and (1b) are proofs, (1c) is the product rule.)
#
#   (1e) [bonus consistency] one instance with deg e1 = 3 (K = K(x') NON-
#        constant): the Schur identity (i) of ledger Q2 holds verbatim
#        (constancy of K is NOT needed for (i), only for the pencil step).
#
# All arithmetic exact (QQ / Q(sigma)).  Determinants: method='berkowitz'.

import sys
import itertools
import random

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
import sympy as sp
from hc4lib import hessian

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Xp = (x1, x2, x3)
sigma = sp.Symbol('sigma')


def generic_poly(vars_, degs, prefix):
    """Fully generic polynomial: one fresh symbol per monomial of each degree."""
    p = sp.Integer(0)
    syms = []
    idx = 0
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            c = sp.Symbol(f'{prefix}{idx}')
            idx += 1
            syms.append(c)
            t = c
            for i in m:
                t *= vars_[i]
            p += t
    return p, syms


def rand_poly(vars_, degs, rnd):
    """Dense random polynomial with nonzero rational coefficients."""
    p = sp.Integer(0)
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            c = sp.Rational(rnd.choice([q for q in range(-9, 10) if q != 0]),
                            rnd.randint(1, 4))
            t = c
            for i in m:
                t *= vars_[i]
            p += t
    return sp.expand(p)


# --------------------------------------------------------------------------
# (1a) bordered form, fully symbolic, linear in coefficients  => PROOF
# --------------------------------------------------------------------------
e0g, ce0 = generic_poly(Xp, [0, 1, 2, 3, 4], 'a')   # 35 generic coefficients
e1g, ce1 = generic_poly(Xp, [0, 1, 2], 'b')         # 10 generic coefficients

fg = e0g + x4*e1g + sp.Rational(1, 2)*sigma*x4**2
H = hessian(fg, [x4, x1, x2, x3])                    # coordinates (x4, x')
gvec = [sp.diff(e1g, w) for w in Xp]                 # g = grad' e1
inner = hessian(sp.expand(e0g + x4*e1g), list(Xp))   # Hess3(e0 + x4 e1)

B = sp.zeros(4, 4)
B[0, 0] = sigma
for i in range(3):
    B[0, i+1] = gvec[i]
    B[i+1, 0] = gvec[i]
    for j in range(3):
        B[i+1, j+1] = inner[i, j]

for i in range(4):
    for j in range(4):
        assert sp.expand(H[i, j] - B[i, j]) == 0, f'bordered form FAILED at {(i, j)}'
print('(1a) PROOF: Hess f == [[sigma, g^T],[g, Hess3(e0+x4 e1)]] '
      'for ALL e0 (deg<=4), e1 (deg<=2), symbolic sigma  [linear in coeffs]')

# --------------------------------------------------------------------------
# (1b) Schur identity as polynomial identity in matrix entries  => PROOF
# --------------------------------------------------------------------------
m11, m12, m13, m22, m23, m33 = sp.symbols('m11 m12 m13 m22 m23 m33')
G1, G2, G3 = sp.symbols('G1 G2 G3')
M = sp.Matrix([[m11, m12, m13], [m12, m22, m23], [m13, m23, m33]])
Gv = sp.Matrix([G1, G2, G3])
Big = sp.zeros(4, 4)
Big[0, 0] = sigma
for i in range(3):
    Big[0, i+1] = Gv[i]
    Big[i+1, 0] = Gv[i]
    for j in range(3):
        Big[i+1, j+1] = M[i, j]

detBig = Big.det(method='berkowitz')
Schur = M - (Gv * Gv.T) / sigma
rhs_schur = sigma * Schur.det(method='berkowitz')
# clear denominators (worst power sigma^-2 on the rhs) and expand:
assert sp.expand(sigma**2 * (detBig - rhs_schur)) == 0, 'Schur matrix identity FAILED'
# rank-1 expansion form: det = sigma*det M - g^T adj(M) g
rhs_adj = sigma * M.det(method='berkowitz') - (Gv.T * M.adjugate() * Gv)[0, 0]
assert sp.expand(detBig - rhs_adj) == 0, 'adjugate form of Schur identity FAILED'
print('(1b) PROOF: det[[sigma,g^T],[g,M]] == sigma*det3(M - g g^T/sigma) '
      '== sigma*detM - g^T adjM g   [fully symbolic matrix entries]')

# --------------------------------------------------------------------------
# (1c) block identification (hand argument above; symbolic check + instances)
# --------------------------------------------------------------------------
def block_check(e0, e1, tag):
    et0 = e0 - e1**2 / (2*sigma)
    K = hessian(e1, list(Xp))
    g = [sp.diff(e1, w) for w in Xp]
    L = hessian(sp.expand(et0), list(Xp)) + (x4 + e1/sigma) * K
    R = hessian(sp.expand(e0 + x4*e1), list(Xp)) - sp.Matrix(3, 3, lambda i, j: g[i]*g[j]) / sigma
    for i in range(3):
        for j in range(3):
            d = sp.expand(sigma * (L[i, j] - R[i, j]))
            assert d == 0, f'[{tag}] block identity FAILED at {(i, j)}'

block_check(e0g, e1g, 'symbolic')
print('(1c) CONSISTENCY (quadratic in e1-coeffs; hand argument = product rule '
      'in header): Hess3(etilde0) + (x4+e1/sigma)K == Hess3(e0+x4 e1) - g g^T/sigma '
      'checked with FULLY SYMBOLIC generic coefficients')

rnd = random.Random(20260809)
for k in range(3):
    e0r = rand_poly(Xp, [2, 3, 4], rnd)
    e1r = rand_poly(Xp, [0, 1, 2], rnd)
    block_check(e0r, e1r, f'inst{k+1}')
print('(1c) CONSISTENCY: block identity re-checked on 3 random dense rational instances')

# --------------------------------------------------------------------------
# (1d) full Schur identity on random dense rational instances, symbolic sigma
# --------------------------------------------------------------------------
def schur_check(e0, e1, tag):
    f = sp.expand(e0 + x4*e1 + sp.Rational(1, 2)*sigma*x4**2)
    D = hessian(f, [x4, x1, x2, x3]).det(method='berkowitz')
    et0 = sp.expand(e0 - e1**2/(2*sigma))
    K = hessian(e1, list(Xp))
    N = hessian(et0, list(Xp)) + (x4 + e1/sigma) * K
    R = sigma * N.det(method='berkowitz')
    # rhs has worst denominator sigma^-2 after the sigma prefactor; sigma^2 clears:
    d = sp.expand(sigma**2 * (D - R))
    assert d == 0, f'[{tag}] full Schur identity FAILED'
    print(f'[{tag}] det Hess f == sigma*det3(Hess3 etilde0 + (x4+e1/sigma)K)  OK')

for k in range(4):
    e0r = rand_poly(Xp, [2, 3, 4], rnd)
    e1r = rand_poly(Xp, [0, 1, 2], rnd)
    schur_check(e0r, e1r, f'(1d) instance {k+1}/4')
print('(1d) CONSISTENCY: full Schur identity holds on 4 dense rational instances '
      '(with (1a)+(1b) proofs and (1c) product rule, this closes the identity '
      'for ALL e0 deg<=4, e1 deg<=2, sigma != 0)')

# --------------------------------------------------------------------------
# (1e) bonus: deg e1 = 3, K(x') non-constant -- Schur identity (i) still holds
# --------------------------------------------------------------------------
e1c = rand_poly(Xp, [0, 1, 2, 3], rnd)
e0c = rand_poly(Xp, [2, 3, 4], rnd)
schur_check(e0c, e1c, '(1e) deg e1 = 3')
print('(1e) CONSISTENCY: ledger Q2(i) needs no constancy of K '
      '(constancy of K enters only at the pencil step (ii))')

print('deg4_q2v_bordered_schur: ALL CHECKS PASS')
