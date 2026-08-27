# deg4_q2v_pencil_forcing.py -- branch q2v, TASK Q2V part (3).
#
# Pencil-forcing step of Theorem Q2 (scaffold sec. 1a, line 33-34):
# if det Hess f == c (constant) for
#     f = e0(x') + x4*e1(x') + (sigma/2)*x4^2,  sigma != 0, deg e1 <= 2,
# then with etilde0 := e0 - e1^2/(2 sigma), K := Hess3 e1 (CONSTANT),
#     det3( Hess3 etilde0 + s*K ) == c/sigma   identically in (s, x').   (P)
#
# WHY IT IS FORCED (hand argument, machine-instantiated below): the maps
#     (x', x4) |-> (x', s := x4 + e1(x')/sigma)   and
#     (x', s ) |-> (x', x4 := s - e1(x')/sigma)
# are mutually inverse TRIANGULAR POLYNOMIAL substitutions over Q(sigma)
# (this is the precise form of "u = x4 + e1/sigma sweeps C for fixed x'").
# Substituting x4 = s - e1/sigma into the verified Schur identity
#     det Hess f = sigma * det3(Hess3 etilde0 + (x4 + e1/sigma)*K)
# (certificate deg4_q2v_bordered_schur.py) turns its right side literally into
# sigma * det3(Hess3 etilde0 + s*K); if the left side is the constant c, the
# result is the identity (P).  Constants are preserved by ANY substitution,
# so no genericity is needed -- only sigma != 0.
#
# CHECKS (all fail-closed, exact):
#  (3a) [PROOF] the two substitutions are mutually inverse:
#       (x4 + e1/sigma)|_{x4 -> s - e1/sigma} == s  and
#       (s - e1/sigma)|_{s -> x4 + e1/sigma} == x4, with FULLY SYMBOLIC
#       generic e1 (deg <= 2, 10 coefficient symbols) -- linear in coeffs.
#  (3b) [PROOF, matrix-entry identity] the s^3-coefficient of
#       det3(M + s*K) equals det3(K) for FULLY SYMBOLIC symmetric M, K.
#       COROLLARY: (P) with a CONSTANT right side forces det3 K = 0, i.e.
#       rank K <= 2 in any constant-det instance (checked on instances).
#  (3c) [consistency, 3 random dense rational instances, symbolic sigma]
#       substitution identity  (det Hess f)|_{x4 -> s - e1/sigma}
#            == sigma * det3(Hess3 etilde0 + s*K)  identically in (s, x')
#       (det is not linear in coefficients: instances; generality follows
#       from the Schur identity certificate + (3a)).
#  (3d) [the actual forcing, on 3 CONSTANT-DET instances, K of rank 0, 1, 2,
#       symbolic sigma; random dense rational data inside a triangular family]
#       det Hess f is literally the constant -sigma, and
#       det3(Hess3 etilde0 + s*K), expanded as a polynomial in (s, x1,x2,x3),
#       is LITERALLY the constant -1 == c/sigma  (total degree 0).
#
# Instance family for (3d): etilde0 = x1*x3 + x2^2/2 + phi(x1) with phi a
# random dense quartic; then Hess3 etilde0 + s*K stays "triangular" for the
# chosen K and its det is -1 for ALL (s, x'):  rank 0: e1 affine (K = 0);
# rank 1: e1 = (r/2)x1^2 + affine (K = r*E11); rank 2: e1 = (r/2)x1^2 +
# u*x1*x2 + affine (K = [[r,u,0],[u,0,0],[0,0,0]], det2 = -u^2 != 0).
# e0 := etilde0 + e1^2/(2 sigma) makes f an EXACT member of the Q2 class with
# det Hess f == -sigma (c = -sigma, c/sigma = -1).
#
# All arithmetic exact.  Determinants: method='berkowitz'.  No Groebner needed.

import sys
import itertools
import random

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
import sympy as sp
from hc4lib import hessian

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Xp = (x1, x2, x3)
sigma = sp.Symbol('sigma')
s = sp.Symbol('s')


def generic_poly(vars_, degs, prefix):
    p = sp.Integer(0)
    idx = 0
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            c = sp.Symbol(f'{prefix}{idx}')
            idx += 1
            t = c
            for i in m:
                t *= vars_[i]
            p += t
    return p


def rand_poly(vars_, degs, rnd):
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
# (3a) the substitutions are mutually inverse  [PROOF: linear in e1-coeffs]
# --------------------------------------------------------------------------
e1g = generic_poly(Xp, [0, 1, 2], 'b')
u_fwd = x4 + e1g/sigma
assert sp.expand(sigma*(u_fwd.subs(x4, s - e1g/sigma) - s)) == 0, '(3a) fwd FAILED'
assert sp.expand(sigma*((s - e1g/sigma).subs(s, x4 + e1g/sigma) - x4)) == 0, \
    '(3a) bwd FAILED'
print('(3a) PROOF: x4 <-> s = x4 + e1(x\')/sigma are mutually inverse triangular '
      'polynomial substitutions (generic e1, symbolic sigma)')

# --------------------------------------------------------------------------
# (3b) s^3-coefficient of det3(M + sK) is det3 K  [PROOF: matrix entries]
# --------------------------------------------------------------------------
msym = sp.symbols('m11 m12 m13 m22 m23 m33')
ksym = sp.symbols('k11 k12 k13 k22 k23 k33')
M = sp.Matrix([[msym[0], msym[1], msym[2]],
               [msym[1], msym[3], msym[4]],
               [msym[2], msym[4], msym[5]]])
K = sp.Matrix([[ksym[0], ksym[1], ksym[2]],
               [ksym[1], ksym[3], ksym[4]],
               [ksym[2], ksym[4], ksym[5]]])
pencil_det = sp.expand((M + s*K).det(method='berkowitz'))
c3 = sp.Poly(pencil_det, s).coeff_monomial(s**3)
assert sp.expand(c3 - K.det(method='berkowitz')) == 0, '(3b) FAILED'
print('(3b) PROOF: [s^3] det3(M + sK) == det3 K  (fully symbolic symmetric M, K); '
      'hence (P) with constant value forces det3 K = 0, i.e. rank K <= 2')

# --------------------------------------------------------------------------
# (3c) substitution identity on random (non-constant-det) instances
# --------------------------------------------------------------------------
rnd = random.Random(20260811)

def substitution_check(e0, e1, tag):
    f = sp.expand(e0 + x4*e1 + sp.Rational(1, 2)*sigma*x4**2)
    D = hessian(f, [x4, x1, x2, x3]).det(method='berkowitz')
    et0 = sp.expand(e0 - e1**2/(2*sigma))
    Ki = hessian(e1, list(Xp))
    P = sigma * (hessian(et0, list(Xp)) + s*Ki).det(method='berkowitz')
    d = sp.expand(sigma**2 * (D.subs(x4, s - e1/sigma) - P))
    assert d == 0, f'[{tag}] substitution identity FAILED'
    print(f'[{tag}] (det Hess f)|_(x4 -> s - e1/sigma) == '
          'sigma*det3(Hess3 etilde0 + sK)  OK')

for k in range(3):
    e0r = rand_poly(Xp, [2, 3, 4], rnd)
    e1r = rand_poly(Xp, [0, 1, 2], rnd)
    substitution_check(e0r, e1r, f'(3c) instance {k+1}/3')
print('(3c) CONSISTENCY: substitution identity holds identically in (s, x\') even '
      'for NON-constant det Hess f; constancy of det Hess f therefore transfers '
      'verbatim to the pencil')

# --------------------------------------------------------------------------
# (3d) constant-det instances: the pencil is literally constant in (s, x')
# --------------------------------------------------------------------------
def nz(rnd):
    return sp.Rational(rnd.choice([q for q in range(-9, 10) if q != 0]),
                       rnd.randint(1, 4))

def forcing_check(e1, expected_rank, tag):
    phi = rand_poly((x1,), [0, 1, 2, 3, 4], rnd)          # dense quartic in x1
    et0 = sp.expand(x1*x3 + sp.Rational(1, 2)*x2**2 + phi)
    e0 = sp.expand(et0 + e1**2/(2*sigma))
    f = sp.expand(e0 + x4*e1 + sp.Rational(1, 2)*sigma*x4**2)
    Ki = hessian(e1, list(Xp))
    # K constant (deg e1 <= 2) and of the advertised rank, det3 K = 0:
    assert all(not Ki[i, j].free_symbols & set(Xp)
               for i in range(3) for j in range(3)), f'[{tag}] K not constant'
    assert Ki.rank() == expected_rank, f'[{tag}] rank K != {expected_rank}'
    assert Ki.det(method='berkowitz') == 0, f'[{tag}] det3 K != 0'
    # det Hess f is LITERALLY the constant c = -sigma:
    D = hessian(f, [x4, x1, x2, x3]).det(method='berkowitz')
    d = sp.expand(sigma**4 * (D + sigma))
    assert d == 0, f'[{tag}] det Hess f != -sigma'
    # pencil det, expanded in (s, x'): literally the constant -1 == c/sigma:
    E = sp.expand((hessian(et0, list(Xp)) + s*Ki).det(method='berkowitz'))
    assert E == -1, f'[{tag}] pencil not the constant -1: {E}'
    assert sp.Poly(E, s, x1, x2, x3).total_degree() == 0, f'[{tag}] not constant'
    # and it matches c/sigma exactly:
    assert sp.cancel(E - (-sigma)/sigma) == 0, f'[{tag}] pencil != c/sigma'
    print(f'[{tag}] det Hess f == -sigma (constant) and '
          'det3(Hess3 etilde0 + sK) == -1 == c/sigma LITERALLY in (s, x\')  OK')

# rank K = 0: e1 affine
e1_r0 = sp.expand(nz(rnd)*x1 + nz(rnd)*x2 + nz(rnd)*x3 + nz(rnd))
forcing_check(e1_r0, 0, '(3d) rank K = 0')
# rank K = 1: e1 = (r/2) x1^2 + affine
e1_r1 = sp.expand(nz(rnd)/2*x1**2 + nz(rnd)*x1 + nz(rnd)*x3 + nz(rnd))
forcing_check(e1_r1, 1, '(3d) rank K = 1')
# rank K = 2: e1 = (r/2) x1^2 + u x1 x2 + affine
e1_r2 = sp.expand(nz(rnd)/2*x1**2 + nz(rnd)*x1*x2 + nz(rnd)*x1
                  + nz(rnd)*x2 + nz(rnd))
forcing_check(e1_r2, 2, '(3d) rank K = 2')
print('(3d) rank K = 3 is IMPOSSIBLE in a constant-det instance by (3b) '
      '(the s^3-coefficient det3 K must vanish) -- the rank 0/1/2 instances '
      'above exhaust the possible K-ranks of the Q2 pencil')

print('deg4_q2v_pencil_forcing: ALL CHECKS PASS')
