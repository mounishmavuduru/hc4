# deg4_q2v_ap_identities.py -- branch q2v, TASK Q2V part (4).
#
# Machine-verification of the identity machinery behind scaffold sec. 1b
# (the sigma = 0 affine-pivot class AP4):  f = e0(x') + x4*e1(x'),
# g := grad' e1 (affine entries when deg e1 <= 2), K := Hess3 e1 (constant),
#     det Hess f = - g^T adj3( Hess3 e0 + x4*K ) g   == c            (AP)
# identically in (x', x4), and (AP) is EQUIVALENT to the three identities
#     (c0)  g^T (adj3 Hess3 e0) g == -c
#     (c1)  g^T L(Hess3 e0, K) g == 0,   L(P,K) := adj(P+K) - adjP - adjK
#     (c2)  g^T (adj3 K) g == 0.
#
# CHECKS (all fail-closed, exact arithmetic):
#  (4a) [PROOF, polynomial identity in matrix entries] For FULLY SYMBOLIC
#       GENERAL (not even symmetric) 3x3 matrices P, K and a scalar t:
#           adj3(P + t*K) == adj3(P) + t*L(P,K) + t^2*adj3(K),
#       entrywise.  (adj entries are quadratic forms in the matrix entries;
#       the t-expansion of adj(P+tK) is quadratic in t with mixed term
#       bilinear in (P,K); evaluating at t = 1 identifies the mixed term as
#       adj(P+K) - adjP - adjK.  The expand==0 below IS a proof: it is a
#       polynomial identity in the 18 entry symbols and t.)
#  (4b) [PROOF, polynomial identity in matrix entries] Zero-corner bordered
#       determinant, general 3x3 D and general column vectors b, c:
#           det [[0, b^T], [c, D]] == - b^T adj3(D) c.
#  (4c) [PROOF, matrix entries] For symmetric symbolic P, K and symbolic
#       g, t, c:
#           -g^T adj3(P + tK) g  ==  -g^T adjP g - t*(g^T L g) - t^2*(g^T adjK g),
#       and the coefficients of (as a polynomial in t) of
#       -g^T adj3(P+tK) g - c  are EXACTLY
#           t^0: -g^T adjP g - c,   t^1: -g^T L g,   t^2: -g^T adjK g.
#       EQUIVALENCE (AP) <=> (c0)&(c1)&(c2): in the application P = Hess3 e0,
#       K, g are x4-free and t = x4; a polynomial in x4 with x4-free
#       coefficients vanishes identically iff every coefficient vanishes.
#       This gives both directions of the equivalence at once.
#  (4d) [PROOF, linear in coefficients] For generic e0 (deg <= 4) and e1
#       (deg <= 2) with FULLY SYMBOLIC coefficients: in coordinates (x4, x'),
#           Hess(e0 + x4*e1) == [[0, g^T], [g, Hess3 e0 + x4*K]],
#       K == Hess3 e1 is x-free (constant), and every entry of g has total
#       degree <= 1 in x' (affine entries).  All statements linear in the
#       coefficients => proof for all such e0, e1.
#       CHAIN: (4d) + (4b) (specialized to b = c = g, D = Hess3 e0 + x4 K)
#       prove  det Hess f == -g^T adj3(Hess3 e0 + x4 K) g  for ALL e0, e1:
#       specialization of a proven polynomial identity is rigorous.
#  (4e) [consistency, 3 random dense rational instances] the composed
#       statements: det Hess f == -g^T adj3(Hess3 e0 + x4*K) g, and its
#       x4-coefficients match (-g^T adjP g, -g^T L g, -g^T adjK g)
#       (belt-and-braces for the chain (4d)+(4b)+(4a); det/adj are not
#       linear in the coefficients, hence instance-labelled per rules).
#
# All arithmetic exact.  Determinants: method='berkowitz'.

import sys
import itertools
import random

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
import sympy as sp
from hc4lib import hessian

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Xp = (x1, x2, x3)
t = sp.Symbol('t')
c_sym = sp.Symbol('c')


def generic_poly(vars_, degs, prefix):
    p = sp.Integer(0)
    idx = 0
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            cc = sp.Symbol(f'{prefix}{idx}')
            idx += 1
            tt = cc
            for i in m:
                tt *= vars_[i]
            p += tt
    return p


def rand_poly(vars_, degs, rnd):
    p = sp.Integer(0)
    for d in degs:
        for m in itertools.combinations_with_replacement(range(len(vars_)), d):
            cc = sp.Rational(rnd.choice([q for q in range(-9, 10) if q != 0]),
                             rnd.randint(1, 4))
            tt = cc
            for i in m:
                tt *= vars_[i]
            p += tt
    return sp.expand(p)


def Lmix(P, K):
    return (P + K).adjugate() - P.adjugate() - K.adjugate()


# --------------------------------------------------------------------------
# (4a) adjugate pencil expansion, GENERAL 3x3 P, K  [PROOF]
# --------------------------------------------------------------------------
Pg = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'p{i+1}{j+1}'))
Kg = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'k{i+1}{j+1}'))
lhs = (Pg + t*Kg).adjugate()
rhs = Pg.adjugate() + t*Lmix(Pg, Kg) + t**2*Kg.adjugate()
for i in range(3):
    for j in range(3):
        assert sp.expand(lhs[i, j] - rhs[i, j]) == 0, f'(4a) FAILED at {(i, j)}'
print('(4a) PROOF: adj3(P + tK) == adjP + t*(adj(P+K)-adjP-adjK) + t^2*adjK '
      'for FULLY SYMBOLIC GENERAL 3x3 P, K (18 entry symbols)')

# --------------------------------------------------------------------------
# (4b) zero-corner bordered determinant, general D, b, c  [PROOF]
# --------------------------------------------------------------------------
Dg = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'd{i+1}{j+1}'))
bg = sp.Matrix([sp.Symbol('b1'), sp.Symbol('b2'), sp.Symbol('b3')])
cg = sp.Matrix([sp.Symbol('c1'), sp.Symbol('c2'), sp.Symbol('c3')])
Big = sp.zeros(4, 4)
for i in range(3):
    Big[0, i+1] = bg[i]
    Big[i+1, 0] = cg[i]
    for j in range(3):
        Big[i+1, j+1] = Dg[i, j]
det_b = Big.det(method='berkowitz')
assert sp.expand(det_b + (bg.T * Dg.adjugate() * cg)[0, 0]) == 0, '(4b) FAILED'
print('(4b) PROOF: det[[0, b^T],[c, D]] == -b^T adj3(D) c '
      'for FULLY SYMBOLIC general 3x3 D and vectors b, c')

# --------------------------------------------------------------------------
# (4c) (AP) <=> (c0),(c1),(c2) coefficient decomposition  [PROOF]
# --------------------------------------------------------------------------
ms = sp.symbols('P11 P12 P13 P22 P23 P33')
ks = sp.symbols('K11 K12 K13 K22 K23 K33')
Ps = sp.Matrix([[ms[0], ms[1], ms[2]], [ms[1], ms[3], ms[4]], [ms[2], ms[4], ms[5]]])
Ks = sp.Matrix([[ks[0], ks[1], ks[2]], [ks[1], ks[3], ks[4]], [ks[2], ks[4], ks[5]]])
gs = sp.Matrix([sp.Symbol('g1'), sp.Symbol('g2'), sp.Symbol('g3')])
Q = sp.expand((-(gs.T * (Ps + t*Ks).adjugate() * gs))[0, 0])
c0_expr = sp.expand((gs.T * Ps.adjugate() * gs)[0, 0])       # g^T adjP g
c1_expr = sp.expand((gs.T * Lmix(Ps, Ks) * gs)[0, 0])        # g^T L g
c2_expr = sp.expand((gs.T * Ks.adjugate() * gs)[0, 0])       # g^T adjK g
assert sp.expand(Q - (-c0_expr - t*c1_expr - t**2*c2_expr)) == 0, '(4c) expansion FAILED'
QP = sp.Poly(Q - c_sym, t)
assert sp.expand(QP.coeff_monomial(sp.Integer(1)) - (-c0_expr - c_sym)) == 0, '(4c) t^0 FAILED'
assert sp.expand(QP.coeff_monomial(t) - (-c1_expr)) == 0, '(4c) t^1 FAILED'
assert sp.expand(QP.coeff_monomial(t**2) - (-c2_expr)) == 0, '(4c) t^2 FAILED'
print('(4c) PROOF: -g^T adj3(P+tK) g - c has t-coefficients exactly '
      '(-g^T adjP g - c, -g^T L g, -g^T adjK g)  [symmetric symbolic P, K; '
      'since P, K, g are x4-free and t = x4, (AP) <=> (c0) & (c1) & (c2)]')

# --------------------------------------------------------------------------
# (4d) AP bordered form, generic coefficients  [PROOF: linear in coeffs]
# --------------------------------------------------------------------------
e0g = generic_poly(Xp, [0, 1, 2, 3, 4], 'a')
e1g = generic_poly(Xp, [0, 1, 2], 'b')
fg = e0g + x4*e1g
Kg2 = hessian(e1g, list(Xp))
gv = [sp.diff(e1g, w) for w in Xp]
H = hessian(fg, [x4, x1, x2, x3])
inner = hessian(e0g, list(Xp)) + x4*Kg2
for i in range(3):
    assert sp.expand(H[0, i+1] - gv[i]) == 0, '(4d) border FAILED'
    assert sp.expand(H[i+1, 0] - gv[i]) == 0, '(4d) border FAILED'
    for j in range(3):
        assert sp.expand(H[i+1, j+1] - inner[i, j]) == 0, '(4d) inner FAILED'
assert H[0, 0] == 0, '(4d) corner FAILED'
# K constant, g affine:
assert all(not Kg2[i, j].free_symbols & set(Xp)
           for i in range(3) for j in range(3)), '(4d) K not constant'
for gi in gv:
    dpoly = sp.Poly(gi, *Xp)
    assert dpoly.total_degree() <= 1, '(4d) g not affine'
print('(4d) PROOF: Hess(e0 + x4 e1) == [[0, g^T],[g, Hess3 e0 + x4 K]], '
      'K constant, g affine, for ALL e0 (deg<=4), e1 (deg<=2)  [linear in '
      'coeffs].  With (4b) (b = c = g, D = Hess3 e0 + x4 K) this proves '
      '(AP)\'s determinant formula for the whole class by specialization.')

# --------------------------------------------------------------------------
# (4e) composed statements on random dense rational instances  [consistency]
# --------------------------------------------------------------------------
rnd = random.Random(20260812)
for kk in range(3):
    e0r = rand_poly(Xp, [2, 3, 4], rnd)
    e1r = rand_poly(Xp, [0, 1, 2], rnd)
    fr = sp.expand(e0r + x4*e1r)
    Kr = hessian(e1r, list(Xp))
    gr = sp.Matrix([sp.diff(e1r, w) for w in Xp])
    Pr = hessian(e0r, list(Xp))
    D = hessian(fr, [x4, x1, x2, x3]).det(method='berkowitz')
    R = sp.expand((-(gr.T * (Pr + x4*Kr).adjugate() * gr))[0, 0])
    assert sp.expand(D - R) == 0, f'(4e) inst {kk+1}: det formula FAILED'
    DP = sp.Poly(sp.expand(D), x4)
    a0 = sp.expand((gr.T * Pr.adjugate() * gr)[0, 0])
    a1 = sp.expand((gr.T * Lmix(Pr, Kr) * gr)[0, 0])
    a2 = sp.expand((gr.T * Kr.adjugate() * gr)[0, 0])
    assert sp.expand(DP.coeff_monomial(sp.Integer(1)) + a0) == 0, f'(4e) inst {kk+1}: c0 FAILED'
    assert sp.expand(DP.coeff_monomial(x4) + a1) == 0, f'(4e) inst {kk+1}: c1 FAILED'
    assert sp.expand(DP.coeff_monomial(x4**2) + a2) == 0, f'(4e) inst {kk+1}: c2 FAILED'
    print(f'(4e) instance {kk+1}/3: det Hess(e0 + x4 e1) == -g^T adj3(Hess3 e0 '
          '+ x4 K) g, x4-coefficients == -(c0|c1|c2)-forms  OK')
print('(4e) CONSISTENCY: composed (AP) machinery re-checked on 3 dense rational '
      'instances (general case follows from the (4a)-(4d) proofs by specialization)')

print('deg4_q2v_ap_identities: ALL CHECKS PASS')
