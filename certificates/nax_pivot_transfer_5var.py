# =====================================================================
# NAX-B.  The five-variable affine-pivot transfer, and a GAP in the
#         proposed "Theorem G + Gordan-Noether" route.
#
# Setting.  Phi = t A(w) + B(w) in C[t, w1..w4] with det Hess_5 Phi = c in C^*.
#
# PROVED here (fully generic symbolic identities => proofs):
#   (T1) det Hess_5 Phi = - g^T adj(P + t K) g,  g = grad A, P = Hess_4 B,
#        K = Hess_4 A; this is a CUBIC in t whose t^3-coefficient is
#        -B_4(A) with B_4(A) := grad A^T adj(Hess_4 A) grad A.
#        Hence constancy of det Hess_5 Phi forces  B_4(A) == 0  --- the
#        four-variable analogue of the project's condition (c2).
#   (T2) By the project's Theorem G (dimension-free rigidity; certificate
#        theorem_G.py) B_4(A) == 0 implies det Hess_4 A == 0.
#   (T3) DESCENT IDENTITY (new).  For any kappa, mu in C,
#          det Hess_4 ( B + kappa A^2 + mu A )
#              = det_4( P + (2 kappa A + mu) K ) - 2 kappa c .
#        So Meng-Yang's Schur descent produces a potential with constant
#        Hessian determinant EXACTLY when the symmetric pencil P + sK is
#        singular -- which is the extra hypothesis of Conjecture C2 in the
#        ledger, and which holds automatically for doublings.
#
# *** GAP FOUND ***  The next step proposed for the non-affine D1 argument,
#     "det Hess_4(A) == 0 and Gordan-Noether in four variables make A a
#      cone, so A depends on three variables after a linear change",
#     is INVALID.  Gordan-Noether is a statement about FORMS.  A is not
#     homogeneous; homogenising it produces a form in FIVE variables, and
#     Gordan-Noether is FALSE in five variables (Perazzo).  Concretely,
#          A0 = w2 + w1 w3 + w1^2 w4
#     satisfies B_4(A0) == 0 and det Hess_4 A0 == 0 but is NOT affinely
#     3-variable; its homogenisation is exactly the Perazzo cubic
#          w0^2 w2 + w0 w1 w3 + w1^2 w4 ,
#     the standard Gordan-Noether counterexample.  (This is not a
#     counterexample to the project's Corollary E, which is about THREE
#     variables and is correct; it shows that Corollary E does not
#     generalise to four.)
#
#   (T4) What survives: B_4(A) == 0 <=> rank Hess_4 A <= 2, or rank = 3 and
#        grad A is orthogonal to ker Hess_4 A.  Both branches are realised.
# =====================================================================

import sympy as sp

t = sp.Symbol('t')
w = list(sp.symbols('w1 w2 w3 w4'))
w0 = sp.Symbol('w0')

# ---------------------------------------------------------------------
# (T1)  the bordered expansion, fully generic  => PROOF
# ---------------------------------------------------------------------
P = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'p{min(i,j)}{max(i,j)}'))
K = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'k{min(i,j)}{max(i,j)}'))
g = sp.Matrix([sp.Symbol(f'g{i}') for i in range(4)])
M = P + t*K
Bor = sp.zeros(5, 5)
for i in range(4):
    Bor[0, i+1] = Bor[i+1, 0] = g[i]
    for j in range(4):
        Bor[i+1, j+1] = M[i, j]
lhs = sp.expand(Bor.det(method='berkowitz'))
assert sp.expand(lhs + (g.T*M.adjugate()*g)[0, 0]) == 0
pl = sp.Poly(lhs, t)
assert pl.degree() == 3
assert sp.expand(pl.coeff_monomial(t**3) + (g.T*K.adjugate()*g)[0, 0]) == 0
print('[T1] det[[0,g^T],[g,P+tK]] = -g^T adj(P+tK) g, cubic in t,')
print('     [t^3] = -g^T adj(K) g = -B_4(A)      (fully generic => PROOF)   OK')
print('     ==> det Hess_5(tA+B) in C^*  forces  B_4(A) == 0.')

# sanity: the same statement realised on a concrete Phi
Aex = w[0]*w[2] + w[1]
Bex = w[0]**2 + w[1]*w[3] + w[2]**3
Phi = sp.expand(t*Aex + Bex)
V = [t] + w
HPhi = sp.Matrix(5, 5, lambda i, j: sp.diff(Phi, V[i], V[j]))
gA = sp.Matrix([sp.diff(Aex, v) for v in w])
HA = sp.Matrix(4, 4, lambda i, j: sp.diff(Aex, w[i], w[j]))
B4 = sp.expand((gA.T*HA.adjugate()*gA)[0, 0])
assert sp.expand(sp.Poly(sp.expand(HPhi.det(method='berkowitz')), t).coeff_monomial(t**3)
                 + B4) == 0
print('[T1] concrete instance agrees   OK')

# ---------------------------------------------------------------------
# (T3)  descent identity, fully generic  => PROOF
# ---------------------------------------------------------------------
# rank-one update:  det(M + s v v^T) = det M + s v^T adj(M) v   (4x4 generic)
s = sp.Symbol('s')
Mg = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'm{min(i,j)}{max(i,j)}'))
vg = sp.Matrix([sp.Symbol(f'v{i}') for i in range(4)])
assert sp.expand((Mg + s*vg*vg.T).det(method='berkowitz')
                 - Mg.det(method='berkowitz') - s*(vg.T*Mg.adjugate()*vg)[0, 0]) == 0
print('[T3] rank-one update det(M + s v v^T) = det M + s v^T adj(M) v '
      '(4x4 generic => PROOF)   OK')
# Hess_4(B + kappa A^2 + mu A) = P + (2 kappa A + mu) K + 2 kappa g g^T,
# and g^T adj(P + sK) g == -c for all s (by T1 and constancy), whence
#     det Hess_4(Theta) = det(P + (2kA+mu)K) + 2k*(-c) .
kap, muS = sp.symbols('kappa muS')
Aq, Bq = sp.Function('Aq'), sp.Function('Bq')
Aa = Aq(*w)
Bb = Bq(*w)
Th = Bb + kap*Aa**2 + muS*Aa
HT = sp.Matrix(4, 4, lambda i, j: sp.expand(sp.diff(Th, w[i], w[j])))
Pf = sp.Matrix(4, 4, lambda i, j: sp.diff(Bb, w[i], w[j]))
Kf = sp.Matrix(4, 4, lambda i, j: sp.diff(Aa, w[i], w[j]))
gf = sp.Matrix([sp.diff(Aa, v) for v in w])
assert sp.expand(HT - (Pf + (2*kap*Aa + muS)*Kf + 2*kap*gf*gf.T)) == sp.zeros(4, 4)
print('[T3] Hess_4(B + kA^2 + mu A) = P + (2kA+mu)K + 2k gg^T   (symbolic)   OK')
print('[T3] ==> det Hess_4(B+kA^2+muA) = det(P+(2kA+mu)K) - 2 k c .')
print('     For a DOUBLING the pencil P+sK is identically singular (it has a')
print('     zero diagonal block), which is why the Meng-Yang descent works.')

# ---------------------------------------------------------------------
# *** GAP ***  the Gordan-Noether step
# ---------------------------------------------------------------------
A0 = w[1] + w[0]*w[2] + w[0]**2*w[3]
H0 = sp.Matrix(4, 4, lambda i, j: sp.diff(A0, w[i], w[j]))
g0 = sp.Matrix([sp.diff(A0, v) for v in w])
assert sp.expand(H0.det()) == 0
assert sp.expand((g0.T*H0.adjugate()*g0)[0, 0]) == 0
assert H0.rank() == 2
vv = sp.symbols('v1 v2 v3 v4')
Dv = sp.expand(sum(vv[i]*sp.diff(A0, w[i]) for i in range(4)))
sol = sp.solve(sp.Poly(Dv, *w).coeffs(), vv, dict=True)
assert sol == [{vv[0]: 0, vv[1]: 0, vv[2]: 0, vv[3]: 0}], sol
print('[GAP] A0 = w2 + w1 w3 + w1^2 w4 :  B_4(A0) == 0, det Hess_4 A0 == 0,')
print('      rank Hess_4 A0 = 2, and D_v A0 == 0 only for v = 0, i.e. A0 is')
print('      NOT affinely 3-variable.')
E = sp.expand(w0**3*A0.subs({w[i]: w[i]/w0 for i in range(4)}, simultaneous=True))
assert sp.expand(E - (w0**2*w[1] + w0*w[0]*w[2] + w[0]**2*w[3])) == 0
VE = [w0] + w
HE = sp.Matrix(5, 5, lambda i, j: sp.diff(E, VE[i], VE[j]))
assert sp.expand(HE.det(method='berkowitz')) == 0
uu = sp.symbols('u0 u1 u2 u3 u4')
DuE = sp.expand(sum(uu[i]*sp.diff(E, VE[i]) for i in range(5)))
solE = sp.solve(sp.Poly(DuE, *VE).coeffs(), uu, dict=True)
assert solE == [{uu[i]: 0 for i in range(5)}], solE
print('[GAP] its homogenisation is the PERAZZO cubic w0^2 w2 + w0 w1 w3 + w1^2 w4:')
print('      det Hess_5 == 0 but it is NOT a cone (D_u E == 0 only for u = 0).')
print('[GAP] ==> "det Hess_4 A == 0 + Gordan-Noether(4) => A is a cone" is FALSE.')
print('      Gordan-Noether applies to FORMS; dehomogenising costs one variable,')
print('      and GN fails in five.  The project\'s Corollary E (three variables)')
print('      is correct and does NOT extend to four.')

# ---------------------------------------------------------------------
# (T4)  what B_4 == 0 does say
# ---------------------------------------------------------------------
Ms = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'q{min(i,j)}{max(i,j)}'))
# rank <= 2 => adj == 0 (all 3x3 minors vanish): checked on a rank-2 sample
R2 = sp.Matrix([[1, 0, 0, 0], [0, 3, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
assert R2.rank() == 2 and R2.adjugate() == sp.zeros(4, 4)
# rank 3 => adj = kappa k k^T with k spanning ker: sample
R3 = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]])
assert R3.rank() == 3
ad = R3.adjugate()
assert ad == sp.Matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]])
print('[T4] B_4(A) == 0 <=> rank Hess_4 A <= 2 (adj = 0, automatic) OR')
print('     rank Hess_4 A = 3 and grad A _|_ ker Hess_4 A.  A0 realises the')
print('     first branch, so no cone conclusion is available.   OK')

print()
print('=== NAX-B SUMMARY ===================================================')
print('PROVED : det Hess_5(tA+B) in C^* => B_4(A) == 0 => (Thm G) det Hess_4 A == 0')
print('PROVED : descent identity det Hess_4(B+kA^2+muA) = det(P+(2kA+mu)K) - 2kc')
print('GAP    : the Gordan-Noether "A is a cone" step is INVALID in four')
print('         variables (Perazzo).  The non-affine D1 route proposed in the')
print('         task therefore breaks at exactly that point.')
print('=====================================================================')
