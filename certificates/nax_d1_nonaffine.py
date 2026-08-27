# =====================================================================
# NAX-A.  Removing the "affine only" caveat from D1  --  what is true
#         and what is FALSE.
#
# D1 (existing, `pivot_obstruction_d1.py`) proves: for lam != 0, mu, no
# CONSTANT v != 0 has D_v^2 psi_{lam,mu} == const.  Equivalently: no AFFINE
# change of coordinates of C^5 makes a member of the Meng-Yang HC_5 family
# affine-linear in a variable.
#
# This certificate establishes, in exact arithmetic:
#
#  (A0) REFORMULATION.  For a polynomial automorphism sigma of C^5,
#       Psi o sigma is affine-linear in the coordinate t  <=>  u^2 Psi == 0,
#       where u := sigma_*(d/dt) is the polynomial derivation with
#       u(h) o sigma = d(h o sigma)/dt.  u is a *coordinate* locally
#       nilpotent derivation.  (For sigma affine, u = D_v with v constant:
#       this is exactly D1's hypothesis.)
#
#  (A1) **THE NAIVE STRENGTHENING OF D1 IS FALSE.**  There is an explicit
#       polynomial automorphism tau of C^5 (triangular: linear in y with
#       coefficients in C[x1,x2], det J tau = 1) such that EVERY member
#       Psi_{lam,mu} of the Meng-Yang family becomes affine-linear in the
#       first new coordinate, with pivot coefficient  -x2(1+x1x2).
#       The underlying derivation is  u = 3x1 d/dy1 - (1+x1x2) d/dy2,
#       which satisfies u(A) == 0 outright.
#       So "no affine change of coordinates" CANNOT be upgraded to
#       "no polynomial change of coordinates": the caveat is not removable
#       as stated, and the extra hypothesis det Hess(Psi o sigma) in C^*
#       is *load-bearing*, not cosmetic.
#
#  (A2) The determinant is destroyed by tau, exactly as the corrected
#       theorem demands:  det Hess(Psi_{lam,mu} o tau) = (lam/9) * R^2 with
#       R an explicit NON-constant polynomial (for every lam != 0, mu).
#
#  (A3) MASTER FORMULA (new, proved here, fully generic => a proof).
#       Let a, b in C[x1,x2]^3, F(x1,x2,x3) := a(x1,x2) x3 + b(x1,x2),
#       and  Theta := <b,z> + (lam/2)<a,z>^2 + mu<a,z>  in C[x1,x2,z1,z2,z3].
#       Then
#            det Hess_5 Theta  =  lam * ( Jac F )^2 |_{x3 = lam<a,z> + mu}.
#       Meng-Yang's family is the case a = dF/dx3, b = F|_{x3=0} for
#       Alpoege's F; Jac F = -2 gives det Hess = 4 lam, reproving the
#       family identity structurally.  This *explains* the Schur descent:
#       it is the Keller property of F, transported.
#
#  (A4) FIBER-LINEAR REDUCTION (new).  For any automorphism of the form
#       sigma(x,y) = (phi(x), Q(x) y + r(x)) with phi in Aut(C^2),
#       Q in GL_3(C[x1,x2]), r in C[x1,x2]^3:
#           Psi o sigma is affine-linear in a y-coordinate AND has
#           det Hess in C^*
#       <=>  G := Q^T . (F o (phi x id))  is a KELLER map of C^3 whose
#            first component does not involve x3.
#       Necessary conditions on G = (g(x1,x2), G2, G3), G_i = a'_i x3 + b'_i:
#            (i)  a'_3 D a'_2 = a'_2 D a'_3,   D := J(g, .),
#           (ii)  a'_3 D b'_2 - a'_2 D b'_3 = kappa in C^*,
#       whence D a'_i = h a'_i and, with W := a'_3 b'_2 - a'_2 b'_3,
#           J(g, W) = h W + kappa.
#       In particular grad g vanishes NOWHERE (at a critical point of g,
#       (ii) would read kappa = 0), and when h == 0 the pair (g, W) is a
#       planar KELLER pair, so by Moh's theorem g is a coordinate of
#       C[x1,x2] whenever deg g, deg W <= 100.
#
# Everything below is exact (rational / symbolic).  Fail-closed.
# =====================================================================

import sympy as sp

x1, x2, x3 = sp.symbols('x1 x2 x3')
y1, y2, y3 = sp.symbols('y1 y2 y3')
z1, z2, z3 = sp.symbols('z1 z2 z3')
lam, mu = sp.symbols('lam mu')

U = 1 + x1*x2

# ---------------------------------------------------------------------
# 0.  The data: Alpoege's F (dim 3) and the Meng-Yang family.
# ---------------------------------------------------------------------
a = [U**3, 3*x1*U**2, -x1**3]                      # dF/dx3
b = [x2**2*U*(4 + 3*x1*x2),
     x2 + 3*x1*x2**2*(4 + 3*x1*x2),
     2*x1 - 3*x1**2*x2]                            # F|_{x3=0}
F = [sp.expand(a[i]*x3 + b[i]) for i in range(3)]

JF = sp.Matrix(3, 3, lambda i, j: sp.diff(F[i], [x1, x2, x3][j]))
assert sp.expand(JF.det(method='berkowitz')) == -2, 'Alpoege F is not Keller with Jac -2'
print('[0] Alpoege F: F = a*x3 + b with a,b in C[x1,x2]^3, Jac F == -2   OK')

A = sp.expand(sum(a[i]*[y1, y2, y3][i] for i in range(3)))
B = sp.expand(sum(b[i]*[y1, y2, y3][i] for i in range(3)))
# ledger form:  A = y1 U^3 + 3 x1 y2 U^2 - x1^3 y3,  B as recorded in B0.2
assert sp.expand(A - (y1*U**3 + 3*x1*y2*U**2 - x1**3*y3)) == 0
assert sp.expand(B - (y1*x2**2*U*(4 + 3*x1*x2)
                      + y2*(x2 + 3*x1*x2**2*(4 + 3*x1*x2))
                      + y3*(2*x1 - 3*x1**2*x2))) == 0
Psi = sp.expand(B + sp.Rational(1, 2)*lam*A**2 + mu*A)
print('[0] Meng-Yang family psi_{lam,mu} = B + (lam/2)A^2 + mu A reconstructed   OK')

# ---------------------------------------------------------------------
# 1.  MASTER FORMULA  (A3).  Proof in two fully generic steps.
# ---------------------------------------------------------------------
# Step 1a.  For any 2x2 symmetric X, 2x3 Y, 3-vector c:
#     det [[X, Y],[Y^T, lam c c^T]] = lam * det[c ; Y_row1 ; Y_row2]^2,
# in particular it does NOT depend on X.  Fully generic => proof.
Xb = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'X{min(i,j)}{max(i,j)}'))
Yb = sp.Matrix(2, 3, lambda i, j: sp.Symbol(f'Y{i}{j}'))
cb = sp.Matrix([sp.Symbol(f'c{i}') for i in range(3)])
Hb = sp.Matrix(sp.BlockMatrix([[Xb, Yb], [Yb.T, lam*cb*cb.T]]))
db = sp.expand(Hb.det(method='berkowitz'))
for i in range(2):
    for j in range(2):
        assert sp.diff(db, Xb[i, j]) == 0
Mb = sp.Matrix([[cb[0], cb[1], cb[2]],
                [Yb[0, 0], Yb[0, 1], Yb[0, 2]],
                [Yb[1, 0], Yb[1, 1], Yb[1, 2]]])
assert sp.expand(db - lam*sp.expand(Mb.det())**2) == 0
print('[1a] block identity det[[X,Y],[Y^T,lam cc^T]] = lam det[c;Y1;Y2]^2, X-free '
      '(fully generic => PROOF)   OK')

# Step 1b.  Assemble: for Theta = <b,z> + (lam/2)<a,z>^2 + mu<a,z> + e(x),
# the (z,z)-block of Hess Theta is lam a a^T and the (x,z)-block is
#   Y_{p,j} = d_p b_j + tau d_p a_j + lam a_j <d_p a, z>,  tau = lam<a,z>+mu,
# so det[a ; Y_1 ; Y_2] = det[a ; v_1 ; v_2] = (Jac F)|_{x3=tau}.
# Verified end-to-end on the Meng-Yang data (and with a generic extra e(x)).
Z = [z1, z2, z3]


def master_check(av, bv, extra=0):
    """Uses [1a] (already a proof) to avoid recomputing 5x5 determinants:
    it suffices to check the Hessian block shape and the 3x3 identity
        det[a ; Y_1 ; Y_2] == (Jac F)|_{x3 = lam<a,z>+mu}."""
    S = sp.expand(sum(av[i]*Z[i] for i in range(3)))
    Th = sp.expand(sum(bv[i]*Z[i] for i in range(3))
                   + sp.Rational(1, 2)*lam*S**2 + mu*S + extra)
    V = [x1, x2, z1, z2, z3]
    H = sp.Matrix(5, 5, lambda i, j: sp.expand(sp.diff(Th, V[i], V[j])))
    # (z,z) block must be lam * a a^T
    for i in range(3):
        for j in range(3):
            if sp.expand(H[2+i, 2+j] - lam*av[i]*av[j]) != 0:
                return False
    Y = sp.Matrix(2, 3, lambda p, j: H[p, 2+j])
    Md = sp.Matrix([[av[0], av[1], av[2]],
                    [Y[0, 0], Y[0, 1], Y[0, 2]],
                    [Y[1, 0], Y[1, 1], Y[1, 2]]])
    Fv = [sp.expand(av[i]*x3 + bv[i]) for i in range(3)]
    Jv = sp.Matrix(3, 3, lambda i, j: sp.diff(Fv[i], [x1, x2, x3][j]))
    jac = sp.expand(Jv.det(method='berkowitz'))
    return sp.expand(sp.expand(Md.det()) - jac.subs(x3, lam*S + mu)) == 0


assert master_check(a, b)
print('[1b] MASTER FORMULA verified on the Meng-Yang/Alpoege data   OK')
# a second, unrelated instance (guards against an accidental coincidence)
a_t = [x1*x2 + 1, x2**2, x1 - 3]
b_t = [x1**3 - x2, 5*x1*x2**2, x2**3 + x1**2]
assert master_check(a_t, b_t)
assert master_check(a_t, b_t, extra=x1**4 - 7*x1*x2 + x2**5)
print('[1b] MASTER FORMULA verified on an independent instance, incl. an '
      'x-only additive term e(x)   OK')

# Consequence: det Hess psi_{lam,mu} = lam * (Jac F)^2 = lam * (-2)^2 = 4 lam,
# which is ledger entry B0.2 -- now with a structural reason.  (The direct
# 5x5 determinant is in verify_meng_yang_hc5.py / verify_mengyang_fast.py;
# recomputing it here would only duplicate that work.)
assert sp.expand(lam*(-2)**2 - 4*lam) == 0
print('[1c] master formula => det Hess psi_{lam,mu} = lam (Jac F)^2 = 4 lam '
      '(structural reproof of ledger B0.2)   OK')

# ---------------------------------------------------------------------
# 2.  The derivation u and the explicit non-affine automorphism tau  (A0/A1).
# ---------------------------------------------------------------------
def u(h):
    return sp.expand(3*x1*sp.diff(h, y1) - U*sp.diff(h, y2))


assert sp.expand(u(A)) == 0, 'u must annihilate the pivot coefficient A'
assert sp.expand(u(Psi) - (-x2*U)) == 0
assert sp.expand(u(u(Psi))) == 0
print('[2a] u := 3x1 d/dy1 - (1+x1x2) d/dy2 satisfies u(A)=0, u(psi) = -x2(1+x1x2), '
      'u^2(psi) = 0   OK')
# u is nilpotent as a derivation (u^2 = 0 on each coordinate) hence locally
# nilpotent; it is a COORDINATE derivation because the row (3x1, -U, 0) is
# unimodular over C[x1,x2] and we complete it EXPLICITLY:
assert sp.expand((-x2/3)*(3*x1) + (-1)*(-U)) == 1, 'unimodularity witness'
N = sp.Matrix([[-x2/3, -1, 0], [U, 3*x1, 0], [0, 0, 1]])
Ninv = sp.Matrix([[3*x1, 1, 0], [-U, -x2*sp.Rational(1, 3), 0], [0, 0, 1]])
assert sp.expand(N*Ninv - sp.eye(3)) == sp.zeros(3, 3)
assert sp.expand(Ninv*N - sp.eye(3)) == sp.zeros(3, 3)
assert sp.expand(N.det()) == 1
print('[2b] N in SL_3(C[x1,x2]) with polynomial inverse; rows of N are '
      'u-adapted   OK')

# tau : (x1,x2,z) -> (x1,x2, Ninv(x) z)   is a polynomial automorphism of C^5
zv = sp.Matrix([z1, z2, z3])
ysub = Ninv*zv
sub = {y1: ysub[0], y2: ysub[1], y3: ysub[2]}
# inverse of tau
yv = sp.Matrix([y1, y2, y3])
zback = N*yv
chk = [sp.expand(zback[i].subs(sub, simultaneous=True)) for i in range(3)]
assert chk == [z1, z2, z3], 'tau^{-1} o tau != id'
Jtau = sp.Matrix(5, 5, lambda i, j: sp.diff(
    [x1, x2, ysub[0], ysub[1], ysub[2]][i], [x1, x2, z1, z2, z3][j]))
assert sp.expand(Jtau.det()) == 1
print('[2c] tau is a polynomial automorphism of C^5 with det J tau = 1 '
      '(inverse verified explicitly); tau is NOT affine   OK')
assert sp.total_degree(sp.expand(ysub[0])) > 1
print('[2c] (tau has a component of degree 2: 3*x1*z1 + z2)   OK')

Pt = sp.expand(Psi.subs(sub, simultaneous=True))
assert sp.expand(sp.diff(Pt, z1, 2)) == 0, 'psi o tau must be affine-linear in z1'
piv = sp.expand(sp.diff(Pt, z1))
assert sp.expand(piv - (-x2*U)) == 0
print('[2d] *** psi_{lam,mu} o tau IS AFFINE-LINEAR IN z1 for EVERY lam, mu, ***')
print('     pivot coefficient  d(psi o tau)/dz1 = -x2(1+x1x2).')
print('     ==> the naive strengthening of D1 ("no polynomial change of')
print('         coordinates makes psi affine-linear in a variable") is FALSE.')

At = sp.expand(A.subs(sub, simultaneous=True))
assert sp.expand(At - (U**2*z2 - x1**3*z3)) == 0
print('[2e] A o tau = (1+x1x2)^2 z2 - x1^3 z3  (z1 has disappeared)   OK')

# ---------------------------------------------------------------------
# 3.  ... but tau destroys the Hessian determinant  (A2).
# ---------------------------------------------------------------------
w = [x1, x2, z2, z3]
Atil = piv
Btil = sp.expand(Pt.subs(z1, 0))
assert sp.expand(Pt - z1*Atil - Btil) == 0
gt = sp.Matrix([sp.diff(Atil, v) for v in w])
Kt = sp.Matrix(4, 4, lambda i, j: sp.diff(Atil, w[i], w[j]))
Ptm = sp.Matrix(4, 4, lambda i, j: sp.diff(Btil, w[i], w[j]))
Mm = sp.expand(Ptm + z1*Kt)


def mnr(Mx, rows, cols):
    return Mx.extract(rows, cols).det(method='berkowitz')


adj11 = mnr(Mm, [1, 2, 3], [1, 2, 3])
adj22 = mnr(Mm, [0, 2, 3], [0, 2, 3])
adj12 = -mnr(Mm, [0, 2, 3], [1, 2, 3])
dH = sp.expand(-(gt[0]**2*adj11 + 2*gt[0]*gt[1]*adj12 + gt[1]**2*adj22))
# cross-check against the master formula applied to (alpha', beta')
alp = [sp.expand(sum(Ninv[i, k]*a[i] for i in range(3))) for k in range(3)]
bet = [sp.expand(sum(Ninv[i, k]*b[i] for i in range(3))) for k in range(3)]
assert alp[0] == 0
S2 = sp.expand(sum(alp[i]*Z[i] for i in range(3)))
Fp = [sp.expand(alp[i]*x3 + bet[i]) for i in range(3)]
Jp = sp.Matrix(3, 3, lambda i, j: sp.diff(Fp[i], [x1, x2, x3][j]))
jacp = sp.expand(Jp.det(method='berkowitz'))
assert sp.expand(dH - lam*sp.expand(jacp.subs(x3, lam*S2 + mu))**2) == 0
print('[3a] det Hess(psi o tau) agrees with the master formula applied to '
      'Q = Ninv   OK')

R = sp.expand(jacp.subs(x3, lam*S2 + mu))
assert sp.expand(dH - lam*R**2) == 0
# R is NOT constant for any lam != 0 (it has the monomial 4*lam*x1^8*x2^3*z3):
cf = sp.Poly(R, x1, x2, z1, z2, z3).coeff_monomial(x1**8*x2**3*z3)
print('[3b] coefficient of x1^8 x2^3 z3 in R :', sp.simplify(cf))
assert sp.simplify(cf + 4*lam) == 0   # = -4 lam, nonzero for every lam != 0
print('[3b] det Hess(psi o tau) = lam * R^2 with R non-constant whenever '
      'lam != 0   ==> det Hess is NOT in C^*   OK')
print('     (equivalently: G = Ninv^T . F is NOT a Keller map of C^3.)')

# ---------------------------------------------------------------------
# 4.  FIBER-LINEAR REDUCTION (A4) and its necessary conditions.
# ---------------------------------------------------------------------
# For sigma(x,y) = (x, Q(x) y + r(x)), the master formula gives
#   det Hess(psi o sigma) = lam * ( Jac(Q^T F) )^2 |_{x3 = lam<Q^T a, z> + rho},
# and "affine-linear in z1" <=> (Q^T a)_1 = 0.  Since a is unimodular over
# C[x1,x2] (certified below) and Q is invertible, (Q^T a)_2, (Q^T a)_3 have
# no common zero, so x3 = lam<Q^T a,z>+rho sweeps C over every (x1,x2).
# Hence:  det Hess in C^*  <=>  Jac(Q^T F) is a nonzero constant.
assert sp.groebner([a[0], a[1], a[2]], x1, x2).exprs == [sp.Integer(1)], \
    'a must be unimodular'
print('[4a] a = (U^3, 3x1 U^2, -x1^3) is unimodular over C[x1,x2] '
      '(Groebner basis = {1})   OK')

# necessary conditions on G = (g, a2' x3 + b2', a3' x3 + b3'), g in C[x1,x2]
gg, A2s, A3s, B2s, B3s = sp.symbols('g A2 A3 B2 B3', cls=sp.Function)
gs = gg(x1, x2)
Gs = [gs, A2s(x1, x2)*x3 + B2s(x1, x2), A3s(x1, x2)*x3 + B3s(x1, x2)]
JG = sp.Matrix(3, 3, lambda i, j: sp.diff(Gs[i], [x1, x2, x3][j]))
dJG = sp.expand(JG.det())


def JJ(f, h):
    return sp.diff(f, x1)*sp.diff(h, x2) - sp.diff(f, x2)*sp.diff(h, x1)


pred = sp.expand(x3*(A3s(x1, x2)*JJ(gs, A2s(x1, x2)) - A2s(x1, x2)*JJ(gs, A3s(x1, x2)))
                 + (A3s(x1, x2)*JJ(gs, B2s(x1, x2)) - A2s(x1, x2)*JJ(gs, B3s(x1, x2))))
assert sp.expand(dJG - pred) == 0
print('[4b] Jac G = x3 (a3 J(g,a2) - a2 J(g,a3)) + (a3 J(g,b2) - a2 J(g,b3))  '
      '(symbolic, arbitrary functions => identity)   OK')
# W := a3 b2 - a2 b3 satisfies J(g,W) = h W + kappa when J(g,a_i) = h a_i
hs = sp.Function('h')(x1, x2)
subs_h = {}
W = A3s(x1, x2)*B2s(x1, x2) - A2s(x1, x2)*B3s(x1, x2)
lhs4 = sp.expand(JJ(gs, W))
rhs4 = sp.expand(B2s(x1, x2)*JJ(gs, A3s(x1, x2)) - B3s(x1, x2)*JJ(gs, A2s(x1, x2))
                 + (A3s(x1, x2)*JJ(gs, B2s(x1, x2)) - A2s(x1, x2)*JJ(gs, B3s(x1, x2))))
assert sp.expand(lhs4 - rhs4) == 0
print('[4c] J(g,W) = [b2 J(g,a3) - b3 J(g,a2)] + kappa ; with J(g,a_i) = h a_i '
      'this is J(g,W) = hW + kappa   OK')
print('[4d] at any critical point of g the Jacobian identity forces kappa = 0, '
      'so grad g vanishes nowhere; when h == 0, (g,W) is a planar KELLER pair.')

# The pivot coefficient g is forced to be <q,b> for q a SYZYGY of a:
q1 = sp.Matrix([3*x1, -U, 0])
q2 = sp.Matrix([0, x1**2, 3*U**2])
for q in (q1, q2):
    assert sp.expand(sum(q[i]*a[i] for i in range(3))) == 0
print('[4e] s1 = (3x1, -U, 0) and s2 = (0, x1^2, 3U^2) are syzygies of a   OK')
# and they GENERATE:  p U^3 + 3 q x1 U^2 - r x1^3 = 0  forces (gcd(U,x1)=1)
#   U^2 | r, write r = U^2 rho; then p U + 3 q x1 = rho x1^3 forces x1 | p,
#   write p = x1 pi; then q = (rho x1^2 - pi U)/3.  So (p,q,r) = (pi/3) s1 +
#   (rho/3) s2.  The two divisibility steps are the only inputs:
pi_, rho_ = sp.symbols('pi_ rho_')
gen = sp.Matrix([x1*pi_, (rho_*x1**2 - pi_*U)/3, U**2*rho_])
assert sp.expand(sum(gen[i]*a[i] for i in range(3))) == 0
assert sp.expand(gen - (pi_*q1/3 + rho_*q2/3)) == sp.zeros(3, 1)
print('[4f] the syzygy module of a is FREE on s1, s2 (parametrisation verified)   OK')
g1 = sp.expand(sum(q1[i]*b[i] for i in range(3)))
g2 = sp.expand(sum(q2[i]*b[i] for i in range(3)))
print('[4g] <s1,b> =', sp.factor(g1))
print('[4g] <s2,b> =', sp.factor(g2))
# g = <s1,b> = -x2(1+x1x2) is NOT a coordinate of C[x1,x2]: its fibre over 0
# is reducible.  (A coordinate has every fibre isomorphic to C, in particular
# irreducible.)  So by [4d] the h == 0 branch is closed for q = s1.
assert sp.factor(g1) == sp.factor(-x2*U)
assert len(sp.factor_list(g1)[1]) == 2
print('[4h] <s1,b> = -x2(1+x1x2) factors into 2 distinct irreducibles, so its '
      'zero fibre is REDUCIBLE and it is not a coordinate of C[x1,x2];')
print('     hence (Moh / JC_2 in degree <= 100) the h == 0 branch is impossible '
      'for q = s1.')

print()
print('=== NAX-A SUMMARY ===================================================')
print('PROVED : master formula det Hess_5 Theta = lam (Jac F)^2|_{x3=lam<a,z>+mu}')
print('PROVED : psi_{lam,mu} o tau is affine-linear in a coordinate for an')
print('         explicit NON-affine automorphism tau  ==>  the "affine only"')
print('         caveat in D1 is NOT removable by simply deleting the word.')
print('PROVED : that tau destroys constancy: det Hess(psi o tau) = lam R^2,')
print('         R non-constant.  So the corrected statement (pivot AND')
print('         det Hess in C^*) survives this attack.')
print('PROVED : fiber-linear reduction  <=>  Q^T F Keller with G_1 x3-free,')
print('         + the necessary conditions (i),(ii),(iii: grad g nowhere 0).')
print('OPEN   : the corrected theorem itself, even for fiber-linear sigma.')
print('=====================================================================')
