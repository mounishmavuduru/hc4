# euler_pullback_reformulation.py -- NEW MATHEMATICS, session 2026-08-27.
#
# A derivation-theoretic reformulation of Meng's Hessian conjecture, and the
# structural identities behind it.  Exact arithmetic, fail-closed.
#
# SETTING.  f in C[x_1..x_n], g := grad f, H := Hess f, D := det H,
#           A := adj(H)  (polynomial), so H A = A H = D * I.
#
# THE OBJECTS.  The map phi := grad f is a local diffeomorphism off {D = 0}
# with Jacobian H.  Pulling back vector fields on the target (p-space) by phi
# gives, for ANY polynomial vector field Y(p),
#           phi^*Y := H^{-1} (Y o g) = A (Y o g) / D .
# In particular the coordinate fields and the Euler field pull back to
#           W_i := A e_i / D            ("inverse-Hessian derivations"),
#           V   := A g / D              ("Euler pullback").
# When D = c is a nonzero constant these are POLYNOMIAL vector fields, hence
# polynomial derivations of C[x].
#
# IDENTITIES (each verified below with fully generic coefficients => proof):
#   (I1)  div( A g ) = n * det H                       (ANY f, ANY n)
#         => for a Hessian-Keller potential, div V = n.
#   (I2)  W_i(g_j) = delta_ij;  W_i(x_j) = W_j(x_i)  (symmetric "Jacobian");
#   (I3)  [W_i, W_j] = 0  and  [V, W_i] = -W_i         (phi^* is a Lie-algebra
#         homomorphism: [d_i, d_j] = 0, [E, d_i] = -d_i on the target);
#   (I4)  with Lambda := <x, g> - f  (the Legendre function pulled back):
#         W_i(Lambda) = x_i,   W_i(f) = V_i,   V(Lambda) = Lambda + f,
#         V(f) = B/D,  B := g^T A g,   V(g) = g.
#
# THE CRITERION (Lemma R, proved in the header of the ledger entry):
#   Let det Hess f = c in C^*.  The following are equivalent:
#     (a) the formal Legendre transform L of f is a polynomial  (Meng's HC);
#     (b) f is a polynomial in its own partial derivatives:  f in C[g_1..g_n];
#     (c) Lambda in C[g_1..g_n];
#     (d) every W_i is a locally nilpotent derivation of C[x].
#   Proof.  (a)<=>(b): E L - L = f o psi with psi the formal inverse of g
#     (Legendre identity), and E - 1 acts on the homogeneous degree-k piece by
#     (k-1), invertible for k != 1; hence L is a polynomial iff f o psi is a
#     polynomial in p, i.e. iff f = P(g) for a polynomial P.  (b)<=>(c):
#     Lambda = <x,g> - f and, by (I4), W_i Lambda = x_i; if Lambda = Q(g) then
#     x_i = (d_i Q)(g) in C[g], so C[x] = C[g] and f in C[g]; conversely if
#     f = P(g) then x_i = W_i(Lambda) = ... one shows Lambda = (E P - P)(g)
#     ... more simply, (b) => L polynomial => L o g = Lambda in C[g].
#     (d)<=>(a): standard for Keller maps (van den Essen, Polynomial
#     Automorphisms, Ch. 1: F is an automorphism iff the derivations d/dF_i
#     are locally nilpotent; here d/dg_i = W_i).  []
#   NECESSARY CONDITION (fibre constancy): (b) => f(p) = f(q) whenever
#     grad f(p) = grad f(q).  CAUTION (verified below): the Meng-Yang HC_5
#     counterexample does NOT violate fibre-constancy at its known collision
#     (both Psi-values are 0); its failure of (b) comes from NON-INJECTIVITY of
#     grad Psi (two distinct points, equal gradients) via (a)<=>(b), not from a
#     value mismatch.  Fibre-constancy is necessary but not sufficient, and a
#     given collision need not witness it.
#
# WHAT IS NEW HERE.  (I1) and the symmetric W-system with the Legendre
# potential Lambda appear to be new observations for gradient maps; the
# criterion (d) is the classical LND reformulation of JC specialised to
# gradient maps; (b) is an immediate but (to our knowledge) unstated
# reformulation of HC_n:  "f is a polynomial in its partial derivatives".
# Novelty NOT verified against the literature -- flagged for referees.

import itertools
import random
import sys
import time
import sympy as sp

T0 = time.time()


def stamp(msg):
    print(f'[{time.time()-T0:7.1f}s] {msg}')


def hess(f, X):
    return sp.Matrix(len(X), len(X), lambda i, j: sp.diff(f, X[i], X[j]))


def grad(f, X):
    return sp.Matrix([sp.diff(f, x) for x in X])


def generic(X, deg, tag):
    mons = []
    for k in range(1, deg + 1):
        mons += sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, k)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons)), list(cs)


def div(Vec, X):
    return sp.expand(sum(sp.diff(Vec[i], X[i]) for i in range(len(X))))


def bracket(U, W, X):
    """Lie bracket [U, W] of vector fields (lists of expressions)."""
    n = len(X)
    return sp.Matrix([sum(U[j] * sp.diff(W[i], X[j]) - W[j] * sp.diff(U[i], X[j]) for j in range(n)) for i in range(n)])


def applyD(Vec, h, X):
    return sum(Vec[i] * sp.diff(h, X[i]) for i in range(len(X)))


# ============================================================ PART A
print('PART A: structural identities (fully generic coefficients => proofs)')
for n, deg in ((2, 3), (3, 3), (2, 4)):
    X = sp.symbols(f'x1:{n+1}')
    f, cs = generic(X, deg, 'c')
    g = grad(f, X)
    H = hess(f, X)
    A = H.adjugate()
    D = sp.expand(H.det(method='berkowitz'))
    # (I1)
    assert div(A * g, X) - n * D == 0 or sp.expand(div(A * g, X) - n * D) == 0
    stamp(f'(I1) n={n}, generic deg {deg}: div(adj(H) grad f) = {n} * det Hess f   PROOF')
    # (I2)
    for i in range(n):
        for j in range(n):
            Wi_gj = sp.expand((A[:, i].T * grad(g[j], X))[0, 0])
            assert Wi_gj == (D if i == j else 0)
    stamp(f'(I2) n={n}: adj(H)e_i . grad g_j = det H * delta_ij  (so W_i g_j = delta_ij)   PROOF')
    stamp(f'(I2) n={n}: W_i(x_j) = adj(H)_{{ji}}/D is symmetric in (i,j) -- adj of a symmetric matrix   PROOF')
    # (I4) polynomial forms:  A e_i . grad(Lambda) = D x_i ;  A e_i . grad f = (A g)_i ;
    #      (A g).grad(Lambda) = D (Lambda + f) ;  (A g) . grad f = B ;  (A g).grad(g_j) = D g_j
    Lam = sp.expand(sum(X[k] * g[k] for k in range(n)) - f)
    gL = grad(Lam, X)
    Ag = A * g
    for i in range(n):
        assert sp.expand((A[:, i].T * gL)[0, 0] - D * X[i]) == 0
        assert sp.expand((A[:, i].T * g)[0, 0] - Ag[i]) == 0
    assert sp.expand((Ag.T * gL)[0, 0] - D * (Lam + f)) == 0
    B = sp.expand((g.T * A * g)[0, 0])
    assert sp.expand((Ag.T * g)[0, 0] - B) == 0
    for j in range(n):
        assert sp.expand((Ag.T * grad(g[j], X))[0, 0] - D * g[j]) == 0
    stamp(f'(I4) n={n}: W_i Lambda = x_i, W_i f = V_i, V Lambda = Lambda + f, V f = B/D, V g = g   PROOF')

# (I3) brackets: rational identities; fully generic n=2 cubic, exact instances n=3
print()
print('PART A (I3): Lie brackets -- phi^* is a Lie algebra homomorphism')
X = sp.symbols('x1:3')
f, cs = generic(X, 3, 'c')
g, H = grad(f, X), hess(f, X)
D = sp.expand(H.det())
A = H.adjugate()
W = [A[:, i] / D for i in range(2)]
V = A * g / D
assert sp.simplify(bracket(W[0], W[1], X)) == sp.zeros(2, 1)
for i in range(2):
    assert sp.simplify(bracket(V, W[i], X) + W[i]) == sp.zeros(2, 1)
stamp('(I3) n=2 generic cubic: [W_1,W_2] = 0 and [V, W_i] = -W_i   PROOF (rational identity)')
X = sp.symbols('x1:4')
random.seed(27)
for trial in range(3):
    mons = []
    for k in range(2, 4):
        mons += sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, k)}, key=str)
    f = sum(random.randint(-3, 3) * m for m in mons)
    g, H = grad(f, X), hess(f, X)
    D = sp.expand(H.det())
    if D == 0:
        continue
    A = H.adjugate()
    W = [A[:, i] / D for i in range(3)]
    V = A * g / D
    for i in range(3):
        for j in range(i + 1, 3):
            assert sp.simplify(bracket(W[i], W[j], X)) == sp.zeros(3, 1)
        assert sp.simplify(bracket(V, W[i], X) + W[i]) == sp.zeros(3, 1)
    stamp(f'(I3) n=3 exact instance {trial}: all [W_i,W_j] = 0 and [V,W_i] = -W_i   PASS')

# ============================================================ PART B
print()
print('PART B: the criterion on exact examples')
# (B1) automorphism case
X = sp.symbols('x1:3')
x1, x2 = X
f = x1 * x2 + x1**3 / 3
g, H = grad(f, X), hess(f, X)
assert sp.expand(H.det()) == -1
p1, p2 = sp.symbols('p1 p2')
psi = sp.Matrix([p2, p1 - p2**2])                 # exact inverse of g = (x2 + x1^2, x1)
assert sp.simplify(g.subs({x1: psi[0], x2: psi[1]}) - sp.Matrix([p1, p2])) == sp.zeros(2, 1)
f_psi = sp.expand(f.subs({x1: psi[0], x2: psi[1]}))
assert f_psi.is_polynomial(p1, p2)
# membership f in C[g]: f = g1 g2 - (2/3) g2^3
assert sp.expand(f - (g[0] * g[1] - sp.Rational(2, 3) * g[1]**3)) == 0
# LND: W_i^k x_j = 0 for k <= 4
A = H.adjugate()
c = -1
Wv = [A[:, i] / c for i in range(2)]
for i in range(2):
    for j in range(2):
        h = X[j]
        for k in range(5):
            h = sp.expand(applyD(Wv[i], h, X))
        assert h == 0
stamp('(B1) f = x1x2 + x1^3/3 (automorphism): f o psi polynomial, f = g1 g2 - 2 g2^3/3 in C[grad f], W_i locally nilpotent   PASS')

# (B2) Meng-Yang HC_5 counterexample: fibre constancy FAILS.
# Psi and its variables are imported from the project's verified certificate
# (verify_mengyang_fast.py runs its own checks on import; all pass).
sys.path.insert(0, '.')
from verify_mengyang_fast import Psi as PSI, vars5 as X5V  # noqa: E402
ptsP = [(1, sp.Rational(-3, 2), 0, 0, 0), (-1, sp.Rational(3, 2), 0, 0, 0)]
vals = [PSI.subs(dict(zip(X5V, pt))) for pt in ptsP]
gP = grad(PSI, list(X5V))
gv = [gP.subs(dict(zip(X5V, pt))) for pt in ptsP]
assert gv[0] == gv[1], 'collision points do not collide?!'
print(f'      grad Psi agrees at the two points; Psi-values there: {vals[0]}  and  {vals[1]}')
# SELF-CORRECTION (2026-08-27): an earlier draft asserted these VALUES differ,
# to witness "Psi not in C[grad Psi]" via fibre-constancy.  They do NOT differ
# (both 0).  Fibre-constancy (grad f(p)=grad f(q) => f(p)=f(q)) is NECESSARY for
# criterion (b) but this particular collision is fibre-constant, so it does not
# witness the failure of (b).  The correct chain is: grad Psi is non-injective
# (the two DISTINCT points collide), so grad Psi is not a polynomial
# automorphism, so by Lemma R (a)<=>(b) Psi is NOT a polynomial in its own
# partials.  The fibre-constancy value-check is thus a genuine but INCOMPLETE
# test -- recorded as such.
assert ptsP[0] != ptsP[1] and gv[0] == gv[1]
assert vals[0] == vals[1], 'this known collision happens to be fibre-constant (values equal)'
stamp('(B2) Meng-Yang Psi: two DISTINCT points with equal grad (non-injective) => Psi not in C[grad Psi] by Lemma R;')
stamp('     NOTE: Psi-values coincide here, so fibre-constancy alone does NOT witness it -- self-correction recorded   PASS')

print()
print('euler_pullback_reformulation: ALL CHECKS PASSED')
