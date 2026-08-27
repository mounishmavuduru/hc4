# ============================================================================
# advG_flow_proof.py   -- ADVERSARIAL VERIFICATION of THEOREM G (agent key advG)
#
# Independent second-party audit of
#     THEOREM G.  e in C[x_1..x_n],  B(e) := grad(e)^T adj(Hess e) grad(e).
#                 B(e) == 0  ==>  det Hess(e) == 0.
#
# VERDICT OF THIS AUDIT: the theorem is TRUE.  The project's Legendre/analytic-
# continuation proof is essentially correct but its load-bearing paragraph
# ("this expression IS the analytic continuation of psi along the whole cone")
# is a non-sequitur as written.  This script certifies a REPLACEMENT PROOF that
# is completely elementary, uses only a LOCAL holomorphic ODE (equivalently a
# formal power series), and never continues anything analytically.
#
# ----------------------------------------------------------------------------
# NEW PROOF (advG).  Suppose B(e) == 0 and D := det Hess e =/= 0.
#
# Notation:  g = grad e, H = Hess e, A = adj(H), D = det H, W = A g
#   (polynomial), V = H^{-1} g = W/D  (rational, regular on Omega = {D =/= 0}),
#   M := sum_k (d_k H) V_k, i.e. M_ij = sum_k e_ijk V_k  (symmetric),
#   Mt := sum_k (d_k H) W_k = D * M   (polynomial).
#
# Step 0.  There is x0 with D(x0) =/= 0 and p0 := g(x0) =/= 0.
#          [If g vanished on the nonempty Zariski-open {D =/= 0} it would vanish
#           identically, forcing H == 0, hence D == 0.]
#
# Step 1 (J1).  From V = H^{-1}g,  d_j(H^{-1}) = -H^{-1}(d_jH)H^{-1}  and
#          d_j g = H e_j:
#              d_j V = e_j - H^{-1}(d_jH) V ,      i.e.   Jac(V) = I - H^{-1} M .
#          Polynomial form (P1):  D*Jac(W) - W*(grad D)^T = D^2 * I - A * Mt .
#
# Step 2 (J2).  d_j(g^T V) = (d_j g)^T V + g^T d_j V
#                          = (HV)_j + g_j - V^T M e_j = 2 g_j - (MV)_j , i.e.
#              grad( g^T V ) = 2 g - M V .
#          Polynomial form (P2):  D*grad(B) - B*grad(D) = 2 D^2 g - Mt * W ,
#          because g^T V = B/D.
#          (Both Step 1 and Step 2 are two lines of matrix calculus, valid in
#           every dimension; P1/P2 are machine-verified below for n = 2,3,4 with
#           a FULLY GENERIC 3-jet, and for n = 5,6 on random exact jets.)
#
# Step 3.  B == 0  =>  g^T V == 0 on Omega  =>  (J2)  M V = 2 g
#                   =>  (J1)  Jac(V) V = V - H^{-1} M V = V - 2 H^{-1} g = -V.
#
# Step 4 (local flow).  V is holomorphic near x0.  Let x(s) solve
#          x'(s) = V(x(s)),  x(0) = x0,  on a disc |s| < eps.  Then
#            (a) d/ds [ grad e (x(s)) ] = H x'(s) = H V = g(x(s)),
#                so   grad e (x(s)) = e^s p0 ;
#            (b) x''(s) = Jac(V) x'(s) = Jac(V) V = -V = -x'(s)  (Step 3),
#                so x'' + x' = 0 with x(0)=x0, x'(0)=V(x0):
#                   x(s) = a e^{-s} + b,   a = -V(x0),  b = x0 + V(x0).
#
# Step 5 (polynomiality kills it).  Put P(tau) := grad e(a tau + b) in C[tau]^n
#          and Q(tau) := tau P(tau) - p0 in C[tau]^n.  By 4(a),(b),
#            Q(e^{-s}) = e^{-s} grad e(a e^{-s} + b) - p0 = e^{-s} e^{s} p0 - p0 = 0
#          for all |s| < eps.  Since e^{-s} takes infinitely many values, Q == 0
#          as a polynomial.  But Q(0) = -p0 =/= 0.   CONTRADICTION.          []
#
# ----------------------------------------------------------------------------
# WHY THE PROOF IS SHARP.  Everything through Step 4 holds verbatim for RATIONAL
# e.  Only Step 5 uses that grad e is a polynomial.  And rational counter-
# examples DO exist: every e homogeneous of degree 0 has B(e) == 0 identically
# (Euler: H x = -g and <g,x> = 0, so B = x^T H adj(H) H x = det(H) x^T H x =
# -det(H)<g,x> = 0), while det Hess e can be nonzero -- e.g. e = x2/x1 in n = 2
# has det Hess = -1/x1^4.  For that e the proof runs: V = -x, x(s) = e^{-s}x0,
# grad e(x(s)) = e^s p0, and grad e(a tau) = tau^{-1} p0 -- a NEGATIVE power of
# tau, exactly what a polynomial gradient cannot produce.  This script verifies
# the whole chain on such witnesses: the proof is tight and no step is vacuous.
#
# Cross-check with the project's Identity E (Euler contraction):
#   e homogeneous of degree d >= 2  ==>  B = d/(d-1) * e * det Hess e.
#   The factor d/(d-1) is never 0 or 1 for d >= 2; it is 0 exactly at d = 0 --
#   precisely the rational witnesses above.
#
# ALL CHECKS ARE FAIL-CLOSED (assert / explicit exit) AND EXACT.
# ============================================================================

import itertools
import random
import sys
import sympy as sp

FAILS = []


def check(name, cond):
    print(f'   {"PASS" if cond else "*** FAIL ***"}  {name}')
    if not cond:
        FAILS.append(name)


# ---------------------------------------------------------------------------
# Truncated polynomial ring  R = K[x_1..x_n]/(x)^2  ("multi-dual numbers").
# An element is (c0, [c1_1,...,c1_n]) representing c0 + sum_k c1_k x_k.
# Every quantity in P1/P2 evaluated at x = 0 (value and first derivative) is
# computable in R, since e is taken to be a CUBIC and hence H is linear in x.
# ---------------------------------------------------------------------------
class Tr:
    __slots__ = ('c0', 'c1', 'n')

    def __init__(self, c0, c1):
        self.c0 = sp.expand(c0)
        self.c1 = [sp.expand(t) for t in c1]
        self.n = len(c1)

    @staticmethod
    def const(v, n):
        return Tr(sp.sympify(v), [sp.Integer(0)] * n)

    def __add__(self, o):
        return Tr(self.c0 + o.c0, [a + b for a, b in zip(self.c1, o.c1)])

    def __sub__(self, o):
        return Tr(self.c0 - o.c0, [a - b for a, b in zip(self.c1, o.c1)])

    def __neg__(self):
        return Tr(-self.c0, [-a for a in self.c1])

    def __mul__(self, o):
        return Tr(self.c0 * o.c0,
                  [self.c0 * b + a * o.c0 for a, b in zip(self.c1, o.c1)])

    def is_zero(self):
        return self.c0 == 0 and all(t == 0 for t in self.c1)


def tr_det(Mx, n):
    """Determinant by full permutation expansion, in the ring R."""
    from sympy.combinatorics import Permutation
    tot = Tr.const(0, n) if n else Tr.const(1, 0)
    idx = list(range(len(Mx)))
    for perm in itertools.permutations(idx):
        sgn = Permutation(list(perm)).signature()
        term = Tr.const(sgn, n)
        for i, j in enumerate(perm):
            term = term * Mx[i][j]
        tot = tot + term
    return tot


def tr_adj(Mx, n):
    """Adjugate: adj(M)_{ij} = (-1)^{i+j} * minor_{ji}."""
    m = len(Mx)
    out = [[None] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            sub = [[Mx[r][c] for c in range(m) if c != i]
                   for r in range(m) if r != j]
            d = tr_det(sub, n) if m > 1 else Tr.const(1, n)
            out[i][j] = d if (i + j) % 2 == 0 else -d
    return out


def jet_symbols(n, tag, randomize=False, seed=0):
    """Fully generic (or random exact) 3-jet: a_i, h_ij=h_ji, t_ijk symmetric."""
    rnd = random.Random(seed)

    def mk(name):
        if randomize:
            return sp.Integer(rnd.randint(-9, 9))
        return sp.Symbol(name)

    a = [mk(f'{tag}a{i}') for i in range(n)]
    h = {}
    for i in range(n):
        for j in range(i, n):
            h[(i, j)] = mk(f'{tag}h{i}{j}')
    t = {}
    for i in range(n):
        for j in range(i, n):
            for k in range(j, n):
                t[(i, j, k)] = mk(f'{tag}t{i}{j}{k}')

    def H(i, j):
        return h[(min(i, j), max(i, j))]

    def T(i, j, k):
        s = tuple(sorted((i, j, k)))
        return t[s]
    return a, H, T


def jet_check(n, tag, randomize=False, seed=0):
    """Verify P1 and P2 at x = 0 for a cubic with the given 3-jet."""
    a, Hf, Tf = jet_symbols(n, tag, randomize, seed)
    # H(x) = h + sum_k T[:,:,k] x_k     (exact for a cubic)
    Hx = [[Tr(Hf(i, j), [Tf(i, j, k) for k in range(n)]) for j in range(n)]
          for i in range(n)]
    # g(x) = a + h x + O(x^2)
    gx = [Tr(a[i], [Hf(i, k) for k in range(n)]) for i in range(n)]

    Dx = tr_det(Hx, n)
    Ax = tr_adj(Hx, n)
    Wx = [Tr.const(0, n) for _ in range(n)]
    for i in range(n):
        acc = Tr.const(0, n)
        for k in range(n):
            acc = acc + Ax[i][k] * gx[k]
        Wx[i] = acc
    Bx = Tr.const(0, n)
    for i in range(n):
        Bx = Bx + gx[i] * Wx[i]

    D0 = Dx.c0
    gradD = Dx.c1
    W0 = [w.c0 for w in Wx]
    JW = [[Wx[i].c1[j] for j in range(n)] for i in range(n)]   # d_j W_i at 0
    B0 = Bx.c0
    gradB = Bx.c1
    A0 = [[Ax[i][j].c0 for j in range(n)] for i in range(n)]
    g0 = [gi.c0 for gi in gx]
    Mt = [[sum(Tf(i, j, k) * W0[k] for k in range(n)) for j in range(n)]
          for i in range(n)]

    # P1: D*JW - W*(grad D)^T == D^2 I - A*Mt
    ok1 = True
    for i in range(n):
        for j in range(n):
            lhs = D0 * JW[i][j] - W0[i] * gradD[j]
            rhs = (D0**2 if i == j else 0) - sum(A0[i][k] * Mt[k][j]
                                                 for k in range(n))
            if sp.expand(lhs - rhs) != 0:
                ok1 = False
    # P2: D*grad B - B*grad D == 2 D^2 g - Mt*W
    ok2 = True
    for i in range(n):
        lhs = D0 * gradB[i] - B0 * gradD[i]
        rhs = 2 * D0**2 * g0[i] - sum(Mt[i][k] * W0[k] for k in range(n))
        if sp.expand(lhs - rhs) != 0:
            ok2 = False
    # non-vacuity: on the generic/random jet, B and D are NOT zero
    nonvac = (sp.expand(B0) != 0 and sp.expand(D0) != 0)
    return ok1, ok2, nonvac


# ===========================================================================
print('=' * 76)
print('PART 1.  The two jet identities P1, P2 -- the engine of the new proof.')
print('=' * 76)
print("""
Both sides of P1 and P2 are polynomial functions of the 3-JET (g,H,T) of e at
the evaluation point only (no derivative of e of order > 3 occurs).  Every
3-jet with H, T symmetric is realised at x = 0 by a cubic polynomial.  Hence
verifying P1, P2 for a FULLY GENERIC 3-jet is a COMPLETE PROOF of J1, J2 for
every polynomial (indeed every rational) e in that number of variables.
""")

for n in (2, 3, 4):
    ok1, ok2, nonvac = jet_check(n, f'g{n}_', randomize=False)
    check(f'n={n}: P1 on the FULLY GENERIC 3-jet  =>  J1 for all e', ok1)
    check(f'n={n}: P2 on the FULLY GENERIC 3-jet  =>  J2 for all e', ok2)
    check(f'n={n}: non-vacuous (generic B =/= 0 and det Hess =/= 0)', nonvac)

for n in (5, 6):
    allok = True
    for seed in (1, 2, 3):
        ok1, ok2, nonvac = jet_check(n, f'r{n}_', randomize=True, seed=seed)
        allok = allok and ok1 and ok2 and nonvac
    check(f'n={n}: P1 & P2 on 3 random exact 3-jets (consistency check)', allok)

print("""
Consequence chain (pure algebra):
  B == 0  =>  P2 gives  Mt*W = 2 D^2 g,  i.e.  M V = 2 g  (divide by D^2);
          =>  P1 gives  Jac(V) = I - H^{-1} M, so
              Jac(V) V = V - H^{-1} M V = V - 2 H^{-1} g = V - 2V = -V.
""")


# ===========================================================================
def hess(e, xs):
    n = len(xs)
    return sp.Matrix(n, n, lambda i, j: sp.diff(e, xs[i], xs[j]))


def grad(e, xs):
    return sp.Matrix([sp.diff(e, v) for v in xs])


def jacm(vec, xs):
    return sp.Matrix(len(vec), len(xs),
                     lambda i, j: sp.diff(vec[i], xs[j]))


def data(e, xs):
    n = len(xs)
    g = grad(e, xs)
    H = hess(e, xs)
    D = sp.expand(H.det(method='berkowitz'))
    A = H.adjugate()
    W = sp.expand(A * g)
    B = sp.expand((g.T * W)[0, 0])
    return g, H, D, A, W, B


print('=' * 76)
print('PART 2.  P1/P2 verified on explicit polynomials with B =/= 0 (the')
print('         correction terms are present and match exactly).')
print('=' * 76)

cases = [
    (3, lambda X: X[0]**4 + X[0]*X[1]*X[2] + X[1]**3 - 2*X[2]**2
                  + X[0] + 3*X[1]),
    (4, lambda X: X[0]**3 + X[1]**3 + X[2]**3 + X[3]**3
                  + X[0]*X[1]*X[2] - X[1]*X[3] + X[2]),
]
for n, build in cases:
    xs = list(sp.symbols(f'y1:{n+1}'))
    e = build(xs)
    g, H, D, A, W, B = data(e, xs)
    assert D != 0 and B != 0
    Mt = sp.Matrix(n, n, lambda i, j: sp.expand(
        sum(sp.diff(e, xs[i], xs[j], xs[k]) * W[k] for k in range(n))))
    P1 = sp.expand(D * jacm(W, xs) - W * grad(D, xs).T
                   - (D**2 * sp.eye(n) - A * Mt))
    P2 = sp.expand(D * grad(B, xs) - B * grad(D, xs)
                   - (2 * D**2 * g - Mt * W))
    check(f'n={n} explicit poly (B =/= 0, det Hess =/= 0): P1 at symbolic x',
          P1 == sp.zeros(n, n))
    check(f'n={n} same instance: P2 at symbolic x', P2 == sp.zeros(n, 1))

# ===========================================================================
print()
print('=' * 76)
print('PART 3.  Tightness: RATIONAL witnesses with B == 0 and det Hess =/= 0.')
print('         Steps 0-4 hold verbatim for them; only Step 5 (polynomiality)')
print('         fails -- so no step of the proof is vacuous or removable.')
print('=' * 76)

s = sp.Symbol('s')
tau = sp.Symbol('tau')
u1, u2 = sp.symbols('u1 u2')
v1, v2, v3 = sp.symbols('v1 v2 v3')
q1, q2, q3, q4 = sp.symbols('q1 q2 q3 q4')

witnesses = [
    ('n=2  e = x2/x1', [u1, u2], u2 / u1),
    ('n=3  e = (x2*x3)/x1**2', [v1, v2, v3], v2 * v3 / v1**2),
    ('n=3  e = x2/x1 + x3/x2', [v1, v2, v3], v2 / v1 + v3 / v2),
    ('n=4  e = (x2*x3*x4)/x1**3', [q1, q2, q3, q4], q2 * q3 * q4 / q1**3),
]

for name, xs, e in witnesses:
    n = len(xs)
    g = grad(e, xs)
    H = hess(e, xs)
    D = sp.cancel(sp.together(H.det(method='berkowitz')))
    A = H.adjugate()
    B = sp.cancel(sp.together((g.T * A * g)[0, 0]))
    check(f'{name}: B == 0', sp.simplify(B) == 0)
    check(f'{name}: det Hess =/= 0', sp.simplify(D) != 0)
    V = sp.Matrix([sp.cancel(sp.together(c)) for c in (H.inv() * g)])
    check(f'{name}: V = H^(-1) grad e == -x   (=> x(s) = e^(-s) x0)',
          sp.simplify(V + sp.Matrix(xs)) == sp.zeros(n, 1))
    JV = jacm(V, xs)
    check(f'{name}: Jac(V) V == -V   (Step 3 conclusion holds)',
          sp.simplify(JV * V + V) == sp.zeros(n, 1))
    sub = {xv: sp.exp(-s) * xv for xv in xs}
    check(f'{name}: grad e(x(s)) == e^s p0   (Step 4a holds)',
          sp.simplify(g.subs(sub, simultaneous=True) - sp.exp(s) * g)
          == sp.zeros(n, 1))
    tw = tau * g.subs({xv: tau * xv for xv in xs}, simultaneous=True)
    check(f'{name}: tau*grad e(tau x) is CONSTANT in tau -- possible only'
          ' because grad e is NOT a polynomial',
          sp.simplify(sp.Matrix([sp.diff(sp.cancel(sp.together(c)), tau)
                                 for c in tw])) == sp.zeros(n, 1))

print("""
   => Steps 0-4 are non-vacuous and are SATISFIED by genuine objects; the
      contradiction in Step 5 is produced by polynomiality alone.  Exactly the
      sharpness one wants from a rigidity theorem.
""")

# ===========================================================================
print('=' * 76)
print("PART 4.  Cross-checks: bordered Hessian, Identity E, converse failure.")
print('=' * 76)


def generic_poly(n, deg, tag, homog=False):
    xs = [sp.Symbol(f'{tag}_x{i}') for i in range(1, n + 1)]
    rng = [deg] if homog else range(1, deg + 1)
    mons = []
    for d in rng:
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(xs, d)},
                       key=str)
    cs = [sp.Symbol(f'{tag}_c{i}') for i in range(len(mons))]
    return xs, sp.expand(sum(c * m for c, m in zip(cs, mons)))


for n in (2, 3):
    xs, e = generic_poly(n, 3, f'z{n}')
    g, H, D, A, W, B = data(e, xs)
    Bord = sp.Matrix(sp.BlockMatrix([[H, g], [g.T, sp.zeros(1, 1)]]))
    check(f'n={n}: det[[H,g],[g^T,0]] == -B   (generic cubic, symbolic x)',
          sp.expand(Bord.det(method='berkowitz') + B) == 0)

for n, d in ((2, 3), (2, 4), (3, 2), (3, 3), (3, 4), (4, 3)):
    xs, e = generic_poly(n, d, f'm{n}{d}', homog=True)
    g, H, D, A, W, B = data(e, xs)
    check(f'Identity E, n={n} d={d}: B == d/(d-1)*e*det Hess e (generic form)',
          sp.expand(B - sp.Rational(d, d - 1) * e * D) == 0)

print("""
   d/(d-1) is never 0 or 1 for an integer d >= 2, and equals 0 exactly at d = 0:
   an independent confirmation that the obstruction in Theorem G is precisely
   "Euler-degree 0", which no non-constant polynomial can have.
""")

w1, w2, w3 = sp.symbols('w1 w2 w3')
h = w1 * w2 + w1**2 * w3
_, _, Dh, _, _, Bh = data(h, [w1, w2, w3])
check('converse fails: h = x1x2 + x1^2x3 has det Hess == 0 but B =/= 0',
      Dh == 0 and Bh != 0)

# ===========================================================================
print()
print('=' * 76)
if FAILS:
    print(f'*** {len(FAILS)} CHECK(S) FAILED ***')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('ALL CHECKS PASSED.')
print("""
CONCLUSION OF THIS CERTIFICATE (agent advG):
  Theorem G is PROVED by the elementary argument in this header.  It needs no
  analytic continuation, no Legendre transform, no local inverse psi, and no
  genericity beyond one point x0 with det Hess e(x0) =/= 0, grad e(x0) =/= 0.
  It works over any field of characteristic 0 if the local holomorphic flow is
  replaced by its formal power-series solution x(s) in K[[s]]: the same
  identities hold formally, and e^{-s} is transcendental over K, so
  Q(e^{-s}) = 0 in K[[s]] still forces Q == 0 in K[tau], whence Q(0) = -p0 = 0.
""")
