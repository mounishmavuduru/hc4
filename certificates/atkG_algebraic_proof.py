# ============================================================================
# atkG_algebraic_proof.py -- THIRD independent adversarial audit of THEOREM G
# (agent key atkG).  Fail-closed, exact arithmetic only.
#
#   THEOREM G.  e in C[x_1..x_n],  B(e) := grad(e)^T adj(Hess e) grad(e).
#               B(e) == 0 identically  ==>  det Hess(e) == 0 identically.
#
# VERDICT OF THIS AUDIT: TRUE.  The eigensplit proof of theorem_G_paper.tex
# Sec. 2 was re-derived line by line and found complete (see report).  The
# ORIGINAL sketch (theorem_G.py header, main_result.tex) is imprecise at the
# words "this expression IS the analytic continuation of psi": psi is only a
# local inverse and nothing continues psi AS AN INVERSE; what is continued is
# the Legendre transform L (by the explicit affine-in-t formula, holomorphic
# on the whole punctured cone), and the identity grad e o grad(L_ext) = id is
# then transported by the identity theorem.  That is exactly the paper's
# repair, so the defect is one of exposition, not of substance.
#
# ---------------------------------------------------------------------------
# A PURELY ALGEBRAIC PROOF (atkG) -- no analysis, no Puiseux series, no flow.
# Works over any field k of characteristic 0.
#
# Suppose B == 0 and D := det H != 0 in k[x].  Let K := k(x_1..x_n) and put
#     g := grad e,  H := Hess e,  A := adj H,  W := A g,  V := W/D = H^{-1} g,
#     delta := D_V = sum_j V_j d/dx_j        (a derivation of the field K).
# Matrix calculus (identities (J1),(J2) below, certified here for fully
# generic 3-jets) gives, with M_ij := sum_k e_ijk V_k:
#     (i)   delta g  = H V = g,
#     (ii)  Jac(V)   = I - H^{-1} M,          grad(g^T V) = 2 g - M V .
# Since g^T V = B/D == 0, (ii) yields M V = 2 g and hence
#     (iii) delta V  = Jac(V) V = V - 2 H^{-1} g = -V .
# Put  b := x + V  and  a := V.   Then  x = b - a,  delta b = V - V = 0,
# delta a = -a.   Each g_i is a POLYNOMIAL, so Taylor expansion is finite:
#     g_i(x) = g_i(b - a) = sum_{k=0}^{d-1} c_{i,k},
#     c_{i,k} := sum_{|alpha|=k} (-1)^k (d^alpha g_i)(b) a^alpha / alpha! ,
# and c_{i,k} is homogeneous of degree k in the a_j with delta-constant
# coefficients, so  delta c_{i,k} = -k c_{i,k}.   Now (i) says
#     sum_k (k+1) c_{i,k} = delta g_i - g_i + ... = 0 ,  i.e.
#     sum_k (k+1) c_{i,k} = 0 .
# Eigenvectors of a derivation for pairwise distinct eigenvalues are linearly
# independent over its constants (take a shortest nontrivial relation, apply
# delta, subtract lambda_1 times the relation: a shorter one), so every
# (k+1) c_{i,k} = 0, hence (char 0) every c_{i,k} = 0, hence g_i = 0.  Thus
# g == 0, H == 0, D == 0: contradiction.                                  []
#
# Where polynomiality enters: ONLY in the finiteness of the eigen-expansion
# of g(b - a).  For the rational witness e = x2/x1 (B == 0, det Hess != 0) one
# has V = -x, b = 0, a = -x and g(b - a) = g(x) is homogeneous of degree -1 in
# a: eigenvalue +1 -- exactly what a polynomial map can never produce.  For
# the algebraic sharpness witness F = x1 - sqrt(x1^2 - 2 x2) the same happens
# (checked below): g(b - a*tau) = tau^{-1} grad F(x0).
#
# This is the s-free skeleton of the advG flow proof (x(s) = b - a e^{-s}) and
# of the pzG Puiseux proof (chi = b + a t^{-1}); all three are the same
# mechanism.  It is recorded because it is the shortest and needs no auxiliary
# objects (no ODE, no Puiseux field, no cone chart).
#
# WHAT THIS SCRIPT CERTIFIES
#  (A) (J1) D*Jac(W) - W*(grad D)^T = D^2 I - A*Mt        (Mt := D*M)
#      (J2) D*grad(B) - B*grad(D)   = 2 D^2 g - Mt*W
#      (K)  (D*Jac(W) - W*(grad D)^T) W + D^2 W = A (D grad B - B grad D)
#           [ = D^3 (Jac(V) V + V) = D^3 H^{-1} grad(B/D) ]
#      at the jet level for FULLY GENERIC 3-jets, n = 2, 3 (default run) and
#      n = 4 (34 symbols; separate run `py atkG_algebraic_proof.py n4generic`,
#      recorded in atkG_algebraic_proof_n4generic.out.txt, ~6 min) -- since
#      every quantity at a point depends only on the 3-jet of e there, this is
#      a PROOF of (J1),(J2),(K) for all e in all dimensions <= 4; n = 4, 5, 6
#      also on random exact jets.  (Independent implementation from advG.)
#  (B) (K) as an exact polynomial identity on random rational polynomials,
#      n = 2 (deg 4, 5), n = 3 (deg 3, 4), n = 4 (deg 3); plus delta g = g and
#      delta(x + V) = H^{-1} grad(B/D) as exact rational identities.
#  (C) sharpness / non-vacuity: on the algebraic witness F, delta b = 0,
#      delta a = -a, delta g = g hold exactly and g(b - a tau) = tau^{-1} p0;
#      on the rational witness x2/x1 likewise; on a polynomial with
#      det Hess != 0, delta b != 0 (consistent: there B != 0).
#  (D) the eigen-comparison step as a generic theorem: for fully symbolic
#      a, b and fully generic e (n = 2 deg 5; n = 3 deg 4), g(b - a*tau) is a
#      polynomial in tau of degree <= d-1 -- only eigenvalues 0,-1,..,-(d-1)
#      occur, never +1.
# ============================================================================

import itertools
import random
import sys
import time
import sympy as sp

T0 = time.time()
FAILS = []


def check(name, cond):
    print(f'   {"PASS" if cond else "*** FAIL ***"}  {name}', flush=True)
    if not cond:
        FAILS.append(name)


# ---------------------------------------------------------------- helpers
def trunc(P, k=1):
    """Drop all monomials of total degree > k (work modulo (x)^(k+1))."""
    d = {m: c for m, c in P.as_dict().items() if sum(m) <= k}
    if not d:
        return sp.Poly(0, *P.gens, domain=P.domain)
    return sp.Poly.from_dict(d, *P.gens, domain=P.domain)


def tmul(a, b, k=1):
    return trunc(a * b, k)


def tdet(M, k=1):
    """Determinant by Laplace expansion along the first row, truncated."""
    m = len(M)
    if m == 1:
        return M[0][0]
    tot = None
    for j in range(m):
        sub = [row[:j] + row[j+1:] for row in M[1:]]
        term = tmul(M[0][j], tdet(sub, k), k)
        if j % 2:
            term = -term
        tot = term if tot is None else tot + term
    return tot


def tadj(M, k=1):
    m = len(M)
    out = [[None]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            sub = [[M[r][c] for c in range(m) if c != i]
                   for r in range(m) if r != j]
            d = tdet(sub, k) if m > 1 else sp.Poly(1, *M[0][0].gens, domain=M[0][0].domain)
            out[i][j] = d if (i + j) % 2 == 0 else -d
    return out


def jet_identities(n, randomize=False, seed=0):
    """Verify (J1),(J2),(K) at x = 0 for a cubic with generic/random 3-jet."""
    xs = sp.symbols(f'x1:{n+1}')
    rnd = random.Random(seed)
    names = []
    for i in range(n):
        names.append(f'a{i}')
    for i in range(n):
        for j in range(i, n):
            names.append(f'h{i}{j}')
    for i in range(n):
        for j in range(i, n):
            for k in range(j, n):
                names.append(f't{i}{j}{k}')
    if randomize:
        vals = {nm: sp.Integer(rnd.randint(-7, 7)) for nm in names}
        dom = sp.QQ
        sym = lambda nm: vals[nm]
    else:
        syms = sp.symbols(names)
        smap = dict(zip(names, syms))
        dom = sp.ZZ[tuple(syms)]
        sym = lambda nm: smap[nm]
    # generic cubic (no constant term; constants never enter)
    e = 0
    for i in range(n):
        e += sym(f'a{i}') * xs[i]
    for i in range(n):
        for j in range(i, n):
            e += sym(f'h{i}{j}') * xs[i] * xs[j]
    for i in range(n):
        for j in range(i, n):
            for k in range(j, n):
                e += sym(f't{i}{j}{k}') * xs[i] * xs[j] * xs[k]
    P = lambda ex: sp.Poly(ex, *xs, domain=dom)
    g = [trunc(P(sp.diff(e, v))) for v in xs]
    H = [[trunc(P(sp.diff(e, u, v))) for v in xs] for u in xs]
    # third derivatives (constants)
    Hk = [[[P(sp.diff(e, u, v, w)) for v in xs] for u in xs] for w in xs]
    A = tadj(H)
    D = tdet(H)
    W = []
    for i in range(n):
        acc = None
        for j in range(n):
            t = tmul(A[i][j], g[j])
            acc = t if acc is None else acc + t
        W.append(acc)
    Bp = None
    for i in range(n):
        t = tmul(g[i], W[i])
        Bp = t if Bp is None else Bp + t
    zero = {v: 0 for v in xs}
    ev = lambda p: p.as_expr().subs(zero)
    lin = lambda p, j: p.as_expr().coeff(xs[j], 1).subs(zero)   # d/dx_j at 0
    D0 = ev(D)
    gradD0 = sp.Matrix([lin(D, j) for j in range(n)])
    W0 = sp.Matrix([ev(w) for w in W])
    g0 = sp.Matrix([ev(gi) for gi in g])
    A0 = sp.Matrix([[ev(A[i][j]) for j in range(n)] for i in range(n)])
    JacW0 = sp.Matrix([[lin(W[i], j) for j in range(n)] for i in range(n)])
    B0 = ev(Bp)
    gradB0 = sp.Matrix([lin(Bp, j) for j in range(n)])
    Mt0 = sp.zeros(n, n)
    for k in range(n):
        Mt0 += sp.Matrix([[ev(Hk[k][i][j]) for j in range(n)] for i in range(n)]) * W0[k]
    I = sp.eye(n)
    J1 = (D0*JacW0 - W0*gradD0.T) - (D0**2*I - A0*Mt0)
    J2 = (D0*gradB0 - B0*gradD0) - (2*D0**2*g0 - Mt0*W0)
    K = (D0*JacW0 - W0*gradD0.T)*W0 + D0**2*W0 - A0*(D0*gradB0 - B0*gradD0)
    ok1 = all(sp.expand(t) == 0 for t in J1)
    ok2 = all(sp.expand(t) == 0 for t in J2)
    okK = all(sp.expand(t) == 0 for t in K)
    return ok1, ok2, okK


MODE = sys.argv[1] if len(sys.argv) > 1 else 'fast'
# 'fast'      : everything, with n = 4 on random exact jets (about 5 minutes)
# 'n4generic' : ONLY the n = 4 FULLY GENERIC 3-jet identities (about 6 minutes;
#               recorded in atkG_algebraic_proof_n4generic.out.txt).  Split so
#               that each run stays under the environment's ten-minute kill.
if MODE == 'n4generic':
    t0 = time.time()
    ok1, ok2, okK = jet_identities(4)
    check(f'n=4 fully generic 3-jet (34 symbols): (J1) {ok1}, (J2) {ok2}, (K) {okK}  '
          f'[{time.time()-t0:.1f}s]', ok1 and ok2 and okK)
    if FAILS:
        print('FAILED:', FAILS)
        sys.exit(1)
    print('n=4 generic jet identities PASS')
    sys.exit(0)

print('(A) jet-level identities (J1),(J2),(K)  [generic 3-jet => proof]')
for n in (2, 3):
    t0 = time.time()
    ok1, ok2, okK = jet_identities(n)
    check(f'n={n} fully generic 3-jet: (J1) {ok1}, (J2) {ok2}, (K) {okK}  '
          f'[{time.time()-t0:.1f}s]', ok1 and ok2 and okK)
for n, seed in ((4, 7), (4, 8), (5, 1), (6, 2)):
    t0 = time.time()
    ok1, ok2, okK = jet_identities(n, randomize=True, seed=seed)
    check(f'n={n} random exact 3-jet: (J1) {ok1}, (J2) {ok2}, (K) {okK}  '
          f'[{time.time()-t0:.1f}s]', ok1 and ok2 and okK)


# ---------------------------------------------------------------- (B)
print()
print('(B) (K) and the derivation identities on random exact polynomials')


def rand_poly(xs, d, rnd, lo=-5, hi=5):
    mons = []
    for k in range(1, d+1):
        mons += [sp.prod(c) for c in itertools.combinations_with_replacement(xs, k)]
    return sum(sp.Rational(rnd.randint(lo, hi), rnd.randint(1, 3))*m for m in mons)


def poly_K_check(e, xs):
    n = len(xs)
    dom = sp.QQ
    P = lambda ex: sp.Poly(ex, *xs, domain=dom)
    g = [P(sp.diff(e, v)) for v in xs]
    Hm = sp.Matrix([[sp.diff(e, u, v) for v in xs] for u in xs])
    A = Hm.adjugate()
    Dexpr = Hm.det(method='berkowitz')
    A = [[P(A[i, j]) for j in range(n)] for i in range(n)]
    D = P(Dexpr)
    if D.is_zero:
        return None
    W = [sum((A[i][j]*g[j] for j in range(n)), P(0)) for i in range(n)]
    B = sum((g[i]*W[i] for i in range(n)), P(0))
    gradD = [D.diff(v) for v in xs]
    gradB = [B.diff(v) for v in xs]
    JacW = [[W[i].diff(v) for v in xs] for i in range(n)]
    gDW = sum((gradD[j]*W[j] for j in range(n)), P(0))          # grad D . W
    ok = True
    for i in range(n):
        lhs = sum((D*JacW[i][j]*W[j] for j in range(n)), P(0)) - W[i]*gDW + D*D*W[i]
        rhs = sum((A[i][j]*(D*gradB[j] - B*gradD[j]) for j in range(n)), P(0))
        if not (lhs - rhs).is_zero:
            ok = False
    return ok, B.is_zero


rnd = random.Random(20260824)
for n, degs in ((2, (4, 5)), (3, (3, 4)), (4, (3,))):
    xs = sp.symbols(f'x1:{n+1}')
    for d in degs:
        t0 = time.time()
        for trial in range(2):
            e = rand_poly(xs, d, rnd)
            r = poly_K_check(e, xs)
            if r is None:
                continue
            ok, Bzero = r
            check(f'n={n} deg={d} trial {trial}: (K) exact polynomial identity '
                  f'(B==0: {Bzero})  [{time.time()-t0:.1f}s]', ok)

# rational-function form on a small instance: delta g = g, delta(x+V) = H^{-1} grad(B/D)
xs = sp.symbols('x1 x2 x3')
e = xs[0]**3 + xs[0]*xs[1]*xs[2] + xs[1]**2 - xs[2]**2 + 2*xs[0]*xs[2] + xs[1]
Hm = sp.Matrix([[sp.diff(e, u, v) for v in xs] for u in xs])
g = sp.Matrix([sp.diff(e, v) for v in xs])
D = Hm.det(method='berkowitz')
assert sp.expand(D) != 0
V = (Hm.adjugate()*g)/D
delta = lambda f: sum(V[j]*sp.diff(f, xs[j]) for j in range(3))
Bq = (g.T*Hm.adjugate()*g)[0, 0]
dg = sp.Matrix([delta(g[i]) for i in range(3)])
check('n=3 instance: delta g == g (exact rational)',
      all(sp.cancel(dg[i] - g[i]) == 0 for i in range(3)))
b = sp.Matrix(xs) + V
db = sp.Matrix([delta(b[i]) for i in range(3)])
rhs = Hm.adjugate()*sp.Matrix([sp.diff(Bq/D, v) for v in xs])/D
check('n=3 instance: delta(x+V) == H^{-1} grad(B/D) (exact rational)',
      all(sp.cancel(db[i] - rhs[i]) == 0 for i in range(3)))
check('n=3 instance: B != 0 and delta(x+V) != 0 (positive control)',
      sp.expand(Bq) != 0 and any(sp.cancel(db[i]) != 0 for i in range(3)))


# ---------------------------------------------------------------- (C)
print()
print('(C) sharpness witnesses: the eigen-mechanism on non-polynomial e')
x1, x2 = sp.symbols('x1 x2')
tau = sp.Symbol('tau')


def witness_chain(F, X, x0, label):
    """delta b = 0, delta a = -a, delta g = g and g(b - a tau) = tau^{-1} p0."""
    Hm = sp.Matrix([[sp.diff(F, u, v) for v in X] for u in X])
    g = sp.Matrix([sp.diff(F, v) for v in X])
    D = sp.simplify(Hm.det())
    assert sp.simplify(D) != 0
    Bw = sp.simplify((g.T*Hm.adjugate()*g)[0, 0])
    check(f'{label}: B == 0 and det Hess != 0', Bw == 0)
    V = sp.simplify(Hm.inv()*g)
    delta = lambda f: sum(V[j]*sp.diff(f, X[j]) for j in range(len(X)))
    bb = sp.Matrix(X) + V
    check(f'{label}: delta(x+V) == 0 exactly',
          all(sp.simplify(delta(bb[i])) == 0 for i in range(len(X))))
    check(f'{label}: delta V == -V exactly',
          all(sp.simplify(delta(V[i]) + V[i]) == 0 for i in range(len(X))))
    check(f'{label}: delta g == g exactly',
          all(sp.simplify(delta(g[i]) - g[i]) == 0 for i in range(len(X))))
    sub0 = dict(zip(X, x0))
    a0 = sp.Matrix([sp.nsimplify(sp.simplify(V[i].subs(sub0))) for i in range(len(X))])
    b0 = sp.Matrix([sp.nsimplify(sp.simplify(bb[i].subs(sub0))) for i in range(len(X))])
    p0 = sp.Matrix([sp.nsimplify(sp.simplify(g[i].subs(sub0))) for i in range(len(X))])
    pt = b0 - a0*tau
    return g, pt, p0, a0, b0


# rational witness e = x2/x1 (Euler degree 0)
g, pt, p0, a0, b0 = witness_chain(x2/x1, (x1, x2), (sp.Rational(2), sp.Rational(3)), 'e = x2/x1')
gl = sp.Matrix([sp.simplify(g[i].subs({x1: pt[0], x2: pt[1]})) for i in range(2)])
check('e = x2/x1: b = 0, a = -x0, and g(b - a tau) == tau^{-1} p0 (eigenvalue +1)',
      b0 == sp.zeros(2, 1) and a0 == -sp.Matrix([2, 3]) and
      all(sp.simplify(gl[i] - p0[i]/tau) == 0 for i in range(2)))

# algebraic witness F = x1 - sqrt(x1^2 - 2 x2)
R = sp.sqrt(x1**2 - 2*x2)
F = x1 - R
x0 = (sp.Rational(5, 2), sp.Integer(2))            # R(x0) = 3/2 rational
g, pt, p0, a0, b0 = witness_chain(F, (x1, x2), x0, 'F = x1 - sqrt(x1^2-2x2)')
rad = sp.expand(pt[0]**2 - 2*pt[1])                 # radicand along b - a tau
# polynomial square root via square-free factorisation (exact)
cst, facs = sp.Poly(rad, tau).sqf_list()
cst_sqrt = sp.sqrt(cst)
is_square = cst_sqrt.is_Rational and all(k % 2 == 0 for _, k in facs)
branch = None
if is_square:
    branch = cst_sqrt * sp.prod([f.as_expr()**(k//2) for f, k in facs])
    # choose the sign that agrees with R(x0) = 3/2 at tau = 1
    if sp.simplify(branch.subs(tau, 1) - sp.Rational(3, 2)) != 0:
        branch = -branch
check('F: radicand along b - a tau is a perfect square (branch is single-valued)',
      branch is not None and sp.expand(branch**2 - rad) == 0 and
      branch.subs(tau, 1) == sp.Rational(3, 2))
gF = sp.Matrix([1 - pt[0]/branch, 1/branch])
check('F: g(b - a tau) == tau^{-1} p0 exactly -- eigenvalue +1 component, impossible '
      'for polynomial g', all(sp.simplify(gF[i] - p0[i]/tau) == 0 for i in range(2)))
print(f'      (a = {list(a0)}, b = {list(b0)}, p0 = {list(p0)}, branch = {branch})')


# ---------------------------------------------------------------- (D)
print()
print('(D) eigen-comparison step, generic: g(b - a tau) has tau-degree <= d-1')
for n, dmax in ((2, 5), (3, 4)):
    vs = sp.symbols(f'y1:{n+1}')
    mons = []
    for d in range(1, dmax+1):
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'e{n}_0:{len(mons)}')
    e_gen = sum(c*m for c, m in zip(cs, mons))
    av = sp.symbols(f'a1:{n+1}')
    bv = sp.symbols(f'b1:{n+1}')
    sub = {v: bv[i] - av[i]*tau for i, v in enumerate(vs)}
    maxdeg = 0
    for v in vs:
        comp = sp.expand(sp.diff(e_gen, v).subs(sub))
        maxdeg = max(maxdeg, sp.Poly(comp, tau).degree())
    check(f'n={n}, deg<={dmax}: tau-degree of g(b - a tau) = {maxdeg} <= {dmax-1}; '
          f'no tau^{-1} term, i.e. eigenvalue +1 never occurs', maxdeg <= dmax-1)

print()
print(f'total {time.time()-T0:.1f}s')
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print('atkG: all pillars of the algebraic proof PASS. Theorem G verdict: TRUE.')
sys.exit(0)
