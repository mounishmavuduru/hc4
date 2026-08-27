# ============================================================================
# chkG_algebraic_audit.py -- INDEPENDENT ADVERSARIAL AUDIT of THEOREM G
# (agent key chkG).  Fail-closed, exact arithmetic only.
#
#   THEOREM G.  e in C[x_1..x_n],  B(e) := grad(e)^T adj(Hess e) grad(e).
#               B(e) == 0 identically  ==>  det Hess(e) == 0 identically.
#
# VERDICT: TRUE (SOUND).  Details of the audit of the original write-up:
#
# (1) "V := H^{-1} g is phi-related to the Euler field":  D(phi) V = H V = g
#     = phi, and the Euler field at phi(x) is phi(x).  Correct.
# (2) "E L - L = e o psi":  L = <psi,p> - e(psi).  d_j L = psi_j +
#     sum_i d_j psi_i (p_i - d_i e(psi)) = psi_j, so grad L = psi,
#     E L = <psi,p>, and E L - L = e o psi.  Correct.  Moreover
#     E(e o psi) = <g(psi), Dpsi p> = p^T H(psi)^{-1} p = (B/D)(psi) = 0, so
#     (E^2 - E) L = 0.  Correct.
# (3) "(E^2-E)L = 0 forces L = A + Btilde, homogeneous of degree 0, 1, on a
#     cone neighbourhood over C":  set Btilde := E L, A := L - E L on a ball
#     U0 around p0 (0 not in U0).  Then E A = 0 and E Btilde = Btilde on U0.
#     For p in U0 the set S_p := {s in C : s p in U0} is the preimage of a
#     ball under a real-linear embedding C -> C^n, hence convex, hence
#     connected; on S_p the holomorphic functions s -> A(s p) and
#     s -> Btilde(s p)/s have derivative (E A)(sp)/s = 0 and
#     ((E Btilde) - Btilde)(sp)/s^2 = 0, so A(sp) = A(p), Btilde(sp) = s Btilde(p)
#     for s in S_p.  This is exactly what makes L_ext(t p) := A(p) + t Btilde(p)
#     WELL DEFINED on the punctured cone Gamma := C^* U0 (two representations
#     t p = t' p' have p' = (t/t') p with t/t' in S_p).  Nothing is lost at
#     t = 0 because 0 is not in Gamma; branch points of psi never enter
#     because psi is never continued.  Correct (this is the paper's Step 5).
# (4) THE CONTINUATION STEP.  The original sketch (theorem_G.py header,
#     main_result.tex, lemma_ledger.md) says "this expression IS the analytic
#     continuation of psi along the whole cone ... the identity phi o psi = id
#     continues along the cone too".  As WRITTEN this is a non-sequitur: psi
#     is a local inverse, the ray may meet the branch locus {D = 0} of phi,
#     and a local inverse need not continue as an inverse.  BUT the step is
#     REPAIRABLE and IS repaired in theorem_G_paper.tex Sec. 2 (first proof):
#     psi_ext := grad L_ext is a globally defined holomorphic SECTION-CANDIDATE
#     on the connected open cone Gamma given by an explicit formula; the map
#     q -> grad e(psi_ext(q)) - q is holomorphic on Gamma and vanishes on the
#     open set U0, hence on all of Gamma by the identity theorem.  Only this
#     transported identity is used.  So grad e(t^{-1} a + b) = t p0 for all
#     t in C^*, a Laurent-polynomial identity whose t^1 coefficient gives
#     p0 = 0.  (A side remark: the identity theorem then forces
#     D(psi_ext(q)) != 0 on all of Gamma, so the section automatically avoids
#     the branch locus -- no hypothesis about the ray is needed.)
#     VERDICT on (4): exposition gap in the original sketch, no gap in the
#     repaired proof.  main_result.tex and lemma_ledger.md still carry the
#     original wording and should be aligned with the paper.
#
# (5) PURELY ALGEBRAIC PROOF (re-derived by hand in this audit; it is the
#     derivation form of the advG flow proof / atkG proof, recorded here with
#     every step spelled out, because it needs no analysis at all).
#     Assume B == 0 and D := det H != 0 in k[x], k any field of char 0.
#     Work in the field K := k(x).  V := H^{-1} g in K^n (H invertible over K
#     since D != 0), delta := sum_j V_j d/dx_j, a derivation of K.
#       (i)   delta g_i = sum_j V_j e_ij = (H V)_i = g_i.
#       (ii)  d_j V = d_j(H^{-1}) g + H^{-1} d_j g = -H^{-1}(d_j H) V + e_j,
#             and ((d_j H) V)_i = sum_k e_ikj V_k = M_ij with M symmetric
#             (e_ijk totally symmetric); so Jac(V) = I - H^{-1} M.
#       (iii) d_j <g,V> = <d_j g, V> + <g, d_j V> = (H V)_j + g_j - (V^T M)_j
#             = 2 g_j - (M V)_j.
#     Since <g,V> = g^T H^{-1} g = B/D = 0 in K, (iii) gives M V = 2 g, and
#     then (ii) gives delta V = Jac(V) V = V - 2 H^{-1} g = -V.
#     Put b := x + V, a := V in K^n.  Then delta b = V - V = 0 and
#     delta a = -a.  For any polynomial P in k[x], the finite Taylor formula
#     is an exact identity in K:  P(x) = P(b - a) = sum_alpha (d^alpha P)(b)
#     (-a)^alpha / alpha!.  Group by |alpha| = m:  P(x) = sum_m c_m with
#     delta c_m = -m c_m  (delta kills every polynomial in the b_j, and
#     delta(a^alpha) = -|alpha| a^alpha).  Apply this to P = g_i and use (i):
#         0 = delta g_i - g_i = sum_m (-m - 1) c_m,   i.e.  sum_m (m+1) c_m = 0.
#     LEMMA (eigenvectors of a derivation).  If u_0, ..., u_r in K satisfy
#     delta u_m = mu_m u_m with pairwise distinct mu_m in k, and
#     sum_m u_m = 0, then every u_m = 0.  [Take a relation with the fewest
#     nonzero terms; apply delta and subtract mu_{m0} times the relation:
#     sum_{m != m0} (mu_m - mu_{m0}) u_m = 0 is a shorter nontrivial
#     relation, contradiction.]  With u_m := (m+1) c_m, mu_m := -m, we get
#     (m+1) c_m = 0 for all m, so c_m = 0 (char 0), so g_i = sum_m c_m = 0.
#     Hence g == 0, H == 0, D == 0: contradiction.                        []
#     Polynomiality enters ONLY through the finiteness of the Taylor sum;
#     for the rational witness e = x2/x1 (B == 0, D != 0) one finds V = -x,
#     b = 0, a = -x, and g(b - a) = g(x) is homogeneous of degree -1 in a:
#     an eigenvalue +1 component, which a polynomial can never have.
#
# WHAT THIS SCRIPT CERTIFIES (own implementation, independent of advG/atkG)
#  (A) identities (i),(ii),(iii) and the consequence delta V + V =
#      H^{-1} grad(B/D) as EXACT rational-function identities on random
#      exact polynomials, n = 2 (deg 3,4,5), n = 3 (deg 3,4), n = 4 (deg 3);
#  (B) the same identities FULLY GENERICALLY: n = 2 as a complete polynomial
#      identity in x and the 9 coefficients of a generic cubic (so for every
#      e, at every point, since only the 3-jet enters); n = 3 for the
#      generic cubic (19 coefficients) at the point 0 (all 3-jets);
#  (C) the finite Taylor/eigen-decomposition: for generic e (n = 2 deg 5,
#      n = 3 deg 4, n = 4 deg 3) and symbolic a, b, g(b - a tau) is a
#      polynomial in tau of degree <= d-1: only eigenvalues 0,-1,..,-(d-1);
#  (D) non-vacuity / sharpness: on Euler-degree-0 rational e (n = 2, 3)
#      B == 0, D != 0, delta b == 0, delta a == -a, delta g == g all hold and
#      g(b - a tau) carries tau^{-1}; on a polynomial e with D != 0 the chain
#      breaks exactly at delta b = H^{-1} grad(B/D) != 0;
#  (E) the Legendre identities of the analytic proof on an exact example:
#      grad L = psi, E L - L = e o psi, (E^2 - E) L = (B/D) o psi != 0 when
#      B != 0 (so affineness of L along rays is exactly the B == 0 input);
#  (F) the convexity fact used in (3): for a ball U0 not containing 0 and
#      p in U0, {s in C : s p in U0} is a disc (explicitly computed).
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


def monomials(vs, k):
    return sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, k)}, key=str)


def rand_poly(xs, d, rnd, dens=0.7):
    mons = []
    for k in range(1, d + 1):
        mons += monomials(xs, k)
    return sum(sp.Rational(rnd.randint(-5, 5), rnd.randint(1, 3))*m for m in mons if rnd.random() < dens)


def derivation_identities(e, xs, rational=True):
    """Return dict of booleans for (i),(ii),(iii),(iv) on the polynomial e (exact)."""
    n = len(xs)
    g = sp.Matrix([sp.diff(e, v) for v in xs])
    H = sp.Matrix([[sp.diff(e, u, v) for v in xs] for u in xs])
    D = sp.expand(H.det(method='berkowitz'))
    if D == 0:
        return None
    A = H.adjugate()
    B = sp.expand((g.T*A*g)[0, 0])
    V = (A*g)/D                                  # H^{-1} g, exact rational
    delta = lambda f: sum(V[j]*sp.diff(f, xs[j]) for j in range(n))
    M = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = sum(sp.diff(e, xs[i], xs[j], xs[k])*V[k] for k in range(n))
    Hinv = A/D
    out = {}
    # (i) delta g = g
    out['i'] = all(sp.cancel(delta(g[i]) - g[i]) == 0 for i in range(n))
    # (ii) Jac(V) = I - H^{-1} M
    JacV = sp.Matrix([[sp.diff(V[i], xs[j]) for j in range(n)] for i in range(n)])
    R2 = JacV - (sp.eye(n) - Hinv*M)
    out['ii'] = all(sp.cancel(R2[i, j]) == 0 for i in range(n) for j in range(n))
    # (iii) grad(<g,V>) = 2 g - M V
    gV = (g.T*V)[0, 0]
    R3 = sp.Matrix([sp.diff(gV, v) for v in xs]) - (2*g - M*V)
    out['iii'] = all(sp.cancel(R3[i]) == 0 for i in range(n))
    # (iv) delta V + V = H^{-1} grad(B/D)   (so B == 0 forces delta V = -V)
    R4 = sp.Matrix([delta(V[i]) + V[i] for i in range(n)]) - Hinv*sp.Matrix([sp.diff(B/D, v) for v in xs])
    out['iv'] = all(sp.cancel(R4[i]) == 0 for i in range(n))
    out['B'] = B
    out['gV_is_B_over_D'] = sp.cancel(gV - B/D) == 0
    return out


# ---------------------------------------------------------------- (A)
print('(A) derivation identities (i)-(iv) as exact rational identities, random polynomials')
rnd = random.Random(20260825)
for n, degs in ((2, (3, 4, 5)), (3, (3, 4)), (4, (3,))):
    xs = sp.symbols(f'x1:{n+1}')
    for d in degs:
        t0 = time.time()
        r = None
        while r is None:
            r = derivation_identities(rand_poly(xs, d, rnd), xs)
        check(f'n={n} deg={d}: (i) delta g = g: {r["i"]}, (ii) Jac V = I - H^-1 M: {r["ii"]}, '
              f'(iii) grad<g,V> = 2g - MV: {r["iii"]}, (iv) delta V + V = H^-1 grad(B/D): {r["iv"]}, '
              f'<g,V> = B/D: {r["gV_is_B_over_D"]}  [{time.time()-t0:.1f}s]',
              r['i'] and r['ii'] and r['iii'] and r['iv'] and r['gV_is_B_over_D'])

# ---------------------------------------------------------------- (B)
print()
print('(B) the identities FULLY GENERICALLY')
t0 = time.time()
xs = sp.symbols('x1 x2')
mons = monomials(xs, 1) + monomials(xs, 2) + monomials(xs, 3)
cs = sp.symbols(f'g0:{len(mons)}')
e_gen = sum(c*m for c, m in zip(cs, mons))
r = derivation_identities(e_gen, xs)
check(f'n=2, generic cubic (9 symbolic coefficients), identities (i)-(iv) as polynomial identities '
      f'in x AND the coefficients => proof for all e, n=2  [{time.time()-t0:.1f}s]',
      r is not None and r['i'] and r['ii'] and r['iii'] and r['iv'])


def generic_jet_check_n3():
    """n = 3, generic cubic (19 coefficients): (ii) and (iii) at the point 0 in polynomial form.
    D^2 Jac(V) = D Jac(W) - W gradD^T,  so (ii) <=> D Jac(W) - W gradD^T = D^2 I - A Mt,
    Mt := sum_k (d_k H) W_k;  (iii) <=> D grad B - B grad D = 2 D^2 g - Mt W.
    Everything is evaluated at x = 0 after differentiating; only the 3-jet enters."""
    xs = sp.symbols('y1 y2 y3')
    n = 3
    mons = monomials(xs, 1) + monomials(xs, 2) + monomials(xs, 3)
    cs = sp.symbols(f'j0:{len(mons)}')
    dom = sp.QQ[cs]
    e = sum(c*m for c, m in zip(cs, mons))
    P = lambda ex: sp.Poly(ex, *xs, domain=dom)

    def trunc1(p):
        d = {m: c for m, c in p.as_dict().items() if sum(m) <= 1}
        return sp.Poly.from_dict(d, *xs, domain=dom) if d else P(0)

    def mul(p, q):
        return trunc1(p*q)

    g = [trunc1(P(sp.diff(e, v))) for v in xs]
    H = [[trunc1(P(sp.diff(e, u, v))) for v in xs] for u in xs]
    T3 = [[[P(sp.diff(e, u, v, w)) for w in xs] for v in xs] for u in xs]

    def det3(Mx):
        tot = P(0)
        for perm, sgn in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
                          ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
            t = mul(mul(Mx[0][perm[0]], Mx[1][perm[1]]), Mx[2][perm[2]])
            tot = tot + t if sgn == 1 else tot - t
        return tot

    def adj3(Mx):
        out = [[None]*3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                rows = [r for r in range(3) if r != j]
                cols = [c for c in range(3) if c != i]
                mnr = mul(Mx[rows[0]][cols[0]], Mx[rows[1]][cols[1]]) - mul(Mx[rows[0]][cols[1]], Mx[rows[1]][cols[0]])
                out[i][j] = mnr if (i + j) % 2 == 0 else -mnr
        return out

    D = det3(H)
    A = adj3(H)
    W = [sum((mul(A[i][j], g[j]) for j in range(n)), P(0)) for i in range(n)]
    B = sum((mul(g[i], W[i]) for i in range(n)), P(0))
    zero = {v: 0 for v in xs}
    val = lambda p: sp.expand(p.as_expr().subs(zero))
    der = lambda p, j: sp.expand(sp.diff(p.as_expr(), xs[j]).subs(zero))
    D0, B0 = val(D), val(B)
    g0 = sp.Matrix([val(x) for x in g])
    W0 = sp.Matrix([val(x) for x in W])
    A0 = sp.Matrix([[val(A[i][j]) for j in range(n)] for i in range(n)])
    JW = sp.Matrix([[der(W[i], j) for j in range(n)] for i in range(n)])
    gD = sp.Matrix([der(D, j) for j in range(n)])
    gB = sp.Matrix([der(B, j) for j in range(n)])
    Mt = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            Mt[i, j] = sum(val(T3[i][j][k])*W0[k] for k in range(n))
    lhs2 = D0*JW - W0*gD.T - (D0**2*sp.eye(n) - A0*Mt)
    lhs3 = D0*gB - B0*gD - (2*D0**2*g0 - Mt*W0)
    return all(sp.expand(t) == 0 for t in lhs2), all(sp.expand(t) == 0 for t in lhs3)


t0 = time.time()
ok2, ok3 = generic_jet_check_n3()
check(f'n=3, generic cubic (19 symbolic coefficients) at the point 0: (ii) {ok2}, (iii) {ok3} '
      f'=> proof for all e, n=3 (only the 3-jet enters)  [{time.time()-t0:.1f}s]', ok2 and ok3)

# ---------------------------------------------------------------- (C)
print()
print('(C) finite Taylor/eigen-decomposition: g(b - a tau) has tau-degree <= d-1 (generic e)')
tau = sp.Symbol('tau')
for n, dmax in ((2, 5), (3, 4), (4, 3)):
    vs = sp.symbols(f'u1:{n+1}')
    mons = []
    for d in range(1, dmax + 1):
        mons += monomials(vs, d)
    cs = sp.symbols(f'h{n}_0:{len(mons)}')
    e_g = sum(c*m for c, m in zip(cs, mons))
    av = sp.symbols(f'a1:{n+1}')
    bv = sp.symbols(f'b1:{n+1}')
    sub = {v: bv[i] - av[i]*tau for i, v in enumerate(vs)}
    degs = [sp.Poly(sp.expand(sp.diff(e_g, v).subs(sub)), tau).degree() for v in vs]
    check(f'n={n}, generic deg<={dmax} ({len(cs)} coefficients): tau-degrees {degs} <= {dmax-1}; '
          f'eigenvalues of delta on g are among 0,-1,..,-{dmax-1}, never +1', max(degs) <= dmax - 1)

# ---------------------------------------------------------------- (D)
print()
print('(D) non-vacuity and sharpness of the eigen-mechanism')


def chain_on(e, xs, x0, label):
    n = len(xs)
    g = sp.Matrix([sp.diff(e, v) for v in xs])
    H = sp.Matrix([[sp.diff(e, u, v) for v in xs] for u in xs])
    D = sp.simplify(H.det())
    assert D != 0
    A = H.adjugate()
    B = sp.simplify((g.T*A*g)[0, 0])
    V = sp.simplify(A*g/D)
    delta = lambda f: sum(V[j]*sp.diff(f, xs[j]) for j in range(n))
    b = sp.Matrix(xs) + V
    a = V
    db = [sp.simplify(delta(b[i])) for i in range(n)]
    da = [sp.simplify(delta(a[i]) + a[i]) for i in range(n)]
    dg = [sp.simplify(delta(g[i]) - g[i]) for i in range(n)]
    sub0 = dict(zip(xs, x0))
    a0 = sp.Matrix([sp.simplify(a[i].subs(sub0)) for i in range(n)])
    b0 = sp.Matrix([sp.simplify(b[i].subs(sub0)) for i in range(n)])
    p0 = sp.Matrix([sp.simplify(g[i].subs(sub0)) for i in range(n)])
    pt = dict(zip(xs, list(b0 - a0*tau)))
    gl = [sp.simplify(g[i].subs(pt)) for i in range(n)]
    return B, D, db, da, dg, gl, p0


# rational Euler-degree-0 witnesses: B == 0 by Euler, D != 0, chain holds, tau^{-1} appears
x1, x2, x3 = sp.symbols('x1 x2 x3')
for e, xs, x0, label in ((x2/x1, (x1, x2), (2, 3), 'e = x2/x1'),
                         ((x1**2 + x2*x3)/x3**2, (x1, x2, x3), (1, 2, 1), 'e = (x1^2 + x2 x3)/x3^2')):
    B, D, db, da, dg, gl, p0 = chain_on(e, xs, x0, label)
    check(f'{label}: B == 0, D = {D} != 0; delta b == 0: {all(t == 0 for t in db)}, '
          f'delta a == -a: {all(t == 0 for t in da)}, delta g == g: {all(t == 0 for t in dg)}',
          B == 0 and all(t == 0 for t in db) and all(t == 0 for t in da) and all(t == 0 for t in dg))
    check(f'{label}: g(b - a tau) == tau^(-1) p0 exactly (eigenvalue +1 -- impossible for polynomial g)',
          all(sp.simplify(gl[i] - p0[i]/tau) == 0 for i in range(len(xs))))

# a polynomial with D != 0: the chain must break exactly at delta b = H^{-1} grad(B/D) != 0
e_poly = x1**3 + x1*x2*x3 + x2**2 - x3**2 + 2*x1*x3 + x2
xs3 = (x1, x2, x3)
B, D, db, da, dg, gl, p0 = chain_on(e_poly, xs3, (1, 1, 1), 'poly')
check(f'polynomial e with D != 0: delta g == g holds ({all(t == 0 for t in dg)}), B != 0 ({B != 0}), '
      f'and delta b != 0 ({any(t != 0 for t in db)}) -- chain breaks exactly where B != 0',
      all(t == 0 for t in dg) and B != 0 and any(t != 0 for t in db))

# ---------------------------------------------------------------- (E)
print()
print('(E) Legendre identities of the analytic proof on an exact example')
p1, p2 = sp.symbols('p1 p2', positive=True)
e_ex = x1**3/3 + x2**2/2                     # grad = (x1^2, x2); psi = (sqrt p1, p2)
psi = (sp.sqrt(p1), p2)
L = sp.expand(psi[0]*p1 + psi[1]*p2 - e_ex.subs({x1: psi[0], x2: psi[1]}))
gradL = [sp.simplify(sp.diff(L, p) - s) for p, s in zip((p1, p2), psi)]
E = lambda f: p1*sp.diff(f, p1) + p2*sp.diff(f, p2)
e_psi = e_ex.subs({x1: psi[0], x2: psi[1]})
check('grad L == psi', all(t == 0 for t in gradL))
check('E L - L == e o psi', sp.simplify(E(L) - L - e_psi) == 0)
Hx = sp.Matrix([[sp.diff(e_ex, u, v) for v in (x1, x2)] for u in (x1, x2)])
gx = sp.Matrix([sp.diff(e_ex, v) for v in (x1, x2)])
BoverD = ((gx.T*Hx.adjugate()*gx)[0, 0]/Hx.det()).subs({x1: psi[0], x2: psi[1]})
check('(E^2 - E) L == (B/D) o psi  (= p1^(1/2)/2 != 0 here: affineness along rays IS the B == 0 input)',
      sp.simplify(E(E(L)) - E(L) - BoverD) == 0 and sp.simplify(BoverD) != 0)

# ---------------------------------------------------------------- (F)
print()
print('(F) convexity fact: {s : s p in Ball(p0, r)} is a disc when 0 not in the ball')
sr, si = sp.symbols('sr si', real=True)
p0v = sp.Matrix([3, 0])
pv = sp.Matrix([sp.Rational(5, 2), sp.Rational(1, 2)])       # a point of the ball of radius 1 about p0
s = sr + sp.I*si
dist2 = sp.expand(sum((sp.re(s*pv[i]) - p0v[i])**2 + sp.im(s*pv[i])**2 for i in range(2)))
# dist2 = |p|^2 |s - s*|^2 + const  =>  the set {dist2 < r^2} is a disc in the s-plane
pp = (pv.T*pv)[0, 0]
sstar = (pv.T*p0v)[0, 0]/pp                                   # real here since p0 is real
quad = sp.expand(pp*((sr - sstar)**2 + si**2) + (p0v.T*p0v)[0, 0] - sstar**2*pp)
check('|s p - p0|^2 == |p|^2 |s - s*|^2 + const (explicit): sublevel sets are discs, hence connected',
      sp.expand(dist2 - quad) == 0)

print()
print(f'total {time.time()-T0:.1f}s')
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print('chkG: every machine-checkable pillar of the algebraic proof of Theorem G PASSES.')
print('Verdict: Theorem G is TRUE; the original continuation wording is an exposition gap,')
print('repaired in theorem_G_paper.tex; the derivation proof above needs no analysis at all.')
sys.exit(0)
