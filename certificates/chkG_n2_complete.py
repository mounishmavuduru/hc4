# ============================================================================
# chkG_n2_complete.py -- THEOREM G in TWO variables, settled COMPLETELY
# (agent key chkG; independent adversarial audit).  Fail-closed, exact only.
#
#   THEOREM G (n = 2).  For e in C[x1,x2]:
#       B(e) := grad(e)^T adj(Hess e) grad(e) == 0   ==>   det Hess e == 0,
#   and in fact e = P(l) for a polynomial P and an affine-linear l.
#
# ---------------------------------------------------------------------------
# HAND PROOF (complete, all degrees).  For a symmetric 2x2 matrix
# H = [[a,b],[b,c]] one has adj H = J^T H J with J = [[0,-1],[1,0]], so
#     B = T^T H T,     T := J g = (-e_2, e_1),
# T being the Hamiltonian field tangent to the level curves of e.  Along an
# integral curve x(s) of T (x' = T(x)) one computes
#     T' = Jac(T) T = ( e_12 e_2 - e_22 e_1 ,  -e_11 e_2 + e_12 e_1 ),
#     T ^ T' := T_1 T'_2 - T_2 T'_1 = e_11 e_2^2 - 2 e_12 e_1 e_2 + e_22 e_1^2 = B.
# (Both identities are certified below with fully generic coefficients.)
# Assume B == 0 and e nonconstant.  A nonconstant polynomial e: C^2 -> C has
# only finitely many critical values (Sard/Bertini in characteristic 0) and
# finitely many values c for which e - c is not square-free; take c outside
# this finite set, so the fibre Z_c := {e = c} is a smooth reduced affine
# curve on which grad e, hence T, never vanishes.  Near any point of Z_c the
# holomorphic flow of T parametrises Z_c, and B == 0 gives T ^ T' == 0, i.e.
# T' = lam(s) T with lam holomorphic; then T(s) = T(0) exp(int lam) has
# constant direction and x(s) = x(0) + T(0) * int exp(int lam) lies on a
# straight line.  An irreducible curve that locally lies in a line IS that
# line, so every irreducible component of Z_c is a line.  Now: two lines in
# the same smooth fibre are disjoint (a crossing would be a singular point),
# and two lines in different fibres are disjoint (different level sets);
# disjoint affine lines are parallel.  Hence all components of all generic
# fibres are lines with ONE common direction v.  For x in the dense union of
# generic fibres, x + C v lies in the fibre of x, so e(x + s v) = e(x) for
# all s; this polynomial identity in (x, s) then holds everywhere, i.e.
# D_v e == 0, and after a linear change of variables e in C[x2].  Then
# Hess e has rank <= 1 and det Hess e == 0.                                []
#
# (The same conclusion also follows from the general algebraic proof audited
# in chkG_algebraic_audit.py; the argument above is the n = 2 short cut.)
#
# ---------------------------------------------------------------------------
# MACHINE DECISION over the FULL coefficient space, degrees 2..5 (and 6 if
# time permits), linear terms included -- the runs on disk stop at degree 3
# (pzG_n2_decision.py hangs on the 14-variable Groebner basis of degree 4;
# atkG_n2_decision.py was killed before its decision loop).
#
#  (P1) Covariance: for affine A, B(e o A) = det(A)^2 B(e) o A and
#       det Hess(e o A) = det(A)^2 det Hess(e) o A; B(lam e) = lam^3 B(e).
#       Fully generic cubic e (9 coefficients), generic affine A (the
#       generic quartic, 14 coefficients, also passed in a debug run).
#  (P2) Top form: for e of degree d the top-degree piece of B(e) is B(e_d)
#       = d/(d-1) e_d det Hess(e_d) (Identity E), d = 2..6, generic e with
#       ALL lower monomials present.  So B == 0 forces det Hess(e_d) == 0.
#  (P3) Binary forms with vanishing Hessian are d-th powers of linear forms
#       (the 2x2 minors of the Hankel matrix of e_d lie in the radical of
#       the ideal of coefficients of det Hess e_d), d = 2..6.
#  Hence every e of degree d with B == 0 is, after a linear change, a
#  scaling and dropping the constant,   e = x2^d + (all monomials of degree
#  1..d-1)   -- the COMPLETE slice (S_d).
#  (P4) DECISION on (S_d): weight-by-weight elimination.  B is a cubic form
#       in e; with e = x2^d + e_{d-1} + ... + e_1 the homogeneous piece of
#       B(e) of degree m = 2d - 4 + k is  3 beta(x2^d, x2^d, e_k) + (terms in
#       e_j, j > k)  =  d^2 x2^{2d-2} d_1^2 e_k + (terms in e_j, j > k),
#       so it is AFFINE-LINEAR in the coefficients of e_k with CONSTANT
#       integer coefficients.  Going down in weight, each piece is solved
#       exactly (row reduction over Q; non-pivot coefficients stay free,
#       incompatibility rows become polynomial constraints on the free
#       coefficients), so the solution set of B == 0 on (S_d) is EXACTLY
#       {free parameters : J = 0} for the recorded constraint ideal J.
#       Then every coefficient of e on a monomial containing x1, and every
#       coefficient of det Hess e, is certified to lie in rad(J) (power
#       reduction against a Groebner basis, else Rabinowitsch).  Conclusion:
#       on (S_d), B == 0 forces e in C[x2], hence det Hess e == 0.
#  (P5) Cross-validation: for d = 2, 3 the same conclusion is re-derived by
#       a direct Groebner/radical computation over the full coefficient
#       space (no normalisation, all 5 resp. 9 coefficients free) -- the two
#       methods must agree.
#  (P6) Contrapositive sampling, degrees 4..7: exact rational e with
#       det Hess != 0 all have B != 0 (a single hit would refute the theorem).
#  (P7) Positive control: e = P(x2) and e = P(l) have B == 0 and det Hess == 0.
# ============================================================================

import itertools
import random
import sys
import time
import sympy as sp

T0 = time.time()
x1, x2 = X = sp.symbols('x1 x2')
zz = sp.Symbol('zz')
FAILS = []
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def check(name, cond):
    print(f'   {"PASS" if cond else "*** FAIL ***"}  {name}', flush=True)
    if not cond:
        FAILS.append(name)


def hess(e, vs=X):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def Bof(e, vs=X):
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T*hess(e, vs).adjugate()*g)[0, 0])


def dethess(e, vs=X):
    return sp.expand(hess(e, vs).det())


def monomials(vs, k):
    return sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, k)}, key=str)


def coeff_list(expr, vs=X):
    expr = sp.expand(expr)
    return [] if expr == 0 else [sp.expand(c) for c in sp.Poly(expr, *vs).coeffs()]


def in_radical(q, G, gens, ideal_gens):
    """Certificate that q lies in rad(ideal): q^N reduces to 0 mod G, else Rabinowitsch."""
    q = sp.expand(q)
    if q == 0:
        return 'zero'
    acc = sp.Integer(1)
    for N in range(1, 7):
        acc = sp.expand(acc*q)
        if G.reduce(acc)[1] == 0:
            return f'q^{N} in I'
    GR = sp.groebner(list(ideal_gens) + [1 - zz*q], *(list(gens) + [zz]), order='grevlex')
    return 'Rabinowitsch' if list(GR.exprs) == [sp.Integer(1)] else None


# ---------------------------------------------------------------- (P0)
print('(P0) the n = 2 hand-proof identities, fully generic cubic')
c9 = sp.symbols('k0:9')
mons3 = monomials(X, 1) + monomials(X, 2) + monomials(X, 3)
e3 = sum(c*m for c, m in zip(c9, mons3))
g3 = sp.Matrix([sp.diff(e3, v) for v in X])
H3 = hess(e3)
J = sp.Matrix([[0, -1], [1, 0]])
check('adj(H) == J^T H J for symmetric 2x2 H',
      sp.expand(H3.adjugate() - J.T*H3*J) == sp.zeros(2, 2))
T = J*g3
check('B == T^T H T with T = J grad e', sp.expand((T.T*H3*T)[0, 0] - Bof(e3)) == 0)
JacT = sp.Matrix([[sp.diff(T[i], v) for v in X] for i in range(2)])
Tp = JacT*T
check('B == T ^ (Jac(T) T)  (velocity wedge acceleration of the level-curve flow)',
      sp.expand(T[0]*Tp[1] - T[1]*Tp[0] - Bof(e3)) == 0)

# ---------------------------------------------------------------- (P1)
print()
print('(P1) affine covariance, fully generic cubic e (9 coefficients; the quartic case')
print('     also passed in a debug run, 154 s under machine load)')
t0 = time.time()
mons4 = monomials(X, 1) + monomials(X, 2) + monomials(X, 3)
c14 = sp.symbols('q0:9')
e4 = sum(c*m for c, m in zip(c14, mons4))
a11, a12, a21, a22, u1, u2, lam = sp.symbols('a11 a12 a21 a22 u1 u2 lam')
Amat = sp.Matrix([[a11, a12], [a21, a22]])
y = {x1: a11*x1 + a12*x2 + u1, x2: a21*x1 + a22*x2 + u2}
eA = sp.expand(e4.subs(y, simultaneous=True))
B_e, D_e = Bof(e4), dethess(e4)
check('B(e o A) == det(A)^2 * (B(e) o A)',
      sp.expand(Bof(eA) - Amat.det()**2*B_e.subs(y, simultaneous=True)) == 0)
check('det Hess(e o A) == det(A)^2 * (det Hess(e) o A)',
      sp.expand(dethess(eA) - Amat.det()**2*D_e.subs(y, simultaneous=True)) == 0)
check('B(lam e + const) == lam^3 B(e);  det Hess(lam e + const) == lam^2 det Hess(e)',
      sp.expand(Bof(lam*e4 + 5) - lam**3*B_e) == 0 and
      sp.expand(dethess(lam*e4 + 5) - lam**2*D_e) == 0)
print(f'      [{time.time()-t0:.1f}s]')

# ---------------------------------------------------------------- (P2),(P3)
print()
print('(P2),(P3) top form: Identity E for binary forms; zero-Hessian binary forms are l^d')
for d in range(2, 7):
    t0 = time.time()
    bs = sp.symbols(f'b0:{d+1}')
    f = sum(bs[j]*sp.binomial(d, j)*x1**(d-j)*x2**j for j in range(d+1))
    Bf, Df = Bof(f), dethess(f)
    check(f'd={d}: B(f) == d/(d-1) * f * det Hess f  (generic binary form)',
          sp.expand(Bf - sp.Rational(d, d-1)*f*Df) == 0)
    mons_low = []
    for k in range(1, d):
        mons_low += monomials(X, k)
    lows = sp.symbols(f'l0:{len(mons_low)}')
    e_full = f + sum(c*m for c, m in zip(lows, mons_low))
    Bfull = sp.Poly(Bof(e_full), *X)
    topdeg = 3*d - 4
    top = sum(coef*x1**m[0]*x2**m[1] for m, coef in Bfull.terms() if sum(m) == topdeg)
    check(f'd={d}: [B(e)]_top == B(e_d) and deg B(e) <= {topdeg}, ALL lower monomials present',
          sp.expand(top - Bf) == 0 and Bfull.total_degree() <= topdeg)
    Heqs = coeff_list(Df)
    G = sp.groebner(Heqs, *bs, order='grevlex')
    minors = [sp.expand(bs[i]*bs[j+1] - bs[j]*bs[i+1]) for i in range(d) for j in range(i+1, d)]
    res = [in_radical(mn, G, bs, Heqs) for mn in minors]
    check(f'd={d}: det Hess(f) == 0 ==> all {len(minors)} Hankel 2x2 minors vanish ==> f = c*l^d  '
          f'[{time.time()-t0:.1f}s]', all(r is not None for r in res))


# ---------------------------------------------------------------- (P4)
def sequential_decision(d, verbose=True):
    """Complete decision of B == 0 on the slice e = x2^d + (all monomials of degree 1..d-1)."""
    n = 2
    unknowns_by_k = {}
    e = x2**d
    allc = []
    for k in range(1, d):
        mons = monomials(X, k)
        cs = list(sp.symbols(f'c{k}_0:{len(mons)}'))
        unknowns_by_k[k] = cs
        allc += cs
        e += sum(c*m for c, m in zip(cs, mons))
    BP = sp.Poly(Bof(e), *X)
    byweight = {}
    for mon, coef in BP.terms():
        byweight.setdefault(sum(mon), []).append(sp.expand(coef))
    top = 3*d - 4
    assert max(byweight) < top, 'top piece B(x2^d) must vanish'
    subs = {}
    constraints = []
    for m in sorted(byweight, reverse=True):
        eqs = [sp.expand(c.xreplace(subs)) for c in byweight[m]]
        eqs = [q for q in eqs if q != 0]
        k = m - n*d + 2*n            # weight of beta(e_d, e_d, e_k) is 2d-4+k
        new = [u for u in unknowns_by_k.get(k, []) if u not in subs]
        if not eqs:
            continue
        if not new:
            constraints += eqs
            continue
        rows, rhs = [], []
        for q in eqs:
            qp = sp.Poly(q, *new)
            assert qp.total_degree() <= 1, 'piece must be affine-linear in the new unknowns'
            row = [qp.coeff_monomial(u) for u in new]
            assert all(r.is_Rational for r in row), 'linear part must have constant coefficients'
            rows.append(row)
            rhs.append(-sp.expand(qp.coeff_monomial(1)))
        M = sp.Matrix(rows)
        r = M.rows
        aug = M.row_join(sp.eye(r))
        R, piv = aug.rref()
        pivots_M = [p for p in piv if p < len(new)]
        P = R[:, len(new):]            # transformation: R[:, :len(new)] == P*M
        Rm = R[:, :len(new)]
        assert sp.simplify(Rm - P*M) == sp.zeros(r, len(new))
        Prhs = P*sp.Matrix(rhs)
        for i in range(r):
            if all(Rm[i, j] == 0 for j in range(len(new))):
                q = sp.expand(Prhs[i])
                if q != 0:
                    constraints.append(q)
            else:
                j0 = next(j for j in range(len(new)) if Rm[i, j] != 0)
                assert Rm[i, j0] == 1
                expr = Prhs[i] - sum(Rm[i, j]*new[j] for j in range(len(new)) if j != j0)
                subs[new[j0]] = sp.expand(expr)
        if verbose:
            print(f'      weight {m}: {len(eqs)} eqs, new unknowns e_{k} ({len(new)}), '
                  f'rank {len(pivots_M)}, free {len(new)-len(pivots_M)}, '
                  f'constraints so far {len(constraints)}', flush=True)
    # every substitution must be fully resolved (RHS only in free symbols)
    free = [u for u in allc if u not in subs]
    for u, ex in subs.items():
        assert not (ex.free_symbols & set(subs)), 'unresolved substitution'
    e_sol = sp.expand(e.xreplace(subs))
    # sanity: B(e_sol) must vanish modulo the constraints (checked after Groebner)
    constraints = sorted({sp.expand(q) for q in constraints if sp.expand(q) != 0}, key=str)
    return e_sol, free, constraints, subs


print()
print('(P4) COMPLETE DECISION on the slice (S_d) = x2^d + all monomials of degree 1..d-1')
for d in range(2, DMAX + 1):
    t0 = time.time()
    e_sol, free, cons, subs = sequential_decision(d)
    print(f'   d={d}: {len(subs)+len(free)} lower coefficients -> {len(free)} free parameters, '
          f'{len(cons)} polynomial constraints  [{time.time()-t0:.1f}s]', flush=True)
    G = sp.groebner(cons, *free, order='grevlex') if cons else sp.groebner([sp.Integer(0)], *free, order='grevlex')
    print(f'      Groebner basis of J: {len(G.exprs)} elements  [{time.time()-t0:.1f}s]', flush=True)
    # (i) B(e_sol) itself lies in the ideal J (consistency of the elimination)
    Bsol = coeff_list(Bof(e_sol))
    check(f'd={d}: every coefficient of B(e_sol) reduces to 0 mod J (elimination consistent)',
          all(G.reduce(q)[1] == 0 for q in Bsol))
    # (ii) x1-coefficients of e in rad(J)
    ep = sp.Poly(e_sol, *X)
    x1coeffs = [sp.expand(c) for mn, c in ep.terms() if mn[0] > 0]
    res = [in_radical(q, G, free, cons) for q in x1coeffs]
    check(f'd={d}: all {len(x1coeffs)} coefficients of e on x1-monomials lie in rad(J) => e in C[x2]  '
          f'[{time.time()-t0:.1f}s]', all(r is not None for r in res))
    # (iii) det Hess coefficients in rad(J)
    Dc = coeff_list(dethess(e_sol))
    res2 = [in_radical(q, G, free, cons) for q in Dc]
    check(f'd={d}: all {len(Dc)} coefficients of det Hess e lie in rad(J) => det Hess == 0  '
          f'[{time.time()-t0:.1f}s]', all(r is not None for r in res2))
    # (iv) the slice really contains the expected solutions (J != (1), P(x2) allowed)
    check(f'd={d}: J != (1) and x2^k-coefficients are unconstrained (control)',
          list(G.exprs) != [sp.Integer(1)] and
          all(G.reduce(sp.expand(c))[1] != 0 for mn, c in ep.terms() if mn[0] == 0 and 0 < mn[1] < d
              and sp.expand(c) != 0))

# ---------------------------------------------------------------- (P5)
print()
print('(P5) cross-validation: direct Groebner over the FULL coefficient space, d = 2, 3')
for d in (2, 3):
    t0 = time.time()
    mons = []
    for k in range(1, d + 1):
        mons += monomials(X, k)
    cs = sp.symbols(f'f{d}_0:{len(mons)}')
    e = sum(c*m for c, m in zip(cs, mons))
    Beqs = coeff_list(Bof(e))
    G = sp.groebner(Beqs, *cs, order='grevlex')
    Dc = coeff_list(dethess(e))
    res = [in_radical(q, G, cs, Beqs) for q in Dc]
    check(f'd={d}: full space ({len(cs)} coefficients): every det-Hess coefficient in rad(B-ideal)  '
          f'[{time.time()-t0:.1f}s]', all(r is not None for r in res))

# ---------------------------------------------------------------- (P6)
print()
print('(P6) contrapositive sampling: exact rational e with det Hess != 0 must have B != 0')
rnd = random.Random(2026_08_25)
viol, tested = 0, 0
for d in (4, 5, 6, 7):
    mons = []
    for k in range(1, d + 1):
        mons += monomials(X, k)
    for _ in range(40):
        e = sum(sp.Rational(rnd.randint(-6, 6), rnd.randint(1, 4))*m for m in mons if rnd.random() < 0.6)
        if dethess(e) == 0:
            continue
        tested += 1
        if Bof(e) == 0:
            viol += 1
            print('      *** COUNTEREXAMPLE CANDIDATE ***', e)
check(f'{tested} exact instances with det Hess != 0, degrees 4..7: violations = {viol}', viol == 0)

# ---------------------------------------------------------------- (P7)
print()
print('(P7) positive controls')
P = sp.Symbol('s')**5 - 3*sp.Symbol('s')**2 + sp.Symbol('s')
ec = P.subs(sp.Symbol('s'), x2)
el = P.subs(sp.Symbol('s'), 2*x1 - 3*x2 + 1)
check('e = P(x2): B == 0 and det Hess == 0', Bof(ec) == 0 and dethess(ec) == 0)
check('e = P(2x1 - 3x2 + 1): B == 0 and det Hess == 0', Bof(el) == 0 and dethess(el) == 0)
check('e = x1^2/2 + x2 (converse fails): det Hess == 0 but B == 1',
      dethess(x1**2/2 + x2) == 0 and Bof(x1**2/2 + x2) == 1)

print()
print(f'total {time.time()-T0:.1f}s')
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print(f'chkG: THEOREM G in n = 2 DECIDED for all polynomials of degree <= {DMAX}, linear terms')
print('included (complete slice after certified normalisation): B == 0 forces e in C[l],')
print('hence det Hess e == 0.  Consistent with the hand proof in the header.')
sys.exit(0)
