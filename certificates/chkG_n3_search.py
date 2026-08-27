# ============================================================================
# chkG_n3_search.py -- refutation attempt on THEOREM G in THREE variables
# (agent key chkG).  Fail-closed, exact arithmetic only.
#
#   Looking for e in C[w1,w2,w3] with B(e) == 0 and det Hess e != 0
#   (Theorem G says there is none).  Degrees 3 and 4.
#
# METHOD (complete decision per slice, weight-by-weight elimination).
#  * Top form.  The top-degree piece of B(e) is B(e_d) = d/(d-1) e_d det Hess e_d
#    (Identity E, certified for ternary forms d = 3, 4), so B == 0 forces
#    det Hess e_d == 0; by Gordan-Noether (n = 3) e_d is then a cone: after a
#    linear change e_d in C[w1,w2].  (For d = 3 this is re-certified here:
#    the 3x3 minors of the coefficient matrix of (d_1 e_3, d_2 e_3, d_3 e_3)
#    lie in the radical of the ideal of coefficients of det Hess e_3.)
#  * Binary normal forms under GL_2(C) (complete lists):
#      d = 3:  w1^3 + w2^3  (3 distinct roots),  w1^2 w2,  w1^3;
#      d = 4:  w1 w2 (w1+w2)(w1+lam w2), lam not in {0,1} (modulus lam),
#              w1^2 w2 (w1+w2),  w1^2 w2^2,  w1^3 w2,  w1^4.
#    e_d = 0 is the lower-degree case.  Scaling e is harmless (B is a form).
#  * On the slice  e = F + (ALL monomials of degree 1..d-1)  (19 unknowns for
#    d = 4, 9 for d = 3) the homogeneous piece of B(e) of degree
#    m = 3d + k - 6 is  4 beta(F,F,F,e_k) + (terms in e_j, j > k)  -- B is a
#    quartic form in e when n = 3 -- so it is affine-linear in the
#    coefficients of e_k with a CONSTANT rational coefficient matrix.  Going
#    down in weight, each piece is solved exactly by row reduction over Q;
#    non-pivot coefficients stay free, incompatibility rows become polynomial
#    constraints; pieces without new unknowns are pure constraints.  Hence the
#    solution set of B == 0 on the slice is EXACTLY {free params : J = 0}.
#  * Then every coefficient of det Hess e (as a polynomial in the free
#    parameters) is certified to lie in rad(J).  Also checked: B(e_sol)
#    reduces to 0 mod J (consistency of the elimination), and J != (1)
#    (the slice contains the affinely-2-variable solutions).
#  For the p1111 family the modulus lam is SAMPLED (each value is a complete
#  decision for that lam); the lam-generic statement is not decided here.
#
# USAGE:  py -u chkG_n3_search.py deg3
#         py -u chkG_n3_search.py deg4:p4 | deg4:p31 | deg4:p22 | deg4:p211 | deg4:p1111
#         py -u chkG_n3_search.py sample
# ============================================================================

import itertools
import random
import sys
import time
import sympy as sp

T0 = time.time()
w1, w2, w3 = X = sp.symbols('w1 w2 w3')
zz = sp.Symbol('zz')
FAILS = []
MODE = sys.argv[1] if len(sys.argv) > 1 else 'deg3'


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
    return sp.expand(hess(e, vs).det(method='berkowitz'))


def monomials(vs, k):
    return sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, k)}, key=str)


def coeff_list(expr, vs=X):
    expr = sp.expand(expr)
    return [] if expr == 0 else [sp.expand(c) for c in sp.Poly(expr, *vs).coeffs()]


def in_radical(q, G, gens, ideal_gens, maxpow=6):
    q = sp.expand(q)
    if q == 0:
        return 'zero'
    acc = sp.Integer(1)
    for N in range(1, maxpow + 1):
        acc = sp.expand(acc*q)
        if G.reduce(acc)[1] == 0:
            return f'q^{N} in I'
    GR = sp.groebner(list(ideal_gens) + [1 - zz*q], *(list(gens) + [zz]), order='grevlex')
    return 'Rabinowitsch' if list(GR.exprs) == [sp.Integer(1)] else None


def B_poly(e, vs, unknowns):
    """B(e) as a Poly in vs over QQ[unknowns], via cofactors (n = 3)."""
    dom = sp.QQ[unknowns] if unknowns else sp.QQ
    P = lambda ex: sp.Poly(ex, *vs, domain=dom)
    n = len(vs)
    g = [P(sp.diff(e, v)) for v in vs]
    H = [[P(sp.diff(e, u, v)) for v in vs] for u in vs]
    A = [[None]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rows = [r for r in range(n) if r != j]
            cols = [c for c in range(n) if c != i]
            if n == 3:
                mnr = H[rows[0]][cols[0]]*H[rows[1]][cols[1]] - H[rows[0]][cols[1]]*H[rows[1]][cols[0]]
            else:
                mnr = H[rows[0]][cols[0]]
            A[i][j] = mnr if (i + j) % 2 == 0 else -mnr
    W = [sum((A[i][j]*g[j] for j in range(n)), P(0)) for i in range(n)]
    B = sum((g[i]*W[i] for i in range(n)), P(0))
    return B, dom


def sequential_decision(vs, F, d, label, verbose=True):
    n = len(vs)
    unknowns_by_k, allc = {}, []
    e = F
    for k in range(1, d):
        mons = monomials(vs, k)
        cs = list(sp.symbols(f'c{k}_0:{len(mons)}'))
        unknowns_by_k[k] = cs
        allc += cs
        e += sum(c*m for c, m in zip(cs, mons))
    t0 = time.time()
    BP, dom = B_poly(e, vs, allc)
    byweight = {}
    for mon, coef in BP.terms():
        byweight.setdefault(sum(mon), []).append(sp.expand(dom.to_sympy(coef)))
    top = (n + 1)*d - 2*n
    if verbose:
        print(f'   [{label}] {len(allc)} unknowns, {sum(len(v) for v in byweight.values())} equations, '
              f'weights {min(byweight)}..{max(byweight)} (top {top})  [{time.time()-t0:.1f}s]', flush=True)
    assert max(byweight) < top, 'top piece B(F) must vanish (F is a cone)'
    subs, constraints = {}, []
    for m in sorted(byweight, reverse=True):
        eqs = [sp.expand(c.xreplace(subs)) for c in byweight[m]]
        eqs = [q for q in eqs if q != 0]
        k = m - n*d + 2*n
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
        R, piv = M.row_join(sp.eye(r)).rref()
        Pm = R[:, len(new):]
        Rm = R[:, :len(new)]
        assert Rm - Pm*M == sp.zeros(r, len(new))
        Prhs = Pm*sp.Matrix(rhs)
        rank = 0
        for i in range(r):
            if all(Rm[i, j] == 0 for j in range(len(new))):
                q = sp.expand(Prhs[i])
                if q != 0:
                    constraints.append(q)
            else:
                rank += 1
                j0 = next(j for j in range(len(new)) if Rm[i, j] != 0)
                assert Rm[i, j0] == 1
                subs[new[j0]] = sp.expand(Prhs[i] - sum(Rm[i, j]*new[j] for j in range(len(new)) if j != j0))
        if verbose:
            print(f'      weight {m}: {len(eqs)} eqs, new unknowns e_{k} ({len(new)}), rank {rank}, '
                  f'free {len(new)-rank}, constraints so far {len(constraints)}', flush=True)
    for u, ex in subs.items():
        assert not (ex.free_symbols & set(subs)), 'unresolved substitution'
    free = [u for u in allc if u not in subs]
    e_sol = sp.expand(e.xreplace(subs))
    constraints = sorted({sp.expand(q) for q in constraints if sp.expand(q) != 0}, key=str)
    return e_sol, free, constraints


def decide_slice(F, d, label):
    t0 = time.time()
    e_sol, free, cons = sequential_decision(X, F, d, label)
    print(f'   [{label}] -> {len(free)} free parameters, {len(cons)} constraints  [{time.time()-t0:.1f}s]', flush=True)
    G = sp.groebner(cons, *free, order='grevlex') if cons else sp.groebner([sp.Integer(0)], *free, order='grevlex')
    print(f'   [{label}] Groebner basis of J: {len(G.exprs)} elements  [{time.time()-t0:.1f}s]', flush=True)
    check(f'[{label}] J != (1): the slice has solutions (affinely 2-variable ones)',
          list(G.exprs) != [sp.Integer(1)])
    Bsol = coeff_list(Bof(e_sol))
    check(f'[{label}] every coefficient of B(e_sol) reduces to 0 mod J (elimination consistent)  '
          f'[{time.time()-t0:.1f}s]', all(G.reduce(q)[1] == 0 for q in Bsol))
    Dc = coeff_list(dethess(e_sol))
    res = []
    for q in Dc:
        rr = in_radical(q, G, free, cons)
        res.append(rr)
        if rr is None:
            print(f'      *** det-Hess coefficient NOT in rad(J):  {q}')
    check(f'[{label}] all {len(Dc)} coefficients of det Hess e lie in rad(J) => det Hess == 0 on V(B)  '
          f'[{time.time()-t0:.1f}s]', all(rr is not None for rr in res))
    if any(rr is None for rr in res):
        # try to exhibit an explicit solution with det Hess != 0 (a refutation) -- diagnostic only
        print('      (diagnostic) Groebner basis of J:', list(G.exprs)[:10])


# ---------------------------------------------------------------- pillars
print(f'chkG n=3, mode = {MODE}')
print('(pillars) Identity E for ternary forms; cone lemma for ternary cubics')
for d in (3, 4):
    bs = sp.symbols(f'b{d}_0:{len(monomials(X, d))}')
    f = sum(c*m for c, m in zip(bs, monomials(X, d)))
    check(f'd={d}: B(f) == d/(d-1) f det Hess f (generic ternary form, {len(bs)} coefficients)',
          sp.expand(Bof(f) - sp.Rational(d, d-1)*f*dethess(f)) == 0)
    mons_low = []
    for k in range(1, d):
        mons_low += monomials(X, k)
    if d == 3:
        lows = sp.symbols(f'l{d}_0:{len(mons_low)}')
        e_full = f + sum(c*m for c, m in zip(lows, mons_low))
        BP = sp.Poly(Bof(e_full), *X)
        topdeg = 4*d - 6
        top = sum(coef*w1**m[0]*w2**m[1]*w3**m[2] for m, coef in BP.terms() if sum(m) == topdeg)
        check(f'd={d}: [B(e)]_top == B(e_d) for generic e with ALL lower monomials',
              sp.expand(top - Bof(f)) == 0 and BP.total_degree() <= topdeg)

if MODE in ('deg3', 'pillars'):
    t0 = time.time()
    bs = sp.symbols('b3_0:10')
    f = sum(c*m for c, m in zip(bs, monomials(X, 3)))
    Heqs = coeff_list(dethess(f))
    G = sp.groebner(Heqs, *bs, order='grevlex')
    print(f'      Groebner basis of the ternary-cubic zero-Hessian ideal: {len(G.exprs)} elements  '
          f'[{time.time()-t0:.1f}s]', flush=True)
    m2 = monomials(X, 2)
    Cm = sp.Matrix([[sp.Poly(sp.diff(f, v), *X).coeff_monomial(mm) for mm in m2] for v in X])
    minors = [sp.expand(Cm.extract([0, 1, 2], list(cols)).det()) for cols in itertools.combinations(range(6), 3)]
    res = [in_radical(q, G, bs, Heqs) for q in minors]
    check(f'ternary cubic: det Hess == 0 ==> partials linearly dependent ==> cone (all {len(minors)} '
          f'3x3 minors in radical)  [{time.time()-t0:.1f}s]', all(r is not None for r in res))

# ---------------------------------------------------------------- modes
if MODE == 'deg3':
    print()
    print('DEGREE 3: complete decision on the three binary-cubic normal forms')
    for F, label in ((w1**3 + w2**3, 'deg3, e_3 = w1^3 + w2^3'),
                     (w1**2*w2, 'deg3, e_3 = w1^2 w2'),
                     (w1**3, 'deg3, e_3 = w1^3')):
        decide_slice(F, 3, label)

elif MODE.startswith('deg4:'):
    typ = MODE.split(':')[1]
    forms = {
        'p4': [(w1**4, 'deg4, e_4 = w1^4')],
        'p31': [(w1**3*w2, 'deg4, e_4 = w1^3 w2')],
        'p22': [(w1**2*w2**2, 'deg4, e_4 = w1^2 w2^2')],
        'p211': [(w1**2*w2*(w1 + w2), 'deg4, e_4 = w1^2 w2 (w1+w2)')],
        'p1111': [(sp.expand(w1*w2*(w1 + w2)*(w1 + lam*w2)), f'deg4, e_4 = w1 w2 (w1+w2)(w1+{lam} w2)')
                  for lam in (2, -1, sp.Rational(1, 2), 3)],
    }[typ]
    print()
    for F, label in forms:
        decide_slice(F, 4, label)

elif MODE == 'sample':
    print()
    print('contrapositive sampling: exact e with det Hess != 0 must have B != 0')
    rnd = random.Random(20260825)
    viol, tested = 0, 0
    for d in (3, 4, 5):
        mons = []
        for k in range(1, d + 1):
            mons += monomials(X, k)
        for _ in range(25):
            e = sum(sp.Rational(rnd.randint(-5, 5), rnd.randint(1, 3))*m for m in mons if rnd.random() < 0.5)
            if dethess(e) == 0:
                continue
            tested += 1
            if Bof(e) == 0:
                viol += 1
                print('      *** COUNTEREXAMPLE CANDIDATE ***', e)
    check(f'{tested} exact instances (degrees 3..5) with det Hess != 0: violations = {viol}', viol == 0)
    # structured ansatz: e = a(w1,w2) + w3 b(w1,w2) + w3^2 c(w1,w2), random exact, deg <= 4
    viol2, tested2 = 0, 0
    for _ in range(40):
        a = sum(sp.Rational(rnd.randint(-4, 4), 1)*m for k in range(1, 5) for m in monomials((w1, w2), k) if rnd.random() < 0.5)
        b = sum(sp.Rational(rnd.randint(-4, 4), 1)*m for k in range(0, 4) for m in (monomials((w1, w2), k) if k else [sp.Integer(1)]) if rnd.random() < 0.5)
        c = sum(sp.Rational(rnd.randint(-4, 4), 1)*m for k in range(0, 3) for m in (monomials((w1, w2), k) if k else [sp.Integer(1)]) if rnd.random() < 0.5)
        e = sp.expand(a + w3*b + w3**2*c)
        if e == 0 or dethess(e) == 0:
            continue
        tested2 += 1
        if Bof(e) == 0:
            viol2 += 1
            print('      *** COUNTEREXAMPLE CANDIDATE ***', e)
    check(f'{tested2} structured instances a + w3 b + w3^2 c with det Hess != 0: violations = {viol2}', viol2 == 0)
    hdb = w1*w2 + w1**2*w3
    check(f'de Bondt-van den Essen h = w1 w2 + w1^2 w3: det Hess == 0, B = {Bof(hdb)} != 0',
          dethess(hdb) == 0 and Bof(hdb) != 0)
    ec = (2*w1 - w2 + 3*w3 + 1)**4 + (w1 + w3)**2*(2*w1 - w2 + 3*w3 + 1)
    check('control: affinely 2-variable quartic has B == 0 and det Hess == 0',
          Bof(ec) == 0 and dethess(ec) == 0)

print()
print(f'total {time.time()-T0:.1f}s')
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print(f'mode {MODE}: all checks PASS -- no counterexample to Theorem G in this slice.')
sys.exit(0)
