# ============================================================================
# atkG_n3_search.py -- computational REFUTATION ATTEMPTS on THEOREM G in
# THREE variables (agent key atkG).  Exact, fail-closed.
#
#   Theorem G:  B(e) := grad(e)^T adj(Hess e) grad(e) == 0  ==>  det Hess e == 0.
#
# A counterexample is an e with B(e) == 0 and det Hess e != 0.  We search for
# one by exact Groebner elimination over complete coefficient slices.
#
# NORMALISATION used for e of degree d in C[w1,w2,w3] (all steps classical or
# certified here):
#   * [B(e)]_top = B(e_d) = (d/(d-1)) e_d det Hess(e_d)   (Identity E; certified
#     below for ternary forms d = 3, 4), so B == 0 forces det Hess(e_d) == 0;
#   * a ternary form with vanishing Hessian is a cone (Gordan--Noether, n = 3;
#     classical, cited from Watanabe--de Bondt in the project ledger, T2), so
#     after a linear change e_d in C[w1,w2];
#   * (degree 4 only) the binary quartic e_4 is moved by GL_2 acting on
#     (w1,w2) -- which preserves the shape e_4 in C[w1,w2] -- to one of the
#     root-multiplicity normal forms
#        (4)      w1^4              (3,1)    w1^3 w2
#        (2,2)    w1^2 w2^2         (2,1,1)  w1^2 w2 (w1+w2)
#        (1,1,1,1) w1 w2 (w1+w2)(w1 + lam w2),  lam(lam-1) != 0  [lam a symbol]
#     (PGL_2 is 3-transitive on P^1; e_4 = 0 is degree < 4);
#   * translations w -> w + v add D_v e_4 to e_3 (certified), which lets us
#     kill the coefficients of e_3 listed per slice; scaling e is free;
#     the constant term is irrelevant.
# Modes (command-line argument):
#   deg3              complete decision, cubic part in C[w1,w2] (13 unknowns)
#   deg4:<form>       complete decision on the named slice
#                     form in {p4, p31, p22, p211, p1111}
#   pivot             structured ansatz e = a(w1,w2) + w3 b(w1,w2) + gam w3^2/2
#                     of degree <= 4  (affine-pivot shapes, dBvdE family)
#   sample            random exact contrapositive sampling deg 3,4,5
# Each mode certifies, for every exact solution of B == 0 in its slice, that
# det Hess e == 0 (radical membership: power reduction, else Rabinowitsch).
# A slice that does not finish within the environment's ten-minute kill is
# reported by its missing "DONE" line; that is NOT a refutation.
# ============================================================================

import itertools
import random
import sys
import time
import sympy as sp

T0 = time.time()
w1, w2, w3 = W = sp.symbols('w1 w2 w3')
zz = sp.Symbol('zz')
FAILS = []


def check(name, cond):
    print(f'   {"PASS" if cond else "*** FAIL ***"}  {name}', flush=True)
    if not cond:
        FAILS.append(name)


def hess(e, vs=W):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def Bof(e, vs=W):
    g = sp.Matrix([sp.diff(e, a) for a in vs])
    return sp.expand((g.T*hess(e, vs).adjugate()*g)[0, 0])


def dethess(e, vs=W):
    return sp.expand(hess(e, vs).det(method='berkowitz'))


def coeffs(expr, vs=W):
    return [sp.expand(c) for c in sp.Poly(expr, *vs).coeffs()] if sp.expand(expr) != 0 else []


def monomials(vs, degs):
    out = []
    for k in degs:
        out += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, k)}, key=str)
    return out


def decide(e, params, label, extra_eqs=(), nonzero=()):
    """Certify: every coefficient of det Hess e lies in rad(<coeffs of B(e)> + extra),
    saturated by the product of `nonzero` (Rabinowitsch with 1 - zz*prod)."""
    t0 = time.time()
    Beqs = coeffs(Bof(e)) + list(extra_eqs)
    Deqs = [q for q in coeffs(dethess(e)) if q != 0]
    print(f' [{label}] {len(params)} unknowns, {len(Beqs)} equations, '
          f'{len(Deqs)} det-Hess coefficients', flush=True)
    gens = list(params)
    if nonzero:
        Beqs = Beqs + [1 - zz*sp.prod(nonzero)]
        gens = gens + [zz]
    G = sp.groebner(Beqs, *gens, order='grevlex')
    print(f'   Groebner basis: {len(G.exprs)} elements [{time.time()-t0:.1f}s]', flush=True)
    if list(G.exprs) == [sp.Integer(1)]:
        print('   (B == 0 has NO solution in this slice: vacuous, hence true)')
        check(f'{label}: no solutions of B==0 in slice', True)
        return
    zz2 = sp.Symbol('zz2')
    n_direct = n_rab = 0
    for q in Deqs:
        acc = sp.Integer(1)
        done = False
        for N in range(1, 5):
            acc = sp.expand(acc*q)
            if G.reduce(acc)[1] == 0:
                done = True
                n_direct += 1
                break
        if not done:
            GR = sp.groebner(Beqs + [1 - zz2*q], *(gens + [zz2]), order='grevlex')
            if list(GR.exprs) == [sp.Integer(1)]:
                n_rab += 1
                done = True
        if not done:
            check(f'{label}: det-Hess coefficient {q} NOT in radical -- '
                  f'THEOREM G CONTRADICTED in this slice', False)
            return
    check(f'{label}: all {len(Deqs)} det-Hess coefficients in radical '
          f'({n_direct} by powers, {n_rab} Rabinowitsch) [{time.time()-t0:.1f}s]', True)


mode = sys.argv[1] if len(sys.argv) > 1 else 'deg3'
print(f'atkG n=3 search, mode = {mode}')

# ---------------------------------------------------------------- shared pillars
if mode in ('deg3', 'deg4:p4'):
    print('(pillars) Identity E for ternary forms and translation action on e_3')
    for d in (3, 4):
        mons = monomials(W, [d])
        cs = sp.symbols(f'f{d}_0:{len(mons)}')
        f = sum(c*m for c, m in zip(cs, mons))
        check(f'd={d}: B(f) == d/(d-1) f det Hess f (generic ternary form)',
              sp.expand(Bof(f) - sp.Rational(d, d-1)*f*dethess(f)) == 0)
    v = sp.symbols('v1 v2 v3')
    mons4 = monomials(W, [4])
    cs4 = sp.symbols(f'q0:{len(mons4)}')
    f4 = sum(c*m for c, m in zip(cs4, mons4))
    shifted = sp.expand(f4.subs(dict(zip(W, [a+b for a, b in zip(W, v)])), simultaneous=True))
    deg3part = sum(t for t in sp.Add.make_args(shifted) if sp.Poly(t, *W).total_degree() == 3)
    check('translation: [e_4(w+v)]_3 == D_v e_4  (generic quartic)',
          sp.expand(deg3part - sum(v[i]*sp.diff(f4, W[i]) for i in range(3))) == 0)

# ---------------------------------------------------------------- deg3
if mode == 'deg3':
    c3 = sp.symbols('a0:4')
    q3 = sp.symbols('b0:6')
    l3 = sp.symbols('d0:3')
    cub = c3[0]*w1**3 + c3[1]*w1**2*w2 + c3[2]*w1*w2**2 + c3[3]*w2**3
    quad = sum(a*m for a, m in zip(q3, monomials(W, [2])))
    lin = sum(a*m for a, m in zip(l3, W))
    e = cub + quad + lin
    decide(e, list(c3) + list(q3) + list(l3), 'deg<=3, cubic part in C[w1,w2] (complete)')
    # control: an affinely 2-variable cubic with a genuinely 3-variable second form
    ectrl = (w1 + 2*w2 - w3)**3 + (w1 + 2*w2 - w3)*(w2 + 3*w3) + (w2 + 3*w3)**2
    check('control: affinely-2-variable cubic has B == 0 and det Hess == 0',
          Bof(ectrl) == 0 and dethess(ectrl) == 0)
    hdb = w1*w2 + w1**2*w3
    check('control: dBvdE witness w1 w2 + w1^2 w3 has det Hess == 0, B != 0',
          dethess(hdb) == 0 and Bof(hdb) != 0)

# ---------------------------------------------------------------- deg4 slices
if mode.startswith('deg4:'):
    form = mode.split(':')[1]
    lam = sp.Symbol('lam')
    forms = {
        'p4':    (w1**4,                       ['w1**3']),
        'p31':   (w1**3*w2,                    ['w1**3', 'w1**2*w2']),
        'p22':   (w1**2*w2**2,                 ['w1*w2**2', 'w1**2*w2']),
        'p211':  (w1**2*w2*(w1 + w2),          ['w1**3', 'w1*w2**2']),
        'p1111': (w1*w2*(w1 + w2)*(w1 + lam*w2), ['w1**3', 'w2**3']),
    }
    e4, kill = forms[form]
    mons3 = monomials(W, [3])
    mons2 = monomials(W, [2])
    a3 = sp.symbols(f'a0:{len(mons3)}')
    b2 = sp.symbols(f'b0:{len(mons2)}')
    d1 = sp.symbols('d0:3')
    # certify that the killed e_3-coefficients are reachable by translation:
    # the map v -> D_v e_4 must hit those monomials independently.
    v = sp.symbols('v1 v2 v3')
    Dv = sp.expand(sum(v[i]*sp.diff(e4, W[i]) for i in range(3)))
    killm = [sp.sympify(k) for k in kill]
    Mk = sp.Matrix([[sp.Poly(Dv, *W).coeff_monomial(m).coeff(v[i]) for i in range(2)] for m in killm])
    rk = Mk.rank(simplify=True) if form != 'p1111' else Mk.subs(lam, 7).rank()
    check(f'{form}: translations can kill e_3 coefficients {kill} (rank {rk} = {len(kill)})',
          rk == len(kill))
    params = [a for a, m in zip(a3, mons3) if m not in killm]
    e3 = sum(a*m for a, m in zip(a3, mons3) if m not in killm)
    e2 = sum(b*m for b, m in zip(b2, mons2))
    e1 = sum(c*m for c, m in zip(d1, W))
    e = e4 + e3 + e2 + e1
    allp = params + list(b2) + list(d1)
    if form == 'p1111':
        allp = allp + [lam]
        decide(e, allp, f'deg 4, e_4 = {e4} (complete slice)', nonzero=[lam, lam - 1])
    else:
        decide(e, allp, f'deg 4, e_4 = {e4} (complete slice)')
    print('DONE', flush=True)

# ---------------------------------------------------------------- pivot ansatz
if mode == 'pivot':
    Y = (w1, w2)
    gam = sp.Symbol('gam')
    for da, db in ((4, 3), (4, 2), (3, 3), (5, 3)):
        ma = monomials(Y, range(1, da+1))
        mb = monomials(Y, range(0, db+1)) if db > 0 else [sp.Integer(1)]
        ca = sp.symbols(f'A0:{len(ma)}')
        cb = sp.symbols(f'C0:{len(mb)}')
        a = sum(c*m for c, m in zip(ca, ma))
        b = sum(c*m for c, m in zip(cb, mb))
        e = a + w3*b + gam*w3**2/2
        decide(e, list(ca) + list(cb) + [gam], f'pivot ansatz deg a<={da}, deg b<={db}')
    print('DONE', flush=True)

# ---------------------------------------------------------------- sampling
if mode == 'sample':
    rnd = random.Random(2026_08_24)
    viol = tested = 0
    for d in (3, 4, 5):
        mons = monomials(W, range(1, d+1))
        for _ in range(25):
            e = sum(sp.Rational(rnd.randint(-4, 4), rnd.randint(1, 3))*m for m in mons)
            if dethess(e) == 0:
                continue
            tested += 1
            if Bof(e) == 0:
                viol += 1
                print('  REFUTATION:', e)
    check(f'{tested} random exact instances with det Hess != 0 all have B != 0 '
          f'(violations: {viol})', viol == 0)
    # structured families that MUST satisfy B == 0 (positive controls)
    for k in (2, 3, 4):
        L1 = 2*w1 - w2 + 3*w3
        L2 = w1 + w2 - w3
        e = L1**k + L1*L2**(k-1) - 3*L2**2 + L1
        check(f'control: affinely 2-variable degree {k} has B == 0 and det Hess == 0',
              Bof(e) == 0 and dethess(e) == 0)
    print('DONE', flush=True)

print(f'total {time.time()-T0:.1f}s')
if FAILS:
    print('FAILED:', FAILS)
    sys.exit(1)
print('mode', mode, ': all checks PASS')
sys.exit(0)
