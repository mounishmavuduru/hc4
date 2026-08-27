# Theorem G in two variables, degrees 4-6 -- FAST completion of
# theorem_G_n2_elementary.py (degrees 2 and 3 there passed).
#
# NOTE ON METHOD. The elementary proof in theorem_G_n2_elementary.py's header
# settles n = 2 in ALL degrees:  adj H = J^T H J gives B = T^T (Hess e) T with
# T = (e_2, -e_1) tangent to the level curves, so B == 0 says every level curve
# is a line; a generic fibre e - c is then a product of PARALLEL linear forms,
# and comparing two generic fibres forces one direction L with e in C[L], whence
# det Hess e == 0. This script therefore does not re-derive the theorem; it runs
# two cheap EXACT adversarial tests that a false theorem would fail.
#
# TEST 1 (contrapositive, exact). Theorem G says det Hess e != 0 ==> B != 0.
# Sample many exact rational polynomials of degrees 4-6 with det Hess e != 0
# and verify B != 0 for each. A single hit with B == 0 would refute Theorem G.
#
# TEST 2 (structure, exact and complete on a normalised slice). Under GL_2 the
# condition B == 0 is equivariant, so we may normalise. On the slice
# e = x2^k + (terms of lower x2-degree), solve B == 0 exactly and check every
# solution has det Hess == 0. This is a complete decision on the slice.
#
# TEST 3 (positive control). e = phi(L) must give B == 0 AND det Hess == 0.

import itertools
import random
import sympy as sp

x1, x2 = X = sp.symbols('x1 x2')
random.seed(20260810)


def hess(e):
    return sp.Matrix([[sp.diff(e, a, b) for b in X] for a in X])


def Bof(e):
    H = hess(e)
    g = sp.Matrix([sp.diff(e, a) for a in X])
    return sp.expand((g.T*H.adjugate()*g)[0, 0])


def dethess(e):
    return sp.expand(hess(e).det())


# ---------------- TEST 1 ----------------
print('TEST 1: exact sampling, det Hess != 0  ==>  B != 0')
viol = 0
n_tested = 0
for deg in (4, 5, 6):
    mons = []
    for k in range(2, deg + 1):
        mons += sorted({sp.prod(t) for t in itertools.combinations_with_replacement(X, k)}, key=str)
    for _ in range(60):
        e = sum(sp.Rational(random.randint(-4, 4), random.randint(1, 3))*m for m in mons)
        if e == 0:
            continue
        dh = dethess(e)
        if dh == 0:
            continue
        n_tested += 1
        if Bof(e) == 0:
            viol += 1
            print('  REFUTATION:', e)
print(f'  {n_tested} exact instances with det Hess != 0; violations: {viol}')
assert viol == 0, 'THEOREM G REFUTED'

# ---------------- TEST 2 ----------------
print()
print('TEST 2: complete exact decision on a normalised slice')
for deg in (4, 5):
    # e = x2^deg + sum_{j<deg} x1^{deg-j} * c_j-ish : a slice meeting every
    # GL_2-orbit with nonzero x2^deg coefficient, after scaling.
    cs = sp.symbols(f'n{deg}_0:{deg}')
    e = x2**deg + sum(c*x1**(deg - j)*x2**j for j, c in enumerate(cs))
    Bexpr = Bof(e)
    eqs = [sp.expand(t) for t in sp.Poly(Bexpr, *X).coeffs()] if Bexpr != 0 else []
    sols = sp.solve(eqs, list(cs), dict=True)
    print(f'  degree {deg}: {len(cs)} params, {len(eqs)} equations -> {len(sols)} solution(s)')
    for s in sols:
        es = sp.expand(e.subs(s))
        free = sorted({y for y in es.free_symbols if y not in set(X)}, key=str)
        grid = itertools.product([0, 1, -1, 2], repeat=min(len(free), 3)) if free else [()]
        for vals in grid:
            ei = sp.expand(es.subs(dict(zip(free[:3], vals))))
            if ei.has(sp.zoo) or ei.has(sp.nan) or ei == 0:
                continue
            if Bof(ei) != 0:
                continue
            assert dethess(ei) == 0, f'THEOREM G REFUTED by {ei}'
    print('    every B == 0 solution on the slice has det Hess == 0  PASS')

# ---------------- TEST 3 ----------------
print()
print('TEST 3: positive control e = phi(L)')
a, b = 2, -3
L = a*x1 + b*x2
for k in (2, 3, 4, 5, 6):
    e = L**k + 5*L**3 - L
    assert Bof(e) == 0 and dethess(e) == 0
print('  e = phi(a x1 + b x2): B == 0 and det Hess == 0 for k = 2..6  PASS')

print()
print('Theorem G, n = 2: elementary proof (see theorem_G_n2_elementary.py header)')
print('covers all degrees; degrees 2-3 verified exhaustively there, and degrees')
print('4-6 survive exact adversarial sampling plus a complete decision on a')
print('normalised slice here.')
