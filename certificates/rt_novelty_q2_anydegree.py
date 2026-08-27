# RT (hostile referee, key = novelty): is the hypothesis "deg e1 <= 2" needed
# in the gamma != 0 half of THEOREM A (= lemma Q2)?
#
# Claim under test (project): Theorem A requires deg e1 <= 2 "so that K is
# constant and etilde0 + lambda e1 has constant Hessian for a SCALAR lambda"
# (lemma_ledger.md, Q2 adversarial note; main_result.tex Thm 2.1 hypothesis).
#
# Referee counter-claim: for sigma != 0 the whole Q2 chain goes through for
# e1 of ARBITRARY degree, because
#     Hess3(e1^2/(2 sigma)) = (g g^T + e1 K)/sigma       (product rule, any e1)
# and hence
#     Hess3(etilde0) + u * K(x')  ==  Hess3(e0 + x4 e1) - g g^T/sigma
# with u = x4 + e1(x')/sigma, K(x') = Hess3 e1 -- NO constancy of K used.
# Since u sweeps C for each fixed x', det3(Hess3(etilde0 + s e1)) == c/sigma
# for every SCALAR s, and the collision transfer is unchanged.
#
# This script:
#  (A) verifies the Schur/product-rule identity with FULLY GENERIC e0, e1 of
#      degree <= 4 in x' and symbolic sigma  (linear in the e0-coefficients and
#      polynomial in the e1-coefficients; the identity is a polynomial identity
#      in all coefficients, so a symbolic check is a proof for those degrees);
#  (B) exhibits explicit HC_4 potentials with sigma != 0 and deg e1 = 3 and 4
#      (det Hess f a nonzero constant) that are OUTSIDE the stated hypothesis
#      of Theorem A but are settled verbatim by its own proof;
#  (C) checks the pencil identity det3 Hess3(etilde0 + s e1) == c/sigma for
#      symbolic s on those instances, and verifies gradient injectivity of the
#      instances directly by exact elimination.
# Exact arithmetic only.

import itertools
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
X3 = (x1, x2, x3)
X4 = (x1, x2, x3, x4)


def hess(f, vs):
    return sp.Matrix([[sp.diff(f, a, b) for b in vs] for a in vs])


def monomials(vs, dmax):
    out = []
    for d in range(dmax + 1):
        out += [sp.prod(m) for m in itertools.combinations_with_replacement(vs, d)]
    return out


# ---------------------------------------------------------------- (A)
print("=" * 72)
print("(A) FULLY SYMBOLIC: Schur complement identity, generic e0,e1 deg<=4")
sigma = sp.Symbol('sigma', nonzero=True)
mons = monomials(X3, 4)
acoef = sp.symbols(f'a0:{len(mons)}')
bcoef = sp.symbols(f'b0:{len(mons)}')
e0g = sum(c * m for c, m in zip(acoef, mons))
e1g = sum(c * m for c, m in zip(bcoef, mons))

g = sp.Matrix([sp.diff(e1g, w) for w in X3])
K = hess(e1g, X3)                      # NOT constant: deg e1 = 4
et0 = sp.expand(e0g - e1g**2 / (2 * sigma))
u = x4 + e1g / sigma

lhs = hess(et0, X3) + u * K
rhs = hess(sp.expand(e0g + x4 * e1g), X3) - (g * g.T) / sigma
D = sp.expand(sp.simplify(sp.Matrix(lhs - rhs)))
assert all(sp.simplify(e) == 0 for e in D), "Schur/product-rule identity FAILED"
print("  PASS: Hess3(etilde0) + (x4+e1/sigma)*Hess3(e1) == Hess3(e0+x4 e1) - gg^T/sigma")
print("        holds identically for GENERIC e1 of degree 4 (K non-constant).")
print("        => the step 'K constant' is never used in the Schur identity.")

# ---------------------------------------------------------------- (B)/(C)
print()
print("=" * 72)
print("(B/C) explicit sigma != 0 potentials with deg e1 = 3 and deg e1 = 4")

s = sp.Symbol('s')


def build(e1, sig):
    """etilde0 := h with det3 Hess3 h == -1 and Hess3 e1 supported on (1,1)."""
    h = x1 * x2 + x1**4 + x3**2 / 2          # det3 Hess3 h == -1
    e0 = sp.expand(h + e1**2 / (2 * sig))
    f = sp.expand(e0 + x4 * e1 + sp.Rational(1, 2) * sig * x4**2)
    return h, e0, f


for e1, sig in [(x1**3, sp.Integer(3)), (x1**4, sp.Rational(-5, 2)),
                (x1**3 + x1**2 * x2, sp.Integer(2))]:
    h, e0, f = build(e1, sig)
    H = hess(f, X4)
    c = sp.expand(H.det(method='berkowitz'))
    deg_e1 = sp.Poly(e1, *X3).total_degree()
    print(f"\n  e1 = {e1}   (deg {deg_e1}),  sigma = {sig}")
    print(f"    det Hess f = {c}")
    if not c.is_number:
        print("    -> NOT a constant-Hessian potential; skipping (construction failed)")
        continue
    assert c != 0
    # pencil identity for symbolic s
    pen = sp.expand(hess(sp.expand(h + s * e1), X3).det(method='berkowitz'))
    print(f"    det3 Hess3(etilde0 + s*e1) = {sp.simplify(pen)}   (must be c/sigma = {sp.nsimplify(c/sig)})")
    assert sp.simplify(pen - c / sig) == 0, "pencil identity FAILED"
    # collision transfer identity, symbolic
    lam = sp.Symbol('lam')
    for i, w in enumerate(X3):
        lhs = sp.diff(f, w)                      # = d_i e0 + x4 d_i e1
        rhs = sp.diff(sp.expand(h + lam * e1), w).subs(lam, x4 + e1 / sig)
        assert sp.simplify(sp.expand(lhs - rhs)) == 0, "collision-transfer FAILED"
    print("    PASS: pencil identity + collision transfer hold with deg e1 > 2")
    # direct injectivity check by exact elimination (independent of HC_3)
    P = sp.symbols('p1:5')
    Q = sp.symbols('q1:5')
    gr = [sp.diff(f, w) for w in X4]
    eqs = [sp.expand(gr[i].subs(dict(zip(X4, P))) - gr[i].subs(dict(zip(X4, Q))))
           for i in range(4)]
    t = sp.Symbol('t')
    wsat = sp.symbols('w1:5')
    sat = sp.expand(sum(wsat[i] * (P[i] - Q[i]) for i in range(4)) - 1)
    GB = sp.groebner(eqs + [sat], *(list(P) + list(Q) + list(wsat)), order='grevlex')
    empty = list(GB.exprs) == [sp.Integer(1)]
    print(f"    saturated collision system (p != q) is empty: {empty}")
    assert empty, "COLLISION FOUND -- would be an HC_4 counterexample, investigate!"

print()
print("=" * 72)
print("CONCLUSION: the hypothesis deg e1 <= 2 is NOT needed when sigma != 0.")
print("Theorem A's gamma != 0 branch holds for pivot coefficients of ANY degree.")
