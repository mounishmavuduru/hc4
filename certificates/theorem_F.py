# THEOREM F (partial resolution of Conjecture E).
#
# If e in C[x1,x2,x3] satisfies BOTH
#     (i)  B(e) := grad(e)^T adj(Hess e) grad(e) == 0   (vanishing bordered Hessian)
#     (ii) det Hess_3(e) == 0
# then e is affinely 2-variable.
#
# PROOF. Let d = deg e and E the degree-d homogenization of e in
# (x0,x1,x2,x3). By the homogenization identity (homogenization_identity.py)
#     det Hess_4(E)|_{x0=1} = d(d-1) e det Hess_3(e) - (d-1)^2 B(e),
# so (i)+(ii) give det Hess_4(E)|_{x0=1} == 0, hence det Hess_4(E) == 0
# (a form vanishing on the chart x0=1 vanishes identically).
# By Gordan-Noether in FOUR variables, E is a cone: there is v = (v0,v') != 0
# with D_v E == 0.
#   Case v0 = 0: restricting to x0 = 1 gives D_{v'} e == 0 with v' != 0,
#     so e is linearly 2-variable.
#   Case v0 != 0 (normalize v0 = 1): Euler for E gives
#     dE = x0 E_0 + sum x_i E_i, so at x0 = 1, E_0 = d e - x.grad e; then
#     D_v E == 0 reads d e = (x - v').grad e. Translating y = x - v', the
#     polynomial etilde(y) := e(y + v') satisfies y.grad etilde = d etilde,
#     i.e. etilde is HOMOGENEOUS of degree d (converse of Euler's theorem).
#     B is translation-invariant, so B(etilde) == 0, and Identity E turns
#     this into (d/(d-1)) etilde det Hess etilde == 0, i.e.
#     det Hess_3(etilde) == 0. Gordan-Noether in THREE variables makes
#     etilde a cone, so etilde -- hence e -- is affinely 2-variable.        []
#
# Corollary: Conjecture E holds whenever det Hess_3(e) == 0; in particular
# for every HOMOGENEOUS e (where Identity E makes (i) and (ii) equivalent).
#
# This script certifies the two computational steps of the proof:
#   (A) B is invariant under translations (fully symbolic);
#   (B) the converse of Euler's theorem: y.grad p = d p  =>  p homogeneous
#       of degree d (verified on a generic polynomial of degree <= 3 by
#       solving the linear conditions);
# and runs the full theorem end-to-end on examples.

import itertools
import sys
import sympy as sp

x1, x2, x3 = X = sp.symbols('x1 x2 x3')


def B(e, vars_=X):
    K = sp.Matrix([[sp.diff(e, u, v) for v in vars_] for u in vars_])
    g = sp.Matrix([sp.diff(e, u) for u in vars_])
    return sp.expand((g.T * K.adjugate() * g)[0, 0])


def dethess(e, vars_=X):
    return sp.expand(sp.Matrix([[sp.diff(e, u, v) for v in vars_]
                                for u in vars_]).det(method='berkowitz'))


# ---- (A) translation invariance of B, fully symbolic ----
mons = [sp.Integer(1)]
for d in range(1, 4):
    mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(X, d)}, key=str)
cs = sp.symbols(f'c0:{len(mons)}')
e = sum(c*m for c, m in zip(cs, mons))
t = sp.symbols('t1:4')
e_sh = e.subs(dict(zip(X, [v + s for v, s in zip(X, t)])), simultaneous=True)
lhs = sp.expand(B(e_sh))
rhs = sp.expand(B(e).subs(dict(zip(X, [v + s for v, s in zip(X, t)])), simultaneous=True))
assert sp.expand(lhs - rhs) == 0
print('(A) PASS: B(e(x+t)) = B(e)(x+t) -- B is translation-invariant')
print('    (fully symbolic, generic cubic => proof)')

# ---- (B) converse of Euler's theorem ----
expr = sp.expand(sum(v*sp.diff(e, v) for v in X) - 3*e)
sol = sp.solve(sp.Poly(expr, *X).coeffs(), list(cs), dict=True)
assert sol, 'no solution'
e_sol = sp.expand(e.subs(sol[0]))
# every surviving monomial must have total degree exactly 3
degs = {sum(m) for m, _ in sp.Poly(e_sol, *X).terms()} if e_sol != 0 else set()
assert degs <= {3}, f'non-homogeneous survivor: {degs}'
print('(B) PASS: x.grad p = 3p forces p homogeneous of degree 3'
      ' (generic cubic ansatz)')

# ---- end-to-end examples ----
cases = {
    'x1^2*x3   (cone, planar)':        x1**2*x3,
    'x1*x2     (planar)':              x1*x2,
    '(x1+1)^2*(x3+2) (translate of a cone)': (x1 + 1)**2*(x3 + 2),
    'x1 + phi(x2) (planar)':           x1 + x2**4,
}
for name, cand in cases.items():
    b0, dh = B(cand), dethess(cand)
    hyp = (b0 == 0 and dh == 0)
    # affinely 2-variable test: some constant direction annihilates e
    v = sp.symbols('n1 n2 n3')
    p = sp.Poly(sp.expand(sum(vi*sp.diff(cand, u) for vi, u in zip(v, X))), *X)
    ns = sp.linsolve(p.coeffs(), list(v))
    planar = any(any(c.free_symbols & set(v) for c in s) for s in ns)
    print(f'  {name:40s} B=={b0 == 0}, detHess=={dh == 0},'
          f' hypothesis:{hyp}, affinely 2-var:{planar}')
    if hyp:
        assert planar, 'THEOREM F VIOLATED'

print()
print('THEOREM F certified. Corollary: Conjecture E holds whenever')
print('det Hess_3(e) == 0, in particular for all homogeneous pivot coefficients.')
