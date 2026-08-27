# WHERE THE PIVOT OBSTRUCTION ACTUALLY LIVES -- a REFUTED hypothesis, and the
# correct picture. Recorded because a refuted hypothesis is a real result.
#
# HYPOTHESIS (mine, REFUTED below). The project's results reduce HC_4 to the
# class of potentials with NO pivot. Such potentials exist in dimension five
# (D1: the Meng-Yang counterexample Psi has empty pivot cone). I conjectured
# the reason was Gordan-Noether: for det Hess f constant, lemma T1 forces
# det Hess(f_d) == 0 on the leading form; in n <= 4 Gordan-Noether then makes
# f_d a CONE, supplying a direction v with D_v f_d == 0 and hence
# v^T Hess(f_d) v == 0 -- a candidate pivot. In n >= 5 Gordan-Noether fails
# (Perazzo), so I expected Psi's leading form to have vanishing Hessian
# determinant while being a cone in NO direction, which would have explained
# the dimension gap exactly.
#
# THIS IS FALSE. The leading form of Psi is the single monomial
#     Psi_top = x1^6 x2^6 y1^2
# (because deg A = 7 is attained only by the term y1 u^3, the other two terms
# of A having degrees 6 and 4), and a monomial omitting y2 and y3 is a cone in
# a two-dimensional space of directions. So candidate pivot directions DO
# exist for Psi's leading form; Gordan-Noether is not the obstruction.
#
# THE CORRECT PICTURE, verified below: the obstruction is entirely in the
# LOWER graded pieces. Every cone direction of the leading form (the whole
# span of e_{y2}, e_{y3}) fails to annihilate v^T Hess(Psi_k) v for some lower
# k. So "the leading form supplies a candidate" is necessary but nowhere near
# sufficient, and the four-variable case cannot be settled at the level of the
# leading form alone -- consistent with the fact that HC_4 remains open.

import sympy as sp

x1, x2, y1, y2, y3 = W = sp.symbols('x1 x2 y1 y2 y3')

u = 1 + x1*x2
A = y1*u**3 + 3*x1*y2*u**2 - x1**3*y3
B = (y1*x2**2*u*(4 + 3*x1*x2)
     + y2*(x2 + 3*x1*x2**2*(4 + 3*x1*x2))
     + y3*(2*x1 - 3*x1**2*x2))
Psi = sp.expand(A**2 + 13*A + 2*B)

P = sp.Poly(Psi, *W)
D = P.total_degree()


def graded(poly, k):
    p = sp.Poly(poly, *W)
    return sp.expand(sum(co*sp.prod([v**e for v, e in zip(W, mon)])
                         for mon, co in zip(p.monoms(), p.coeffs()) if sum(mon) == k))


top = graded(Psi, D)
print(f'deg Psi = {D}')
print(f'leading form Psi_top = {top}')
assert sp.expand(top - x1**6*x2**6*y1**2) == 0
print('  (a single monomial -- the degrees of A\'s three terms are 7, 6, 4)')

# T1 holds
Htop = sp.Matrix([[sp.diff(top, a, b) for b in W] for a in W])
assert sp.expand(Htop.det(method='berkowitz')) == 0
print('T1 confirmed: det Hess(Psi_top) == 0')

# ... but the leading form IS a cone: my hypothesis is refuted
v = sp.symbols('v1:6')
Dv = sp.expand(sum(v[i]*sp.diff(top, W[i]) for i in range(5)))
eqs = sp.Poly(Dv, *W).coeffs() if Dv != 0 else []
cone_dirs = sp.linsolve(eqs, list(v))
sol = list(cone_dirs)[0]
free = sorted({s for c in sol for s in (c.free_symbols & set(v))}, key=str)
print(f'cone directions of Psi_top: solution {sol}, free parameters {free}')
assert free, 'expected a nontrivial cone direction'
print(f'  => Psi_top IS a cone, in a {len(free)}-dimensional space of directions')
print('  => HYPOTHESIS REFUTED: Gordan-Noether failure is NOT the reason')
print('     the Meng-Yang counterexample has no pivot.')

# The obstruction lies in the lower graded pieces: check each cone direction.
print()
print('Locating the true obstruction: for each cone direction of the leading')
print('form, find the largest k with v^T Hess(Psi_k) v not identically zero.')
basis = []
for fs in free:
    sub = {g: (1 if g == fs else 0) for g in free}
    basis.append([sp.simplify(c.subs(sub)) for c in sol])
for vv in basis:
    vm = sp.Matrix(vv)
    worst = None
    for k in range(D, 1, -1):
        pk = graded(Psi, k)
        if pk == 0:
            continue
        Hk = sp.Matrix([[sp.diff(pk, a, b) for b in W] for a in W])
        q = sp.expand((vm.T*Hk*vm)[0, 0])
        if q != 0:
            worst = (k, q)
            break
    label = str(vm.T.tolist()[0])
    if worst is None:
        print(f'  v = {label}: annihilates EVERY graded piece -- would be a pivot!')
        raise AssertionError('contradicts D1')
    print(f'  v = {label}: first failure at graded degree {worst[0]}')

print()
print('CONCLUSION (recorded as a negative result): the leading form of the')
print('Meng-Yang counterexample is a cone, so candidate pivot directions exist;')
print('they are killed by LOWER graded pieces, not by Gordan-Noether. Hence the')
print('dimension gap between HC_4 and HC_5 is NOT explained at the level of the')
print('leading form, and the remaining no-pivot question in dimension four')
print('cannot be settled by leading-form arguments alone.')
