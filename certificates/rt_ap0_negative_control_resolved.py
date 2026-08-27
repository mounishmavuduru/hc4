# Resolution of the failed assertion in rt_instances_AP0.py (agent-authored).
#
# That script's "S3 negative control" asserted that
#     f = x1*x2 + x2^2*x3 + x3*x4
# has NON-constant Hessian determinant, expecting "x2-degree >= 2 fiber
# coupling breaks constancy". The assertion FAILED. This script establishes
# why: the expectation was wrong, and f is in fact a perfectly valid member
# of the project's AP0 class -- so the intended negative control was actually
# a positive instance, and Theorem A predicts it is injective.
#
#   Hess f = [[0,1,0,0],[1,2x3,2x2,0],[0,2x2,0,1],[0,0,1,0]],  det == 1.
#
# In AP0 terms: pivot v = e4, e1 = x3 (AFFINE), e0 = x1x2 + x2^2 x3, and the
# AP0 condition det2 Hess_(x1,x2) e0 == -c reads 0*2x3 - 1^2 = -1 = -c with
# c = 1. So Theorem A (case gamma = 0, deg e1 <= 1) applies and grad f must be
# injective. Verified below by saturated collision Groebner.
#
# Conclusion: no defect in the project's results; the agent's control was
# mis-designed. Recorded because a failed assert must be explained, not
# silently dropped.

import sympy as sp

x1, x2, x3, x4 = X = sp.symbols('x1 x2 x3 x4')

f = sp.expand(x1*x2 + x2**2*x3 + x3*x4)
H = sp.hessian(f, X)
d = sp.expand(H.det(method='berkowitz'))
print('f =', f)
print('det Hess f =', d)
assert d == 1, 'expected the constant 1'
print('=> CONSTANT, so the agent\'s "negative control" is a Keller potential.')

# it is an AP0 member: f = e0(x1,x2,x3) + x4*e1(x1,x2,x3) with e1 affine
e1 = sp.diff(f, x4)
e0 = sp.expand(f - x4*e1)
print('e1 =', e1, ' (degree', sp.Poly(e1, x1, x2, x3).total_degree(), '-> affine: AP0)')
print('e0 =', e0)
assert sp.expand(sp.diff(e1, x4)) == 0 and sp.Poly(e1, x1, x2, x3).total_degree() <= 1

# AP0 condition: det2 Hess_(x1,x2) e0 == -c
det2 = sp.expand(sp.diff(e0, x1, 2)*sp.diff(e0, x2, 2) - sp.diff(e0, x1, x2)**2)
print('det2 Hess_(x1,x2) e0 =', det2, ' (must equal -c = -1)')
assert det2 == -1

# Theorem A therefore predicts injectivity; verify it independently.
p = sp.symbols('pn1:5'); q = sp.symbols('qn1:5'); w = sp.symbols('wn1:5')
gr = [sp.diff(f, v) for v in X]
eqs = [sp.expand(g.subs(dict(zip(X, p)), simultaneous=True)
                 - g.subs(dict(zip(X, q)), simultaneous=True)) for g in gr]
sat = sp.expand(sum(w[i]*(p[i]-q[i]) for i in range(4)) - 1)
G = sp.groebner(eqs + [sat], *(list(p)+list(q)+list(w)), order='grevlex')
print('saturated collision Groebner basis =', list(G.exprs))
assert list(G.exprs) == [sp.Integer(1)]
print()
print('RESOLVED: the failed assertion was a mis-designed control, not a defect.')
print('The potential is an AP0 instance and its gradient is injective, exactly')
print('as Theorem A predicts.')
