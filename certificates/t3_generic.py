# T3 fully generic: for f = f4(x1,x2,x3) + f3(x1..x4) + f2(x1..x4) with ALL
# coefficients symbolic, the degree-7 graded piece of det Hess f equals
# det3(Hess3 f4) * (f3)_{x4x4}.  The identity is multilinear in the graded
# pieces of f, so the fully symbolic verification below is a PROOF for all f
# of degree 4 (with GN-normalized leading form).

import itertools
import sys
import sympy as sp

sys.path.insert(0, r'C:\Users\mouni\hc4\src')
from hc4lib import hessian, graded_part

x1, x2, x3, x4 = X = sp.symbols('x1 x2 x3 x4')

def generic_form(vars_, deg, tag):
    mons = list(itertools.combinations_with_replacement(range(len(vars_)), deg))
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    f = sp.Integer(0)
    for c, m in zip(cs, mons):
        term = c
        for i in m:
            term *= vars_[i]
        f += term
    return f

f4 = generic_form([x1, x2, x3], 4, 'a')   # GN-normalized leading form, generic
f3 = generic_form(list(X), 3, 'b')
f2 = generic_form(list(X), 2, 'c')
f = f4 + f3 + f2

D = sp.expand(hessian(f, list(X)).det(method='berkowitz'))
lhs = graded_part(D, list(X), 7)
rhs = sp.expand(hessian(f4, [x1, x2, x3]).det(method='berkowitz') * sp.diff(f3, x4, x4))
assert sp.expand(lhs - rhs) == 0
print('T3 GENERIC PASS: [det Hess f]_7 = det3(Hess3 f4) * (f3)_{x4x4}')
print('holds with fully symbolic coefficients (proof for all degree-4 f).')
