# Certificate check: the Alpoge 3D counterexample to the Jacobian conjecture
# (as reported in press/search results, to be cross-checked against arXiv:2608.00222).
# Claim: F = (P,Q,R) has det JF == -2 identically, yet three distinct points share one image.
# Everything below is exact rational arithmetic in sympy.

import sympy as sp

x, y, z = sp.symbols('x y z')

P = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
Q = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
R = 2*x - 3*x**2*y - x**3*z

F = sp.Matrix([P, Q, R])
J = F.jacobian([x, y, z])
detJ = sp.expand(J.det(method='berkowitz'))
print("det JF =", detJ)
assert detJ == -2, "Jacobian determinant is not identically -2"

pts = [
    (sp.Integer(0), sp.Integer(0), sp.Rational(-1, 4)),
    (sp.Integer(1), sp.Rational(-3, 2), sp.Rational(13, 2)),
    (sp.Integer(-1), sp.Rational(3, 2), sp.Rational(13, 2)),
]
images = [tuple(sp.simplify(c.subs({x: p[0], y: p[1], z: p[2]})) for c in F) for p in pts]
for p, im in zip(pts, images):
    print(f"F{p} = {im}")

assert len(set(pts)) == 3, "points are not distinct"
assert images[0] == images[1] == images[2], "images differ"
print("VERIFIED: det JF = -2 identically; 3 distinct rational points share the image", images[0])
