# RED-TEAM (referee) independent re-derivation of the 2026 landscape inputs.
# Written from the PRIMARY TeX statements only (lit/meng-yang/jc_hc_status_note.tex,
# lit/gao/Jacobian_CE.tex), NOT from the project's own certificate scripts.
# Exact rational arithmetic; fail-closed asserts.

import sympy as sp

x1, x2, x3 = sp.symbols('x1 x2 x3')
y1, y2, y3 = sp.symbols('y1 y2 y3')

# ---------------------------------------------------------------- Alpoge map
u = 1 + x1*x2
F = sp.Matrix([
    u**3*x3 + x2**2*u*(4 + 3*x1*x2),
    x2 + 3*x1*u**2*x3 + 3*x1*x2**2*(4 + 3*x1*x2),
    2*x1 - 3*x1**2*x2 - x1**3*x3,
])
V = [x1, x2, x3]
J = F.jacobian(V)
jac = sp.expand(J.det(method='berkowitz'))
print('Jac F =', jac)
assert jac == -2

degs = [sp.Poly(F[i], *V).total_degree() for i in range(3)]
print('component degrees =', degs)
assert degs == [7, 6, 4]

pts = [(0, 0, sp.Rational(-1, 4)),
       (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
       (-1, sp.Rational(3, 2), sp.Rational(13, 2))]
imgs = [tuple(sp.simplify(F.subs(dict(zip(V, p)))))for p in pts]
print('images:', imgs)
for im in imgs:
    assert im == (sp.Rational(-1, 4), 0, 0), im
assert len(set(pts)) == 3

# ------------------------------------------------------- doubling  (HC_6)
W6 = [x1, x2, x3, y1, y2, y3]
phi = y1*F[0] + y2*F[1] + y3*F[2]
H6 = sp.hessian(phi, W6)
d6 = sp.expand(H6.det(method='berkowitz'))
print('det Hess phi (6 vars) =', d6)
assert d6 == -4
assert sp.Poly(phi, *W6).total_degree() == 8

# gradient collision of the doubled potential at (P,0)
g6 = sp.Matrix([sp.diff(phi, v) for v in W6])
vals = [tuple(sp.simplify(g6.subs(dict(zip(W6, list(p) + [0, 0, 0]))))) for p in pts]
print('grad phi at (P,0):', vals)
for v in vals:
    assert v == (0, 0, 0, sp.Rational(-1, 4), 0, 0), v

# --------------------------------------------------- Schur descent (HC_5)
# phi = x3*A + B   (split by x3-degree; phi is affine-linear in x3)
pol = sp.Poly(sp.expand(phi), x3)
assert pol.degree() == 1
A = sp.expand(pol.coeff_monomial(x3))
B = sp.expand(pol.coeff_monomial(1))
assert sp.expand(phi - (x3*A + B)) == 0

# paper's displayed A and B
A_paper = y1*u**3 + 3*x1*y2*u**2 - x1**3*y3
B_paper = (y1*x2**2*u*(4 + 3*x1*x2)
           + y2*(x2 + 3*x1*x2**2*(4 + 3*x1*x2))
           + y3*(2*x1 - 3*x1**2*x2))
assert sp.expand(A - A_paper) == 0
assert sp.expand(B - B_paper) == 0
print('A, B match the displayed formulas in Meng-Yang Sec. 6.')

W5 = [x1, x2, y1, y2, y3]
Psi = sp.expand(A**2 + 13*A + 2*B)
pP = sp.Poly(Psi, *W5)
print('deg Psi =', pP.total_degree(), ' #monomials =', len(pP.terms()))
assert pP.total_degree() == 14
assert len(pP.terms()) == 42
assert all(c.is_integer for c in pP.coeffs())

H5 = sp.hessian(Psi, W5)
d5 = sp.expand(H5.det(method='berkowitz'))
print('det Hess Psi =', d5)
assert d5 == 128

Pp = {x1: 1, x2: sp.Rational(-3, 2), y1: 0, y2: 0, y3: 0}
Pm = {x1: -1, x2: sp.Rational(3, 2), y1: 0, y2: 0, y3: 0}
g5 = sp.Matrix([sp.diff(Psi, v) for v in W5])
vp = tuple(sp.simplify(g5.subs(Pp)))
vm = tuple(sp.simplify(g5.subs(Pm)))
print('grad Psi(P+) =', vp, ' grad Psi(P-) =', vm)
assert vp == vm == (0, 0, sp.Rational(-1, 2), 0, 0)

# ------------------------------------------- two-parameter family det = 4*lam
lam, mu = sp.symbols('lam mu')
psi = B + sp.Rational(1, 2)*lam*A**2 + mu*A
Hf = sp.hessian(sp.expand(psi), W5)
df = sp.expand(Hf.det(method='berkowitz'))
df = sp.simplify(df)
print('det Hess psi_{lam,mu} =', df)
assert sp.expand(df - 4*lam) == 0

print()
print('ALL LANDSCAPE CHECKS PASSED (referee, independent).')
