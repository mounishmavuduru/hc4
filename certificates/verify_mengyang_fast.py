# Fast exact verification (domain-ge over QQ(x)) of the two Meng-Yang determinant identities.
import sympy as sp

x1, x2, x3, y1, y2, y3 = sp.symbols('x1 x2 x3 y1 y2 y3')
u = 1 + x1*x2
F1 = u**3*x3 + x2**2*u*(4 + 3*x1*x2)
F2 = x2 + 3*x1*u**2*x3 + 3*x1*x2**2*(4 + 3*x1*x2)
F3 = 2*x1 - 3*x1**2*x2 - x1**3*x3

phi = y1*F1 + y2*F2 + y3*F3
vars6 = [x1, x2, x3, y1, y2, y3]
H6 = sp.hessian(phi, vars6)
det6 = sp.factor(H6.det(method='domain-ge'))
print("det Hess phi =", det6)
assert det6 == -4

A = y1*u**3 + 3*x1*y2*u**2 - x1**3*y3
B = (y1*x2**2*u*(4 + 3*x1*x2)
     + y2*(x2 + 3*x1*x2**2*(4 + 3*x1*x2))
     + y3*(2*x1 - 3*x1**2*x2))
assert sp.expand(phi - (x3*A + B)) == 0

lam, mu = sp.symbols('lam mu')
vars5 = [x1, x2, y1, y2, y3]
psi_fam = B + lam/2*A**2 + mu*A
H5f = sp.hessian(psi_fam, vars5)
det5f = sp.factor(H5f.det(method='domain-ge'))
print("det Hess psi_{lam,mu} =", det5f)
assert sp.expand(det5f - 4*lam) == 0

Psi = sp.expand(A**2 + 13*A + 2*B)
print("deg =", sp.Poly(Psi, *vars5).total_degree(), " monomials =", len(sp.Poly(Psi, *vars5).terms()))
grad5 = [sp.diff(Psi, v) for v in vars5]
Pp = {x1: 1, x2: sp.Rational(-3, 2), y1: 0, y2: 0, y3: 0}
Pm = {x1: -1, x2: sp.Rational(3, 2), y1: 0, y2: 0, y3: 0}
gp = tuple(g.subs(Pp) for g in grad5)
gm = tuple(g.subs(Pm) for g in grad5)
print("grad Psi(P+) =", gp, " grad Psi(P-) =", gm)
assert gp == gm == (0, 0, sp.Rational(-1, 2), 0, 0)
# det Hess Psi = 2^5 * det Hess psi_{1, 13/2}  (Psi = 2*psi_{1,13/2}), i.e. 32*4 = 128
det5 = sp.factor(sp.hessian(Psi, vars5).det(method='domain-ge'))
print("det Hess Psi =", det5)
assert det5 == 128
print("ALL FAST CHECKS PASSED")
