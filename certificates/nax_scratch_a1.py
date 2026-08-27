# scratch: does a NON-affine automorphism make Meng-Yang psi affine-linear in a variable?
import sympy as sp

x1, x2, y1, y2, y3 = sp.symbols('x1 x2 y1 y2 y3')
lam, mu = sp.symbols('lam mu')
U = 1 + x1*x2
A = y1*U**3 + 3*x1*y2*U**2 - x1**3*y3
b1 = x2**2*U*(4 + 3*x1*x2)
b2 = x2 + 3*x1*x2**2*(4 + 3*x1*x2)
b3 = 2*x1 - 3*x1**2*x2
B = y1*b1 + y2*b2 + y3*b3
Psi = sp.expand(B + sp.Rational(1, 2)*lam*A**2 + mu*A)

# derivation u = 3x1 d/dy1 - U d/dy2
def u(h):
    return sp.expand(3*x1*sp.diff(h, y1) - U*sp.diff(h, y2))

print("u(A) =", sp.simplify(u(A)))
print("u(Psi) =", sp.factor(sp.simplify(u(Psi))))
print("u2(Psi) =", sp.simplify(u(u(Psi))))
