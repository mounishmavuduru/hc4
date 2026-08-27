# scratch: explicit automorphism tau, Psi o tau affine in z1?  is det Hess still constant?
import sympy as sp

x1, x2, y1, y2, y3 = sp.symbols('x1 x2 y1 y2 y3')
z1, z2, z3 = sp.symbols('z1 z2 z3')
lam, mu = sp.symbols('lam mu')
U = 1 + x1*x2
A = y1*U**3 + 3*x1*y2*U**2 - x1**3*y3
b1 = x2**2*U*(4 + 3*x1*x2)
b2 = x2 + 3*x1*x2**2*(4 + 3*x1*x2)
b3 = 2*x1 - 3*x1**2*x2
B = y1*b1 + y2*b2 + y3*b3
Psi = sp.expand(B + sp.Rational(1, 2)*lam*A**2 + mu*A)

N = sp.Matrix([[-x2/3, -1, 0], [U, 3*x1, 0], [0, 0, 1]])
Ninv = sp.Matrix([[3*x1, 1, 0], [-U, -x2*sp.Rational(1, 3), 0], [0, 0, 1]])
assert sp.expand(N*Ninv - sp.eye(3)) == sp.zeros(3, 3)
assert sp.expand(N.det()) == 1
print("N invertible over C[x1,x2], det N = 1  OK")

zv = sp.Matrix([z1, z2, z3])
ysub = Ninv*zv
sub = {y1: ysub[0], y2: ysub[1], y3: ysub[2]}
Pt = sp.expand(Psi.subs(sub, simultaneous=True))
At = sp.expand(A.subs(sub, simultaneous=True))
print("A o tau =", sp.factor(At))
print("d2 Psit / dz1^2 =", sp.simplify(sp.diff(Pt, z1, 2)))
print("d Psit / dz1 =", sp.factor(sp.simplify(sp.diff(Pt, z1))))

# bordered determinant:  Psit = z1*Atil(w) + Btil(w),  w = (x1,x2,z2,z3)
w = [x1, x2, z2, z3]
Atil = sp.expand(sp.diff(Pt, z1))
Btil = sp.expand(Pt.subs(z1, 0))
assert sp.expand(Pt - z1*Atil - Btil) == 0
g = sp.Matrix([sp.diff(Atil, v) for v in w])
K = sp.Matrix(4, 4, lambda i, j: sp.diff(Atil, w[i], w[j]))
P = sp.Matrix(4, 4, lambda i, j: sp.diff(Btil, w[i], w[j]))
M = sp.expand(P + z1*K)
print("g =", list(g))

# det Hess_5 Psit = - g^T adj(M) g ; g supported in coords 0,1
def minor(Mx, rows, cols):
    return Mx.extract(rows, cols).det(method='berkowitz')

adj11 = minor(M, [1, 2, 3], [1, 2, 3])
adj22 = minor(M, [0, 2, 3], [0, 2, 3])
adj12 = -minor(M, [0, 2, 3], [1, 2, 3])
val = sp.expand(g[0]**2*adj11 + 2*g[0]*g[1]*adj12 + g[1]**2*adj22)
dh = sp.expand(-val)
print("det Hess(Psi o tau) =", sp.factor(dh))
