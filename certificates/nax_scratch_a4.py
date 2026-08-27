# scratch: (1) 5-var affine-pivot expansion -> (c3) is B_4(A)=0
#          (2) the "Gordan-Noether in 4 variables => A is a cone" step: TEST IT
import sympy as sp

# (1) fully generic bordered identity in 4+1 variables
t = sp.Symbol('t')
P = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'p{min(i,j)}{max(i,j)}'))
K = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'k{min(i,j)}{max(i,j)}'))
g = sp.Matrix([sp.Symbol(f'g{i}') for i in range(4)])
M = P + t*K
Bor = sp.zeros(5, 5)
Bor[0, 0] = 0
for i in range(4):
    Bor[0, i+1] = g[i]
    Bor[i+1, 0] = g[i]
    for j in range(4):
        Bor[i+1, j+1] = M[i, j]
lhs = sp.expand(Bor.det(method='berkowitz'))
rhs = sp.expand(-(g.T*M.adjugate()*g)[0, 0])
print('bordered identity (5x5, fully generic):', sp.expand(lhs - rhs) == 0)
pl = sp.Poly(lhs, t)
print('degree in t:', pl.degree())
top = sp.expand(pl.coeff_monomial(t**3))
predicted = sp.expand(-(g.T*K.adjugate()*g)[0, 0])
print('t^3 coefficient == -g^T adj(K) g :', sp.expand(top - predicted) == 0)

# (2) candidate counterexample to "B_4(A)=0 => A affinely 3-variable"
w = sp.symbols('w1 w2 w3 w4')
A0 = w[1] + w[0]*w[2] + w[0]**2*w[3]
H = sp.Matrix(4, 4, lambda i, j: sp.diff(A0, w[i], w[j]))
gA = sp.Matrix([sp.diff(A0, v) for v in w])
print('A0 =', A0)
print('det Hess_4 A0 =', sp.expand(H.det()))
print('B_4(A0) = ', sp.expand((gA.T*H.adjugate()*gA)[0, 0]))
print('rank Hess_4 A0 =', H.rank())
# affinely 3-variable?  D_v A0 == 0 for some v != 0 ?
vv = sp.symbols('v1 v2 v3 v4')
Dv = sp.expand(sum(vv[i]*sp.diff(A0, w[i]) for i in range(4)))
sol = sp.solve(sp.Poly(Dv, *w).coeffs(), vv, dict=True)
print('solutions of D_v A0 == 0 :', sol)
# homogenization: is it the Perazzo cubic (GN counterexample in 5 vars)?
w0 = sp.Symbol('w0')
E = sp.expand(w0**3*A0.subs({w[i]: w[i]/w0 for i in range(4)}))
print('homogenization E =', sp.simplify(E))
V = [w0] + list(w)
HE = sp.Matrix(5, 5, lambda i, j: sp.diff(E, V[i], V[j]))
print('det Hess_5 E =', sp.expand(HE.det(method='berkowitz')))
DvE = sp.expand(sum(sp.Symbol(f'u{i}')*sp.diff(E, V[i]) for i in range(5)))
solE = sp.solve(sp.Poly(DvE, *V).coeffs(), [sp.Symbol(f'u{i}') for i in range(5)], dict=True)
print('E a cone? D_u E == 0 solutions:', solE)
