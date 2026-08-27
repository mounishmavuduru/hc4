import sympy as sp
lam = sp.Symbol('lam')
X = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f'X{min(i,j)}{max(i,j)}'))
Y = sp.Matrix(2, 3, lambda i, j: sp.Symbol(f'Y{i}{j}'))
a = sp.Matrix([sp.Symbol(f'a{i}') for i in range(3)])
Z = lam*a*a.T
H = sp.Matrix(sp.BlockMatrix([[X, Y], [Y.T, Z]]))
d = sp.expand(H.det(method='berkowitz'))
print('independent of X block:', all(sp.diff(d, X[i, j]) == 0 for i in range(2) for j in range(2)))
M = sp.Matrix([[a[0], a[1], a[2]], [Y[0, 0], Y[0, 1], Y[0, 2]], [Y[1, 0], Y[1, 1], Y[1, 2]]])
print('det[[X,Y],[Y^T,lam a a^T]] == lam*det[a;Y1;Y2]^2 :',
      sp.expand(d - lam*sp.expand(M.det())**2) == 0)
