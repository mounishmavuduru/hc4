# scratch: master formula  det Hess_5(<b,z> + (lam/2)<a,z>^2 + mu<a,z>) = lam * (Jac F)^2|_{x3 = lam<a,z>+mu}
# with F(x1,x2,x3) = a(x1,x2) x3 + b(x1,x2).   Fully generic a,b => proof.
import sympy as sp

x1, x2 = sp.symbols('x1 x2')
z = sp.symbols('z1 z2 z3')
lam, mu = sp.symbols('lam mu')
X = [x1, x2]

def genpoly(name, deg):
    n = (deg+1)*(deg+2)//2
    c = [sp.Symbol(f'{name}_{j}') for j in range(n)]
    t = 0
    k = 0
    for d in range(deg+1):
        for i in range(d+1):
            t += c[k]*x1**i*x2**(d-i)
            k += 1
    return sp.expand(t)

for deg in (1, 2):
    a = [genpoly(f'a{i}', deg) for i in range(3)]
    b = [genpoly(f'b{i}', deg) for i in range(3)]
    az = sum(a[i]*z[i] for i in range(3))
    bz = sum(b[i]*z[i] for i in range(3))
    Th = sp.expand(bz + sp.Rational(1, 2)*lam*az**2 + mu*az)
    W = [x1, x2, z[0], z[1], z[2]]
    H = sp.Matrix(5, 5, lambda i, j: sp.diff(Th, W[i], W[j]))
    lhs = sp.expand(H.det(method='berkowitz'))
    x3 = sp.Symbol('x3')
    F = [sp.expand(a[i]*x3 + b[i]) for i in range(3)]
    JF = sp.Matrix(3, 3, lambda i, j: sp.diff(F[i], [x1, x2, x3][j]))
    jac = sp.expand(JF.det(method='berkowitz'))
    rhs = sp.expand(lam*(jac.subs(x3, lam*az + mu))**2)
    print(f'deg {deg}: master formula holds?', sp.expand(lhs - rhs) == 0)
