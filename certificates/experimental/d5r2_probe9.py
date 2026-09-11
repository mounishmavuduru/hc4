# probe: degree-5 rank-2, the [det Hess f]_9 graded piece under N4 = l*u*u^T.
import itertools, random
import sympy as sp
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
def hess(g): return sp.Matrix(4,4,lambda i,j: sp.diff(g,X[i],X[j]))
def topdeg(poly,d):
    P=sp.Poly(sp.expand(poly),*X)
    return sp.expand(sum(c*sp.prod(v**e for v,e in zip(X,m)) for m,c in P.terms() if sum(m)==d))
def det2(g,vs): return sp.expand(sp.Matrix(2,2,lambda i,j: sp.diff(g,vs[i],vs[j])).det())
rng=random.Random(11)
def rpoly(vs,d,lo=-3,hi=3):
    mons=[sp.prod(t) for t in itertools.combinations_with_replacement(vs,d)]
    return sum(rng.randint(lo,hi)*m for m in mons)

mons5=[x1**(5-i)*x2**i for i in range(6)]
for trial in range(2):
    f5=sum(rng.randint(-4,4)*m for m in mons5)
    f3=rpoly(X,3)
    # ---- sub-case (a): u = e3 constant: f4 = A(x1,x2,x3) + x4*B(x1,x2)
    A=rpoly((x1,x2,x3),4); B=rpoly((x1,x2),3)
    f4a=A + x4*B
    f=f5+f4a+f3
    d9=topdeg(hess(f).det(),9)
    ell=sp.expand(sp.diff(A,x3,2))                 # l = d^2 A / dx3^2 (quadratic)
    f3_44=sp.expand(sp.diff(f3,x4,2))               # d^2 f3 / dx4^2 (linear)
    claim=sp.expand(det2(f5,(x1,x2))*ell*f3_44)
    print(f'(a) trial {trial}: [det]_9 == det2(Hess2 f5)*(A_x3x3)*(f3_x4x4)?  {sp.expand(d9-claim)==0}')
    if sp.expand(d9-claim)!=0:
        resid=sp.expand(d9-claim)
        print('    residual factor:', sp.factor(resid))
    # ---- sub-case (b): u = (L1,L2) LINEAR in (x1,x2): f4 = (al/2)(L1 x3 + L2 x4)^2 + C x3 + D x4 + E
    L1=rng.randint(1,3)*x1+rng.randint(-3,3)*x2; L2=rng.randint(-3,3)*x1+rng.randint(1,3)*x2
    al=rng.randint(1,3)
    C=rpoly((x1,x2),3); D=rpoly((x1,x2),3); E=rpoly((x1,x2),4)
    f4b=sp.Rational(al,2)*(L1*x3+L2*x4)**2 + C*x3 + D*x4 + E
    # sanity: N4 = al*(L1,L2)(L1,L2)^T
    N4=sp.Matrix(2,2,lambda i,j: sp.diff(f4b,(x3,x4)[i],(x3,x4)[j]))
    assert sp.expand(N4[0,0]-al*L1**2)==0 and sp.expand(N4[0,1]-al*L1*L2)==0 and sp.expand(N4[1,1]-al*L2**2)==0
    fb=f5+f4b+f3
    d9b=topdeg(hess(fb).det(),9)
    print(f'(b) trial {trial}: [det]_9 degree-9 piece, #terms={len(sp.Poly(d9b,*X).terms()) if d9b!=0 else 0}')
    print('    factor([det]_9) =', sp.factor(d9b))
