# verify the DISCOVERED structure of [det]_9 in rank-2 sub-case (a) (u=e3 const):
#   [det]_9 = l * ( det2(Hess2 f5) * f3_x4x4 + S7 ),   l = A_x3x3,
# where S7 in C[x1,x2] (degree 7) is the Laplace cross-term from the (1,4),(2,4)
# column pairs, built from Hess2 f5 and B_1, B_2 (B = the x4-coefficient of f4).
import itertools, random
import sympy as sp
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
def hess(g,vs): return sp.Matrix(len(vs),len(vs),lambda i,j: sp.diff(g,vs[i],vs[j]))
def topdeg(poly,d):
    P=sp.Poly(sp.expand(poly),*X)
    return sp.expand(sum(c*sp.prod(v**e for v,e in zip(X,m)) for m,c in P.terms() if sum(m)==d))
rng=random.Random(11)
def rpoly(vs,d,lo=-3,hi=3):
    mons=[sp.prod(t) for t in itertools.combinations_with_replacement(vs,d)]
    return sum(rng.randint(lo,hi)*m for m in mons)
ok=True
for trial in range(2):
    f5=sum(rng.randint(-4,4)*x1**(5-i)*x2**i for i in range(6))
    f3=rpoly(X,3); A=rpoly((x1,x2,x3),4); B=rpoly((x1,x2),3)
    f=f5+A+x4*B+f3
    d9=topdeg(hess(f,X).det(),9)
    ell=sp.expand(sp.diff(A,x3,2))
    q,r=sp.div(sp.Poly(d9,*X), sp.Poly(ell,*X))
    div=(r.is_zero)
    print(f'trial {trial}: [det]_9 divisible by l=A_x3x3? {div}')
    if div:
        Qp=sp.expand(q.as_expr())
        core=sp.expand(Qp - sp.expand(hess(f5,(x1,x2)).det())*sp.diff(f3,x4,2))
        # the remaining S7 must be in C[x1,x2] only, degree 7
        onlyx12 = not (x3 in core.free_symbols or x4 in core.free_symbols)
        deg = sp.Poly(core,*X).total_degree() if core!=0 else None
        print(f'   S7 := [det]_9/l - det2Hess2f5*f3_x4x4 : in C[x1,x2] only? {onlyx12}, degree {deg}')
        # S7 should be built from Hess2 f5 and B_1,B_2 only: check it is invariant when f3, A are changed (with same f5,B)
        f3b=rpoly(X,3); Ab=rpoly((x1,x2,x3),4)
        fb=f5+Ab+x4*B+f3b
        ellb=sp.expand(sp.diff(Ab,x3,2))
        qb,rb=sp.div(sp.Poly(topdeg(hess(fb,X).det(),9),*X), sp.Poly(ellb,*X))
        coreb=sp.expand(qb.as_expr() - sp.expand(hess(f5,(x1,x2)).det())*sp.diff(f3b,x4,2))
        print(f'   S7 depends only on (f5,B) [same for different A,f3]? {sp.expand(core-coreb)==0}')
        ok&=div and onlyx12 and sp.expand(core-coreb)==0
    else: ok=False
print('STRUCTURE CONFIRMED:',ok)
