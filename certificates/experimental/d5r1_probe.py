# probe: degree-5 RANK-1 branch (f5 = x1^5). First graded reduction.
# Expect [det Hess f]_9 = 20*x1^3 * det3(Hess_{x2,x3,x4} f4)  (only the (1,1)
# entry 20x1^3 times the {2,3,4} minor of the f4-block reaches degree 9).
import itertools, random
import sympy as sp
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
def hess(g,vs): return sp.Matrix(len(vs),len(vs),lambda i,j: sp.diff(g,vs[i],vs[j]))
def topdeg(poly,d):
    P=sp.Poly(sp.expand(poly),*X)
    return sp.expand(sum(c*sp.prod(v**e for v,e in zip(X,m)) for m,c in P.terms() if sum(m)==d))
rng=random.Random(5)
def rpoly(vs,d,lo=-3,hi=3):
    mons=[sp.prod(t) for t in itertools.combinations_with_replacement(vs,d)]
    return sum(rng.randint(lo,hi)*m for m in mons)
ok=True
for trial in range(3):
    f5=x1**5
    f4=rpoly(X,4); f3=rpoly(X,3)
    f=f5+f4+f3
    D=sp.expand(hess(f,X).det())
    d9=topdeg(D,9)
    claim=sp.expand(20*x1**3*hess(f4,(x2,x3,x4)).det())
    m=sp.expand(d9-claim)==0
    ok&=m
    print(f'r=1 trial {trial}: [det]_9 == 20*x1^3*det3(Hess_(x2,x3,x4) f4)?  {m}')
    if not m: print('   residual:', sp.factor(sp.expand(d9-claim)))
    # also: is the degree-10 piece zero (no higher piece)?
    print(f'   [det]_10 == 0? {topdeg(D,10)==0}')
print('ALL MATCH:',ok)
