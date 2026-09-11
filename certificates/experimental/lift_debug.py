# sanity: is the extracted point a genuine survivor (det Hess F == 0 mod p)? and
# are the gauge directions in the weight-9 null space (direct trace computation)?
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
from _d5_close import Y, a5, as_, E_components, coeffs_in_y
y1,y2,y3 = Y; x4 = sp.Symbol('x4'); V4=(y1,y2,y3,x4)
HERE=os.path.dirname(os.path.abspath(__file__))
J=json.load(open(os.path.join(HERE,'survpt','survivor_point.json')))
p=J['p']; C0,C1,C2=J['c']; vals={sp.Symbol(k):int(v)%p for k,v in J['point'].items()}
b3=C2*y1**3+C0*y1**2*y2+C1*y1**2*y3
a5p=sp.expand(a5.subs(vals))
F=a5p+x4*b3+sp.Rational(1,2)*x4**2*y1
def hess(g): return sp.Matrix(4,4,lambda i,j: sp.diff(g,V4[i],V4[j]))
def modp(e):
    P=sp.Poly(sp.expand(e),*V4)
    return sp.Poly({m:(int(c.p)*pow(int(c.q),p-2,p))%p for m,c in P.terms() if (int(c.p)*pow(int(c.q),p-2,p))%p},*V4).as_expr() if P.terms() else 0
HF=hess(F)
d=modp(HF.det(method='berkowitz'))
print('det Hess F mod p == 0 ?', d==0, '' if d==0 else f'(#terms {len(sp.Poly(d,*V4).terms())})')
# E-tower at the point via the library
(E3,E2,E1,E0),A=E_components(a5,b3)
Jg=[q for e in (E3,E2,E1,E0) for q in coeffs_in_y(sp.expand(e)) if q!=0]
print('J at point all zero mod p ?', all(int(sp.expand(q.subs(vals)))%p==0 for q in Jg))
adj=HF.adjugate()
for lab,gexpr in [('d/dy1 F',sp.diff(F,y1)),('y1*d/dx4 F',y1*sp.diff(F,x4)),('y2*d/dx4 F',y2*sp.diff(F,x4))]:
    T=modp(sp.trace(adj*hess(gexpr)))
    print(f'tr(adj(HF) Hess({lab})) mod p == 0 ?', T==0)
# is det Hess F zero over Q at the point? (it is only required mod p)
print('E4 = adj(Hess b3)[0,0] :', sp.expand(sp.Matrix(3,3,lambda i,j: sp.diff(b3,Y[i],Y[j])).adjugate()[0,0]))
