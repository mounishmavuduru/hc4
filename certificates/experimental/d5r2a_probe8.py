# EXPERIMENTAL probe: degree-5 rank-2 sub-case (a), the [det]_8 stage.
# f = f5(x1,x2) + A(x1,x2,x3) + x4 B(x1,x2) + f3 + f2, with the [det]_9 constraint
# already imposed (f3_x4x4 = lambda in C[x1,x2], D*lambda = grad(B)^T adj(H5) grad(B)).
# Question: does [det]_8 = 0 (165 equations, 41 unknowns) admit solutions when
# lambda =/= 0?  Two witnesses: lambda = 0 (f5=x1^4x2, B=x1^3) and lambda =/= 0
# (f5=x1^5+x2^5, B=x1^3+x2^3, lambda=(9/20)(x1+x2)).  Solved mod p in Singular.
import itertools, os, random, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
import _d5_close as base
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
eps = sp.Symbol('eps')
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'probe8'); os.makedirs(ND, exist_ok=True)
p = 32003
def hess(g, vs=X): return sp.Matrix(len(vs),len(vs),lambda i,j: sp.diff(g,vs[i],vs[j]))
def graded(poly,d):
    P=sp.Poly(sp.expand(poly),*X)
    return sp.expand(sum(c*sp.prod(v**e for v,e in zip(X,m)) for m,c in P.terms() if sum(m)==d))
def gen_form(vs,d,name):
    mons=[sp.prod(c) for c in itertools.combinations_with_replacement(vs,d)]
    cs=sp.symbols(f'{name}0:{len(mons)}'); return sum(c*m for c,m in zip(cs,mons)), list(cs)
def c8_formula(f5,f4,f3,f2):
    H5=hess(f5,(x1,x2)); H4=hess(f4); H3=hess(f3); H2=hess(f2)
    D=H5.det(); N4=H4[2:,2:]; R3=H3[2:,2:]; R2=H2[2:,2:]
    H4e=H4+eps*H3
    S8=sum(H5[i,j]*sp.diff(H4e.cofactor(i,j),eps).subs(eps,0) for i in range(2) for j in range(2))
    return sp.expand(D*(R3.det()+(N4.adjugate()*R2).trace())+S8+H4.det())
def s_poly(e):
    return base.s_poly(e) if hasattr(base,'s_poly') else str(e).replace('**','^')
def run(name,t,tmo=500):
    fn=os.path.join(ND,name); open(fn,'w',newline='\n').write(t)
    w='/mnt/c/'+os.path.abspath(fn)[3:].replace('\\','/')
    return subprocess.run(base.WSL+['bash','-c',f'timeout {tmo} Singular -q < "{w}" 2>&1'],capture_output=True,text=True).stdout

A,cA=gen_form((x1,x2,x3),4,'A'); f3,cf3=gen_form(X,3,'k'); f2,cf2=gen_form(X,2,'s')
# guard: formula == direct graded piece on a random instance
rng=random.Random(5)
f5g=sum(rng.randint(-3,3)*x1**(5-i)*x2**i for i in range(6)); Bg=sum(rng.randint(-3,3)*x1**(3-i)*x2**i for i in range(4))
inst={s:rng.randint(-3,3) for s in cA+cf3+cf2}
fg=(f5g+A+x4*Bg+f3+f2).subs(inst)
print('guard: c8 formula == [det Hess f]_8 ?', sp.expand(c8_formula(f5g,(A+x4*Bg).subs(inst),f3.subs(inst),f2.subs(inst))-graded(hess(fg).det(),8))==0)

kdict={sp.Poly(f3,*X).as_expr():None}
P3=sp.Poly(f3,*X)
def coef_of(mon):
    for m,c in P3.terms():
        if sp.prod(v**e for v,e in zip(X,m))==mon: return c
for label,f5w,Bw,lam in [('lambda=0', x1**4*x2, x1**3, 0),
                         ('lambda!=0', x1**5+x2**5, x1**3+x2**3, sp.Rational(9,20)*(x1+x2))]:
    # impose f3_x4x4 = lam: coefficients of x3x4^2, x4^3 -> 0; x1x4^2, x2x4^2 -> lam-coeff/2
    fix={coef_of(x3*x4**2):0, coef_of(x4**3):0,
         coef_of(x1*x4**2):sp.Rational(1,2)*lam.coeff(x1) if lam!=0 else 0,
         coef_of(x2*x4**2):sp.Rational(1,2)*lam.coeff(x2) if lam!=0 else 0}
    f3w=f3.subs(fix)
    assert sp.expand(sp.diff(f3w,x4,2)-lam)==0
    Hw=hess(f5w,(x1,x2)); gw=sp.Matrix([sp.diff(Bw,x1),sp.diff(Bw,x2)])
    assert sp.expand(Hw.det()*lam-(gw.T*Hw.adjugate()*gw)[0])==0
    c8=c8_formula(f5w,A+x4*Bw,f3w,f2)
    unk=[s for s in cA+cf3+cf2 if s in c8.free_symbols]
    eqs=[co for co in sp.Poly(c8,*X).coeffs() if co!=0]
    P8=sp.Poly(c8,x4); print(f'{label}: x4-degree of [det]_8 = {P8.degree()}; #eqs={len(eqs)}, #unknowns={len(unk)}')
    # also the x4-top coefficient structure
    txt=[f'ring R={p},({",".join(map(str,unk))}),dp;',
         'ideal I='+',\n'.join(s_poly(e) for e in eqs)+';',
         'ideal S=std(I); "DIM",dim(S); "VDIM",vdim(S);',
         'LIB "primdec.lib"; if(dim(S)>=0){ list L=minAssGTZ(S); "NCOMP",size(L); int i; for(i=1;i<=size(L);i++){ "COMP",i,"dim",dim(std(L[i])); print(L[i]); "ENDCOMP"; } }',
         'quit;']
    out=run(f'{label.replace("!","n")}.sing','\n'.join(txt)+'\n')
    import re
    md=re.search(r'DIM\s+(-?\d+)',out); mn=re.search(r'NCOMP\s+(\d+)',out)
    print(f'   [det]_8 = 0 system mod {p}: dim = {md.group(1) if md else "?"}; ncomp = {mn.group(1) if mn else "?"}')
    for cm in re.findall(r'COMP\s+(\d+)\s+dim\s+(\d+)\s*\n(.*?)ENDCOMP',out,re.S):
        gens=[g.strip().rstrip(',') for g in cm[2].strip().splitlines() if g.strip()]
        print(f'   comp {cm[0]} (dim {cm[1]}): {len(gens)} gens; first: {gens[:6]}')
    if not md: print(out[:1200])
