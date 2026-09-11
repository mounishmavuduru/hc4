# P5.A1, mod p, FIXED c but the survivor point SYMBOLIC: is
#   J(a) + tower(a, F4, F3, F2) + (w*t(a) - 1) + (v*c0 - 1)   the unit ideal?
# (t = random combination of the coefficients of det Hess3 a5 -> rank-3
# localisation; c0 = weight-0 coefficient of det Hess f.)  If (1): NO rank-3
# survivor with this b3 lifts to a constant-Hessian quintic mod p -- a mod-p
# theorem for this c rather than a point sample.  Heavy: 21 + 41 + 2 unknowns.
import os, re, sys, json, subprocess, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, as_, form, E_components, coeffs_in_y, AV
y1,y2,y3 = Y; x4 = sp.Symbol('x4')
HERE = os.path.dirname(os.path.abspath(__file__)); ND = os.path.join(HERE,'survpt')
seed = int(sys.argv[2]) if len(sys.argv)>2 else 21
p = 32003; rng = random.Random(seed)
C0,C1,C2 = (rng.randrange(1,p) for _ in range(3))
half = (p+1)//2
b3 = C2*y1**3 + C0*y1**2*y2 + C1*y1**2*y3
(E3,E2,E1,E0),A = E_components(a5,b3)
Jg = [q for e in (E3,E2,E1,E0) for q in coeffs_in_y(sp.expand(e)) if q != 0]
detA = sp.expand(A.det(method='berkowitz'))
T = sp.Poly(detA,*Y).coeffs()
tcomb = ' + '.join(f'{rng.randrange(1,p)}*({base.s_poly(g)})' for g in T)
def s(e): return base.s_poly(sp.expand(e))
a4,P = form(Y,4,'p'); b2,Q = form(Y,2,'q'); t0 = sp.Symbol('t0')
a3,Rr = form(Y,3,'r'); b1,S = form(Y,1,'s'); a2,U = form(Y,2,'u')
unk = list(P)+list(Q)+[t0]+list(Rr)+list(S)+list(U)
F  = a5 + x4*b3 + half*x4**2*y1
f  = F + a4 + x4*b2 + half*x4**2*t0 + a3 + x4*b1 + a2
txt = [f'ring R={p},({AV},{",".join(map(str,unk))},w,v,y1,y2,y3,x4),dp;',
       'poly f=' + s(f) + ';',
       'matrix H[4][4]; int i; int j; list V=y1,y2,y3,x4;',
       'for(i=1;i<=4;i++){ for(j=1;j<=4;j++){ H[i,j]=diff(diff(f,V[i]),V[j]); } }',
       'poly d=det(H); "DETTERMS",size(d);',
       'matrix C=coef(d,y1*y2*y3*x4); int n=ncols(C); "NMONO",n;',
       'ideal I; poly cc=0; int k;',
       'for(k=1;k<=n;k++){ if(C[1,k]==1){ cc=C[2,k]; } else { I=I,C[2,k]; } }',
       'I=I,v*cc-1;',
       'ideal J=' + ',\n'.join(base.s_poly(g) for g in Jg) + ';',
       f'poly t={tcomb};',
       'I=I,J,w*t-1;',
       f'ring T2={p},({AV},{",".join(map(str,unk))},w,v),dp;',
       'ideal I=imap(R,I); I=simplify(I,2); "NEQ",size(I);',
       'ideal G=std(I); "DIM",dim(G);',
       'quit;']
fn = os.path.join(ND,f'tower_sym_{seed}.sing'); open(fn,'w',newline='\n').write('\n'.join(txt)+'\n')
w = '/mnt/c/'+os.path.abspath(fn)[3:].replace('\\','/')
tmo = int(sys.argv[1]) if len(sys.argv)>1 else 570
out = subprocess.run(base.WSL+['bash','-c',f'timeout {tmo} Singular -q < "{w}" 2>&1'],capture_output=True,text=True).stdout
print(f'c=({C0},{C1},{C2}) mod {p}; symbolic survivor (21 a-vars) + 41 lower unknowns + 2 Rabinowitsch')
for key in ('DETTERMS','NMONO','NEQ','DIM'):
    m = re.search(key+r'\s+(-?\d+)', out); print(f'  {key} = {m.group(1) if m else "?"}')
if 'DIM' not in out: print('  (no DIM: timeout or error)'); print(out[-600:])
