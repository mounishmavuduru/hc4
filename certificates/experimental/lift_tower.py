# P5.A1 (the bottleneck), EXPERIMENTAL, mod p, at ONE explicit rank-3 survivor
# point F = a5 + x4 b3 + (1/2) x4^2 y1 (survpt/survivor_point.json):
# does F lift to a full quintic f = F + F4 + F3 + F2 with det Hess f = c =/= 0 ?
# Unknowns (41): F4 = a4(15) + x4 b2(6) + (1/2) x4^2 t0(1); F3 = a3(10) + x4 b1(3);
# F2 = a2(6) (x4*b0 has zero Hessian, F1/F0 too).  Weighted grading wt y = 1,
# wt x4 = 2: every positive-weight coefficient of det Hess f must vanish and the
# weight-0 coefficient c must be a unit (Rabinowitsch w*c - 1).  Singular does
# the determinant and the Groebner basis.  A CONSISTENT system would give an
# explicit constant-Hessian quintic mod p to examine; an INCONSISTENT one is
# one-point, one-prime evidence that this survivor does not lift.
import os, re, sys, json, subprocess, itertools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, as_, form
y1,y2,y3 = Y; x4 = sp.Symbol('x4')
HERE = os.path.dirname(os.path.abspath(__file__)); ND = os.path.join(HERE,'survpt')
TAG = f'_{sys.argv[2]}' if len(sys.argv)>2 else ''
J = json.load(open(os.path.join(ND,f'survivor_point{TAG}.json')))
p = J['p']; C0,C1,C2 = J['c']
vals = {sp.Symbol(k): int(v) % p for k,v in J['point'].items()}
half = (p+1)//2
b3 = C2*y1**3 + C0*y1**2*y2 + C1*y1**2*y3
a5p = sp.expand(a5.subs(vals))
def s(e): return base.s_poly(sp.expand(e))
a4,P = form(Y,4,'p'); b2,Q = form(Y,2,'q'); t0 = sp.Symbol('t0')
a3,Rr = form(Y,3,'r'); b1,S = form(Y,1,'s'); a2,U = form(Y,2,'u')
unk = list(P)+list(Q)+[t0]+list(Rr)+list(S)+list(U)
F  = a5p + x4*b3 + half*x4**2*y1
f  = F + a4 + x4*b2 + half*x4**2*t0 + a3 + x4*b1 + a2
V4 = [y1,y2,y3,x4]
txt = [f'ring R={p},({",".join(map(str,unk))},w,y1,y2,y3,x4),dp;',
       'poly f=' + s(f) + ';',
       'matrix H[4][4];',
       'int i; int j;',
       'list V=y1,y2,y3,x4;',
       'for(i=1;i<=4;i++){ for(j=1;j<=4;j++){ H[i,j]=diff(diff(f,V[i]),V[j]); } }',
       'poly d=det(H);',
       '"DETTERMS",size(d);',
       'matrix C=coef(d,y1*y2*y3*x4);',
       'int n=ncols(C); "NMONO",n;',
       'ideal I; poly cc=0; int k;',
       'for(k=1;k<=n;k++){ if(C[1,k]==1){ cc=C[2,k]; } else { I=I,C[2,k]; } }',
       'I=I,w*cc-1;',
       f'ring T={p},({",".join(map(str,unk))},w),dp;',
       'ideal I=imap(R,I); I=simplify(I,2); "NEQ",size(I);',
       'ideal G=std(I); "DIM",dim(G);',
       'if(dim(G)==0){ "VDIM",vdim(G); }',
       'if(dim(G)>=0){ LIB "primdec.lib"; list L=minAssGTZ(G); "NCOMP",size(L); int m; for(m=1;m<=size(L);m++){ "COMP",m,"dim",dim(std(L[m])); } }',
       'quit;']
fn = os.path.join(ND,'tower.sing'); open(fn,'w',newline='\n').write('\n'.join(txt)+'\n')
w = '/mnt/c/'+os.path.abspath(fn)[3:].replace('\\','/')
tmo = int(sys.argv[1]) if len(sys.argv)>1 else 570
out = subprocess.run(base.WSL+['bash','-c',f'timeout {tmo} Singular -q < "{w}" 2>&1'],capture_output=True,text=True).stdout
print(f'survivor point c=({C0},{C1},{C2}) mod {p}; 41 unknowns')
for key in ('DETTERMS','NMONO','NEQ','DIM','VDIM','NCOMP'):
    m = re.search(key+r'\s+(-?\d+)', out); print(f'  {key} = {m.group(1) if m else "?"}')
for cm in re.findall(r'COMP\s+(\d+)\s+dim\s+(-?\d+)', out): print(f'  component {cm[0]}: dim {cm[1]}')
if 'DIM' not in out: print('  (no DIM: timeout or error)'); print(out[-800:])
