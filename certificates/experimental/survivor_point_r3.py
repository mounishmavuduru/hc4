# P5.A1 step 0 (rank-3 component): the generic 6-plane slice of V(J) lands on a
# rank<3 point (det Hess3 a5 = 0 there), so the top-dimensional part of V(J) is
# rank-deficient and the rank-3 solutions (M1: detA not in rad J) live on
# LOWER-dimensional components.  Localize at rank 3 with Rabinowitsch:
#   I = J + (w*t - 1), t = random combination of the coefficients of detA;
# dim(I) =: d, slice by d random affine hyperplanes -> 0-dim -> minAssGTZ ->
# a rational component (if any) -> the point.
import os, re, random, subprocess, sys, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, as_, E_components, coeffs_in_y, AV, s_poly
y1,y2,y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'survpt'); os.makedirs(ND, exist_ok=True)
p = 32003
rng = random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 21)
C0,C1,C2 = (rng.randrange(1,p) for _ in range(3))
b3 = C2*y1**3 + C0*y1**2*y2 + C1*y1**2*y3
(E3,E2,E1,E0),A = E_components(a5,b3)
J = [q for e in (E3,E2,E1,E0) for q in coeffs_in_y(sp.expand(e)) if q != 0]
detA = sp.expand(A.det(method='berkowitz'))
T = [c for c in sp.Poly(detA,*Y).coeffs()]
tcomb = ' + '.join(f'{rng.randrange(1,p)}*({s_poly(g)})' for g in T)
def run(name,t,tmo=560):
    fn=os.path.join(ND,name); open(fn,'w',newline='\n').write(t)
    w='/mnt/c/'+os.path.abspath(fn)[3:].replace('\\','/')
    return subprocess.run(base.WSL+['bash','-c',f'timeout {tmo} Singular -q < "{w}" 2>&1'],capture_output=True,text=True).stdout
head=[f'ring R={p},({AV},w),dp;', 'ideal J='+',\n'.join(s_poly(g) for g in J)+';', f'poly t={tcomb};']
# step 1: dimension of the rank-3 locus
out=run('r3dim.sing','\n'.join(head+['ideal I=J,w*t-1; ideal S=std(I); "DIM",dim(S); quit;'])+'\n')
md=re.search(r'DIM\s+(-?\d+)',out)
if not md: print(out[:1000]); sys.exit(1)
d=int(md.group(1)); print(f'c=({C0},{C1},{C2}); rank-3 locus (Rabinowitsch) dim = {d}')
if d<0: print('rank-3 locus EMPTY for this c (contradicts M1?)'); sys.exit(1)
for attempt in range(3):
    lins=[' + '.join(f'{rng.randrange(1,p)}*{a}' for a in as_)+f' + {rng.randrange(1,p)}' for _ in range(d)]
    txt=head+['LIB "primdec.lib";', 'ideal I=J,w*t-1'+(','+','.join(lins) if lins else '')+';',
              'ideal S=std(I); "DIM",dim(S); "VDIM",vdim(S);',
              'list L=minAssGTZ(S); "NCOMP",size(L);',
              'option(redSB); int i; for(i=1;i<=size(L);i++){ ideal P=std(L[i]); "COMP",i,size(P); int k; for(k=1;k<=size(P);k++){ "GBELEM"; print(P[k]); "ENDELEM"; } "ENDCOMP"; kill P; kill k; }',
              'quit;']
    out=run('r3pt.sing','\n'.join(txt)+'\n')
    md=re.search(r'DIM\s+(-?\d+)',out); mv=re.search(r'VDIM\s+(-?\d+)',out); mn=re.search(r'NCOMP\s+(\d+)',out)
    print(f'attempt {attempt}: dim={md.group(1) if md else "?"} vdim={mv.group(1) if mv else "?"} ncomp={mn.group(1) if mn else "?"}')
    if not mn: print(out[:1000]); continue
    syms={str(a):a for a in as_}; wsym=sp.Symbol('w'); syms['w']=wsym; allv=list(as_)+[wsym]
    for comp in re.findall(r'COMP\s+\d+\s+\d+\s*\n(.*?)ENDCOMP',out,re.S):
        gens=[g.strip().replace('\n','') for g in re.findall(r'GBELEM\s*\n(.*?)ENDELEM',comp,re.S)]
        polys=[sp.sympify(g.replace('^','**'),locals=syms) for g in gens if g]
        if not polys or not all(sp.Poly(q,*allv).total_degree()==1 for q in polys):
            print(f'   component with {len(polys)} gens, degrees {sorted(set(sp.Poly(q,*allv).total_degree() for q in polys))[:4]} -> not rational'); continue
        M=[[int(sp.Poly(q,*allv).coeff_monomial(a))%p for a in allv]+[int(sp.Poly(q,*allv).coeff_monomial(1))%p] for q in polys]
        nv=len(allv); r=0; piv=[]
        for c in range(nv):
            k=next((i for i in range(r,len(M)) if M[i][c]),None)
            if k is None: continue
            M[r],M[k]=M[k],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[(v*inv)%p for v in M[r]]
            for i in range(len(M)):
                if i!=r and M[i][c]: f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
            piv.append(c); r+=1
        if r<nv: print(f'   rational component but rank {r}<{nv}'); continue
        vals={allv[c]:(-M[i][nv])%p for i,c in enumerate(piv)}
        avals={a:vals[a] for a in as_}
        dA=sp.expand(detA.subs(avals)); nz=dA!=0 and any(int(co)%p for co in sp.Poly(dA,*Y).coeffs())
        bad=[i for i,q in enumerate(J) if int(sp.expand(q.subs(avals)))%p]
        print(f'   RATIONAL rank-3 point: det Hess3(a5) nonzero mod p? {nz}; J vanishes? {not bad}')
        if nz and not bad:
            tag = f'_{sys.argv[1]}' if len(sys.argv)>1 else ''
            json.dump({'p':p,'c':[C0,C1,C2],'point':{str(k):int(v) for k,v in avals.items()}},
                      open(os.path.join(ND,f'survivor_point{tag}.json'),'w'))
            print(f'saved survpt/survivor_point{tag}.json'); sys.exit(0)
print('no rational rank-3 point found'); sys.exit(1)
