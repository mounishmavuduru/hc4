# EXPERIMENTAL (mod p): the full ordinary-degree tower for a degree-5 potential
# with a prescribed leading structure, solved for the lower-degree parts.
#   f = f5 + f4 + f3 + f2 (+ affine, irrelevant);  det Hess f = c  <=>
#   every positive-degree coefficient of det Hess f vanishes and the constant
#   coefficient is a unit (Rabinowitsch w*c - 1).
# Modes (all rank <= 2 leading forms; see ledger R-D5-R2):
#   r2a   sub-case (a): f5, B fixed (witness), lambda = f3_x4x4 imposed;
#         unknowns A(15) + f3(16) + f2(10).  --lam0 uses the lambda = 0 witness.
#   r2c   N4 = 0: f4 = x3 C + x4 D' + E random; unknowns f3(20) + f2(10).
#   r1    f5 = x1^5, f4 random in C[x1,x2,x3] (so det3 Hess_{234} f4 == 0);
#         unknowns f3(20) + f2(10).
# A CONSISTENT system is examined: dimension, a rational point if one exists,
# independent sympy re-check of det Hess f = c, and the PIVOT CONE
# {v : D_v^2 f constant} (dim <= 0 <=> pivot-free).
import os, re, sys, json, random, itertools, subprocess, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
import _d5_close as base
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
ap = argparse.ArgumentParser()
ap.add_argument('--mode', default='r2a'); ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--p', type=int, default=32003); ap.add_argument('--timeout', type=int, default=560)
ap.add_argument('--lam0', action='store_true')
ap.add_argument('--incremental', action='store_true', help='impose [det]_k = 0 for k >= m, m = 8..1, report the first m with dim = -1')
ap.add_argument('--x4incr', action='store_true', help='stage by x4-degree: impose the x4^j coefficients for j >= m, m = top..1, then everything + unit')
args = ap.parse_args()
p, rng = args.p, random.Random(args.seed)
HERE = os.path.dirname(os.path.abspath(__file__)); ND = os.path.join(HERE, 'tower'); os.makedirs(ND, exist_ok=True)
def form(vs, d, name):
    mons = [sp.prod(c) for c in itertools.combinations_with_replacement(vs, d)]
    cs = sp.symbols(f'{name}0:{len(mons)}'); return sum(c*m for c, m in zip(cs, mons)), list(cs)
def rform(vs, d): return sum(rng.randrange(1, p)*sp.prod(c) for c in itertools.combinations_with_replacement(vs, d))
def s(e): return base.s_poly(sp.expand(e))
def hess(g): return sp.Matrix(4, 4, lambda i, j: sp.diff(g, X[i], X[j]))

f3, K = form(X, 3, 'k'); f2, S2 = form(X, 2, 's')
fixed = {}
if args.mode == 'r2a':
    if args.lam0: f5, B, lam = x1**4*x2, x1**3, 0
    else:         f5, B, lam = x1**5 + x2**5, x1**3 + x2**3, sp.Rational(9, 20)*(x1 + x2)
    A, CA = form((x1, x2, x3), 4, 'A')
    f4 = A + x4*B
    P3 = sp.Poly(f3, *X)
    def coef_of(mon):
        for m, c in P3.terms():
            if sp.prod(v**e for v, e in zip(X, m)) == mon: return c
    fixed = {coef_of(x3*x4**2): 0, coef_of(x4**3): 0,
             coef_of(x1*x4**2): sp.Rational(1, 2)*lam.coeff(x1) if lam != 0 else 0,
             coef_of(x2*x4**2): sp.Rational(1, 2)*lam.coeff(x2) if lam != 0 else 0}
    f3 = f3.subs(fixed); unk = CA + [k for k in K if k not in fixed] + S2
    Hw = sp.Matrix(2, 2, lambda i, j: sp.diff(f5, (x1, x2)[i], (x1, x2)[j])); gw = sp.Matrix([sp.diff(B, x1), sp.diff(B, x2)])
    assert sp.expand(Hw.det()*lam - (gw.T*Hw.adjugate()*gw)[0]) == 0
    desc = f'r2a witness f5={f5}, B={B}, lambda={lam}'
elif args.mode == 'r2c':
    f5 = rform((x1, x2), 5); f4 = x3*rform((x1, x2), 3) + x4*rform((x1, x2), 3) + rform((x1, x2), 4)
    unk = K + S2; desc = 'r2c: N4 = 0, random f5, C, D\', E'
elif args.mode == 'r1':
    f5 = x1**5; f4 = rform((x1, x2, x3), 4)
    assert sp.expand(sp.Matrix(3, 3, lambda i, j: sp.diff(f4, (x2, x3, x4)[i], (x2, x3, x4)[j])).det()) == 0
    unk = K + S2; desc = 'r1: f5 = x1^5, random f4 in C[x1,x2,x3]'
elif args.mode == 'r2c_np':
    # N4 = 0 with f4 = x3 C + x4 D' + E UNKNOWN (f5 random) and NO pivot in span(e3,e4):
    # D_v^2 f = v^T (R3 + R2) v is constant iff v is a common isotropic vector of the
    # constant matrices S_k (R3 = sum_k x_k S_k); generically none iff the binary
    # quadratics q_k = v^T S_k v span all of Sym2 -> impose a random 3x3 minor != 0.
    f5 = rform((x1, x2), 5)
    Cc, CC = form((x1, x2), 3, 'c'); Dd, CD = form((x1, x2), 3, 'd'); Ee, CE = form((x1, x2), 4, 'e')
    f4 = x3*Cc + x4*Dd + Ee
    unk = CC + CD + CE + K + S2; desc = 'r2c_np: N4 = 0, f5 random, C/D\'/E/f3/f2 unknown, no pivot in span(e3,e4)'
elif args.mode in ('r2c_np1', 'r2c_np0'):
    # N4 = 0, f4 = x3 C + x4 D' + E unknown, f5 random.  The tower forces the pure
    # (x3,x4)-part of f3 to be a cube (gamma/6)(x3 + u2 x4)^3 (chart u = (1,u2));
    # then the ONLY candidate pivot direction in span(e3,e4) is w = (u2,-1) and
    # pivot-free <=> mu := w^T R3' w != 0 (R3' = x1 S1 + x2 S2)   [mode r2c_np1, gamma != 0]
    # or, if the cube is absent (S3 = S4 = 0), pivot-free <=> Res(q1,q2) != 0,
    # q_k(v) = v^T S_k v                                          [mode r2c_np0].
    f5 = rform((x1, x2), 5)
    Cc, CC = form((x1, x2), 3, 'c'); Dd, CD = form((x1, x2), 3, 'd'); Ee, CE = form((x1, x2), 4, 'e')
    f4 = x3*Cc + x4*Dd + Ee
    P3 = sp.Poly(f3, *X)
    def coef_of(mon):
        for m, c in P3.terms():
            if sp.prod(v**e for v, e in zip(X, m)) == mon: return c
    pure = [coef_of(x3**3), coef_of(x3**2*x4), coef_of(x3*x4**2), coef_of(x4**3)]
    if args.mode == 'r2c_np1':
        gam, u2 = sp.symbols('gam u2')
        cube = sp.expand(gam/6*(x3 + u2*x4)**3)
        sub = {c: sp.Poly(cube, *X).coeff_monomial(m) for c, m in zip(pure, [x3**3, x3**2*x4, x3*x4**2, x4**3])}
        f3 = f3.subs(sub); unk = CC + CD + CE + [gam, u2] + [k for k in K if k not in sub] + S2
        desc = 'r2c_np1: N4 = 0, f5 random, f4/f3/f2 unknown, f3\'\' = (gam/6)(x3+u2 x4)^3, w = (u2,-1) NOT a pivot (mu != 0)'
    else:
        sub = {c: 0 for c in pure}
        f3 = f3.subs(sub); unk = CC + CD + CE + [k for k in K if k not in sub] + S2
        desc = 'r2c_np0: N4 = 0, f5 random, f4/f3/f2 unknown, f3\'\' = 0, Res(q1,q2) != 0 (no isotropic v)'
elif args.mode == 'r1_np':
    # f5 = x1^5, f4 in C[x1,x2,x3] UNKNOWN, f3, f2 unknown; pivots need v1 = 0 and
    # (generically) v = e4, so "no pivot" = f3_x4x4 != 0 (Rabinowitsch on a random
    # combination of its coefficients).
    f5 = x1**5; A, CA = form((x1, x2, x3), 4, 'A'); f4 = A
    unk = CA + K + S2; desc = 'r1_np: f5 = x1^5, f4 in C[x1,x2,x3] unknown, f3/f2 unknown, e4 not a pivot'
else: raise SystemExit('mode?')
extra_rab = None      # polynomial in the unknowns that must be a UNIT (Rabinowitsch variable w2)
if args.mode == 'r2c_np':
    R3 = sp.Matrix(2, 2, lambda i, j: sp.diff(f3, (x3, x4)[i], (x3, x4)[j]))
    Sk = [R3.applyfunc(lambda e: sp.Poly(e, *X).coeff_monomial(xk)) for xk in X]      # constant 2x2 symmetric matrices
    rows = [[S[0, 0], 2*S[0, 1], S[1, 1]] for S in Sk]                                # q_k in the basis v3^2, v3v4, v4^2
    pick = rng.sample(range(4), 3)
    extra_rab = sp.expand(sp.Matrix([rows[i] for i in pick]).det())
if args.mode == 'r2c_np1':
    R3 = sp.Matrix(2, 2, lambda i, j: sp.diff(f3, (x3, x4)[i], (x3, x4)[j]))
    wv = sp.Matrix([u2, -1]); mu = sp.expand((wv.T*R3*wv)[0])
    assert not (x3 in mu.free_symbols or x4 in mu.free_symbols), 'w not isotropic for the cube part'
    extra_rab = sp.expand(gam*sum(rng.randrange(1, p)*c for c in sp.Poly(mu, *X).coeffs()))
if args.mode == 'r2c_np0':
    R3 = sp.Matrix(2, 2, lambda i, j: sp.diff(f3, (x3, x4)[i], (x3, x4)[j]))
    v3, v4 = sp.symbols('v3 v4')
    q = [sp.expand((sp.Matrix([v3, v4]).T*R3.applyfunc(lambda e: sp.Poly(e, *X).coeff_monomial(xk))*sp.Matrix([v3, v4]))[0]) for xk in (x1, x2)]
    extra_rab = sp.expand(sp.resultant(q[0].subs(v4, 1), q[1].subs(v4, 1), v3))
if args.mode == 'r1_np':
    f344 = sp.diff(f3, x4, 2)
    extra_rab = sp.expand(sum(rng.randrange(1, p)*c for c in sp.Poly(f344, *X).coeffs()))
f = f5 + f4 + f3 + f2
# rationals in f (1/2 from lambda) -> integers mod p
def s_modp(e):
    P = sp.Poly(sp.expand(e), *X, *unk)
    terms = []
    for m, c in P.terms():
        r = sp.Rational(c); cc = (int(r.p)*pow(int(r.q), p-2, p)) % p
        if cc: terms.append(f'{cc}*' + '*'.join(f'{v}^{e}' for v, e in zip(list(X)+unk, m) if e) if any(m) else str(cc))
    return ' + '.join(terms) if terms else '0'
txt = [f'ring R={p},({",".join(map(str,unk))},w,w2,x1,x2,x3,x4),dp;',
       'poly f=' + s_modp(f) + ';',
       ('poly rab2=' + s_modp(extra_rab) + ';') if extra_rab is not None else 'poly rab2=1;',
       'matrix H[4][4]; int i; int j; list V=x1,x2,x3,x4;',
       'for(i=1;i<=4;i++){ for(j=1;j<=4;j++){ H[i,j]=diff(diff(f,V[i]),V[j]); } }',
       'poly d=det(H); "DETTERMS",size(d);',
       'matrix C=coef(d,x1*x2*x3*x4); int n=ncols(C); "NMONO",n;',
       'ideal I; poly cc=0; int k;',
       'for(k=1;k<=n;k++){ if(C[1,k]==1){ cc=C[2,k]; } else { I=I,C[2,k]; } }',
       'I=I,w*cc-1,w2*rab2-1;']
if args.x4incr:
    txt += ['matrix C4=coef(d,x4); int n4=ncols(C4); "X4TOP",deg(C4[1,1]);',
            'int jj; ideal Ik; ideal Gk; int m4; matrix Cj;',
            'for(m4=deg(C4[1,1]); m4>=1; m4--){ Ik=0; for(jj=1;jj<=n4;jj++){ if(deg(C4[1,jj])>=m4){ Cj=coef(C4[2,jj],x1*x2*x3); for(k=1;k<=ncols(Cj);k++){ Ik=Ik,Cj[2,k]; } } } Ik=simplify(Ik,2); Gk=std(Ik); "X4DEG",m4,"NEQK",size(Ik),"DIMK",dim(Gk); if(dim(Gk)==-1){ break; } }',
            'if(dim(Gk)!=-1){ Ik=Ik,I; Ik=simplify(Ik,2); Gk=std(Ik); "X4DEG",0,"NEQK",size(Ik),"DIMK",dim(Gk); }',
            'quit;']
    fn = os.path.join(ND, f'{args.mode}_{args.seed}_{p}{"_lam0" if args.lam0 else ""}_x4incr.sing'); open(fn, 'w', newline='\n').write('\n'.join(txt)+'\n')
    wp = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    out = subprocess.run(base.WSL+['bash', '-c', f'timeout {args.timeout} Singular -q < "{wp}" 2>&1'], capture_output=True, text=True).stdout
    print(f'{desc}; p={p}, seed={args.seed}; unknowns={len(unk)}; x4-DEGREE staged scan')
    for mm in re.findall(r'X4DEG\s+(\d+)\s+NEQK\s+(\d+)\s+DIMK\s+(-?\d+)', out):
        print(f'  x4^j coefficients = 0 for j >= {mm[0]}' + (' (+ all, + unit)' if mm[0] == '0' else '') + f': {mm[1]} equations, dim = {mm[2]}' + ('  <== first inconsistency' if mm[2] == '-1' else ''))
    if 'X4DEG' not in out: print(out[-500:])
    sys.exit(0)
if args.incremental:
    txt += ['int mind; ideal Ik; ideal Gk;',
            'for(mind=8;mind>=1;mind--){ Ik=0; for(k=1;k<=n;k++){ if(deg(C[1,k])>=mind){ Ik=Ik,C[2,k]; } } Ik=simplify(Ik,2); Gk=std(Ik); "MINDEG",mind,"NEQK",size(Ik),"DIMK",dim(Gk); if(dim(Gk)==-1){ break; } }',
            'quit;']
    fn = os.path.join(ND, f'{args.mode}_{args.seed}_{p}{"_lam0" if args.lam0 else ""}_incr.sing'); open(fn, 'w', newline='\n').write('\n'.join(txt)+'\n')
    wp = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    out = subprocess.run(base.WSL+['bash', '-c', f'timeout {args.timeout} Singular -q < "{wp}" 2>&1'], capture_output=True, text=True).stdout
    print(f'{desc}; p={p}, seed={args.seed}; unknowns={len(unk)}; INCREMENTAL depth scan')
    for mm in re.findall(r'MINDEG\s+(\d+)\s+NEQK\s+(\d+)\s+DIMK\s+(-?\d+)', out):
        print(f'  [det]_k = 0 for k >= {mm[0]}: {mm[1]} equations, dim = {mm[2]}' + ('  <== first inconsistency' if mm[2] == '-1' else ''))
    if 'MINDEG' not in out: print(out[-500:])
    sys.exit(0)
txt += [
       f'ring T={p},({",".join(map(str,unk))},w,w2),dp;',
       'ideal I=imap(R,I); I=simplify(I,2); "NEQ",size(I);',
       'ideal G=std(I); "DIM",dim(G);',
       'if(dim(G)==0){ "VDIM",vdim(G); }',
       'if(dim(G)>0){ int dd=dim(G); ideal L; int q; for(q=1;q<=dd;q++){ L=L,' + '+'.join(f'random(1,{p-1})*{u}' for u in unk) + '+random(1,%d); } I=I,L; G=std(I); "SDIM",dim(G); if(dim(G)==0){ "SVDIM",vdim(G); } }' % (p-1),
       'if(dim(G)==0){ LIB "primdec.lib"; list M=minAssGTZ(G); "NCOMP",size(M); option(redSB); int m; for(m=1;m<=size(M);m++){ ideal Pm=std(M[m]); "COMP",m,size(Pm); int t; for(t=1;t<=size(Pm);t++){ "GBELEM"; print(Pm[t]); "ENDELEM"; } "ENDCOMP"; kill Pm; kill t; } }',
       'quit;']
fn = os.path.join(ND, f'{args.mode}_{args.seed}_{p}{"_lam0" if args.lam0 else ""}.sing'); open(fn, 'w', newline='\n').write('\n'.join(txt)+'\n')
wp = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
out = subprocess.run(base.WSL+['bash', '-c', f'timeout {args.timeout} Singular -q < "{wp}" 2>&1'], capture_output=True, text=True).stdout
print(f'{desc}; p={p}, seed={args.seed}; unknowns={len(unk)}')
for key in ('DETTERMS', 'NMONO', 'NEQ', 'DIM', 'VDIM', 'SDIM', 'SVDIM', 'NCOMP'):
    m = re.search(key+r'\s+(-?\d+)', out); print(f'  {key} = {m.group(1) if m else "?"}')
if 'DIM' not in out: print('  (no DIM: timeout/error)'); print(out[-500:]); sys.exit(2)
dim = int(re.search(r'DIM\s+(-?\d+)', out).group(1))
if dim < 0: print('  => INCONSISTENT: no constant-Hessian completion mod p for this leading structure'); sys.exit(0)
# consistent: try to extract a rational point
syms = {str(u): u for u in unk}; wsym = sp.Symbol('w'); syms['w'] = wsym; w2sym = sp.Symbol('w2'); syms['w2'] = w2sym; allv = unk + [wsym, w2sym]
found = None
for comp in re.findall(r'COMP\s+\d+\s+\d+\s*\n(.*?)ENDCOMP', out, re.S):
    gens = [g.strip().replace('\n', '') for g in re.findall(r'GBELEM\s*\n(.*?)ENDELEM', comp, re.S)]
    polys = [sp.sympify(g.replace('^', '**'), locals=syms) for g in gens if g]
    if not polys or not all(sp.Poly(q, *allv).total_degree() == 1 for q in polys): continue
    M = [[int(sp.Poly(q, *allv).coeff_monomial(a)) % p for a in allv] + [int(sp.Poly(q, *allv).coeff_monomial(1)) % p] for q in polys]
    nv = len(allv); r = 0; piv = []
    for c in range(nv):
        kk = next((i for i in range(r, len(M)) if M[i][c]), None)
        if kk is None: continue
        M[r], M[kk] = M[kk], M[r]; inv = pow(M[r][c], p-2, p); M[r] = [(v*inv) % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]: fct = M[i][c]; M[i] = [(a-fct*b) % p for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    if r == nv: found = {allv[c]: (-M[i][nv]) % p for i, c in enumerate(piv)}; break
if not found: print('  => CONSISTENT but no rational point extracted (components not rational / not 0-dim)'); sys.exit(0)
vals = {u: found[u] for u in unk}
fsol = sp.expand(f.subs(vals))
def modp_poly(e):
    P = sp.Poly(sp.expand(e), *X); d = {}
    for m, c in P.terms():
        r = sp.Rational(c); cc = (int(r.p)*pow(int(r.q), p-2, p)) % p
        if cc: d[m] = cc
    return d
dh = modp_poly(hess(fsol).det(method='berkowitz'))
const = dh.get((0, 0, 0, 0), 0); nonconst = {m: c for m, c in dh.items() if any(m)}
print(f'  explicit f found; independent sympy check: det Hess f mod p == {const} (unit? {const != 0}); nonconstant terms: {len(nonconst)}')
# pivot cone: v with D_v^2 f constant
v = sp.symbols('v1:5'); D2 = sp.expand(sum(v[i]*v[j]*sp.diff(fsol, X[i], X[j]) for i in range(4) for j in range(4)))
Pd = sp.Poly(D2, *X); eqs = [c for m, c in Pd.terms() if any(m)]
ptxt = [f'ring R={p},(v1,v2,v3,v4),dp;', 'ideal I=' + ',\n'.join(s_modp_v for s_modp_v in [base.s_poly(sp.expand(sp.Poly(e, *v).as_expr())) for e in eqs]) + ';',
        'ideal G=std(I); "PDIM",dim(G); quit;']
fn2 = os.path.join(ND, f'{args.mode}_{args.seed}_{p}_pivot.sing'); open(fn2, 'w', newline='\n').write('\n'.join(ptxt)+'\n')
wp2 = '/mnt/c/' + os.path.abspath(fn2)[3:].replace('\\', '/')
out2 = subprocess.run(base.WSL+['bash', '-c', f'timeout 120 Singular -q < "{wp2}" 2>&1'], capture_output=True, text=True).stdout
m = re.search(r'PDIM\s+(-?\d+)', out2); pdim = int(m.group(1)) if m else None
print(f'  pivot cone {{v : D_v^2 f const}} has dimension {pdim} (<= 0 means PIVOT-FREE)')
json.dump({'mode': args.mode, 'seed': args.seed, 'p': p, 'f': str(fsol), 'const': const, 'pivot_dim': pdim},
          open(os.path.join(ND, f'{args.mode}_{args.seed}_{p}{"_lam0" if args.lam0 else ""}_solution.json'), 'w'))
print('  solution saved')
