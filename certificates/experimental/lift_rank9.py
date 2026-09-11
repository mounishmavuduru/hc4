# P5.A1 step 1: weight-9 lifting piece at an explicit survivor point (mod p).
# Weights wt(y)=1, wt(x4)=2. F = a5 + x4 b3 + (1/2) x4^2 y1 (weight 5, the
# survivor leading form). F' = a4(y) + x4 b2(y) + (1/2) x4^2 tau0 (weight 4).
# Each Hess(F') entry is exactly one weight below its Hess(F) counterpart, so
#   [det Hess(F+F'+lower)]_{wt 9} = tr( adj(Hess F) . Hess F' )   (Jacobi),
# which is LINEAR in the 22 weight-4 unknowns. Its coefficient matrix' rank over
# F_p (at the point) gives the dimension of admissible weight-4 perturbations.
import os, re, sys, ast, itertools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sympy as sp
from _d5_close import Y, a5, as_, form
y1,y2,y3 = Y
x4 = sp.Symbol('x4')
V4 = (y1,y2,y3,x4)
HERE = os.path.dirname(os.path.abspath(__file__))

# ---- the survivor point (rational mod p) from survpt/survivor_point.json ----
import json
jf = os.path.join(HERE,'survpt','survivor_point.json')
if not os.path.exists(jf):
    print('no survivor point json (run survivor_point_ext.py first)'); sys.exit(2)
J = json.load(open(jf))
p = J['p']; C0,C1,C2 = J['c']
vals = {sp.Symbol(k): int(v) % p for k,v in J['point'].items()}
print(f'point: c=({C0},{C1},{C2}), {len(vals)} a-values mod {p}')

# ---- F at the point, F' with symbolic weight-4 unknowns ----
b3 = C2*y1**3 + C0*y1**2*y2 + C1*y1**2*y3
a5p = sp.expand(a5.subs(vals))
F = a5p + x4*b3 + sp.Rational(1,2)*x4**2*y1
a4, A4 = form(Y,4,'p')          # 15 unknowns
b2, B2 = form(Y,2,'q')          # 6 unknowns
tau0 = sp.Symbol('t0')          # 1 unknown
U = list(A4)+list(B2)+[tau0]
Fp = a4 + x4*b2 + sp.Rational(1,2)*x4**2*tau0

def hess(g): return sp.Matrix(4,4,lambda i,j: sp.diff(g,V4[i],V4[j]))
HF = hess(F); HFp = hess(Fp)
T = sp.expand(sp.trace(HF.adjugate()*HFp))      # weight-9 piece, linear in U
# collect (y,x4)-monomial coefficients -> rows of the linear system in U
P = sp.Poly(T, *V4)
rows = []
for mon, coeff in P.terms():
    cp = sp.Poly(sp.expand(coeff), *U)
    row = [0]*len(U)
    for m2, c2 in cp.terms():
        idx = [i for i,e in enumerate(m2) if e]
        assert len(idx)==1 and m2[idx[0]]==1, 'not linear in unknowns'
        row[idx[0]] = int(c2) % p
    if any(row): rows.append(row)
print(f'weight-9 system: {len(rows)} equations in {len(U)} unknowns (mod {p})')

# ---- rank mod p (Gaussian elimination) ----
def rank_modp(M, p):
    M=[r[:] for r in M]; r=0; n=len(M[0]) if M else 0
    for c in range(n):
        piv=next((i for i in range(r,len(M)) if M[i][c]%p),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        inv=pow(M[r][c],p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
        r+=1
        if r==len(M): break
    return r
rk = rank_modp(rows, p)
print(f'rank = {rk};  null-space dim = {len(U)-rk}  (admissible weight-4 perturbations F\')')
# gauge sanity: F' = d_yi F (y-translations) and y_j * d_x4 F (shears x4 -> x4 + y_j)
# must lie in the null space (det Hess is invariant under these unimodular changes).
def as_vec(expr):
    Pp = sp.Poly(sp.expand(expr), *V4)
    # match against the basis monomials of a4, b2, x4^2/2
    target = sp.Poly(sp.expand(Fp), *V4)
    vec = [0]*len(U)
    for mon, cval in Pp.terms():
        cu = target.coeff_monomial(sp.prod(v**e for v,e in zip(V4,mon)))
        cp = sp.Poly(cu, *U)
        idx = [i for i,e in enumerate(cp.monoms()[0]) if e]
        assert len(idx)==1, f'monomial {mon} not in the weight-4 ansatz'
        rat = lambda q: (int(sp.Rational(q).p) * pow(int(sp.Rational(q).q), p-2, p)) % p   # rationals mod p (1/2 !)
        vec[idx[0]] = (rat(cval) * pow(rat(cp.coeffs()[0]), p-2, p)) % p
    return vec
gauge = [sp.diff(F, v) for v in (y1,y2,y3)] + [v*sp.diff(F, x4) for v in (y1,y2,y3)]
ok = True
for gexpr in gauge:
    vec = as_vec(gexpr)
    ok &= all(sum(r[i]*vec[i] for i in range(len(U))) % p == 0 for r in rows)
print(f'gauge directions (3 translations + 3 shears) in the null space? {ok}  => non-gauge null dim = {len(U)-rk-6}')
print('NOTE: weight-9 is homogeneous-linear, so F\'=0 always works; the obstruction (if any) '
      'sits in the lower weights where the nonzero constant must appear. This rank is the '
      'first structural datum of the lifting tower at a genuine survivor point.')
