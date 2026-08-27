# LEMMA D1 (new result, this project): pivot obstruction for the Meng-Yang family.
#
# Claim: for every lambda != 0 and every mu, the 5-variable potential
#     psi_{lam,mu} = B + (lam/2) A^2 + mu A     (Meng-Yang, arXiv:2607.22198)
# admits NO direction v in C^5 \ {0} with  D_v^2 psi_{lam,mu} == 0  identically.
#
# Consequently no AFFINE change of coordinates on C^5 makes psi_{lam,mu}
# affine-linear in any coordinate, so the Schur-descent pivot required to
# descend from 5 to 4 variables does not exist for any member of the known
# HC_5 counterexample family.  (This makes Remark 4.2 of arXiv:2607.22198
# rigorous for the family produced there.)
#
# Proof strategy (exact, machine-checked below):
#   D_v^2 psi = D_v^2 B + mu D_v^2 A + lam (A D_v^2 A + (D_v A)^2).
#   deg_w D_v^2 B <= 5 and deg_w D_v^2 A <= 5, while A D_v^2 A + (D_v A)^2
#   contributes up to degree 12.  Hence for every w-monomial of total degree
#   >= 6 the coefficient of D_v^2 psi equals lam * (that coefficient in
#   C(v) := A D_v^2 A + (D_v A)^2), a mu-free quadratic form in v.
#   [Checked programmatically below.]
#   It therefore suffices to show that the quadratic forms
#       { [w^alpha] C(v) : |alpha| >= 6 }
#   have no common zero v != 0 over C.  Since all generators are homogeneous
#   (quadratic) in v, their zero set is a cone; it equals {0} iff the leading
#   term ideal of a Groebner basis contains a pure power of each variable
#   (standard finiteness criterion; a finite cone is the origin).
#
# Output: PASS/FAIL for each step.

import sympy as sp

x1, x2, y1, y2, y3 = sp.symbols('x1 x2 y1 y2 y3')
w = [x1, x2, y1, y2, y3]
u = 1 + x1*x2
A = y1*u**3 + 3*x1*y2*u**2 - x1**3*y3
B = (y1*x2**2*u*(4 + 3*x1*x2)
     + y2*(x2 + 3*x1*x2**2*(4 + 3*x1*x2))
     + y3*(2*x1 - 3*x1**2*x2))

v = list(sp.symbols('v1:6'))

def Dv(h):
    return sp.expand(sum(v[i]*sp.diff(h, w[i]) for i in range(5)))

DA = Dv(A)          # D_v A
D2A = Dv(DA)        # D_v^2 A
D2B = Dv(Dv(B))     # D_v^2 B

# --- Step 0: degree bookkeeping used in the proof ---
assert sp.Poly(D2B, *w).total_degree() <= 5
assert sp.Poly(D2A, *w).total_degree() <= 5
print("Step 0 PASS: deg D_v^2 B <= 5, deg D_v^2 A <= 5")

# --- Step 1: the mu-free, lam-scaled high-degree subsystem ---
C = sp.expand(A*D2A + DA**2)
PC = sp.Poly(C, *w)
gens_high = []
for monom, coeff in PC.terms():
    if sum(monom) >= 6:
        gens_high.append(sp.expand(coeff))
print(f"Step 1: {len(gens_high)} coefficient equations of w-degree >= 6")

# each generator must be a homogeneous quadratic in v
for q in gens_high:
    pq = sp.Poly(q, *v)
    assert pq.is_homogeneous and pq.total_degree() == 2, q
print("Step 1 PASS: all generators are homogeneous quadratics in v")

# --- Step 2: Groebner basis; variety = {0} criterion ---
G = sp.groebner(gens_high, *v, order='grevlex')
lts = [sp.Poly(g, *v).LM(order='grevlex') for g in G.exprs]
pure = {i: False for i in range(5)}
for lm in lts:
    exps = lm.exponents if hasattr(lm, 'exponents') else None
for g in G.exprs:
    p = sp.Poly(g, *v)
    m = p.monoms(order='grevlex')[0]
    nz = [i for i, e in enumerate(m) if e > 0]
    if len(nz) == 1:
        pure[nz[0]] = True
print("pure-power leading terms found for variables:",
      [str(v[i]) for i in range(5) if pure[i]])
assert all(pure.values()), "criterion not satisfied by degree>=6 subsystem alone"
print("Step 2 PASS: leading-term ideal contains a pure power of each v_i")
print()
print("LEMMA D1 VERIFIED: for all lam != 0, mu, the only v with")
print("D_v^2 psi_{lam,mu} == 0 is v = 0. No Schur pivot exists on the")
print("Meng-Yang HC_5 counterexample family, for any affine coordinates.")
