# Certificate check: Meng-Yang five-variable counterexample to HC_5 (arXiv:2607.22198 v2)
# and their six-variable warm-up counterexample to HC_6.
# All computations exact over Q.

import sympy as sp

x1, x2, x3, y1, y2, y3 = sp.symbols('x1 x2 x3 y1 y2 y3')
u = 1 + x1*x2

# Alpoge's map F (verified separately in verify_alpoge_3d.py)
F1 = u**3*x3 + x2**2*u*(4 + 3*x1*x2)
F2 = x2 + 3*x1*u**2*x3 + 3*x1*x2**2*(4 + 3*x1*x2)
F3 = 2*x1 - 3*x1**2*x2 - x1**3*x3

# ---- HC_6 counterexample: phi = <y, F(x)>, claimed det Hess = -4 ----
phi = y1*F1 + y2*F2 + y3*F3
vars6 = [x1, x2, x3, y1, y2, y3]
H6 = sp.hessian(phi, vars6)
det6 = sp.expand(H6.det(method='berkowitz'))
print("det Hess phi =", det6)
assert det6 == -4

grad6 = [sp.diff(phi, v) for v in vars6]
pts3 = [(0, 0, sp.Rational(-1, 4)), (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
        (-1, sp.Rational(3, 2), sp.Rational(13, 2))]
g_at = []
for p in pts3:
    sub = {x1: p[0], x2: p[1], x3: p[2], y1: 0, y2: 0, y3: 0}
    g_at.append(tuple(sp.simplify(g.subs(sub)) for g in grad6))
print("grad phi at (P,0):", g_at)
assert g_at[0] == g_at[1] == g_at[2] == (0, 0, 0, sp.Rational(-1, 4), 0, 0)

# ---- HC_5 counterexample: Psi = A^2 + 13A + 2B in (x1,x2,y1,y2,y3) ----
# F = F0 + x3*F1vec  (F affine-linear in x3); A = <y,F1vec>, B = <y,F0>
A = y1*u**3 + 3*x1*y2*u**2 - x1**3*y3
B = (y1*x2**2*u*(4 + 3*x1*x2)
     + y2*(x2 + 3*x1*x2**2*(4 + 3*x1*x2))
     + y3*(2*x1 - 3*x1**2*x2))
# consistency: phi == x3*A + B
assert sp.expand(phi - (x3*A + B)) == 0

Psi = sp.expand(A**2 + 13*A + 2*B)
vars5 = [x1, x2, y1, y2, y3]
assert sp.Poly(Psi, *vars5).total_degree() == 14
print("degree 14 OK; monomial count:", len(sp.Poly(Psi, *vars5).terms()))

H5 = sp.hessian(Psi, vars5)
det5 = sp.expand(H5.det(method='berkowitz'))
print("det Hess Psi =", det5)
assert det5 == 128

grad5 = [sp.diff(Psi, v) for v in vars5]
Pp = {x1: 1, x2: sp.Rational(-3, 2), y1: 0, y2: 0, y3: 0}
Pm = {x1: -1, x2: sp.Rational(3, 2), y1: 0, y2: 0, y3: 0}
gp = tuple(sp.simplify(g.subs(Pp)) for g in grad5)
gm = tuple(sp.simplify(g.subs(Pm)) for g in grad5)
print("grad Psi(P+) =", gp)
print("grad Psi(P-) =", gm)
assert gp == gm == (0, 0, sp.Rational(-1, 2), 0, 0)

# ---- the whole two-parameter family: det Hess (B + (lam/2)A^2 + mu*A) == -lam*c, c=-4 ----
lam, mu = sp.symbols('lam mu')
psi_fam = B + lam/2*A**2 + mu*A
H5f = sp.hessian(psi_fam, vars5)
det5f = sp.expand(H5f.det(method='berkowitz'))
print("det Hess psi_{lam,mu} =", det5f)
assert det5f == 4*lam    # -lam*c with c = det Hess phi ... c here is -4, so -lam*(-4)=4*lam

print("ALL VERIFIED: HC_6 counterexample (det=-4, 3-point gradient collision),")
print("HC_5 counterexample Psi (deg 14, det Hess=128, 2-point gradient collision),")
print("and the full Schur-descent family det Hess psi_{lam,mu} = 4*lam.")
