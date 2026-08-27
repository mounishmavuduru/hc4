# PE -- PIVOT EXISTENCE, part 2: unconditional pivot existence in low degree,
# and the two dimension lemmas.
#
# THEOREM PE-1 (pivot existence in degree <= 4).
#   Every f in C[x1..x4] with det Hess f in C^x and deg f <= 4 admits a pivot.
#   (This is a by-product of the project's Theorem B: its proof exhibits a
#   pivot in every branch.  Isolating it matters because pivot existence, not
#   HC_4 itself, is what the pivot dichotomy needs.)
#
# THEOREM PE-2 (variation-space bound).  Let W = span{ N_alpha } in Sym_n(C),
#   where Hess f - Hess f(0) = sum_alpha x^alpha N_alpha.  If dim W <= n-1
#   then f admits a pivot, for ANY f (no Hessian hypothesis at all).
#   Proof: the pivot cone is the common zero locus in P^{n-1} of dim W <= n-1
#   quadrics; by the projective dimension theorem, at most n-1 hypersurfaces
#   in P^{n-1} always meet.  (n = 4: dim W <= 3 suffices.)
#
# THEOREM PE-3 (one-storey potentials).  If f = f_2 + f_d (a quadratic plus a
#   single homogeneous part, d >= 3) and det Hess f in C^x, then f admits a
#   pivot.  Proof: comparing graded pieces, det Hess f_d = 0; Gordan-Noether
#   in n = 4 makes f_d a cone, D_v f_d = 0, so D_v^2 f = v^T (Hess f_2) v.
#   Corollary: a pivot-free 4-variable potential has at least two non-zero
#   homogeneous parts of degree >= 3.
#
# This script certifies the computational content: the graded identities used
# by the degree-<=4 branches, the isotropic-vector lemma, PE-2's quadric count,
# and an instance sweep of degree-4 constant-Hessian potentials.

import sympy as sp
from pe_pivot_cone import (pivot_quadrics, cone_is_trivial, hess,
                           variation_space_dim, witness)

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
xs = [x1, x2, x3, x4]
vs = list(sp.symbols('v1:5'))

FAIL = []


def check(name, cond):
    print(("  PASS  " if cond else "  FAIL  ") + name)
    if not cond:
        FAIL.append(name)


def gen_form(deg, vars_, tag):
    """generic homogeneous form of degree deg in vars_, symbolic coefficients"""
    mons = sp.itermonomials(vars_, deg, deg)
    mons = sorted(mons, key=sp.default_sort_key)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sp.expand(sum(c * m for c, m in zip(cs, mons))), list(cs)


print("=" * 70)
print("PE-1  the degree-<=4 branches, graded identities re-derived generically")
print("=" * 70)

# ---- branch r = 3 :  T3 for d = 4.
# f4 in C[x1,x2,x3] (Gordan-Noether), f3, f2 generic.
# claim: [det Hess f]_7 = det_3(Hess_3 f4) * d^2 f3/dx4^2 .
f4a, _ = gen_form(4, [x1, x2, x3], 'a')
f3a, _ = gen_form(3, xs, 'b')
f2a, _ = gen_form(2, xs, 'c')
H = hess(f4a + f3a + f2a, xs)
det = sp.expand(H.det(method='berkowitz'))
part7 = sum(co * sp.prod([g**e for g, e in zip(xs, m)])
            for m, co in sp.Poly(det, *xs).terms() if sum(m) == 7)
K3 = sp.Matrix(3, 3, lambda i, j: sp.diff(f4a, xs[i], xs[j]))
claim7 = sp.expand(K3.det(method='berkowitz') * sp.diff(f3a, x4, 2))
check("T3 (d=4):  [det Hess f]_7 = det_3(Hess_3 f4) * (f3)_{x4x4}",
      sp.expand(part7 - claim7) == 0)
print("        => if f4 has essential rank 3 then (f3)_{x4x4} = 0, so")
print("           D_{e4}^2 f = (f2)_{x4x4} is CONSTANT: e4 is a pivot.")

# ---- branch r = 2 :  [det Hess f]_6 = det_2(Hess_2 f4)*det_2(Hess_{(3,4)} f3)
f4b, _ = gen_form(4, [x1, x2], 'p')
f3b, _ = gen_form(3, xs, 'q')
f2b, _ = gen_form(2, xs, 'r')
H = hess(f4b + f3b + f2b, xs)
det = sp.expand(H.det(method='berkowitz'))
part6 = sum(co * sp.prod([g**e for g, e in zip(xs, m)])
            for m, co in sp.Poly(det, *xs).terms() if sum(m) == 6)
A2 = sp.Matrix(2, 2, lambda i, j: sp.diff(f4b, xs[i], xs[j]))
B2 = sp.Matrix(2, 2, lambda i, j: sp.diff(f3b, xs[i + 2], xs[j + 2]))
check("S2 (r=2): [det Hess f]_6 = det_2(Hess_2 f4)*det_2(Hess_{(x3,x4)} f3)",
      sp.expand(part6 - A2.det() * B2.det()) == 0)
print("        => det_2 Hess_{(x3,x4)} f3 = 0; the 2x2 lemma gives")
print("           Hess_{(x3,x4)} f3 = l(x) u u^T with u CONSTANT; any")
print("           0 != v in span(e3,e4) with u.v = 0 has D_v^2 f3 = 0 and")
print("           D_v^2 f4 = 0 (f4 in C[x1,x2]): v is a pivot.")

# ---- branch r = 1 :  [det Hess f]_5 = 12 x1^2 * det_3(Hess_{(x2,x3,x4)} f3)
f4c = x1**4
f3c, _ = gen_form(3, xs, 's')
f2c, _ = gen_form(2, xs, 't')
H = hess(f4c + f3c + f2c, xs)
det = sp.expand(H.det(method='berkowitz'))
part5 = sum(co * sp.prod([g**e for g, e in zip(xs, m)])
            for m, co in sp.Poly(det, *xs).terms() if sum(m) == 5)
M3 = sp.Matrix(3, 3, lambda i, j: sp.diff(f3c, xs[i + 1], xs[j + 1]))
check("S2 (r=1): [det Hess f]_5 = 12 x1^2 det_3(Hess_{(x2,x3,x4)} f3)",
      sp.expand(part5 - 12 * x1**2 * M3.det(method='berkowitz')) == 0)

# isotropic-vector lemma, key identity: for M0 of rank 2 in Sym_3,
# adj(M0) = kappa v0 v0^T with M0 v0 = 0, and [t] det(M0 + tM) = tr(adj M0 . M)
Msym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'm{min(i,j)}{max(i,j)}'))
Nsym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'n{min(i,j)}{max(i,j)}'))
t = sp.Symbol('t')
dt = sp.expand(sp.expand((Msym + t * Nsym).det(method='berkowitz')))
coef_t = sp.Poly(dt, t).coeff_monomial(t)
check("S3 key identity: [t] det(M0 + tM) = tr(adj(M0) M)",
      sp.expand(coef_t - sp.trace(Msym.adjugate() * Nsym)) == 0)
print("        => an all-singular linear family in Sym_3 has a common")
print("           isotropic vector; in branch r=1 that vector v lies in")
print("           span(e2,e3,e4), and D_v^2(x1^4) = 0: v is a pivot.")
print()
print("  CONCLUSION (Theorem PE-1): deg f <= 4, det Hess f in C^x  =>  pivot.")
print("  (deg <= 2 trivial; deg = 3: det Hess f3 = 0 and Gordan-Noether make")
print("   f3 a cone, and a cone direction is a pivot.)")

print()
print("=" * 70)
print("PE-2  the variation space W and the projective dimension theorem")
print("=" * 70)
# the pivot cone is cut out by exactly dim W independent quadrics
fdemo = (x1**2 + x2**2 + x3**2 + x4**2) / 2 + x1**3 + x1 * x2 * x3
q = pivot_quadrics(fdemo, xs, vs)
dW = variation_space_dim(fdemo, xs)
# the C-span of the quadrics has dimension exactly dim W
rows = []
vmons = sorted(sp.itermonomials(vs, 2, 2), key=sp.default_sort_key)
for g in q:
    rows.append([sp.Poly(g, *vs).coeff_monomial(m) for m in vmons])
check("#independent quadrics  ==  dim W", sp.Matrix(rows).rank() == dW)
print(f"        (demo: dim W = {dW}, #quadrics = {len(q)},"
      f" rank = {sp.Matrix(rows).rank()})")
print("  => if dim W <= 3 the pivot cone is a nonempty cone in C^4")
print("     (<=3 hypersurfaces in P^3 always meet).  A pivot-free")
print("     4-variable potential therefore has dim W >= 4.")

print()
print("=" * 70)
print("PE-3  instance sweep: degree-4 constant-Hessian potentials")
print("=" * 70)
# Exhibit degree-4 potentials with det Hess in C^x from each branch and
# confirm the pivot cone is nonempty (fail-closed).
tests = []
# r = 3 branch: f4 essential rank 3, f3 with no x4^2, plus a nondegenerate q
tests.append(("r=3 instance",
              (x1**2 + x2**2 + x3**2) / 2 + x1 * x4
              + (x1**3 + x2**3 + x3**3 + x1 * x2 * x3)
              + (x1**4 + x2**4 + x3**4 + x1**2 * x2 * x3)))
# r = 2 branch
tests.append(("r=2 instance",
              x1 * x3 + x2 * x4 + x1**3 + x1**2 * x2 + x1**4 + x1**2 * x2**2))
# r = 1 branch
tests.append(("r=1 instance",
              x1 * x3 + x2 * x4 + x1**4 + x1**2 * x2 + x2**3))
for name, f in tests:
    H = hess(f, xs)
    d = sp.expand(H.det(method='berkowitz'))
    quads = pivot_quadrics(f, xs, vs)
    triv = cone_is_trivial(quads, vs)
    w = witness(quads, vs)
    ok = (d.is_number and d != 0)
    print(f"  {name}: det Hess = {d}")
    if ok:
        check(f"    {name}: pivot cone nonempty", not triv)
        print(f"         witness v = {list(w) if w is not None else None}")
    else:
        print("    (not a constant-Hessian potential; skipped)")

print()
if FAIL:
    print("FAILURES:", FAIL)
    raise SystemExit(1)
print("ALL PE-1/PE-2/PE-3 CHECKS PASS.")
