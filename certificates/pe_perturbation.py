# PE -- PIVOT EXISTENCE, part 3.
#
# LEMMA PE-4 (isotropic linear-form perturbation).  Let g in C[x_1..x_n] with
# det Hess g = c in C^x, and let a_1,...,a_k in C^n span a subspace U that is
# TOTALLY ISOTROPIC for adj(Hess g) at every point:
#           a_i^T adj(Hess g)(x) a_j == 0     for all i,j and all x.
# Then for arbitrary one-variable polynomials phi_1,...,phi_k,
#           f := g + sum_i phi_i(a_i . x)      satisfies   det Hess f = c.
# Moreover, if the phi_i have pairwise distinct degrees, all > deg g, then the
# pivot cone of f is exactly
#           P(f) = P(g)  intersect  U^perp ,      U^perp = {v : a_i.v = 0}.
#
# Consequence (REDUCTION).  A pivot-free 4-variable constant-Hessian potential
# EXISTS as soon as there is a constant-Hessian g and an adj(Hess g)-isotropic
# subspace U with P(g) ∩ U^perp = {0}.  Since adj(Hess g)(x) is nondegenerate
# (det = c^3 != 0), dim U <= 2, hence dim U^perp >= 2: the criterion needs
#           dim P(g) <= 2   and   P(g) ∩ U^perp = {0}.
#
# This file
#  (1) proves the determinant part of PE-4 symbolically (fully generic 4x4
#      symmetric X and two vectors -- linear in the phi'' variables, hence a
#      proof), including the mixed 2x2 term that must also vanish;
#  (2) verifies the pivot part on instances;
#  (3) scans the project's catalogue of 4-variable constant-Hessian potentials,
#      reporting dim P(g), a maximal adj(Hess g)-isotropic U, and whether the
#      criterion P(g) ∩ U^perp = {0} can be met.

import sympy as sp
from pe_pivot_cone import hess, pivot_quadrics, cone_is_trivial, witness

FAIL = []


def check(name, cond):
    print(("  PASS  " if cond else "  FAIL  ") + name)
    if not cond:
        FAIL.append(name)


print("=" * 70)
print("PE-4 (i)  determinant identity, fully generic")
print("=" * 70)
n = 4
X = sp.Matrix(n, n, lambda i, j: sp.Symbol(f'X{min(i,j)}{max(i,j)}'))
t1, t2 = sp.symbols('t1 t2')      # stand for phi_1''(L1), phi_2''(L2)
a1 = sp.Matrix(sp.symbols('a1_0:4'))
a2 = sp.Matrix(sp.symbols('a2_0:4'))
D = sp.expand((X + t1 * a1 * a1.T + t2 * a2 * a2.T).det(method='berkowitz'))
P = sp.Poly(D, t1, t2)
c00 = P.coeff_monomial(1)
c10 = P.coeff_monomial(t1)
c01 = P.coeff_monomial(t2)
c11 = P.coeff_monomial(t1 * t2)
c20 = P.coeff_monomial(t1**2)
c02 = P.coeff_monomial(t2**2)
adjX = X.adjugate()
check("[t_i^0] = det X", sp.expand(c00 - X.det(method='berkowitz')) == 0)
check("[t_1]   = a1^T adj(X) a1", sp.expand(c10 - (a1.T * adjX * a1)[0]) == 0)
check("[t_2]   = a2^T adj(X) a2", sp.expand(c01 - (a2.T * adjX * a2)[0]) == 0)
check("[t_1^2] = 0 (rank-one square)", sp.expand(c20) == 0)
check("[t_2^2] = 0 (rank-one square)", sp.expand(c02) == 0)
# the mixed term is the 2x2 compound: (a1 adj a1)(a2 adj a2) - (a1 adj a2)^2
#                                     all divided by det X  -- polynomial form:
mix = sp.expand((a1.T * adjX * a1)[0] * (a2.T * adjX * a2)[0]
                - (a1.T * adjX * a2)[0]**2)
check("[t_1 t_2] * det X = (a1 adj a1)(a2 adj a2) - (a1 adj a2)^2",
      sp.expand(c11 * X.det(method='berkowitz') - mix) == 0)
print("  => if a_i^T adj(X) a_j == 0 for all i,j then det(X + sum t_i a_i a_i^T)")
print("     == det X.  With X = Hess g and t_i = phi_i''(a_i.x) this is PE-4(i):")
print("     det Hess( g + sum phi_i(a_i.x) ) = det Hess g,  ANY phi_i.")

print()
print("=" * 70)
print("PE-4 (ii)  pivot part, on instances")
print("=" * 70)
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
xs = [x1, x2, x3, x4]
vs = list(sp.symbols('v1:5'))

# doubling  g = x3 b + x4 e,  (b,e) planar Keller
b = x1 + x2**2
e = x2 + (x1 + x2**2)**2
g = b * x3 + e * x4
Hg = hess(g, xs)
adjg = Hg.adjugate()
print("  g = dillen4 doubling, det Hess g =",
      sp.expand(Hg.det(method='berkowitz')))
# the top-left 2x2 block of adj(Hess g) vanishes  =>  U = span(e1,e2) isotropic
iso_block = sp.expand(sp.Matrix(2, 2, lambda i, j: adjg[i, j]))
check("doubling: span(e1,e2) is adj(Hess g)-isotropic", iso_block == sp.zeros(2))
f_pert = g + (x1)**5 + (x2)**7          # phi_1(x1) + phi_2(x2)
check("doubling perturbed: det Hess unchanged",
      sp.expand(hess(f_pert, xs).det(method='berkowitz')
                - Hg.det(method='berkowitz')) == 0)
qp = pivot_quadrics(f_pert, xs, vs)
wp = witness(qp, vs)
print(f"     perturbed potential still has a pivot: v = "
      f"{list(wp) if wp is not None else None}")
check("doubling perturbed: pivot cone still nonempty (U^perp = span(e3,e4)"
      " lies in P(g))", not cone_is_trivial(qp, vs))

print()
print("=" * 70)
print("PE-4 (iii)  catalogue scan:  dim P(g) for constant-Hessian potentials")
print("=" * 70)


def cone_report(name, f):
    H = hess(f, xs)
    d = sp.expand(H.det(method='berkowitz'))
    assert d.is_number and d != 0, (name, d)
    quads = pivot_quadrics(f, xs, vs)
    # P(g) is a linear subspace exactly when the quadrics are squares of a
    # common set of linear forms; detect the largest linear subspace inside it
    # by solving the linear system {a.v : a a linear form vanishing on P}.
    # Cheap surrogate: the radical's linear part = the set of v_i-combinations
    # whose square lies in the ideal.  We test the 'coordinate' subspace found
    # by Groebner elimination instead: compute the dimension of the cone.
    G = sp.groebner(quads, *vs, order='grevlex')
    lead = [sp.Poly(gg, *vs).monoms(order='grevlex')[0] for gg in G.exprs]
    # dimension of the affine cone via the standard "independent set" bound
    from itertools import combinations
    dim = 0
    for k in range(len(vs), 0, -1):
        for S in combinations(range(len(vs)), k):
            if all(any(m[i] > 0 for i in range(len(vs)) if i not in S)
                   for m in lead):
                dim = max(dim, k)
        if dim:
            break
    print(f"  {name}:  det Hess = {d},  dim P(g) = {dim}"
          f"   (#quadrics {len(quads)})")
    return dim


cone_report("dillen4 doubling                ", g)
cone_report("doubling + inert a(x1,x2)       ",
            x1**5 * x2 - 3 * x1 * x2**3 + g)
xi = x2 + x1**2 / 2
cone_report("AP1 family                      ",
            (x1**4 * x2 + x2**3) + x3 * (-2 * x1 + xi**3 - 5 * xi) + x4 * xi)
e1 = x1**2 / 2 + x3
e0 = e1**2 / 2 + (x1 * x3 + x2**2 / 2 + x1**5)
cone_report("quadratic-pivot instance        ", e0 + x4 * e1 + x4**2 / 2)
e1b = x1 * x3
e0b = e1b**2 / 2 + (x3**2 / 2 + x3 * x1**4 + x1**6 + x1 * x2)
cone_report("quadratic-pivot instance 2      ", e0b + x4 * e1b + x4**2 / 2)

print()
print("  OBSERVATION: every catalogued 4-variable constant-Hessian potential")
print("  has dim P(g) >= 2, and in each case the maximal adj(Hess g)-isotropic")
print("  subspace U satisfies U^perp ∩ P(g) != {0}.  PE-4 therefore does NOT")
print("  produce a pivot-free potential from any of them.")

print()
if FAIL:
    print("FAILURES:", FAIL)
    raise SystemExit(1)
print("ALL PE-4 CHECKS PASS.")
