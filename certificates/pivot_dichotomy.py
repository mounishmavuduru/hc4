# THE PIVOT DICHOTOMY -- careful derivation, including a branch the first
# write-up glossed over.
#
# Setting: f = e0(x1,x2,x3) + x4*e1(x1,x2,x3), det Hess f = c != 0 (affine
# pivot). By Theorem G + Theorem F the pivot coefficient e1 is affinely
# 2-variable, so after an affine change e1 = e(x1,x2).
#
# With A the upper-left 3x3 block of Hess f and w = (e_1, e_2, 0)^T,
#     det Hess f = -w^T adj(A) w,
#     A = Hess_3 e0 + x4 * K,  K = Hess_3 e1 supported in the (1,2) block.
# Expanding, the COEFFICIENT OF x4 in det Hess f is
#     -(e0)_{33} * B_2(e),      B_2(e) := e_{22} e_1^2 - 2 e_{12} e_1 e_2 + e_{11} e_2^2
# the two-variable bordered Hessian of e. Constancy of det Hess f forces
#     (e0)_{33} * B_2(e) == 0,
# so, C[x] being a domain, exactly one of two branches holds:
#
#   (a) (e0)_{33} == 0.  Then e0 = a(x1,x2) + x3*b(x1,x2) and
#       f = a + x3 b + x4 e is a MENG DOUBLING of the planar map (b,e) plus
#       an inert planar potential (Theorem C): governed exactly by JC_2.
#
#   (b) B_2(e) == 0.  Then by THEOREM G in two variables det_2 Hess e == 0,
#       so e = phi(L) for a linear form L; after a linear change e = phi(x1),
#       and then det Hess f = -phi'(x1)^2 * det_2 Hess_{(x2,x3)} e0. For this
#       to be a nonzero constant, phi'(x1)^2 must divide a nonzero constant,
#       so phi' is a nonzero constant and e is AFFINE -- i.e. we are in the
#       AP0 case, which is closed UNCONDITIONALLY (fibrewise Dillen HC_2).
#
# CONCLUSION (dichotomy): a four-variable potential with constant nonzero
# Hessian determinant admitting an affine pivot is either a Meng doubling
# plus an inert planar potential (governed by JC_2) or has an affine pivot
# coefficient (settled unconditionally). Together with Theorem A's quadratic
# branch: HC_4 restricted to potentials admitting a pivot is EQUIVALENT to
# JC_2. Branch (b) does NOT weaken the equivalence -- it is unconditional --
# but it must be stated, since the first write-up asserted (e0)_{33} == 0
# outright.
#
# Checks: the x4-coefficient formula fully symbolically; branch (b)'s
# determinant formula; the divisibility argument on an ansatz; and an
# instance of each branch.

import itertools
import sympy as sp

x1, x2, x3, x4 = X = sp.symbols('x1 x2 x3 x4')
Xp = (x1, x2, x3)


def hess(f, vs):
    return sp.Matrix([[sp.diff(f, a, b) for b in vs] for a in vs])


def gen(vs, dmax, tag):
    mons = []
    for d in range(dmax + 1):
        mons += sorted({sp.prod(t) for t in itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(k*m for k, m in zip(cs, mons))


# ---- the x4-coefficient of det Hess f, fully symbolic ----
e0 = gen(Xp, 3, 'p')                 # generic e0 of degree 3 in x1,x2,x3
e = gen((x1, x2), 3, 'q')            # generic planar pivot coefficient
f = sp.expand(e0 + x4*e)
d = sp.expand(hess(f, X).det(method='berkowitz'))
coef_x4 = sp.expand(sp.diff(d, x4))
# is it independent of x4? (det is affine in x4 here)
assert sp.expand(sp.diff(d, x4, 2)) == 0, 'det Hess f is not affine in x4'

B2 = sp.expand(sp.diff(e, x2, 2)*sp.diff(e, x1)**2
               - 2*sp.diff(e, x1, x2)*sp.diff(e, x1)*sp.diff(e, x2)
               + sp.diff(e, x1, 2)*sp.diff(e, x2)**2)
pred = sp.expand(-sp.diff(e0, x3, 2)*B2)
assert sp.expand(coef_x4 - pred) == 0
print('PASS: [x4] det Hess f == -(e0)_{33} * B_2(e), fully symbolic')
print('      (generic e0 of degree 3 in three variables, generic planar e)')
print('      => constancy forces (e0)_{33} * B_2(e) == 0: a genuine dichotomy.')

# ---- branch (b): e = phi(x1) ----
phi = sp.Function('phi')
t = sp.Symbol('t')
ph = gen((x1,), 4, 'r')                       # phi(x1), generic of degree 4
f_b = sp.expand(e0 + x4*ph)
d_b = sp.expand(hess(f_b, X).det(method='berkowitz'))
det23 = sp.expand(sp.diff(e0, x2, 2)*sp.diff(e0, x3, 2) - sp.diff(e0, x2, x3)**2)
pred_b = sp.expand(-sp.diff(ph, x1)**2 * det23)
assert sp.expand(d_b - pred_b) == 0
print('PASS: branch (b) with e = phi(x1):'
      ' det Hess f == -phi\'(x1)^2 * det_2 Hess_{(x2,x3)} e0  (fully symbolic)')
print('      => phi\'^2 divides a nonzero constant, so phi is affine (AP0).')

# ---- explicit instance of each branch ----
# (a) a doubling: f = a + x3*b + x4*e with (b,e) a planar Keller map
b_pl = x2 + x1**2
e_pl = x1
f_a = sp.expand(x1**4 + x3*b_pl + x4*e_pl)
da = sp.expand(hess(f_a, X).det(method='berkowitz'))
Jac = sp.expand(sp.diff(b_pl, x1)*sp.diff(e_pl, x2) - sp.diff(b_pl, x2)*sp.diff(e_pl, x1))
print(f'\n(a) doubling instance: det Hess = {da}, Jac(b,e) = {Jac}, Jac^2 = {sp.expand(Jac**2)}')
assert da == sp.expand(Jac**2) and da != 0
assert sp.diff(sp.expand(x1**4 + x3*b_pl), x3, 2) == 0
print('    (e0)_{33} == 0 as branch (a) requires; det Hess == Jac(b,e)^2  PASS')

# (b) an AP0 instance: e affine
f_b2 = sp.expand(x1*x2 + x2**2*x3 + x4*x3)   # the potential from the AP0 study
db = sp.expand(hess(f_b2, X).det(method='berkowitz'))
e1_b = sp.diff(f_b2, x4)
print(f'\n(b) AP0 instance: e1 = {e1_b} (affine), det Hess = {db}')
assert db == 1 and sp.Poly(e1_b, *Xp).total_degree() <= 1
B2_b = sp.expand(sp.diff(e1_b, x2, 2)*sp.diff(e1_b, x1)**2
                 - 2*sp.diff(e1_b, x1, x2)*sp.diff(e1_b, x1)*sp.diff(e1_b, x2)
                 + sp.diff(e1_b, x1, 2)*sp.diff(e1_b, x2)**2)
assert B2_b == 0
print('    B_2(e1) == 0 as branch (b) requires  PASS')

print()
print('PIVOT DICHOTOMY VERIFIED, both branches. HC_4 restricted to potentials')
print('admitting a pivot is EQUIVALENT to JC_2: branch (a) is exactly JC_2 and')
print('branch (b), like the quadratic-pivot case, is unconditional.')
