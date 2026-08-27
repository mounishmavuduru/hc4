# litG_consequences.py  --  key: litG  (Theorem G: consequences, chain audit,
#                                       literature placement, attribution)
#
# Independent re-verification of the four parts of the "consequences of
# Theorem G beyond HC_4" task.  Exact arithmetic only; every check is a
# fail-closed assert.  Jet-level identities (identities in free symbols for
# the 1- and 2-jet entries, on which B and det Hess depend polynomially) are
# COMPLETE PROOFS for all polynomials of the stated shape; random-instance
# checks are labelled as consistency checks.
#
# NOTATION.  B(e) := grad(e)^T adj(Hess e) grad(e)  (bordered Hessian up to
# sign = Reilly's / Fox's U(e)).  E := degree-d homogenization of e in x0.
#
# ======================================================================
# PART 1.  THE LOW-DIMENSIONAL CHAIN, AND THE dBvdE FAMILY
# ======================================================================
#
#   B(e) == 0  =(Thm G)=>  det Hess_3 e == 0
#              =(de Bondt-van den Essen, Singular Hessians, J. Algebra 282
#                (2004) 195-204, Thm 3.3 = de Bondt's thesis Prop. 5.3.1,
#                read here from lit/debondt-thesis/homokema.txt lines
#                6521-6531)=>  after a LINEAR change
#                   (A)  e = a(x1,x2) + lam*x3,  lam in {0,1},   or
#                   (B)  e = a1(x1) + a2(x1)x2 + a3(x1)x3
#              =(S1.3 below)=>  e is affinely 2-variable.
#
# The task asked to verify that the family (B) "always has B != 0".  As
# tgc_debondt_chain.py and gcv_family_wronskian.py already found, that is
# FALSE as stated: the correct, and needed, statement is
#
#   (*)  B( sum_i g_i(l) x_i ) = -( det[ a | g(l) | g'(l) ] )^2,   l = a.x,
#
# so B == 0 exactly on the members that are affinely 2-variable.  Here (*)
# is proved for ALL degrees and ALL a at once (S1.2): the Hessian of the
# family is  a v^T + v a^T  and the matrix identity
#        w^T adj(a v^T + v a^T) w = -(det[a,v,w])^2        (a,v,w in C^3)
# is a polynomial identity in 9 free symbols.  This replaces the degree<=4
# generic checks of the earlier certificates by a degree-free proof.
#
# NEW REFORMULATION (S1.1).  det Hess_{n+1}( y*e(x) ) = -y^{n-1} B(e).
# Hence Theorem G reads:  "if the (n+1)-variable polynomial y*e(x) has
# vanishing Hessian then e has vanishing Hessian" -- and more generally
# (S1.1b, the y^{n-1}-coefficient of det Hess(y e + f) is -B(e)):
#
#   COROLLARY LV (linear-variable descent).  If F = y*e(x) + f(x) in
#   C[x1..xn, y] has det Hess_{n+1} F == 0, then det Hess_n e == 0; and for
#   n <= 3, e is affinely 2-variable (Corollary E).
#
# This is consistent with every zero-Hessian normal form of de Bondt's
# thesis (Thm 5.2.10 (i),(ii) and Prop. 5.3.1, Thm 5.3.3), in all of which
# the coefficient of a linear variable is visibly degenerate; and it is
# STRICTLY stronger than "det Hess e == 0" as a hypothesis on e: the dBvdE
# member e = x1x2 + x1^2x3 has det Hess_3 e == 0 but B(e) = -x1^4, so
# y*e + f has non-vanishing Hessian for EVERY f (S1.1c).
#
# ======================================================================
# PART 2.  GEOMETRIC / PROJECTIVE REFORMULATION AND LITERATURE
# ======================================================================
#
# Dictionary (known): B(e) == 0 <=> every level hypersurface has vanishing
# Gauss-Kronecker curvature [Reilly 1986 Prop. 4; Fox, Math. Nachr. 290
# (2017), U(F) = K |dF|^{n+2}]; for a projective hypersurface X = V(F),
# "developable" (degenerate Gauss map) <=> X is contained in its Hessian
# hypersurface, i.e. F | hess(F) [Dolgachev, CAG, Thm 1.1.20 and the
# discussion of parabolic points; Ciliberto-Russo-Simis, Adv. Math. 218
# (2008); De Poi-Ilardi, Ann. Univ. Paedagog. Crac. Stud. Math. 18 (2019)].
#
# PENCIL FORM OF THEOREM G (S2).  Let X_c := V(E - c x0^d) be the projective
# closure of the level set {e = c}.  Dolgachev's affine equation of the
# Hessian [CAG (1.20)-(1.21)], applied to e - c, gives
#     hess(E - c x0^d)|_{x0=1} = d(d-1)(e-c) det Hess_n e - (d-1)^2 B(e),
# so  X_c is developable  <=>  (e - c) | B(e).   Since the e - c are pairwise
# coprime, "(e-c) | B(e) for infinitely many c" forces B(e) == 0.  Hence:
#
#   THEOREM G (pencil form).  If infinitely many (equivalently: all, or a
#   generic) members of the pencil spanned by V(E) and d*H_inf have
#   degenerate Gauss map, then E has vanishing Hessian; more precisely
#   det Hess_n e == 0 and EVERY member E - c x0^d has vanishing Hessian.
#
# By Zak's theorem on the finiteness of the Gauss map of a smooth projective
# variety [Zak, Tangents and Secants, Cor. I.2.8], a developable X_c must be
# singular; for generic c its affine part is smooth, so B(e) == 0 forces
# every X_c to be singular at infinity: the leading forms satisfy
#     { e_d = 0, grad e_d = 0, e_{d-1} = 0 } != {0}.          (S2.4)
# This cheap necessary condition is checked on the B == 0 witnesses.
#
# LITERATURE VERDICT (search log in the structured report).  No statement
# of Theorem G, in either the affine (B == 0 => det Hess == 0) or the pencil
# form, was located.  The closest published items are: Fox 2017, Example
# `gnexample' (lit/fox2/f2.tex l.628-641), which exhibits U(P) == 0 WITH
# H(P) == 0 for the Gordan-Noether quintic-type P = a x3 + b x4 + c x5 as an
# "it can occur" remark, and Example `zamcsection' computing U(F) =
# -(ab' - a'b)^2 for F = a(u)x + b(u)y + Q(u) (the family (B) above); and
# Fassarella, "Logarithmic Hesse's problem", arXiv:1102.1659, Prop. P:Hesse,
# where linearity of the fibres of the polar map is derived from vanishing
# (logarithmic) Hessian -- the vanishing-Hessian analogue of the line-fibre
# mechanism (x'' + x' = 0) in the flow proof of Theorem G.  Neither states
# or implies Theorem G.
#
# ======================================================================
# PART 3.  CONSEQUENCES FOR THE VANISHING-HESSIAN LITERATURE
# ======================================================================
#
#   (i)   Corollary LV above (new, inhomogeneous; on FORMS it collapses to
#         Identity E, see S3.1).
#   (ii)  n = 3: B(e) == 0 <=> e affinely 2-variable (Corollary E) -- Hesse's
#         claim is TRUE for the invariant B in <= 3 variables and FALSE in
#         >= 4 (dehomogenized Perazzo, tgc_beyond_hc4.py R3).
#   (iii) Composite polynomials: B(phi(e)) = phi'(e)^{n+1} B(e) and
#         det Hess(phi(e)) = phi'^n det Hess e + phi'' phi'^{n-1} B(e)  (S3.2),
#         so Theorem G is compatible with composition and the level-set
#         statement is insensitive to composite polynomials.
#   (iv)  No new proof of Gordan-Noether or of the dBvdE classification is
#         obtained; on cones Theorem G is one line from Identity E.
#
# ======================================================================
# PART 4.  ATTRIBUTION OF IDENTITY E (Nagaoka-Yazawa)
# ======================================================================
#
# Source-verified from lit/nagaoka-yazawa/*.tex: theorem-style environments
# share one counter numbered by section (\newtheorem{thm}{..}[section]);
# Section 2 is "Homogeneous polynomials" (l.326); its numbered environments
# are Lemma 2.1 (lem:hom1), Lemma 2.2 (lem:hom2, Euler, credited to
# Anari-Gharan-Vinzant Cor. 4.3) and Proposition 2.3 (identity1):
#     det(-F H_F + s (grad F)^T grad F)
#         = (-1)^{n-1} (r/(r-1)) (s - (r-1)/r) F^n det H_F.
# The acknowledgement thanks H. Ochiai "for offering general facts on
# Proposition identity1" -- the authors themselves treat it as a general
# fact, not as their discovery.  S4 re-proves:
#   (a) NY <=> Identity E (both directions; rank-one expansion),
#   (b) the two CLASSICAL specialisations that are each equivalent to it:
#         log-Hessian:  det Hess(log F) = -(1/(d-1)) F^{-n} det Hess F,
#         power:        det Hess(F^k)   = k^n (kd-1)/(d-1) F^{n(k-1)} det Hess F,
#       (s = 1, resp. s = (k-1)/k... via Hess(F^k) = kF^{k-2}(F H + (k-1) g g^T)),
#   so any pre-2019 source of either specialisation is prior art for NY's
#   Proposition as well.  The bibliographic data in main_result.tex
#   (\bibitem{NY}: J. Algebra 577 (2021) 175-202, Prop. 2.3) is consistent
#   with the TeX source numbering; the journal volume/pages could not be
#   re-verified from a fetchable primary page in this session (ScienceDirect
#   403) -- see the report.

import itertools
import random
import sympy as sp

random.seed(20260824)
x0 = sp.Symbol('x0')
y = sp.Symbol('y')
c = sp.Symbol('c')
s = sp.Symbol('s')


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def grad(e, vs):
    return sp.Matrix([sp.diff(e, a) for a in vs])


def Bof(e, vs):
    H = hess(e, vs)
    g = grad(e, vs)
    return sp.expand((g.T * H.adjugate() * g)[0, 0])


def dethess(e, vs):
    return sp.expand(hess(e, vs).det(method='berkowitz'))


def ess_vars(e, vs):
    n = len(vs)
    cs = sp.symbols(f'vv0:{n}')
    Dv = sp.expand(sum(ci * sp.diff(e, xi) for ci, xi in zip(cs, vs)))
    if Dv == 0:
        return 0
    eqs = sp.Poly(Dv, *vs).coeffs()
    M, _ = sp.linear_eq_to_matrix(eqs, cs)
    return n - len(M.nullspace())


def homogenize(e, vs, d):
    P = sp.Poly(e, *vs)
    out = 0
    for mono, co in zip(P.monoms(), P.coeffs()):
        out += co * sp.prod([v**k for v, k in zip(vs, mono)]) * x0**(d - sum(mono))
    return sp.expand(out)


def sym_matrix(n, tag):
    M = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            M[i, j] = M[j, i] = sp.Symbol(f'{tag}{i}{j}')
    return M


def gen_poly(vs, deg, tag, homogeneous=False):
    degs = [deg] if homogeneous else range(0, deg + 1)
    mons = [m for k in degs for m in sorted({sp.prod(cc) for cc in
            itertools.combinations_with_replacement(vs, k)}, key=str)]
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(ci * m for ci, m in zip(cs, mons))


X2 = sp.symbols('x1:3')
X3 = sp.symbols('x1:4')
X4 = sp.symbols('x1:5')
X5 = sp.symbols('x1:6')
x1, x2, x3, x4, x5 = X5

# ======================================================================
print('=' * 74)
print('S1.1  det Hess_{n+1}( y e(x) ) = -y^{n-1} B(e)      (jet-level PROOF)')
print('=' * 74)
for n in (2, 3, 4):
    H = sym_matrix(n, f'h{n}_')
    g = sp.Matrix(sp.symbols(f'g{n}_0:{n}'))
    M = sp.zeros(n + 1, n + 1)
    M[:n, :n] = y * H
    M[:n, n] = g
    M[n, :n] = g.T
    lhs = sp.expand(M.det(method='berkowitz'))
    rhs = sp.expand(-y**(n - 1) * (g.T * H.adjugate() * g)[0, 0])
    assert sp.expand(lhs - rhs) == 0, f'S1.1 failed n={n}'
    print(f'  n={n}: Hess(y e) = [[yH, g],[g^T, 0]] has det = -y^(n-1) g^T adj(H) g   PROVED')
print('  => B(e) == 0  <=>  det Hess_{n+1}(y e(x)) == 0: Theorem G says that the')
print('     (n+1)-variable polynomial y*e(x) has vanishing Hessian only if e has.')

print()
print('S1.1b  y^{n-1}-coefficient of det Hess_{n+1}(y e + f) is -B(e)   (PROOF)')
for n in (2, 3):
    He = sym_matrix(n, f'p{n}_')
    Hf = sym_matrix(n, f'q{n}_')
    g = sp.Matrix(sp.symbols(f'r{n}_0:{n}'))
    M = sp.zeros(n + 1, n + 1)
    M[:n, :n] = y * He + Hf
    M[:n, n] = g
    M[n, :n] = g.T
    D = sp.expand(M.det(method='berkowitz'))
    assert sp.expand(D + (g.T * (y * He + Hf).adjugate() * g)[0, 0]) == 0
    top = sp.Poly(D, y)
    assert top.degree() <= n - 1
    assert sp.expand(top.coeff_monomial(y**(n - 1)) + (g.T * He.adjugate() * g)[0, 0]) == 0
    print(f'  n={n}: det Hess(y e + f) = -g^T adj(yH_e + H_f) g, degree <= n-1 in y,')
    print(f'         leading coefficient -B(e)                                    PROVED')
print('  => COROLLARY LV: det Hess_{n+1}(y e + f) == 0  =>  B(e) == 0  =(G)=>  det Hess_n e == 0.')

print()
print('S1.1c  the hypothesis of LV is strictly stronger than det Hess e == 0')
e_w = x1 * x2 + x1**2 * x3
f_gen = gen_poly(X3, 3, 'ff')
F = y * e_w + f_gen
D = dethess(F, X3 + (y,))
Dy = sp.Poly(D, y)
assert dethess(e_w, X3) == 0
assert sp.expand(Dy.coeff_monomial(y**2) - x1**4) == 0
print('  e = x1x2 + x1^2x3: det Hess_3 e == 0, but [y^2] det Hess_4(y e + f) = x1^4')
print('  for FULLY GENERIC cubic f  =>  y e + f never has vanishing Hessian.   PROVED')

print()
print('S1.1d  LV against the zero-Hessian normal forms of de Bondt\'s thesis')
# n=3 branch (B): linear variable x2, coefficient a2(x1)  (1 essential var)
A1, A2, A3 = (sum(sp.Integer(random.randint(-3, 3)) * x1**k for k in range(4)) for _ in range(3))
Fb = A1 + A2 * x2 + A3 * x3
assert dethess(Fb, X3) == 0 and dethess(A2, (x1, x3)) == 0
# n=4 Thm 5.3.3 (ii): a1(x1,x2) + a2 x3 + a3 x4 with a2, a3 in C[p]
p = x1**2 + x2
a1 = x1**3 * x2 + x2**2
a2 = p**2 + 1
a3 = p**3 - p
F4 = a1 + a2 * x3 + a3 * x4
assert dethess(F4, X4) == 0
assert dethess(a3, (x1, x2, x3)) == 0 and Bof(a3, (x1, x2, x3)) == 0
# n=4 Thm 5.3.3 (iii) with deg_t f = 1: linear in x4
f0, f1 = 1 + x1**2, x1**3 - 2
b1, b2, b3 = x1, x1**2 + 1, x1**3
al1, al2, al3 = x1, x1**2, x1 + 1
F4b = f0 + f1 * (al1 * x2 + al2 * x3 + al3 * x4) + b1 * x2 + b2 * x3 + b3 * x4
assert dethess(F4b, X4) == 0
e4b = sp.expand(f1 * al3 + b3)
assert dethess(e4b, (x1, x2, x3)) == 0
# Perazzo cubic, linear in x3 (n = 4 remaining variables)
per = x1**2 * x3 + x1 * x2 * x4 + x2**2 * x5
assert dethess(per, X5) == 0
assert dethess(x1**2, (x1, x2, x4, x5)) == 0
print('  thesis Prop 5.3.1(B), Thm 5.3.3(ii),(iii, deg_t f = 1), Perazzo cubic:')
print('  in every case the coefficient of the linear variable has det Hess == 0   PASS')

# ----------------------------------------------------------------------
print()
print('=' * 74)
print('S1.2  THE dBvdE FAMILY h = sum g_i(l) x_i, l = a.x: ALL DEGREES, ALL a')
print('=' * 74)
a = sp.Matrix(sp.symbols('a1c a2c a3c'))
v = sp.Matrix(sp.symbols('v1c v2c v3c'))
w = sp.Matrix(sp.symbols('w1c w2c w3c'))
M2 = a * v.T + v * a.T
lhs = sp.expand((w.T * M2.adjugate() * w)[0, 0])
rhs = sp.expand(-(sp.Matrix.hstack(a, v, w).det())**2)
assert sp.expand(lhs - rhs) == 0
print('  (i)  w^T adj(a v^T + v a^T) w = -(det[a,v,w])^2   for free a,v,w in C^3  PROVED')
# structural facts with symbolic functions g_i(l): grad and Hess of the family
L = sp.Symbol('L')
gfun = [sp.Function(f'g{i}')(L) for i in (1, 2, 3)]
lin = (a.T * sp.Matrix(X3))[0, 0]
h_sym = sum(gi.subs(L, lin) * xi for gi, xi in zip(gfun, X3))
g_l = sp.Matrix([gi.subs(L, lin) for gi in gfun])
gp_l = sp.Matrix([sp.diff(gi, L).subs(L, lin) for gi in gfun])
gpp_l = sp.Matrix([sp.diff(gi, L, 2).subs(L, lin) for gi in gfun])
xv = sp.Matrix(X3)
grad_pred = g_l + (gp_l.T * xv)[0, 0] * a
Hess_pred = a * gp_l.T + gp_l * a.T + (gpp_l.T * xv)[0, 0] * a * a.T
assert sp.simplify(grad(h_sym, X3) - grad_pred) == sp.zeros(3, 1)
assert sp.simplify(hess(h_sym, X3) - Hess_pred) == sp.zeros(3, 3)
print('  (ii) grad h = g + (g\'.x) a,  Hess h = a g\'^T + g\' a^T + (g\'\'.x) a a^T')
print('       for SYMBOLIC functions g_i and symbolic a (chain rule)           PROVED')
# Hess h = a v^T + v a^T with v = g' + (g''.x)/2 a; det[a, v, grad] = det[a, g', g]
vv = gp_l + sp.Rational(1, 2) * (gpp_l.T * xv)[0, 0] * a
assert sp.simplify(a * vv.T + vv * a.T - Hess_pred) == sp.zeros(3, 3)
d1 = sp.Matrix.hstack(a, vv, grad_pred).det()
d2 = sp.Matrix.hstack(a, gp_l, g_l).det()
assert sp.simplify(sp.expand(d1 - d2)) == 0
print('  (iii) Hess h = a v^T + v a^T, and det[a, v, grad h] = det[a, g\', g]    PROVED')
print('  => (*)  B(h) = -( det[a | g(l) | g\'(l)] )^2   for every degree and every a.')
# direct symbolic-function confirmation of (*)
Bsym = sp.expand((grad_pred.T * Hess_pred.adjugate() * grad_pred)[0, 0])
assert sp.simplify(sp.expand(Bsym + d2**2)) == 0
print('       direct symbolic-function confirmation of (*)                    PROVED')

print()
print('S1.3  B(h) == 0  <=>  h affinely 2-variable  (chain closure)')
# weight-2 covariance (matrix level) => WLOG a = e1; then det[e1,g,g'] = W(g2,g3)
T = sp.Matrix(3, 3, sp.symbols('t0:9'))
Mg = sym_matrix(3, 'm')
wv = sp.Matrix(sp.symbols('u0:3'))
lhs = sp.expand(((T.T * wv).T * (T.T * Mg * T).adjugate() * (T.T * wv))[0, 0])
rhs = sp.expand(T.det()**2 * (wv.T * Mg.adjugate() * wv)[0, 0])
assert sp.expand(lhs - rhs) == 0
print('  B(e o T) = det(T)^2 (B(e) o T): weight-2 covariance (matrix level)     PROVED')
e1 = sp.Matrix([1, 0, 0])
gg = sp.Matrix(sp.symbols('G1:4'))
gg_p = sp.Matrix(sp.symbols('Gp1:4'))
assert sp.expand(sp.Matrix.hstack(e1, gg, gg_p).det() - (gg[1] * gg_p[2] - gg[2] * gg_p[1])) == 0
print('  adapted a = e1:  det[e1 | g | g\'] = g2 g3\' - g3 g2\' = W(g2, g3)         PROVED')
print('  Wronskian criterion (char 0): W(g2,g3) == 0 <=> g2, g3 C-dependent')
print('  <=> some (0,alpha,beta) != 0 kills h  <=> h affinely 2-variable.')
print('  Exact instance sweep, degrees up to 8, random a and random translations:')
tL = sp.Symbol('tL')
n_deg = n_nondeg = 0
for k in range(24):
    avec = [sp.Integer(random.randint(-3, 3)) for _ in range(3)]
    if all(vv_ == 0 for vv_ in avec):
        avec[random.randint(0, 2)] = sp.Integer(1)
    lform = sum(vv_ * xx for vv_, xx in zip(avec, X3))
    degs = [random.randint(0, 8) for _ in range(3)]
    gsL = [sum(sp.Integer(random.randint(-3, 3)) * tL**j for j in range(dd + 1)) for dd in degs]
    if k % 4 == 0:      # force a degenerate member: g3 = const * g2
        gsL[2] = sp.Integer(random.choice([-2, 3])) * gsL[1]
    gsl = [sp.expand(GL.subs(tL, lform)) for GL in gsL]
    gpl = [sp.expand(sp.diff(GL, tL).subs(tL, lform)) for GL in gsL]
    h_i = sp.expand(sum(gi * xi for gi, xi in zip(gsl, X3)))
    # random translation (B, det Hess, ess_vars all translation invariant)
    tr = {xx: xx + sp.Integer(random.randint(-2, 2)) for xx in X3}
    h_i = sp.expand(h_i.subs(tr, simultaneous=True))
    Mdet = sp.expand(sp.Matrix([avec, gsl, gpl]).det().subs(tr, simultaneous=True))
    B = Bof(h_i, X3)
    assert dethess(h_i, X3) == 0, 'family member outside the zero-Hessian class'
    assert sp.expand(B + Mdet**2) == 0, 'closed form (*) failed on an instance'
    nv = ess_vars(h_i, X3)
    assert (B == 0) == (nv <= 2), f'CHAIN REFUTED on {h_i}'
    n_deg += (B == 0)
    n_nondeg += (B != 0)
print(f'    24 members (max deg 9): det Hess == 0 for all; B == 0 <=> ess.vars <= 2;')
print(f'    {n_deg} degenerate (B == 0), {n_nondeg} with 3 essential variables (B != 0)   PASS')

print()
print('  branch (A): e = a(x1,x2) + lam x3, B = lam^2 det Hess_2(a)  (jet-level)')
aj = sp.symbols('a1 a2 a11 a12 a22 lam')
Hj = sp.Matrix([[aj[2], aj[3], 0], [aj[3], aj[4], 0], [0, 0, 0]])
gj = sp.Matrix([aj[0], aj[1], aj[5]])
assert sp.expand((gj.T * Hj.adjugate() * gj)[0, 0] - aj[5]**2 * (aj[2] * aj[4] - aj[3]**2)) == 0
print('    PROVED.  With lam = 1, B == 0 <=> det Hess_2 a == 0 <=> (thesis Thm 5.2.10(i),')
print('    n = 2) a ~ a1(x1) + mu x2, so e = a1(x1) + mu x2 + x3 is killed by (0,1,-mu).')
for a_ in [x1**3 + 2 * x1, (x1 + 2 * x2)**4 + (x1 + 2 * x2), x1 * x2, x1**2 + x2**2]:
    e_ = a_ + x3
    B = Bof(e_, X3)
    assert (B == 0) == (ess_vars(e_, X3) <= 2)
print('    instances of branch (A): B == 0 <=> affinely 2-variable                 PASS')

# ======================================================================
print()
print('=' * 74)
print('S2.  PENCIL FORM OF THEOREM G')
print('=' * 74)
print('S2.1  Dolgachev\'s affine Hessian equation for e - c, symbolic c   (PROOF)')
for (n, d, vs) in ((2, 3, X2), (3, 2, X3), (3, 3, X3)):
    e = gen_poly(vs, d, f'hh{n}{d}_')
    E = homogenize(e, vs, d)
    lhs = sp.expand(dethess(sp.expand(E - c * x0**d), (x0,) + tuple(vs)).subs(x0, 1))
    rhs = sp.expand(d * (d - 1) * (e - c) * dethess(e, vs) - (d - 1)**2 * Bof(e, vs))
    assert sp.expand(lhs - rhs) == 0, (n, d)
    print(f'  n={n}, d={d}: hess(E - c x0^d)|_(x0=1) = d(d-1)(e-c) det Hess e - (d-1)^2 B(e)   PROVED')
print('  => X_c = V(E - c x0^d) is contained in its Hessian  <=>  (e - c) | B(e).')
print('     The e - c (c in C) are pairwise coprime, so infinitely many such c force')
print('     B(e) == 0; conversely B == 0 makes every X_c developable.')

print()
print('S2.2  B != 0 witness: NO member of its pencil is developable')
e = x1 * x2 + x1**2 * x3
E = homogenize(e, X3, 3)
Fc = sp.expand(E - c * x0**3)
hessFc = dethess(Fc, (x0,) + X3)
_, rem = sp.reduced(hessFc, [Fc], x0, x1, x2, x3, c, order='lex')
assert sp.expand(rem) != 0
# and no special c: (e - c) | -x1^4 would need e - c = const * x1^k
assert sp.Poly(e - c, x1, x2, x3).coeff_monomial(x1 * x2) == 1
print('  e = x1x2 + x1^2x3 (det Hess_3 e == 0, B = -x1^4): (E - c x0^3) does not divide')
print('  hess(E - c x0^3) for symbolic c, and no special c works (e - c has the')
print('  monomial x1x2, so it cannot divide x1^4)                                PASS')

print()
print('S2.3  B == 0 witnesses: EVERY member of the pencil has vanishing Hessian')
for e, vs, d in [(x2 + x1 * x3 + x1**2 * x4, X4, 3), (x1 * (x2 + x3), X3, 2),
                 (x1 * x2 + x1**2 * x3 + x4, X4, 3), (x1**3 + 3 * x1, X2, 3),  # FIX 2026-08-27: the agent's witness x1^3+3x1+x2 has B = 6 x1 != 0 (its level curves are cubics, not lines); x1^3+3x1 = phi(x1) is the intended B == 0 witness
                 ((x1 + x2 + 2 * x3)**4 + x1 + x2 + 2 * x3, X3, 4)]:
    assert Bof(e, vs) == 0
    assert dethess(e, vs) == 0                       # Theorem G
    E = homogenize(e, vs, d)
    assert dethess(sp.expand(E - c * x0**d), (x0,) + tuple(vs)) == 0
    print(f'    e = {str(e):32s}: B == 0, det Hess == 0, hess(E - c x0^{d}) == 0 (symbolic c)  PASS')

print()
print('S2.4  Zak-type necessary condition: singular points at infinity')
# For B == 0 the developable X_c is singular (Zak), the affine part of a generic
# X_c is smooth, so Sing X_c meets {x0 = 0}:  {e_d = 0, grad e_d = 0, e_{d-1} = 0} != {0}.
for e, vs, d in [(x2 + x1 * x3 + x1**2 * x4, X4, 3), (x1 * (x2 + x3), X3, 2),
                 (x1 * x2 + x1**2 * x3 + x4, X4, 3), (x1**3 + 3 * x1 + x2, X2, 3)]:
    P = sp.Poly(e, *vs)
    ed = sum(co * sp.prod([vv_**k for vv_, k in zip(vs, m)]) for m, co in zip(P.monoms(), P.coeffs()) if sum(m) == d)
    ed1 = sum(co * sp.prod([vv_**k for vv_, k in zip(vs, m)]) for m, co in zip(P.monoms(), P.coeffs()) if sum(m) == d - 1)
    eqs = [ed, ed1] + [sp.diff(ed, vv_) for vv_ in vs]
    G = sp.groebner([q for q in eqs if q != 0], *vs, order='grevlex')
    # the ideal is homogeneous; its zero set is {0} iff every variable has a pure power in LT
    lts = [sp.Poly(q, *vs).monoms()[0] for q in G.exprs]
    pure = all(any(sum(m) == m[i] and m[i] > 0 for m in lts) for i in range(len(vs)))
    assert not pure, f'no singular point at infinity for {e}'
    print(f'    e = {str(e):28s}: {{e_d = grad e_d = e_(d-1) = 0}} is a nontrivial cone   PASS')

# ======================================================================
print()
print('=' * 74)
print('S3.  PLACEMENT IN THE VANISHING-HESSIAN LITERATURE')
print('=' * 74)
print('S3.1  on FORMS, Corollary LV collapses to Identity E')
for (n, d) in ((2, 3), (3, 2), (3, 3)):
    vs = X2 if n == 2 else X3
    Fh = gen_poly(vs, d, f'kk{n}{d}_', homogeneous=True)
    assert sp.expand(Bof(Fh, vs) - sp.Rational(d, d - 1) * Fh * dethess(Fh, vs)) == 0
    assert sp.expand(dethess(y * Fh, tuple(vs) + (y,)) + y**(n - 1) * sp.Rational(d, d - 1) * Fh * dethess(Fh, vs)) == 0
    print(f'  n={n}, d={d}: det Hess(y F) = -y^(n-1) (d/(d-1)) F det Hess F  (F a form)   PROVED')
print('  => for forms, "det Hess(yF) == 0 <=> det Hess F == 0" is classical (Identity E);')
print('     the content of Theorem G is the inhomogeneous case only.')

print()
print('S3.2  composition: B(phi(e)) = phi\'^(n+1) B(e),')
print('      det Hess(phi(e)) = phi\'^n det Hess e + phi\'\' phi\'^(n-1) B(e)   (jet PROOF)')
for n in (2, 3):
    H = sym_matrix(n, f'cH{n}_')
    g = sp.Matrix(sp.symbols(f'cg{n}_0:{n}'))
    p1, p2 = sp.symbols('phi1 phi2')          # phi'(e), phi''(e) at the point
    Hc = p1 * H + p2 * g * g.T
    gc = p1 * g
    Bc = sp.expand((gc.T * Hc.adjugate() * gc)[0, 0])
    assert sp.expand(Bc - p1**(n + 1) * (g.T * H.adjugate() * g)[0, 0]) == 0
    Dc = sp.expand(Hc.det(method='berkowitz'))
    assert sp.expand(Dc - p1**n * H.det(method='berkowitz') - p2 * p1**(n - 1) * (g.T * H.adjugate() * g)[0, 0]) == 0
    print(f'  n={n}: both identities                                                  PROVED')
print('  => B(phi(e)) == 0 <=> B(e) == 0 (phi non-constant): developability of level')
print('     sets is a property of the fibration, and Theorem G is compatible with it.')

# ======================================================================
print()
print('=' * 74)
print('S4.  ATTRIBUTION: Nagaoka-Yazawa Prop. 2.3 and its classical avatars')
print('=' * 74)
for (n, d, vs) in ((2, 2, X2), (2, 3, X2), (2, 4, X2), (3, 2, X3), (3, 3, X3)):
    Fh = gen_poly(vs, d, f'ny{n}{d}_', homogeneous=True)
    H = hess(Fh, vs)
    g = grad(Fh, vs)
    detH = H.det(method='berkowitz')
    NY_l = sp.expand((-Fh * H + s * g * g.T).det(method='berkowitz'))
    NY_r = sp.expand((-1)**(n - 1) * sp.Rational(d, d - 1) * (s - sp.Rational(d - 1, d)) * Fh**n * detH)
    assert sp.expand(NY_l - NY_r) == 0
    # log-Hessian:  Hess(log F) = (F H - g g^T)/F^2
    logH = (Fh * H - g * g.T) / Fh**2
    assert sp.simplify(sp.expand(logH.det(method='berkowitz') * Fh**n * (d - 1) + detH)) == 0
    # power:  Hess(F^k) = k F^(k-2) (F H + (k-1) g g^T)
    k = sp.Symbol('k')
    Hk = k * (Fh * H + (k - 1) * g * g.T)          # times F^(k-2), stripped
    assert sp.expand(Hk.det(method='berkowitz') - k**n * Fh**(n - 1) * sp.Rational(1, d - 1) * (k * d - 1) * Fh * detH) == 0
    print(f'  n={n}, d={d}: NY identity; det Hess(log F) = -(1/(d-1)) F^-n det Hess F;')
    print(f'            det Hess(F^k) = k^n (kd-1)/(d-1) F^(n(k-1)) det Hess F          PROVED')
print('  All three are the SAME affine-in-s identity read at s = generic, s = 1,')
print('  s = (k-1)/k... (rank-one expansion: det(M + s u v^T) is affine in s), so')
print('  each is formally equivalent to Identity E and to NY Prop. 2.3.')

print()
print('ALL CHECKS PASSED  (litG).')
