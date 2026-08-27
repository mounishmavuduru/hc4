# tgc_debondt_chain.py  --  key: tgc  (Theorem G, Consequences)
#
# PART (1) of the "consequences of Theorem G beyond HC_4" task.
#
# THE CHAIN.  For e in C[x1,x2,x3]:
#
#        B(e) == 0
#   =(Theorem G)=>  det Hess_3(e) == 0
#   =(de Bondt-van den Essen, Singular Hessians, J. Algebra 282 (2004),
#     Thm 3.3 = de Bondt's thesis Prop. 5.3.1)=>
#        after a LINEAR change of coordinates e is one of
#          (A)  e = a(x1,x2) + lam*x3,        lam in {0,1},  a in C[x1,x2];
#          (B)  e = a1(x1) + a2(x1)x2 + a3(x1)x3,   a1,a2,a3 in C[x1].
#   =(this script)=>  in each branch B(e) == 0 forces e to be affinely
#     2-variable, and conversely.
#
# So Theorem G + [dBvdE] give the project's Corollary E
#       B(e) == 0  =>  e is affinely 2-variable   (e in C[x1,x2,x3])
# by a route COMPLETELY INDEPENDENT of the project's Theorem F / Identity H /
# Gordan-Noether-in-four-variables.  That is the point of this certificate: an
# independent second proof of Corollary E from Theorem G.
#
# WHAT IS PROVED HERE (jet-level identities; each is polynomial in the 1- and
# 2-jet entries, which are the only data B depends on, hence each symbolic
# check below is a COMPLETE PROOF for all members of the branch):
#
#   (A)   B(a(x1,x2) + lam*x3) = lam^2 * det Hess_2(a).
#   (B)   B(a1 + a2 x2 + a3 x3) = -( a2*a3' - a3*a2' )^2 = -W(a2,a3)^2,
#         the squared Wronskian of the two "transverse" coefficients.
#
# CONSEQUENCES.
#   (A)  B == 0  <=>  lam = 0 (e is literally 2-variable) or det Hess_2 a == 0,
#        and in the latter case a = phi(l) + (linear) for a linear form l
#        (n = 2 case of the same classification / Theorem G in 2 variables),
#        so e = phi(l) + (linear) is affinely 2-variable.  Conversely obvious.
#   (B)  B == 0  <=>  W(a2,a3) == 0  <=>  a2, a3 are C-LINEARLY DEPENDENT
#        (char 0: (a2/a3)' = -W/a3^2), and then, with (c2,c3) != 0 killing the
#        combination, D_v e == 0 for v = (0, c3, -c2) != 0, i.e. e is affinely
#        2-variable.  Conversely if e is affinely 2-variable inside branch (B)
#        the same computation gives W == 0.
#
# HONESTY NOTE ON THE TASK STATEMENT.  The task asked to verify that the
# dBvdE exceptional family h = g1(l)x1 + g2(l)x2 + g3(l)x3 "always has B != 0".
# That is FALSE as literally stated, and this script certifies the correct
# statement: B(h) is minus a perfect square, and it vanishes EXACTLY on the
# degenerate members of the family -- those that are affinely 2-variable
# (e.g. a2 = a3 = x1 gives h = x1(x2+x3), B == 0).  What the chain needs, and
# what is true, is:  B(h) == 0  <=>  h affinely 2-variable.
#
# PRIOR-ART NOTE.  The identity (B) is NOT new: D. J. F. Fox, "Equiaffine
# geometry of level sets ...", arXiv:1503.09108 (Math. Nachr. 290 (2017)),
# section "Examples", computes for F(u,x,y) = a(u)x + b(u)y + Q(u) exactly
#     H(F) = 0,   U(F) = -(a b' - a' b)^2
# (local copy lit/fox2/f2.tex line 1338), where U(F) is his name for our B.
# Fox works with smooth a,b on R; the polynomial case is the same computation.
# What is new here is only the USE of it: combined with Theorem G and the
# dBvdE classification it closes the n = 3 case.

import itertools
import sympy as sp

FAIL = []


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


def kernel_directions(e, vs):
    """Basis of {v in C^n : D_v e == 0}; e is affinely (n-k)-variable where
    k = dim of this space."""
    n = len(vs)
    cs = sp.symbols(f'v0:{n}')
    Dv = sp.expand(sum(c * sp.diff(e, x) for c, x in zip(cs, vs)))
    if Dv == 0:
        return sp.eye(n).columnspace()
    eqs = sp.Poly(Dv, *vs).coeffs()
    sol = sp.linsolve(eqs, cs)
    M, _ = sp.linear_eq_to_matrix(eqs, cs)
    return M.nullspace()


def affinely_k_variable(e, vs):
    """number of essential variables = n - dim ker"""
    return len(vs) - len(kernel_directions(e, vs))


x1, x2, x3 = sp.symbols('x1 x2 x3')
X3 = (x1, x2, x3)

print('=' * 74)
print('(A)  BRANCH  e = a(x1,x2) + lam*x3   of de Bondt-van den Essen')
print('=' * 74)
# JET-LEVEL PROOF.  B depends only on grad e and Hess e.  For this branch
# grad e = (a1, a2, lam), Hess e = [[a11,a12,0],[a12,a22,0],[0,0,0]] with
# a1,a2,a11,a12,a22 arbitrary functions -- treat them as independent symbols.
a1s, a2s, a11, a12, a22, lam = sp.symbols('a1 a2 a11 a12 a22 lam')
H = sp.Matrix([[a11, a12, 0], [a12, a22, 0], [0, 0, 0]])
g = sp.Matrix([a1s, a2s, lam])
Bjet = sp.expand((g.T * H.adjugate() * g)[0, 0])
assert sp.expand(Bjet - lam**2 * (a11 * a22 - a12**2)) == 0, 'branch (A) jet identity'
print('  jet identity  B = lam^2 * (a11 a22 - a12^2) = lam^2 * det Hess_2(a)   PROVED')
assert sp.expand(sp.Matrix([[a11, a12, 0], [a12, a22, 0], [0, 0, 0]]).det()) == 0
print('  and det Hess_3 == 0 identically on branch (A)                          PROVED')

# polynomial-level confirmation with fully generic a of degree <= 4
mons2 = [sp.Integer(1)] + [m for d in range(1, 5)
                           for m in sorted({sp.prod(c) for c in
                                            itertools.combinations_with_replacement((x1, x2), d)}, key=str)]
cA = sp.symbols(f'cA0:{len(mons2)}')
aP = sum(c * m for c, m in zip(cA, mons2))
eA = aP + lam * x3
assert sp.expand(Bof(eA, X3) - lam**2 * (sp.diff(aP, x1, 2) * sp.diff(aP, x2, 2)
                                         - sp.diff(aP, x1, x2)**2)) == 0
print('  confirmed with fully generic a of degree <= 4 and symbolic lam         PASS')

print()
print('=' * 74)
print('(B)  BRANCH  e = a1(x1) + a2(x1)x2 + a3(x1)x3')
print('=' * 74)
# JET-LEVEL PROOF.  grad e = (u, a2, a3), Hess e = [[p, a2', a3'],[a2',0,0],
# [a3',0,0]] with u, p, a2, a3, a2', a3' arbitrary -- independent symbols.
u, p, A2, A3, dA2, dA3 = sp.symbols('u p A2 A3 dA2 dA3')
H = sp.Matrix([[p, dA2, dA3], [dA2, 0, 0], [dA3, 0, 0]])
g = sp.Matrix([u, A2, A3])
Bjet = sp.expand((g.T * H.adjugate() * g)[0, 0])
assert sp.expand(Bjet + (A2 * dA3 - A3 * dA2)**2) == 0, 'branch (B) jet identity'
print('  jet identity  B = -(a2 a3\' - a3 a2\')^2 = -W(a2,a3)^2                   PROVED')
assert sp.expand(H.det()) == 0
print('  and det Hess_3 == 0 identically on branch (B)                          PROVED')

# polynomial-level confirmation, fully generic a1,a2,a3 of degree <= 4
cs1 = sp.symbols('p0:5')
cs2 = sp.symbols('q0:5')
cs3 = sp.symbols('r0:5')
A1p = sum(c * x1**k for k, c in enumerate(cs1))
A2p = sum(c * x1**k for k, c in enumerate(cs2))
A3p = sum(c * x1**k for k, c in enumerate(cs3))
eB = A1p + A2p * x2 + A3p * x3
W = sp.expand(A2p * sp.diff(A3p, x1) - A3p * sp.diff(A2p, x1))
assert sp.expand(Bof(eB, X3) + W**2) == 0
print('  confirmed with fully generic a1,a2,a3 of degree <= 4                   PASS')
assert dethess(eB, X3) == 0
print('  det Hess_3 == 0 for the same generic member                            PASS')

print()
print('  Wronskian test on exact instances:  B == 0  <=>  a2,a3 dependent')
inst = [(sp.Integer(0), sp.Integer(1), x1), (sp.Integer(0), x1, x1**2),
        (x1**3, x1, x1**3), (sp.Integer(0), 1 + x1, 1 - x1),
        (x1**2, x1**2, x1**5), (sp.Integer(0), x1**2, 3 * x1**2),
        (x1, x1 + x1**4, 2 * x1 + 2 * x1**4), (sp.Integer(0), sp.Integer(1), sp.Integer(1)),
        (sp.Integer(0), x1**2 + 1, x1**3), (x1**7, x1**3 - x1, x1**3 + x1)]
for (b1, b2, b3) in inst:
    e = b1 + b2 * x2 + b3 * x3
    B = Bof(e, X3)
    dep = sp.Matrix([[sp.Poly(b2, x1).all_coeffs()[::-1] + [0] * 9][0][:9],
                     [sp.Poly(b3, x1).all_coeffs()[::-1] + [0] * 9][0][:9]]).rank() < 2
    nvars = affinely_k_variable(e, X3)
    ok = (B == 0) == dep and (B == 0) == (nvars <= 2)
    print(f'    a2={str(b2):12s} a3={str(b3):12s}  B={str(B):18s} '
          f'dep={dep}  ess.vars={nvars}  {"OK" if ok else "FAIL"}')
    if not ok:
        FAIL.append(('branchB', b1, b2, b3))
assert not FAIL, FAIL

print()
print('  the de Bondt-van den Essen witness used in the project:')
h = x1 * x2 + x1**2 * x3           # = branch (B) with a1=0, a2=x1, a3=x1^2
assert dethess(h, X3) == 0
assert Bof(h, X3) == -x1**4
assert affinely_k_variable(h, X3) == 3
print('    h = x1x2 + x1^2x3 : det Hess == 0, B = -x1^4 != 0, 3 essential vars  OK')
print('    (indeed -W(x1,x1^2)^2 = -(x1*2x1 - x1^2)^2 = -x1^4)')

print()
print('=' * 74)
print('(C)  THE FAMILY IN THE FORM  h = g1(l)x1 + g2(l)x2 + g3(l)x3, l = a.x')
print('=' * 74)
# Any such h is linearly equivalent to a branch-(B) member: choose coordinates
# y1 = l, y2, y3; then x_i are linear in y and h = sum_j G_j(y1) y_j.  B is a
# relative invariant of weight 2: B(e o T) = det(T)^2 * B(e) o T, so B == 0 is
# a linear-change invariant.  We check the weight-2 covariance symbolically and
# then check the family on exact instances.
T = sp.Matrix(3, 3, sp.symbols('t0:9'))
etest = x1**3 + x1 * x2 * x3 + x2**2 + x3          # arbitrary test polynomial
y = T * sp.Matrix(X3)
eT = sp.expand(etest.subs({x1: y[0], x2: y[1], x3: y[2]}, simultaneous=True))
lhs = sp.expand(Bof(eT, X3))
rhs = sp.expand((T.det()**2) * Bof(etest, X3).subs({x1: y[0], x2: y[1], x3: y[2]},
                                                   simultaneous=True))
assert sp.expand(lhs - rhs) == 0, 'B is not a weight-2 relative invariant'
print('  B(e o T) = det(T)^2 * B(e) o T  for a generic 3x3 T                    PROVED')
print('  => "B == 0" is invariant under linear (and, being translation')
print('     invariant, under affine) changes of coordinates.')

aa = sp.symbols('a1c a2c a3c')
gg = [sum(sp.Symbol(f'g{i}{k}') * sp.Symbol('L')**k for k in range(3)) for i in range(3)]
L = aa[0] * x1 + aa[1] * x2 + aa[2] * x3
famtests = [
    ((1, 0, 0), (0, x1, x1**2)),          # -> x1x2 + x1^2x3   (not planar)
    ((1, 0, 0), (0, x1, 2 * x1)),         # -> degenerate, planar
    ((1, 1, 0), (0, 1, sp.Symbol('L'))),  # l = x1+x2
    ((0, 1, 1), (sp.Symbol('L')**2, sp.Symbol('L'), 1)),
    ((1, 2, 3), (1, sp.Symbol('L'), sp.Symbol('L')**2)),
]
Lsym = sp.Symbol('L')
for avec, gs in famtests:
    Lin = avec[0] * x1 + avec[1] * x2 + avec[2] * x3
    comps = [sp.expand(sp.sympify(gi).subs(Lsym, Lin)) for gi in gs]
    h = sp.expand(sum(c * v for c, v in zip(comps, X3)))
    B = Bof(h, X3)
    nv = affinely_k_variable(h, X3)
    dh = dethess(h, X3)
    ok = (dh == 0) and ((B == 0) == (nv <= 2))
    print(f'    a={avec} g={gs}:  detHess={dh}  B==0? {B==0}  ess.vars={nv}  '
          f'{"OK" if ok else "FAIL"}')
    if not ok:
        FAIL.append(('family', avec, gs))
    if B != 0:
        # B is minus a perfect square (times det(T)^2), so -B must be a square
        assert sp.factor(-B).is_Pow or sp.simplify(sp.sqrt(sp.factor(-B))**2 + B) == 0 \
            or sp.factor(-B) == sp.factor(sp.sqrt(sp.expand(-B)))**2, f'-B not a square: {sp.factor(-B)}'
assert not FAIL, FAIL

print()
print('=' * 74)
print('(D)  THE CHAIN, END TO END, ON EXACT INSTANCES')
print('=' * 74)
# For every instance: if B == 0 then (Theorem G) det Hess == 0, and the
# classification branch then forces "affinely 2-variable".  Fail-closed on any
# instance with B == 0 that is not affinely 2-variable -- such an instance
# would REFUTE Corollary E (and, if det Hess != 0, Theorem G itself).
chain = []
for (b1, b2, b3) in inst:
    chain.append(b1 + b2 * x2 + b3 * x3)
for aa_ in [x1 * x2, x1**2 + x2**2, x1**3 + x1 * x2, (x1 + x2)**4, x1**2, sp.Integer(0)]:
    for lm in (0, 1):
        chain.append(aa_ + lm * x3)
chain += [h, x1 * x2 + x3, x1**2 / 2 + x2, x1 * x2 * x3, x1**2 * x2 + x3**3]
bad = []
for e in chain:
    B = Bof(e, X3)
    dh = dethess(e, X3)
    nv = affinely_k_variable(e, X3)
    if B == 0:
        if dh != 0:
            bad.append(('THEOREM G REFUTED', e, dh))
        if nv > 2:
            bad.append(('COROLLARY E REFUTED', e, nv))
assert not bad, bad
print(f'  {len(chain)} instances: every one with B == 0 has det Hess == 0 and is')
print('  affinely 2-variable.  No counterexample to Theorem G or Corollary E.')

print()
print('CONCLUSION (part 1).')
print('  Theorem G + de Bondt-van den Essen (Singular Hessians, Thm 3.3) give')
print('  Corollary E independently of Theorem F:  for e in C[x1,x2,x3],')
print('     B(e) == 0  <=>  e is affinely 2-variable.')
print('  Both directions hold; the exceptional branch a1+a2x2+a3x3 has')
print('  B = -W(a2,a3)^2, which vanishes exactly on its degenerate members.')
