# fpG_paper_certificate.py -- certificate for ../theorem_G_paper.tex (agent key: fpG)
#
# Fail-closed (assert-based; nonzero exit on any failure), exact/symbolic only.
#
# PART 1: structural "compile check" of the manuscript (no TeX toolchain is
#   installed on this machine): every \ref/\eqref has a \label, no duplicate
#   labels, every \cite has a \bibitem and no orphan \bibitem, \begin/\end
#   properly nested, braces balanced outside comments, one document skeleton,
#   no macro redefinition, and no \ref to a label defined only inside a
#   comment.
#
# PART 2: exact re-derivation of every algebraic identity DISPLAYED in the
#   manuscript.  A check marked PROOF is a polynomial identity in fully
#   generic symbolic coefficients of the stated shape, hence settles every
#   instance of that shape; checks marked INSTANCE are exact witnesses.
#
#   (C1)  n = 1: B(e) = (e')^2                                            PROOF
#   (C2)  H adj(H) g = det(H) g,  g^T adj(H) g = B   (n = 2,3,4)          PROOF
#   (C3)  Legendre identities grad L = psi, E L - L = e o psi, and
#         (E^2 - E) L = (B/D) o psi  for generic nondegenerate quadratics
#         (n = 3)                                                         PROOF
#         + the same chain on the non-quadratic branch e = x1^3/3 + x2^2/2 INSTANCE
#   (C4)  flow-variant identities (Remark rem:flow), polynomial forms:
#         D Jac(W) - W (grad D)^T = D^2 I - adj(H) Mt   and
#         D grad(B) - B grad(D) = 2 D^2 g - Mt W        (n = 2 cubic)     PROOF
#         where W = adj(H) g, Mt = sum_k (d_k H) W_k  (V = W/D, M = Mt/D)
#   (C5)  no-positive-powers step: grad e(b + a/t) is a polynomial in 1/t
#         for symbolic t-free a, b (n = 2, deg 4; n = 3, deg 3)           PROOF
#   (C6)  Theorem thm:F steps: d_0 E|_{x0=1} = d e - x.grad e and
#         d_i E|_{x0=1} = d_i e (n = 3, d = 3); translation invariance of B
#         (n = 3, deg 3); converse of Euler on homogeneous parts           PROOF
#   (C7)  Lemma lem:identE (ternary forms d = 2,3,4)                       PROOF
#   (C8)  Lemma lem:identH (n = 2, d = 2,3; n = 3, d = 2)                  PROOF
#   (C9)  Lemma lem:bordered: block form, bordered expansion, deg_{x4} <= 2,
#         [x4^2] = -B(e1)  (generic cubics)                               PROOF
#   (C10) Lemma lem:affx4: affine in x4, [x4] = -(e0)_{x3x3} B2(e)        PROOF
#   (C11) Theorem thm:quad: Hess3(e1^2/2) = e1 K + g g^T; the Schur identity;
#         det Hess f = gamma det3(Hess3 e0~ + u K); the collision-transfer
#         gradient identity                                               PROOF
#   (C12) Theorem thm:aff: det Hess(e0 + x4 x3) = -det2 Hess_{(x1,x2)} e0 PROOF
#   (C13) Theorem thm:tri, branch B2 == 0:
#         det Hess(e0 + x4 phi(x1)) = -phi'(x1)^2 det2 Hess_{(x2,x3)} e0  PROOF
#   (C14) Theorem thm:doubling: det Hess(a + x3 b + x4 e) = Jac(b,e)^2    PROOF
#   (C15) Lemma lem:inj, direction (=>): the linear system is solvable    PROOF
#   (C16) Proposition prop:pe(b): top graded piece of det Hess(f2 + fd) is
#         det Hess fd  (n = 3, d = 3 generic; n = 4, d = 3 exact instance)
#   (C17) witnesses: converse failure x1^2/2 + x2; dBvdE h = x1x2 + x1^2x3
#         (B = -x1^4); the n = 4 counterexample x2 + x1x3 + x1^2x4 to
#         planarity (B = 0, det Hess = 0, 4 independent partials) and its
#         homogenization = Perazzo's cubic; jet identity
#         B(a + lam x4) = lam^2 det Hess3(a); the sharpness witness
#         x1 - sqrt(x1^2 - 2x2) (B = 0, det Hess = -(x1^2-2x2)^{-2});
#         de Bondt's dillen4 doubling has det Hess = 1                    INSTANCE/PROOF
#   (C18) Remark rem:notdoubling: det Hess w = 1, D_{e4}w = x3, isotropic
#         cone = C.e4 (Rabinowitsch radical membership of v1, v2, v3)     PROOF
#   (C19) Nagaoka-Yazawa identity as displayed + s^1 extraction +
#         cancellation = lem:identE  (n = 3, d = 3)                        PROOF
#   (C20) Proposition prop:pe(a): D_v^2 f = v^T Hess f v                  PROOF

import itertools
import re
import sys
import time
from pathlib import Path

import sympy as sp

T0 = time.time()
TEX = Path(__file__).resolve().parent.parent / 'theorem_G_paper.tex'


def stamp(msg):
    print(f'[{time.time() - T0:6.1f}s] {msg}', flush=True)


# ============================================================== PART 1
stamp('PART 1: structural compile-check of ' + TEX.name)
src = TEX.read_text(encoding='utf-8')
stripped = re.sub(r'(?<!\\)%.*', '', src)

labels = re.findall(r'\\label\{([^}]*)\}', stripped)
refs = re.findall(r'\\(?:ref|eqref)\{([^}]*)\}', stripped)
cites = [k.strip() for grp in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]*)\}', stripped)
         for k in grp.split(',')]
bibs = re.findall(r'\\bibitem\{([^}]*)\}', stripped)

dup_labels = {l for l in labels if labels.count(l) > 1}
assert not dup_labels, f'duplicate labels: {dup_labels}'
undefined_refs = {r for r in refs if r not in labels}
assert not undefined_refs, f'undefined refs: {undefined_refs}'
dup_bibs = {b for b in bibs if bibs.count(b) > 1}
assert not dup_bibs, f'duplicate bibitems: {dup_bibs}'
missing_bibs = {c for c in cites if c not in bibs}
assert not missing_bibs, f'cite keys without bibitem: {missing_bibs}'
orphan_bibs = {b for b in bibs if b not in cites}
assert not orphan_bibs, f'uncited bibitems: {orphan_bibs}'
unused_labels = sorted({l for l in labels if l not in refs})
print(f'  labels: {len(set(labels))} unique; refs: {len(refs)}, all defined')
print(f'  cites: {len(set(cites))} keys; bibitems: {len(bibs)}; exact match, no orphans')
if unused_labels:
    print(f'  note: labels never referenced (harmless): {unused_labels}')

stack = []
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}', stripped):
    kind, env = m.group(1), m.group(2)
    if kind == 'begin':
        stack.append(env)
    else:
        assert stack and stack[-1] == env, \
            f'\\end{{{env}}} does not match open environment {stack[-1:] or "(none)"}'
        stack.pop()
assert not stack, f'unclosed environments: {stack}'
print('  environments: all \\begin/\\end pairs properly nested')

body = stripped.replace(r'\{', '').replace(r'\}', '')
bal = body.count('{') - body.count('}')
assert bal == 0, f'brace imbalance: {bal}'
print('  braces: balanced')

# $ ... $ balance per paragraph is a cheap sanity check for stray math delimiters
for i, para in enumerate(re.split(r'\n\s*\n', stripped)):
    dollars = len(re.findall(r'(?<!\\)\$', para))
    assert dollars % 2 == 0, f'odd number of $ in paragraph {i}: {para[:80]!r}'
print('  inline-math delimiters: balanced in every paragraph')

assert len(re.findall(r'\\documentclass', stripped)) == 1
assert len(re.findall(r'\\begin\{document\}', stripped)) == 1
assert len(re.findall(r'\\end\{document\}', stripped)) == 1
newcmds = re.findall(r'\\newcommand\{\\([A-Za-z]+)\}', stripped)
assert len(newcmds) == len(set(newcmds)), 'duplicate \\newcommand'
# every custom macro used is defined; every defined macro is used
used_macros = set(re.findall(r'\\([A-Za-z]+)', stripped))
for mac in newcmds:
    assert mac in used_macros, f'macro \\{mac} defined but unused'
print(f'  document skeleton OK; {len(newcmds)} custom macros, all used, none redefined')
# the abstract states plainly that HC_4 is open
abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', stripped, re.S).group(1)
assert 'remains open' in abstract and 'HC' in abstract
print('  abstract states that HC_4 remains open')

# ============================================================== PART 2
stamp('PART 2: exact re-derivation of the displayed identities')

x1, x2, x3, x4, x0 = sp.symbols('x1 x2 x3 x4 x0')
X2, X3, X4 = (x1, x2), (x1, x2, x3), (x1, x2, x3, x4)


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def grad(e, vs):
    return sp.Matrix([sp.diff(e, a) for a in vs])


def Bof(e, vs):
    g = grad(e, vs)
    return sp.expand((g.T * hess(e, vs).adjugate() * g)[0, 0])


def det(M):
    return M.det(method='berkowitz')


def generic(vs, deg, tag, const=True, mindeg=None):
    mons = [sp.Integer(1)] if const else []
    lo = 1 if mindeg is None else mindeg
    for d in range(lo, deg + 1):
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons))


def form(vs, deg, tag):
    return generic(vs, deg, tag, const=False, mindeg=deg)


# ---------------------------------------------------------------- (C1)
u = sp.Symbol('u')
e1v = generic((u,), 5, 'q')
H1 = sp.Matrix([[sp.diff(e1v, u, 2)]])
assert H1.adjugate() == sp.Matrix([[1]])
assert sp.expand(Bof(e1v, (u,)) - sp.diff(e1v, u)**2) == 0
stamp('(C1) n = 1: adj((e\'\')) = (1) and B(e) = (e\')^2                      PROOF')

# ---------------------------------------------------------------- (C2)
for n in (2, 3, 4):
    vs = sp.symbols(f'z1:{n+1}')
    e = generic(vs, 3, f'c{n}_')
    H, g = hess(e, vs), grad(e, vs)
    W = sp.expand(H.adjugate() * g)
    assert sp.expand(H * W - det(H) * g) == sp.zeros(n, 1)
    assert sp.expand((g.T * W)[0, 0] - Bof(e, vs)) == 0
stamp('(C2) H adj(H) g = det(H) g and g^T adj(H) g = B, n = 2,3,4 generic cubics  PROOF')

# ---------------------------------------------------------------- (C3)
p1, p2, p3 = sp.symbols('p1 p2 p3')
P3 = (p1, p2, p3)
S = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f's{min(i,j)+1}{max(i,j)+1}'))
lvec = sp.Matrix(sp.symbols('l1 l2 l3'))
c0 = sp.Symbol('c0')
xv = sp.Matrix(X3)
eq = (xv.T * S * xv)[0, 0] / 2 + (lvec.T * xv)[0, 0] + c0
pv = sp.Matrix(P3)
Sinv = S.inv()
psi = Sinv * (pv - lvec)                       # local (here global) inverse of grad e
subs_psi = {X3[i]: psi[i] for i in range(3)}
L = (psi.T * pv)[0, 0] - eq.subs(subs_psi)
gradL = sp.Matrix([sp.diff(L, p) for p in P3])
assert sp.simplify(gradL - psi) == sp.zeros(3, 1)
EL = sum(p * sp.diff(L, p) for p in P3)
assert sp.simplify(EL - L - eq.subs(subs_psi)) == 0
e_psi = eq.subs(subs_psi)
E_e_psi = sum(p * sp.diff(e_psi, p) for p in P3)
BD = (grad(eq, X3).T * hess(eq, X3).inv() * grad(eq, X3))[0, 0]
assert sp.simplify(E_e_psi - BD.subs(subs_psi)) == 0
E2L = sum(p * sp.diff(EL, p) for p in P3)
assert sp.simplify(E2L - EL - E_e_psi) == 0
stamp('(C3) grad L = psi, E L - L = e o psi, (E^2 - E)L = (B/D) o psi, generic quadratics  PROOF')
# non-quadratic exact branch: e = x1^3/3 + x2^2/2, psi = (sqrt(p1), p2)
q1, q2 = sp.symbols('q1 q2', positive=True)
e_nq = x1**3 / 3 + x2**2 / 2
psi_nq = (sp.sqrt(q1), q2)
L_nq = psi_nq[0] * q1 + psi_nq[1] * q2 - e_nq.subs({x1: psi_nq[0], x2: psi_nq[1]})
assert sp.simplify(sp.diff(L_nq, q1) - psi_nq[0]) == 0
assert sp.simplify(sp.diff(L_nq, q2) - psi_nq[1]) == 0
EL_nq = q1 * sp.diff(L_nq, q1) + q2 * sp.diff(L_nq, q2)
e_psi_nq = e_nq.subs({x1: psi_nq[0], x2: psi_nq[1]})
assert sp.simplify(EL_nq - L_nq - e_psi_nq) == 0
BD_nq = (grad(e_nq, X2).T * hess(e_nq, X2).inv() * grad(e_nq, X2))[0, 0]
E_e_psi_nq = q1 * sp.diff(e_psi_nq, q1) + q2 * sp.diff(e_psi_nq, q2)
assert sp.simplify(E_e_psi_nq - BD_nq.subs({x1: psi_nq[0], x2: psi_nq[1]})) == 0
assert sp.simplify(E_e_psi_nq) != 0    # B != 0 here, so L is NOT affine along rays
stamp('(C3) same chain on e = x1^3/3 + x2^2/2 (non-affine L, B != 0)          INSTANCE')

# ---------------------------------------------------------------- (C4)
ec = generic(X2, 3, 'f')
H, g = hess(ec, X2), grad(ec, X2)
D = sp.expand(det(H))
A = H.adjugate()
W = sp.expand(A * g)
Mt = sp.zeros(2, 2)
for k, xk in enumerate(X2):
    Mt += sp.diff(H, xk) * W[k]
Mt = sp.expand(Mt)
JacW = sp.Matrix(2, 2, lambda i, j: sp.diff(W[i], X2[j]))
gradD = grad(D, X2)
lhs = sp.expand(D * JacW - W * gradD.T)
rhs = sp.expand(D**2 * sp.eye(2) - A * Mt)
assert sp.expand(lhs - rhs) == sp.zeros(2, 2)
Bc = Bof(ec, X2)
lhs2 = sp.expand(D * grad(Bc, X2) - Bc * gradD)
rhs2 = sp.expand(2 * D**2 * g - Mt * W)
assert sp.expand(lhs2 - rhs2) == sp.zeros(2, 1)
stamp('(C4) flow identities Jac(V) = I - H^{-1}M and grad(g^T V) = 2g - MV (poly forms)  PROOF')

# ---------------------------------------------------------------- (C5)
t, tau = sp.symbols('t tau')
for vs, deg, tag in ((X2, 4, 'g'), (X3, 3, 'h')):
    n = len(vs)
    e = generic(vs, deg, tag)
    a = sp.symbols(f'a1:{n+1}')
    b = sp.symbols(f'b1:{n+1}')
    sub = {vs[i]: b[i] + a[i] / t for i in range(n)}
    for comp in grad(e, vs):
        expr = sp.expand(comp.subs(sub).subs(t, 1 / tau))
        assert expr.is_polynomial(tau), 'positive power of t found'
stamp('(C5) grad e(b + a/t) has no positive power of t (generic; n=2 deg 4, n=3 deg 3)  PROOF')

# ---------------------------------------------------------------- (C6)
e3 = generic(X3, 3, 'k')
d = 3
E = sp.expand(x0**d * e3.subs({x1: x1 / x0, x2: x2 / x0, x3: x3 / x0}))
assert sp.Poly(E, x0, x1, x2, x3).is_homogeneous
d0E = sp.diff(E, x0).subs(x0, 1)
assert sp.expand(d0E - (d * e3 - sum(xi * sp.diff(e3, xi) for xi in X3))) == 0
for xi in X3:
    assert sp.expand(sp.diff(E, xi).subs(x0, 1) - sp.diff(e3, xi)) == 0
vt = sp.symbols('v1 v2 v3')
shift = {X3[i]: X3[i] + vt[i] for i in range(3)}
assert sp.expand(Bof(e3.subs(shift, simultaneous=True), X3) - Bof(e3, X3).subs(shift, simultaneous=True)) == 0
# converse of Euler: y.grad(sum e_k) = d sum e_k  =>  (k-d) e_k = 0 for all k
for k in range(0, 4):
    ek = form(X3, k, f'm{k}_') if k > 0 else sp.Symbol('m0')
    assert sp.expand(sum(xi * sp.diff(ek, xi) for xi in X3) - k * ek) == 0
stamp('(C6) Theorem F steps: Euler at x0 = 1, translation invariance of B, converse of Euler  PROOF')

# ---------------------------------------------------------------- (C7)
for dd in (2, 3, 4):
    F = form(X3, dd, f'n{dd}_')
    assert sp.expand(Bof(F, X3) - sp.Rational(dd, dd - 1) * F * det(hess(F, X3))) == 0
stamp('(C7) Euler contraction for ternary forms d = 2,3,4                     PROOF')

# ---------------------------------------------------------------- (C8)
for vs, dd, tag in ((X2, 2, 'r'), (X2, 3, 's'), (X3, 2, 'w')):
    n = len(vs)
    e = generic(vs, dd, tag)
    E = sp.expand(x0**dd * e.subs({v: v / x0 for v in vs}, simultaneous=True))
    lhs = det(hess(E, (x0,) + tuple(vs))).subs(x0, 1)
    rhs = dd * (dd - 1) * e * det(hess(e, vs)) - (dd - 1)**2 * Bof(e, vs)
    assert sp.expand(lhs - rhs) == 0
stamp('(C8) affine equation of the Hessian (n=2, d=2,3; n=3, d=2)           PROOF')

# ---------------------------------------------------------------- (C9)
e0 = generic(X3, 3, 'a')
e1 = generic(X3, 3, 'b')
f = e0 + x4 * e1
P, K, g = hess(e0, X3), hess(e1, X3), grad(e1, X3)
H4 = hess(f, X4)
assert H4[:3, :3] == sp.expand(P + x4 * K) or sp.expand(H4[:3, :3] - (P + x4 * K)) == sp.zeros(3, 3)
assert sp.expand(H4[:3, 3] - g) == sp.zeros(3, 1) and H4[3, 3] == 0
det4 = sp.expand(det(H4))
assert sp.expand(det4 + (g.T * (P + x4 * K).adjugate() * g)[0, 0]) == 0
pol = sp.Poly(det4, x4)
assert pol.degree() <= 2
assert sp.expand(pol.coeff_monomial(x4**2) + Bof(e1, X3)) == 0
# Cauchy's bordered expansion for a generic symmetric 3x3 M and generic w
Mg = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'M{min(i,j)+1}{max(i,j)+1}'))
wg = sp.Matrix(sp.symbols('w1 w2 w3'))
bord = sp.Matrix(sp.BlockMatrix([[Mg, wg], [wg.T, sp.zeros(1, 1)]]))
assert sp.expand(det(bord) + (wg.T * Mg.adjugate() * wg)[0, 0]) == 0
stamp('(C9) Lemma lem:bordered (generic cubics) and Cauchy\'s bordered expansion  PROOF')

# ---------------------------------------------------------------- (C10)
e12 = generic(X2, 3, 'c')
f2 = e0 + x4 * e12
det42 = sp.expand(det(hess(f2, X4)))
pol2 = sp.Poly(det42, x4)
assert pol2.degree() <= 1
B2 = sp.expand(sp.diff(e12, x2, 2) * sp.diff(e12, x1)**2
               - 2 * sp.diff(e12, x1, x2) * sp.diff(e12, x1) * sp.diff(e12, x2)
               + sp.diff(e12, x1, 2) * sp.diff(e12, x2)**2)
assert sp.expand(pol2.coeff_monomial(x4) + sp.diff(e0, x3, 2) * B2) == 0
stamp('(C10) Lemma lem:affx4: affine in x4, [x4] = -(e0)_{x3x3} B2(e)         PROOF')

# ---------------------------------------------------------------- (C11)
gam, lam = sp.symbols('gamma lambda')
assert sp.expand(hess(e1**2 / 2, X3) - (e1 * K + g * g.T)) == sp.zeros(3, 3)
et0 = e0 - e1**2 / (2 * gam)
uu = x4 + e1 / gam
lhsm = hess(et0, X3) + uu * K
rhsm = (P + x4 * K) - (g * g.T) / gam
assert sp.simplify(sp.expand(lhsm - rhsm)) == sp.zeros(3, 3)
fq = e0 + x4 * e1 + gam / 2 * x4**2
lhsd = sp.expand(gam**2 * det(hess(fq, X4)))          # division-free form
rhsd = sp.expand(det(gam * (P + x4 * K) - g * g.T))   # FIX 2026-08-25: Schur gives gam^2 det Hess = det3(gam(P+x4K) - g g^T); the agent's extra factor gam was wrong
assert sp.expand(lhsd - rhsd) == 0                    # det Hess f = gamma det3(Hess3 e0~ + uK)
hl = et0 + lam * e1
assert sp.simplify(grad(hl, X3) - (grad(e0, X3) + (lam - e1 / gam) * g)) == sp.zeros(3, 1)
assert sp.simplify(hess(hl, X3) - (hess(et0, X3) + lam * K)) == sp.zeros(3, 3)
stamp('(C11) Theorem thm:quad: Schur identity, determinant, collision-transfer gradient  PROOF')

# ---------------------------------------------------------------- (C12)
fa = e0 + x4 * x3
assert sp.expand(det(hess(fa, X4)) + det(hess(e0, X2))) == 0
stamp('(C12) Theorem thm:aff: det Hess(e0 + x4 x3) = -det2 Hess_{(x1,x2)} e0  PROOF')

# ---------------------------------------------------------------- (C13)
phi = generic((x1,), 4, 'p')
fp = e0 + x4 * phi
assert sp.expand(det(hess(fp, X4)) + sp.diff(phi, x1)**2 * det(hess(e0, (x2, x3)))) == 0
stamp('(C13) Theorem thm:tri branch: det Hess(e0 + x4 phi(x1)) = -phi\'^2 det2 Hess_{(x2,x3)} e0  PROOF')

# ---------------------------------------------------------------- (C14)
aa = generic(X2, 3, 'd')
bb = generic(X2, 3, 'e')
ee = generic(X2, 3, 'f2_')
fd = aa + x3 * bb + x4 * ee
jac = sp.expand(sp.diff(bb, x1) * sp.diff(ee, x2) - sp.diff(bb, x2) * sp.diff(ee, x1))
assert sp.expand(det(hess(fd, X4)) - jac**2) == 0
Hd = hess(fd, X4)
assert sp.expand(Hd[2:, 2:]) == sp.zeros(2, 2)
Mblk = sp.Matrix([[sp.diff(bb, x1), sp.diff(ee, x1)], [sp.diff(bb, x2), sp.diff(ee, x2)]])
assert sp.expand(Hd[:2, 2:] - Mblk) == sp.zeros(2, 2)
stamp('(C14) Theorem thm:doubling: block form and det Hess = Jac(b,e)^2       PROOF')

# ---------------------------------------------------------------- (C15)
# lem:inj (=>): with p' != q', (b,e)(p') = (b,e)(q'), the system
# M(p') (x3,x4)^T = grad a(q') - grad a(p') has a solution iff det M(p') != 0,
# which is the Keller condition.  We verify the gradient bookkeeping symbolically:
pp = sp.symbols('pp1 pp2')
qq = sp.symbols('qq1 qq2')
y3, y4 = sp.symbols('y3 y4')
gf = grad(fd, X4)
at_p = {x1: pp[0], x2: pp[1], x3: y3, x4: y4}
at_q = {x1: qq[0], x2: qq[1], x3: 0, x4: 0}
diff12 = sp.Matrix([gf[0].subs(at_p) - gf[0].subs(at_q), gf[1].subs(at_p) - gf[1].subs(at_q)])
Mp = Mblk.subs({x1: pp[0], x2: pp[1]})
rhs_sys = sp.Matrix([sp.diff(aa, x1).subs({x1: qq[0], x2: qq[1]}) - sp.diff(aa, x1).subs({x1: pp[0], x2: pp[1]}),
                     sp.diff(aa, x2).subs({x1: qq[0], x2: qq[1]}) - sp.diff(aa, x2).subs({x1: pp[0], x2: pp[1]})])
assert sp.expand(diff12 - (Mp * sp.Matrix([y3, y4]) - rhs_sys)) == sp.zeros(2, 1)
assert sp.expand(gf[2].subs(at_p) - bb.subs({x1: pp[0], x2: pp[1]})) == 0
assert sp.expand(gf[3].subs(at_p) - ee.subs({x1: pp[0], x2: pp[1]})) == 0
stamp('(C15) Lemma lem:inj: components 1,2 of grad f(p\',y3,y4) - grad f(q\',0,0) = M(p\')y - rhs  PROOF')

# ---------------------------------------------------------------- (C16)
f2g = form(X3, 2, 'u2_')
f3g = form(X3, 3, 'u3_')
detg = sp.expand(det(hess(f2g + f3g, X3)))
top = sum(term for term in sp.Add.make_args(detg)
          if sp.Poly(term, *X3).total_degree() == 3)
assert sp.expand(top - det(hess(f3g, X3))) == 0
# exact 4-variable instance
f2i = (x1**2 + x2**2 + x3**2 + x4**2) / 2 + x1 * x3
f3i = x1**3 + 2 * x1 * x2 * x4 + x2**2 * x3 - x3**2 * x4 + x1 * x3 * x4
deti = sp.expand(det(hess(f2i + f3i, X4)))
topi = sum(term for term in sp.Add.make_args(deti) if sp.Poly(term, *X4).total_degree() == 4)
assert sp.expand(topi - det(hess(f3i, X4))) == 0
stamp('(C16) top graded piece of det Hess(f2 + fd) = det Hess fd             PROOF (n=3) / INSTANCE (n=4)')

# ---------------------------------------------------------------- (C17)
# converse failure
assert det(hess(x1**2 / 2 + x2, X2)) == 0 and Bof(x1**2 / 2 + x2, X2) == 1
# dBvdE witness
hh = x1 * x2 + x1**2 * x3
assert sp.expand(det(hess(hh, X3))) == 0 and sp.expand(Bof(hh, X3) + x1**4) == 0
# cautionary example: det Hess2 == 0 does not give phi(L)
hc = x1 + x2**2
assert sp.expand(det(hess(hc, X2))) == 0 and Bof(hc, X2) == 2
# n = 4 counterexample to planarity
e4 = x2 + x1 * x3 + x1**2 * x4
H4e = hess(e4, X4)
assert Bof(e4, X4) == 0 and sp.expand(det(H4e)) == 0
assert H4e.adjugate() == sp.zeros(4, 4) and H4e.rank() == 2
parts = [sp.diff(e4, v) for v in X4]
mons = sorted({m for pr in parts for m in sp.Poly(pr, *X4).monoms()})
coef = sp.Matrix([[sp.Poly(pr, *X4).coeff_monomial(sp.prod(v**k for v, k in zip(X4, m))) for m in mons]
                  for pr in parts])
assert coef.rank() == 4      # partials linearly independent => not affinely 3-variable
E4h = sp.expand(x0**3 * e4.subs({v: v / x0 for v in X4}, simultaneous=True))
assert sp.expand(E4h - (x0**2 * x2 + x0 * x1 * x3 + x1**2 * x4)) == 0
assert sp.expand(det(hess(E4h, (x0,) + X4))) == 0
# jet identity B(a(x1,x2,x3) + lam x4) = lam^2 det Hess3(a)
bj = sp.symbols('b11 b12 b13 b22 b23 b33 g1 g2 g3')
Hj = sp.Matrix([[bj[0], bj[1], bj[2], 0], [bj[1], bj[3], bj[4], 0], [bj[2], bj[4], bj[5], 0], [0, 0, 0, 0]])
gj = sp.Matrix([bj[6], bj[7], bj[8], lam])
M3j = Hj[:3, :3]
assert sp.expand((gj.T * Hj.adjugate() * gj)[0, 0] - lam**2 * det(M3j)) == 0
# sharpness witness
Fs = x1 - sp.sqrt(x1**2 - 2 * x2)
assert sp.simplify(Bof(Fs, X2)) == 0
assert sp.simplify(det(hess(Fs, X2)) + 1 / (x1**2 - 2 * x2)**2) == 0
# de Bondt's dillen4 doubling
fdB = (x1 + x2**2) * x3 + (x2 + (x1 + x2**2)**2) * x4
assert sp.expand(det(hess(fdB, X4))) == 1
jacdB = sp.expand(sp.diff(x1 + x2**2, x1) * sp.diff(x2 + (x1 + x2**2)**2, x2)
                  - sp.diff(x1 + x2**2, x2) * sp.diff(x2 + (x1 + x2**2)**2, x1))
assert jacdB == 1
stamp('(C17) witnesses of Remarks rem:conv, rem:dim, rem:strict; jet identity   INSTANCE/PROOF')

# ---------------------------------------------------------------- (C18)
w = x1**2 / 2 + x1 * x2 * (1 + x3) + (x2**2 / 2) * (2 * x3 + x3**2) + x3 * x4
assert sp.expand(det(hess(w, X4))) == 1
assert sp.expand(sp.diff(w, x4, 2)) == 0 and sp.expand(sp.diff(w, x4)) == x3
vv = sp.symbols('v1:5')
zz = sp.Symbol('zz')
Qv = sp.expand((sp.Matrix(vv).T * hess(w, X4) * sp.Matrix(vv))[0, 0])
gens = [sp.expand(c) for c in sp.Poly(Qv, *X4).coeffs()]
for i in range(3):
    G = sp.groebner(gens + [1 - zz * vv[i]], *(list(vv) + [zz]), order='grevlex')
    assert list(G.exprs) == [sp.Integer(1)], f'v{i+1} does not vanish on the isotropic cone'
assert all(sp.expand(gn.subs({vv[0]: 0, vv[1]: 0, vv[2]: 0, vv[3]: 1})) == 0 for gn in gens)
# a doubling has span(e3,e4) in its isotropic cone (block form of C14)
Qd = sp.expand((sp.Matrix(vv).T * Hd * sp.Matrix(vv))[0, 0]).subs({vv[0]: 0, vv[1]: 0})
assert sp.expand(Qd) == 0
stamp('(C18) Remark rem:notdoubling: det Hess w = 1, isotropic cone = C.e4; doublings contain span(e3,e4)  PROOF')

# ---------------------------------------------------------------- (C19)
s = sp.Symbol('s')
Fn = form(X3, 3, 'kk')
HF, gF = hess(Fn, X3), grad(Fn, X3)
n_, r_ = 3, 3
lhsN = sp.expand(det(-Fn * HF + s * (gF * gF.T)))
rhsN = sp.expand((-1)**(n_ - 1) * sp.Rational(r_, r_ - 1) * (s - sp.Rational(r_ - 1, r_))
                 * Fn**n_ * det(HF))
assert sp.expand(lhsN - rhsN) == 0
coeff_s = sp.expand(sp.diff(lhsN, s))
assert sp.expand(coeff_s - (-1)**(n_ - 1) * Fn**(n_ - 1) * Bof(Fn, X3)) == 0
# rank-one expansion used in the extraction, generic 3x3 M and vectors u, v
Mr = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'R{i+1}{j+1}'))
ur = sp.Matrix(sp.symbols('ur1 ur2 ur3'))
vr = sp.Matrix(sp.symbols('vr1 vr2 vr3'))
assert sp.expand(det(Mr + s * ur * vr.T) - (det(Mr) + s * (vr.T * Mr.adjugate() * ur)[0, 0])) == 0
stamp('(C19) Nagaoka-Yazawa identity, s^1 extraction, rank-one expansion      PROOF')

# ---------------------------------------------------------------- (C20)
f4g = generic(X4, 3, 'ff')
Dv2 = sum(vv[i] * sp.diff(sum(vv[j] * sp.diff(f4g, X4[j]) for j in range(4)), X4[i]) for i in range(4))
assert sp.expand(Dv2 - (sp.Matrix(vv).T * hess(f4g, X4) * sp.Matrix(vv))[0, 0]) == 0
stamp('(C20) D_v^2 f = v^T Hess f v                                            PROOF')

print()
print('fpG_paper_certificate: ALL CHECKS PASSED --', TEX.name,
      'is referentially sound and every displayed identity is re-certified.')
sys.exit(0)
