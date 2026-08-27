# refG_paper_integrity.py -- referee-pass certificate for ../theorem_G_paper.tex
#
# Agent key: refG (write-up referee pass, 2026-08-20).
#
# PART 1 (structural "compile check", no LaTeX toolchain available on this
# machine): fail-closed verification that the manuscript has
#   - no undefined \ref / \eqref (every referenced key carries a \label),
#   - no duplicate \label keys,
#   - every \cite key matched by a \bibitem (and no orphan \bibitem),
#   - properly nested \begin{env} ... \end{env} pairs,
#   - balanced braces outside comments,
#   - exactly one \documentclass / \begin{document} / \end{document},
#   - no \newcommand redefinition clashes among the custom macros used.
#
# PART 2 (mathematical spot-audit, exact symbolic arithmetic): independent
# re-derivation of the load-bearing NEW identities exactly as stated in the
# manuscript (all with FULLY GENERIC coefficients of the stated shape, so each
# expand()==0 is a proof for all polynomials of that shape):
#   (I)   Lemma [lem:bordered]: for f = e0(x') + x4*e1(x'),
#         det Hess f = -g^T adj(P + x4 K) g,  deg_{x4} <= 2, and
#         [x4^2] det Hess f = -B(e1).
#   (II)  Lemma [lem:affx4]: for e1 = e(x1,x2),  det Hess f is AFFINE in x4
#         and [x4] det Hess f = -(e0)_{x3x3} * B2(e).
#   (III) Theorem [thm:quad] Schur identity:
#         Hess3(e0 - e1^2/(2*gamma)) + (x4 + e1/gamma) K = (P + x4 K) - g g^T/gamma.
#   (IV)  Theorem [thm:doubling]: det Hess(a + x3 b + x4 e) = Jac(b,e)^2.
#   (V)   Theorem [thm:G] Steps 2-3 structural identities:
#         H adj(H) g = det(H) g  and  g^T adj(H) g = B.
#   (VI)  Attribution [sec:attribution]: the Nagaoka-Yazawa identity as
#         DISPLAYED in the manuscript, and that Lemma [lem:identE] is its
#         s^1-coefficient after cancelling (-1)^(n-1) F^(n-1)  (n = 3, d = 3).
#
# Fail-closed: any assertion failure exits nonzero.

import itertools
import re
import sys
from pathlib import Path

import sympy as sp

TEX = Path(__file__).resolve().parent.parent / 'theorem_G_paper.tex'
src = TEX.read_text(encoding='utf-8')

# ------------------------------------------------------------------ PART 1
print('PART 1: structural compile-check of', TEX.name)

# strip comments (a % not preceded by backslash kills the rest of the line)
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
unused_labels = {l for l in labels if l not in refs}
print(f'  labels: {len(set(labels))} (all unique), refs: {len(refs)} (all defined)')
print(f'  cites: {len(set(cites))} keys, bibitems: {len(bibs)} (exact match, no orphans)')
if unused_labels:
    print(f'  note: labels never referenced (harmless): {sorted(unused_labels)}')

# environment nesting
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

# brace balance (ignore escaped braces)
body = stripped.replace(r'\{', '').replace(r'\}', '')
bal = body.count('{') - body.count('}')
assert bal == 0, f'brace imbalance: {bal}'
print('  braces: balanced')

assert len(re.findall(r'\\documentclass', stripped)) == 1
assert len(re.findall(r'\\begin\{document\}', stripped)) == 1
assert len(re.findall(r'\\end\{document\}', stripped)) == 1
newcmds = re.findall(r'\\newcommand\{\\([A-Za-z]+)\}', stripped)
assert len(newcmds) == len(set(newcmds)), 'duplicate \\newcommand'
print(f'  document skeleton OK; {len(newcmds)} custom macros, no redefinition')

# ------------------------------------------------------------------ PART 2
print()
print('PART 2: exact symbolic re-audit of the load-bearing identities')

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
XP = (x1, x2, x3)


def generic(vs, deg, tag):
    mons = [sp.Integer(1)]
    for d in range(1, deg + 1):
        mons += sorted({sp.prod(c) for c in
                        itertools.combinations_with_replacement(vs, d)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons))


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def grad(e, vs):
    return sp.Matrix([sp.diff(e, a) for a in vs])


def Bof(e, vs):
    return sp.expand((grad(e, vs).T * hess(e, vs).adjugate() * grad(e, vs))[0, 0])


# (I) Lemma [lem:bordered], fully generic cubic e0, e1 in three variables
e0 = generic(XP, 3, 'a')
e1 = generic(XP, 3, 'b')
f = e0 + x4 * e1
P, K, g = hess(e0, XP), hess(e1, XP), grad(e1, XP)
H4 = hess(f, (x1, x2, x3, x4))
det4 = sp.expand(H4.det(method='berkowitz'))
bordered = sp.expand(-(g.T * (P + x4 * K).adjugate() * g)[0, 0])
assert sp.expand(det4 - bordered) == 0
pol = sp.Poly(det4, x4)
assert pol.degree() <= 2
assert sp.expand(pol.coeff_monomial(x4**2) + Bof(e1, XP)) == 0
print('  (I)  lem:bordered: det Hess f = -g^T adj(P+x4 K) g, deg_{x4} <= 2,')
print('       [x4^2] = -B(e1)   (fully generic cubics)   PASS')

# (II) Lemma [lem:affx4], e = generic cubic in (x1,x2)
e = generic((x1, x2), 3, 'c')
f2 = e0 + x4 * e
det42 = sp.expand(hess(f2, (x1, x2, x3, x4)).det(method='berkowitz'))
pol2 = sp.Poly(det42, x4)
assert pol2.degree() <= 1
B2 = sp.expand(sp.diff(e, x2, x2) * sp.diff(e, x1)**2
               - 2 * sp.diff(e, x1, x2) * sp.diff(e, x1) * sp.diff(e, x2)
               + sp.diff(e, x1, x1) * sp.diff(e, x2)**2)
lhs = pol2.coeff_monomial(x4)
assert sp.expand(lhs + sp.diff(e0, x3, x3) * B2) == 0
print('  (II) lem:affx4: det Hess f affine in x4, [x4] = -(e0)_{x3x3} B2(e)')
print('       (fully generic)   PASS')

# (III) Theorem [thm:quad] Schur identity, symbolic gamma
gam = sp.Symbol('gamma')
et0 = sp.expand(e0 - e1**2 / (2 * gam))
u = x4 + e1 / gam
lhsm = hess(et0, XP) + u * K
rhsm = (P + x4 * K) - (g * g.T) / gam
assert sp.simplify(sp.expand(lhsm - rhsm)) == sp.zeros(3, 3)
print('  (III) thm:quad Schur identity: Hess3(e0 - e1^2/2g) + uK = (P+x4K) - gg^T/g')
print('       (fully generic cubics, symbolic gamma)   PASS')

# (IV) Theorem [thm:doubling], fully generic cubics a, b, e in (x1,x2)
a2 = generic((x1, x2), 3, 'd')
b2 = generic((x1, x2), 3, 'e')
e2 = generic((x1, x2), 3, 'f')
fd = a2 + x3 * b2 + x4 * e2
detd = sp.expand(hess(fd, (x1, x2, x3, x4)).det(method='berkowitz'))
jac = sp.expand(sp.diff(b2, x1) * sp.diff(e2, x2) - sp.diff(b2, x2) * sp.diff(e2, x1))
assert sp.expand(detd - jac**2) == 0
print('  (IV) thm:doubling: det Hess(a + x3 b + x4 e) = Jac(b,e)^2 (fully generic)  PASS')

# (V) Theorem [thm:G] structural identities, n = 3 generic cubic
Hh, gh = hess(e0, XP), grad(e0, XP)
W = sp.expand(Hh.adjugate() * gh)
assert sp.expand(Hh * W - Hh.det(method='berkowitz') * gh) == sp.zeros(3, 1)
assert sp.expand((gh.T * W)[0, 0] - Bof(e0, XP)) == 0
print('  (V)  thm:G Steps 2-3: H adj(H) g = det(H) g  and  g^T adj(H) g = B   PASS')

# (VI) Nagaoka-Yazawa identity exactly as displayed in sec:attribution (n=3, d=3)
s = sp.Symbol('s')
n, d = 3, 3
mons = sorted({sp.prod(c) for c in
               itertools.combinations_with_replacement(XP, d)}, key=str)
F = sum(c * m for c, m in zip(sp.symbols(f'k0:{len(mons)}'), mons))
HF, gF = hess(F, XP), grad(F, XP)
lhsN = sp.expand((-F * HF + s * (gF * gF.T)).det(method='berkowitz'))
rhsN = sp.expand((-1)**(n - 1) * sp.Rational(d, d - 1) * (s - sp.Rational(d - 1, d))
                 * F**n * HF.det(method='berkowitz'))
assert sp.expand(lhsN - rhsN) == 0
# and the manuscript's extraction claim: s^1-coefficient = (-1)^(n-1) F^(n-1) B(F)
coeff_s = sp.expand(sp.diff(lhsN, s))
assert sp.expand(coeff_s - (-1)**(n - 1) * F**(n - 1) * Bof(F, XP)) == 0
# cancelling gives Lemma lem:identE
assert sp.expand(Bof(F, XP) - sp.Rational(d, d - 1) * F
                 * HF.det(method='berkowitz')) == 0
print('  (VI) NY identity as displayed + s^1-extraction + cancellation = lem:identE')
print('       (fully generic ternary cubic)   PASS')

print()
print('refG_paper_integrity: ALL CHECKS PASSED --', TEX.name,
      'is referentially sound and its new identities are re-certified.')
sys.exit(0)
