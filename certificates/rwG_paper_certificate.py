# rwG_paper_certificate.py -- certificate for ../theorem_G_paper.tex, revision
# of 2026-08-25 (agent key: rwG = rigorous write-up of Theorem G).
#
# Fail-closed (assert-based; nonzero exit on any failure); exact/symbolic
# arithmetic only (sympy 1.14).  Run:  py -u rwG_paper_certificate.py
#
# PART 1 -- structural "compile check" of the manuscript.  No TeX toolchain is
#   installed on this machine, so the following are checked textually:
#   (S1) every \ref / \eqref has a \label; labels unique; no label only in a
#        comment;
#   (S2) every \cite key has a \bibitem and no orphan \bibitem;
#   (S3) \begin/\end properly nested; braces balanced; $ balanced in every
#        paragraph; \[ \] balanced; \left / \right balanced in every display;
#   (S4) exactly one document skeleton; custom macros defined once and used;
#   (S5) UNDEFINED-MACRO SCAN: every control sequence \name used in the source
#        is either a custom macro or a known LaTeX / amsmath / amssymb /
#        amsthm / hyperref / geometry command (whitelist below);
#   (S6) the abstract says plainly that HC_4 remains open, and the headline
#        corollary (localization) and the no-pivot problem are present;
#   (S7) the three proofs of Theorem G are present with the labels the
#        appendix refers to.
#
# PART 2 -- exact re-derivation of the SECOND PROOF of Theorem G (the
#   derivation-theoretic proof), independently of advG/atkG:
#   (I1)  eq:dV in polynomial form.  With W := adj(H) g = D V and
#         Mt_j := (d_j H) W  (so that Mt = D * (d_j H) V):
#            D * d_j W - W * d_j D = D^2 e_j - adj(H) (d_j H) W .
#         Checked for FULLY GENERIC 3-jets in n = 2 and n = 3 (value and first
#         derivative at x = 0 of a generic cubic); since both sides at a point
#         depend only on the 3-jet of e there, this is a PROOF for all e in
#         those dimensions.  Random exact 3-jet for n = 4.
#   (I2)  eq:MV in polynomial form:  D grad(B) - B grad(D) = 2 D^2 g - M~ W
#         with M~ := sum_l (d_l H) W_l ;  same generic-jet proof.
#   (I3)  On an exact rational instance (n = 2, 3) the field identities of the
#         proof: delta g = g;  delta V = V - H^{-1} M V;
#         grad(g^T V) = 2 g - M V;  and  delta(x + V) = H^{-1} grad(B/D),
#         so that B == 0 forces delta b = 0 with b = x + V.
#   (I4)  Step 2 of the proof (eigen-decomposition) as a formal identity: for
#         a generic polynomial P and the derivation delta = -sum a_j d/da_j on
#         C(a, b) (so delta b = 0, delta a = -a),
#            P(b - a) = sum_k c_k,   c_k = sum_{|al|=k} (-1)^k (d^al P)(b) a^al / al!,
#         and  delta c_k = -k c_k.   Plus invertibility of the Vandermonde
#         matrix ((-k)^m) for d <= 10 (Step 3).
#   (I5)  Reading of Remark rem:flow: for e = x2/x1 (rational, B == 0,
#         det Hess = -x1^{-4}) one has V = -x, b = 0, a = -x, delta g = g,
#         and g(b - a tau) = tau^{-1} g(x): eigenvalue +1; for a generic
#         POLYNOMIAL e (n = 2, deg 5) g(b - a tau) is a polynomial in tau of
#         degree <= d - 1 (eigenvalues 0, -1, ..., -(d-1) only).
#   (I6)  The sharpness witness of Remark rem:conv: F = x1 - sqrt(x1^2 - 2 x2)
#         has B(F) == 0 and det Hess F = -(x1^2 - 2 x2)^{-2}.
#   (I7)  The doubling formula quoted in Section 2 in the form used:
#         det Hess <y, F(x)> = (-1)^n (Jac F)^2 for a fully generic quadratic
#         F in n = 2 (proof for that shape) and an exact instance, n = 3.

import itertools
import random
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

# (S1) labels / refs
labels = re.findall(r'\\label\{([^}]*)\}', stripped)
refs = re.findall(r'\\(?:ref|eqref)\{([^}]*)\}', stripped)
dup = {l for l in labels if labels.count(l) > 1}
assert not dup, f'duplicate labels: {dup}'
undefined = {r for r in refs if r not in labels}
assert not undefined, f'undefined refs: {undefined}'
raw_labels = re.findall(r'\\label\{([^}]*)\}', src)
assert set(raw_labels) == set(labels), 'a label appears only inside a comment'
print(f'  (S1) {len(set(labels))} labels, {len(refs)} refs, all defined, none duplicated')

# (S2) cites / bibitems
cites = [k.strip() for grp in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]*)\}', stripped)
         for k in grp.split(',')]
bibs = re.findall(r'\\bibitem\{([^}]*)\}', stripped)
assert len(bibs) == len(set(bibs)), 'duplicate bibitem'
missing = {c for c in cites if c not in bibs}
assert not missing, f'cite keys without bibitem: {missing}'
orphans = {b for b in bibs if b not in cites}
assert not orphans, f'uncited bibitems: {orphans}'
print(f'  (S2) {len(set(cites))} cite keys <-> {len(bibs)} bibitems, exact match')

# (S3) nesting, braces, $, \[ \], \left \right
stack = []
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}', stripped):
    kind, env = m.group(1), m.group(2)
    if kind == 'begin':
        stack.append(env)
    else:
        assert stack and stack[-1] == env, f'\\end{{{env}}} mismatches {stack[-1:]}'
        stack.pop()
assert not stack, f'unclosed environments: {stack}'
body = stripped.replace(r'\{', '').replace(r'\}', '')
assert body.count('{') == body.count('}'), 'brace imbalance'
for i, para in enumerate(re.split(r'\n\s*\n', stripped)):
    assert len(re.findall(r'(?<!\\)\$', para)) % 2 == 0, f'odd $ count in paragraph {i}'
assert stripped.count('\\[') == stripped.count('\\]'), r'\[ \] imbalance'
# \left / \right inside each display or inline math chunk
math_chunks = re.findall(r'\\\[(.*?)\\\]', stripped, re.S) + \
    re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}', stripped, re.S) + \
    re.findall(r'(?<!\\)\$(.+?)(?<!\\)\$', stripped, re.S)
for ch in math_chunks:
    assert len(re.findall(r'\\left\b', ch)) == len(re.findall(r'\\right\b', ch)), \
        f'\\left/\\right imbalance in: {ch[:60]!r}'
print('  (S3) environments nested, braces / $ / \\[ \\] / \\left \\right balanced')

# (S4) skeleton and macros
assert len(re.findall(r'\\documentclass', stripped)) == 1
assert len(re.findall(r'\\begin\{document\}', stripped)) == 1
assert len(re.findall(r'\\end\{document\}', stripped)) == 1
newcmds = re.findall(r'\\newcommand\{\\([A-Za-z]+)\}', stripped)
assert len(newcmds) == len(set(newcmds)), 'macro defined twice'
used = re.findall(r'\\([A-Za-z]+)', stripped)
for mac in newcmds:
    assert used.count(mac) >= 2, f'macro \\{mac} defined but never used'
print(f'  (S4) one document skeleton; {len(newcmds)} custom macros, each defined once and used')

# (S5) undefined-macro scan
WHITELIST = set('''
documentclass usepackage theoremstyle newtheorem newcommand title author thanks
date begin end maketitle section subsection subsubsection appendix label ref
eqref cite emph textbf texttt textit textsc item bibitem thebibliography
mathbb mathcal mathrm mathit mathbf mathsf mathtt operatorname frac tfrac dfrac
sqrt nabla top det langle rangle cdot cdots ldots dots dotsc dotsb vdots ddots
ne neq le leq ge geq ll gg in notin ni subseteq subset supseteq setminus times
to longmapsto mapsto longrightarrow rightarrow leftarrow Leftarrow Rightarrow
Leftrightarrow leftrightarrow partial sum prod bigcup bigcap cup cap circ equiv
not infty bigl bigr Bigl Bigr biggl biggr left right quad qquad text hbox mbox
widetilde widehat hat tilde bar overline underline vec delta gamma lambda mu
sigma tau alpha beta varphi phi psi chi kappa varepsilon epsilon eta theta
zeta rho xi omega nu pi iota upsilon vartheta varpi varrho varsigma
Delta Omega Lambda Phi Psi Sigma Gamma Theta Xi Pi Upsilon ell
pmatrix smallmatrix vmatrix bmatrix matrix cases ast min max deg dim ker
emptyset iff colon mid star qedhere proof rm bf it sc sf tt small large Large
footnotesize scriptsize normalsize displaystyle textstyle scriptstyle
enumerate itemize description footnote vspace hspace noindent centering
allowbreak hfill newline linebreak pagebreak smallskip medskip bigskip par
phantom hphantom vphantom strut mathstrut lvert rvert lVert rVert lfloor
rfloor lceil rceil boldsymbol forall exists neg lnot land lor vee wedge otimes
oplus cong simeq sim approx propto perp parallel dagger prime nolimits limits
big Big bigg Bigg overset underset stackrel xrightarrow substack tag notag
nonumber intertext mathfrak Bbbk varnothing nmid url href hypersetup l o aa ae
ss S P textsuperscript textsubscript newpage clearpage abstract document
equation align aligned gather split array tabular center flushleft flushright
verbatim quote quotation ensuremath relax kern hskip vskip
'''.split())
unknown = sorted({u for u in set(used) if u not in WHITELIST and u not in newcmds})
assert not unknown, f'control sequences not defined and not standard: {unknown}'
print(f'  (S5) undefined-macro scan: all {len(set(used))} distinct control sequences known')

# (S6) abstract and headline content
abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', stripped, re.S).group(1)
assert 'remains open' in abstract and 'HC' in abstract, 'abstract must state HC_4 open'
assert 'cor:headline' in labels and 'prob:nopivot' in labels and 'thm:tri' in labels
assert 'thm:G' in labels and 'cor:planar' in labels and 'thm:F' in labels
print('  (S6) abstract states that HC_4 remains open; headline corollary, trichotomy,')
print('       planarity corollary and the no-pivot problem are all present')

# (S7) three proofs of Theorem G and the labels used by the appendix
for key in ('First proof of Theorem~\\ref{thm:G}', 'Second proof of Theorem~\\ref{thm:G}',
            'Third proof of Theorem~\\ref{thm:G}'):
    assert key in stripped, f'missing: {key}'
for lab in ('eq:dV', 'eq:dg', 'eq:MV', 'eq:dVV', 'eq:legendre', 'eq:eulerode',
            'eq:section', 'rem:flow', 'rem:cont', 'rem:skeleton', 'rem:hyp'):
    assert lab in labels, f'label {lab} missing'
print('  (S7) three proofs of Theorem G present; proof labels present')
stamp('PART 1 PASSED')

# ============================================================== PART 2
stamp('PART 2: exact re-derivation of the second proof of Theorem G')


def hess(e, vs):
    return sp.Matrix([[sp.diff(e, a, b) for b in vs] for a in vs])


def grad(e, vs):
    return sp.Matrix([sp.diff(e, a) for a in vs])


def generic_cubic(vs, tag, randomize=False, rnd=None):
    """Generic (or random exact) polynomial of degree <= 3 in vs; the 3-jet
    at 0 is then fully general."""
    mons = [sp.Integer(1)]
    for d in range(1, 4):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(vs, d)},
                       key=str)
    if randomize:
        cs = [sp.Integer(rnd.randint(-7, 7)) for _ in mons]
    else:
        cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sum(c * m for c, m in zip(cs, mons))


def jet_identities(n, e, vs):
    """(I1),(I2) at x = 0: value and first derivatives of every quantity are
    computed by differentiating UNEXPANDED expressions and substituting
    x = 0, then expanding polynomials in the jet coefficients only."""
    zero = {v: 0 for v in vs}
    H = hess(e, vs)
    g = grad(e, vs)
    A = H.adjugate()
    D = H.det(method='berkowitz')
    W = A * g
    B = (g.T * A * g)[0, 0]
    H0 = H.subs(zero)
    A0 = A.subs(zero)
    D0 = sp.expand(D.subs(zero))
    W0 = sp.expand(W.subs(zero))
    g0 = g.subs(zero)
    B0 = sp.expand(B.subs(zero))
    gradD0 = [sp.expand(sp.diff(D, v).subs(zero)) for v in vs]
    gradB0 = [sp.expand(sp.diff(B, v).subs(zero)) for v in vs]
    JW0 = [[sp.expand(sp.diff(W[i], vs[j]).subs(zero)) for j in range(n)] for i in range(n)]
    dH0 = [sp.diff(H, v).subs(zero) for v in vs]          # constant matrices
    # (I1): D dW_j - W dD_j = D^2 e_j - A (d_j H) W
    for j in range(n):
        rhs = sp.expand(A0 * dH0[j] * W0)
        for i in range(n):
            lhs = sp.expand(D0 * JW0[i][j] - W0[i] * gradD0[j])
            r = sp.expand((D0**2 if i == j else 0) - rhs[i])
            assert sp.expand(lhs - r) == 0, f'(I1) fails n={n} i={i} j={j}'
    # (I2): D grad B - B grad D = 2 D^2 g - Mt W, Mt = sum_l (d_l H) W_l
    Mt = sp.zeros(n, n)
    for l in range(n):
        Mt += dH0[l] * W0[l]
    MtW = sp.expand(Mt * W0)
    for j in range(n):
        lhs = sp.expand(D0 * gradB0[j] - B0 * gradD0[j])
        rhs = sp.expand(2 * D0**2 * g0[j] - MtW[j])
        assert sp.expand(lhs - rhs) == 0, f'(I2) fails n={n} j={j}'


for n in (2, 3):
    vs = sp.symbols(f'x1:{n+1}')
    e = generic_cubic(vs, f'c{n}_')
    jet_identities(n, e, vs)
    stamp(f'(I1),(I2) n={n}: fully generic 3-jet ({len(e.free_symbols) - n} symbols)  PROOF')
rnd = random.Random(20260825)
for n in (4,):
    vs = sp.symbols(f'x1:{n+1}')
    e = generic_cubic(vs, '', randomize=True, rnd=rnd)
    jet_identities(n, e, vs)
    stamp(f'(I1),(I2) n={n}: random exact 3-jet                       INSTANCE')

# (I3) field identities on exact rational instances
x1, x2, x3 = sp.symbols('x1 x2 x3')


def is_zero_vec(v):
    return all(sp.cancel(sp.together(c)) == 0 for c in v)


for vs, e in (((x1, x2), x1**3/3 + x2**2/2 + x1*x2 + 2*x1),
              ((x1, x2, x3), x1**3/3 + x2**2/2 + x3**2/2 + x1*x2*x3)):
    n = len(vs)
    H = hess(e, vs)
    g = grad(e, vs)
    D = sp.cancel(H.det())
    assert D != 0
    Hinv = H.adjugate() / D
    V = (Hinv * g).applyfunc(sp.cancel)
    B = sp.expand((g.T * H.adjugate() * g)[0, 0])

    def delta(F):
        return sp.cancel(sum(V[j] * sp.diff(F, vs[j]) for j in range(n)))

    M = sp.zeros(n, n)
    for l in range(n):
        M += sp.diff(H, vs[l]) * V[l]
    # delta g = g
    for i in range(n):
        assert sp.cancel(delta(g[i]) - g[i]) == 0
    # delta V = V - H^{-1} M V
    dV = sp.Matrix([delta(V[i]) for i in range(n)])
    assert is_zero_vec(dV - (V - Hinv * M * V))
    # grad(g^T V) = 2 g - M V
    gTV = sp.cancel((g.T * V)[0, 0])
    assert sp.cancel(gTV - B / D) == 0
    gr = sp.Matrix([sp.diff(gTV, v) for v in vs])
    assert is_zero_vec(gr - (2 * g - M * V))
    # delta(x + V) = H^{-1} grad(B/D)
    b = sp.Matrix(vs) + V
    db = sp.Matrix([delta(b[i]) for i in range(n)])
    assert is_zero_vec(db - Hinv * sp.Matrix([sp.diff(B / D, v) for v in vs]))
    assert not is_zero_vec(db)      # here B != 0, so b is not delta-constant
    stamp(f'(I3) n={n}: delta g = g, delta V = V - H^-1 M V, grad(g^T V) = 2g - MV,'
          f' delta(x+V) = H^-1 grad(B/D)  EXACT')

# (I4) eigen-decomposition (Step 2) and Vandermonde (Step 3)
for n, d in ((2, 4), (3, 3)):
    a = sp.symbols(f'a1:{n+1}')
    b = sp.symbols(f'b1:{n+1}')
    y = sp.symbols(f'y1:{n+1}')
    mons = [sp.Integer(1)]
    for k in range(1, d + 1):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(y, k)},
                       key=str)
    cs = sp.symbols(f'p{n}_0:{len(mons)}')
    P = sum(c * m for c, m in zip(cs, mons))
    subs_b = {y[i]: b[i] for i in range(n)}
    c_k = []
    for k in range(d + 1):
        term = 0
        for al in itertools.product(range(k + 1), repeat=n):
            if sum(al) != k:
                continue
            dP = P
            fact = 1
            for i in range(n):
                dP = sp.diff(dP, y[i], al[i])
                fact *= sp.factorial(al[i])
            term += (-1)**k * dP.subs(subs_b) * sp.prod([a[i]**al[i] for i in range(n)]) / fact
        c_k.append(sp.expand(term))
    total = sp.expand(sum(c_k))
    target = sp.expand(P.subs({y[i]: b[i] - a[i] for i in range(n)}))
    assert sp.expand(total - target) == 0, '(I4) Taylor expansion fails'

    def delta_free(F):     # delta = -sum a_j d/da_j  (delta b = 0, delta a = -a)
        return sp.expand(-sum(a[j] * sp.diff(F, a[j]) for j in range(n)))

    for k in range(d + 1):
        assert sp.expand(delta_free(c_k[k]) + k * c_k[k]) == 0, f'(I4) eigenvalue fails k={k}'
    stamp(f'(I4) n={n}, deg {d}, generic P: P(b-a) = sum c_k with delta c_k = -k c_k  PROOF')
for d in range(2, 11):
    Vand = sp.Matrix([[sp.Integer(-k)**m for k in range(d)] for m in range(d)])
    assert Vand.det() != 0
stamp('(I4) Vandermonde ((-k)^m), nodes 0,-1,..,-(d-1), invertible for d <= 10  PROOF')

# (I5) reading of Remark rem:flow
tau = sp.Symbol('tau')
vs = (x1, x2)
e_rat = x2 / x1
H = hess(e_rat, vs)
g = grad(e_rat, vs)
D = sp.cancel(H.det())
assert sp.cancel(D + x1**-4) == 0
assert sp.expand(sp.cancel((g.T * H.adjugate() * g)[0, 0])) == 0
V = sp.simplify(H.inv() * g)
assert sp.simplify(V + sp.Matrix(vs)) == sp.zeros(2, 1)          # V = -x
bvec = sp.Matrix(vs) + V
assert sp.simplify(bvec) == sp.zeros(2, 1)                        # b = 0
avec = V                                                          # a = -x
gsub = g.subs({x1: bvec[0] - avec[0] * tau, x2: bvec[1] - avec[1] * tau})
assert sp.simplify(gsub - g / tau) == sp.zeros(2, 1)              # g(b - a tau) = tau^-1 g(x)
for i in range(2):
    dg = sp.cancel(sum(V[j] * sp.diff(g[i], vs[j]) for j in range(2)))
    assert sp.cancel(dg - g[i]) == 0                              # delta g = g
stamp('(I5) e = x2/x1: B = 0, det Hess = -x1^-4, V = -x, b = 0, delta g = g,'
      ' g(b - a tau) = tau^-1 g(x)  EXACT')
y = sp.symbols('y1 y2')
mons = []
for k in range(1, 6):
    mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(y, k)}, key=str)
cs = sp.symbols(f'q0:{len(mons)}')
e_gen = sum(c * m for c, m in zip(cs, mons))
av = sp.symbols('A1 A2')
bv = sp.symbols('B1 B2')
for i in range(2):
    comp = sp.expand(sp.diff(e_gen, y[i]).subs({y[0]: bv[0] - av[0]*tau, y[1]: bv[1] - av[1]*tau}))
    pol = sp.Poly(comp, tau)
    assert pol.degree() <= 4
stamp('(I5) generic polynomial e (n=2, deg 5): g(b - a tau) has tau-degree <= d-1 = 4,'
      ' no tau^-1 term  PROOF')

# (I6) sharpness witness of Remark rem:conv
F = x1 - sp.sqrt(x1**2 - 2*x2)
HF = hess(F, vs)
gF = grad(F, vs)
assert sp.simplify(HF.det() + (x1**2 - 2*x2)**-2) == 0
assert sp.simplify((gF.T * HF.adjugate() * gF)[0, 0]) == 0
stamp('(I6) F = x1 - sqrt(x1^2 - 2x2): B(F) = 0 and det Hess F = -(x1^2 - 2x2)^-2  EXACT')

# (I7) doubling formula det Hess <y, F(x)> = (-1)^n (Jac F)^2
for n, generic in ((2, True), (3, False)):
    xs = sp.symbols(f'u1:{n+1}')
    ys = sp.symbols(f'v1:{n+1}')
    mons = [sp.Integer(1)]
    for k in range(1, 3):
        mons += sorted({sp.prod(c) for c in itertools.combinations_with_replacement(xs, k)},
                       key=str)
    Fc = []
    for i in range(n):
        if generic:
            cs = sp.symbols(f'f{i}_0:{len(mons)}')
        else:
            cs = [sp.Integer(rnd.randint(-5, 5)) for _ in mons]
        Fc.append(sum(c * m for c, m in zip(cs, mons)))
    phi = sum(ys[i] * Fc[i] for i in range(n))
    allv = list(xs) + list(ys)
    lhs = sp.expand(hess(phi, allv).det(method='berkowitz'))
    J = sp.Matrix([[sp.diff(Fc[i], xs[j]) for j in range(n)] for i in range(n)])
    rhs = sp.expand((-1)**n * J.det(method='berkowitz')**2)
    assert sp.expand(lhs - rhs) == 0
    stamp(f'(I7) n={n}: det Hess <y,F(x)> = (-1)^n (Jac F)^2  '
          + ('PROOF (generic quadratic F)' if generic else 'INSTANCE (exact quadratic F)'))

stamp('rwG_paper_certificate: ALL CHECKS PASSED')
print('theorem_G_paper.tex (revision 2026-08-25) is referentially sound, uses no')
print('undefined control sequences, and every identity of its second proof of')
print('Theorem G is re-derived exactly (generic 3-jets in n = 2, 3 => proof).')
sys.exit(0)
