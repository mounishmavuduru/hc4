# Gradient-equivalence test for the known dimension-3/4 Jacobian counterexamples.
#
# FACT (linear reduction): for a polynomial Keller map F: C^n -> C^n, there exist
# affine L1, L2 with L1 o F o L2 a gradient map  iff  there exists W in GL_n(C)
# with W * J_F(x) symmetric for all x.   [If J_G(x) = S J_F(L2 x) T is symmetric,
# then W := (T^T)^{-1} S works at all points y = L2 x; conversely W J_F symmetric
# means W o F has symmetric Jacobian, i.e. the 1-form sum_i (WF)_i dx_i is closed,
# hence exact with a POLYNOMIAL potential g (Poincare lemma, degree-by-degree
# integration), so W o F = grad g with det Hess g = det(W) Jac F constant.]
#
# Therefore: if some invertible W symmetrizes J_{F4} or J_{F5}, then HC_4 is
# FALSE (g as above is a 4-variable counterexample).  If the symmetrizer space
# contains no invertible matrix, the known dim-4 counterexamples do NOT refute
# HC_4 by affine equivalence -- a negative certificate.
#
# Consistency check: for Alpoge's 3-dim F, an invertible symmetrizer would
# contradict de Bondt's theorem (HC_3 true). Expect none.
#
# The maps F4, F5 are PARSED from Gao's LaTeX source (arXiv:2608.00222) and the
# transcription is checksummed by verifying the constant Jacobian determinant
# stated in Theorems thm:F4 and thm:F5 (-44/9 and 160/29).

import re
import sympy as sp

x, y, z, t = sp.symbols('x y z t')
LOCALS = {'x': x, 'y': y, 'z': z, 't': t}
TEX = open(r'C:\Users\mouni\hc4\lit\gao\Jacobian_CE.tex', encoding='utf-8').read()


def align_block_after(marker):
    i = TEX.index(marker)
    j = TEX.index(r'\begin{align*}', i)
    k = TEX.index(r'\end{align*}', j)
    return TEX[j + len(r'\begin{align*}'):k]


def parse_components(block, prefix):
    # split on component headers  F_{d,k} ={}&
    parts = re.split(r'F_\{\d,\d\}\s*=\{\}&', block)
    assert parts[0].strip() in ('', None) or parts[0].strip() == '', parts[0]
    comps = []
    for chunk in parts[1:]:
        s = chunk
        s = re.sub(r'\\\\(\[\d+pt\])?', '', s)   # line breaks
        s = s.replace('&', '').replace('\n', '').replace(' ', '')
        s = re.sub(r'\\tfrac\{(-?\d+)\}\{(\d+)\}', r'(\1/\2)*', s)
        s = re.sub(r'\^\{(\d+)\}', r'**\1', s)
        s = re.sub(r'\^(\d)', r'**\1', s)
        s = re.sub(r'(?<![\*/\d])(\d+)(?=[xyzt])', r'\1*', s)  # 14x -> 14*x
        s = re.sub(r'(?<=[xyzt\)])(?=[xyzt])', '*', s)          # xy -> x*y
        s = re.sub(r'(\*\*\d+)(?=[xyzt])', r'\1*', s)           # x**6y -> x**6*y
        comps.append(sp.expand(sp.sympify(s, locals=LOCALS)))
    return comps


def symmetrizer_space(F, vars_):
    n = len(F)
    J = sp.Matrix([[sp.diff(F[i], v) for v in vars_] for i in range(n)])
    Jpolys = [[sp.Poly(J[i, j], *vars_) for j in range(n)] for i in range(n)]
    monoms = set()
    for i in range(n):
        for j in range(n):
            monoms.update(Jpolys[i][j].as_dict().keys())
    monoms = sorted(monoms)
    midx = {m: r for r, m in enumerate(monoms)}
    dicts = [[Jpolys[i][j].as_dict() for j in range(n)] for i in range(n)]
    # unknowns W[i][k] -> column index n*i + k
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            # sum_k W[i,k] J[k,j] - W[j,k] J[k,i]  == 0
            eqrows = {}
            for k in range(n):
                for m, c in dicts[k][j].items():
                    eqrows.setdefault(m, [sp.Integer(0)] * (n * n))[n * i + k] += c
                for m, c in dicts[k][i].items():
                    eqrows.setdefault(m, [sp.Integer(0)] * (n * n))[n * j + k] -= c
            rows.extend(eqrows.values())
    Msys = sp.Matrix(rows)
    ns = Msys.nullspace()
    return ns, J


def report(name, F, vars_, expected_jac):
    print(f'--- {name} ---')
    J = sp.Matrix([[sp.diff(F[i], v) for v in vars_] for i in range(len(F))])
    d = sp.factor(J.det(method='domain-ge'))
    print(f'det J = {d} (expected {expected_jac})')
    assert sp.simplify(d - expected_jac) == 0, 'CHECKSUM FAILED - transcription error'
    ns, _ = symmetrizer_space(F, vars_)
    n = len(F)
    print(f'symmetrizer space dimension: {len(ns)}')
    if not ns:
        print(f'RESULT[{name}]: only W = 0; NOT affinely equivalent to any gradient map.')
        return
    cs = sp.symbols(f'c0:{len(ns)}')
    Wgen = sp.zeros(n, n)
    for c, vec in zip(cs, ns):
        Wgen += c * sp.Matrix(n, n, list(vec))
    detW = sp.expand(Wgen.det())
    print('det of generic symmetrizer =', detW)
    if detW == 0:
        print(f'RESULT[{name}]: symmetrizer space nontrivial (dim {len(ns)}) but contains')
        print('NO invertible matrix; NOT affinely equivalent to any gradient map.')
    else:
        print(f'RESULT[{name}]: INVERTIBLE SYMMETRIZER EXISTS -> gradient-equivalent!')
        print('!!! If this is F4/F5 this refutes HC_4 - investigate immediately !!!')


# ---- Alpoge's 3-dim map (consistency check against de Bondt HC_3) ----
u = 1 + x*y
Falp = [u**3*z + y**2*u*(4 + 3*x*y),
        y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
        2*x - 3*x**2*y - x**3*z]
report('Alpoge F (dim 3)', Falp, [x, y, z], -2)

# ---- Gao F4 ----
F4 = parse_components(align_block_after(r'$F_4=(F_{4,1}'), 'F4')
assert len(F4) == 4
degs = [sp.Poly(f, x, y, z, t).total_degree() for f in F4]
print('F4 component degrees:', degs)
assert degs == [4, 11, 12, 21], degs
report('Gao F4 (dim 4, geometric degree 5)', F4, [x, y, z, t], sp.Rational(-44, 9))

# ---- Gao F5 ----
F5 = parse_components(align_block_after(r'$F_5=(F_{5,1}'), 'F5')
assert len(F5) == 4
degs5 = [sp.Poly(f, x, y, z, t).total_degree() for f in F5]
print('F5 component degrees:', degs5)
assert degs5 == [3, 12, 14, 16], degs5
report('Gao F5 (dim 4, geometric degree 10)', F5, [x, y, z, t], sp.Rational(160, 29))

print()
print('DONE.')
