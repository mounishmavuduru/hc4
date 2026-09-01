# _d5_param.py -- PARAMETRIC decision: treat the b3 coefficients (c's) as
# parameters (coefficient field), and compute the Groebner basis of the
# residual ideal J = <coeffs E0,E1,E2,E3> in the a-variables ONLY (the 21
# coefficients of the quintic a5).  The fiber runs showed J=(1) for generic
# numeric b3; if std(J) = (1) over the function field k(c), that PROVES: for
# every b3 outside a proper closed c-locus, NO quintic a5 satisfies the
# residual -> the branch is EMPTY off that locus.  The special c-locus (where
# leading coefficients of the GB vanish) is then a strictly lower-dimensional
# tail to finish separately.
#
# We run over k = F_p (fast, k = F_p(c)) and, if that says (1), over k = Q
# (k = Q(c), a genuine char-0 proof of generic-c inconsistency).
#
#   py _d5_param.py [prime|0]     0 => rationals Q(c)
# writes _lean/param/param_b{k}_{field}.sing  (fed to Singular via stdin).

import itertools, os, sys
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.getcwd(), '_lean', 'param')
os.makedirs(OUT, exist_ok=True)

FIELD = sys.argv[1] if len(sys.argv) > 1 else '32003'   # '0' => Q
y1, y2, y3, s = sp.symbols('y1 y2 y3 s')
Y = (y1, y2, y3)


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def grad(g, vs):
    return sp.Matrix([sp.diff(g, v) for v in vs])


def form(vs, deg, tag):
    mons = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(vs, deg)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sp.expand(sum(c * m for c, m in zip(cs, mons))), list(cs)


def coeffs_in_y(expr):
    e = sp.expand(expr)
    return [] if e == 0 else [sp.expand(c) for c in sp.Poly(e, *Y).coeffs()]


def E_components(a5, b3):
    A, B, q = hess(a5, Y), hess(b3, Y), grad(b3, Y)
    adjA, adjB = A.adjugate(), B.adjugate()
    L = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff((A + s * B).adjugate()[i, j], s).subs(s, 0)))
    E2 = y1 * sp.trace(A * adjB) - (q.T * adjB * q)[0, 0] - 2 * (L * q)[0] - adjA[0, 0]
    E1 = y1 * sp.trace(adjA * B) - (q.T * L * q)[0, 0] - 2 * (adjA * q)[0]
    E0 = y1 * A.det(method='berkowitz') - (q.T * adjA * q)[0, 0]
    E3 = y1 * B.det(method='berkowitz') - 2 * (adjB * q)[0] - L[0, 0]
    return [sp.expand(e) for e in (E3, E2, E1, E0)], A


def s_poly(p):
    return str(sp.expand(p)).replace('**', '^')


def emit(k, avars, cvars, Jpolys, Tpolys):
    fld = 'Q' if FIELD == '0' else FIELD
    par = ','.join(str(c) for c in cvars)
    av = ','.join(str(a) for a in avars)
    coeff = f'(0,{par})' if FIELD == '0' else f'({FIELD},{par})'
    path = os.path.join(OUT, f'param_b{k}_{fld}.sing')
    Jp = [p for p in Jpolys if p != 0]
    Tp = [p for p in Tpolys if p != 0]
    with open(path, 'w', newline='\n') as f:
        f.write(f'// PARAMETRIC branch {k}: k(c) = {coeff}, vars = a only\n')
        f.write(f'ring R = {coeff},({av}),dp;\n')
        f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jp) + ';\n')
        f.write('int t0 = timer;\n')
        f.write('ideal G = std(J);\n')
        f.write(f'"branch {k}, field {coeff}";\n')
        f.write('"  std time 1/100s:", timer - t0, "  GBsize:", size(G);\n')
        f.write('if (size(G)==1 && G[1]==1) {\n')
        f.write('  "  RESULT: J=(1) over k(c) -> NO a5 for generic b3 -> branch EMPTY off a proper c-locus (PROVED for generic c)";\n')
        f.write('} else {\n')
        f.write('  "  dim (over k(c)):", dim(G);\n')
        f.write('  ideal T = ' + ',\n    '.join(s_poly(p) for p in Tp) + ';\n')
        f.write('  ideal red = reduce(T, G);\n')
        f.write('  "  T reduced mod G, nonzero count:", size(simplify(red,2));\n')
        f.write('  "  (if 0: det A in J already; else need radical test)";\n')
        f.write('}\n')
        f.write('quit;\n')
    return path, len(Jp), len(Tp)


b3g, cs = form(Y, 3, 'c')
eqs4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
sols_b = sp.solve(eqs4, cs, dict=True)
a5, as_ = form(Y, 5, 'a')

for k, sb in enumerate(sols_b):
    b3k = sp.expand(b3g.subs(sb))
    if b3k == 0:
        print(f'branch {k}: b3==0 -> EMPTY by domain argument; skip.')
        continue
    Es, A = E_components(a5, b3k)
    Jsym = []
    for e in Es:
        Jsym += coeffs_in_y(e)
    Tsym = coeffs_in_y(A.det(method='berkowitz'))
    # clear denominators generator-wise (keep polynomial in a and c)
    Jp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jsym]
    Tp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Tsym]
    cvars = sorted({v for p in Jsym for v in p.free_symbols if str(v).startswith('c')}, key=str)
    path, nj, nt = emit(k, as_, cvars, Jp, Tp)
    print(f'branch {k}: params {[str(c) for c in cvars]}, {len(as_)} a-vars, {nj} J-gens, {nt} T-gens -> {os.path.basename(path)}')
