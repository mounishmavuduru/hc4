# _d5_lean.py -- build a LEAN (low-degree) Singular decision for the degree-5
# rank-3 pivot-free branches, mod p.
#
# Difference from d5_web_decision.py: we do NOT solve/substitute E3 and do NOT
# clear c5 denominators.  The residual ideal is just the raw y-coefficients of
# E0,E1,E2,E3 -- each of low total degree in the unknowns.  This keeps the
# Groebner input small-degree (the earlier export's degree-28 generators were
# an artifact of substitution + denominator clearing).  It is also strictly
# MORE complete: it keeps the c5 = 0 locus, which the saturated export dropped.
#
# Decision: branch EMPTY  <=>  every y-coefficient of det Hess_3(a5) (ideal T)
# lies in sqrt(J).  We compute G = std(J) once; if G = (1) the residual is
# inconsistent and the branch is EMPTY outright; otherwise test T subset sqrt(J)
# by Rabinowitsch on the reduced basis G.
#
#   py _d5_lean.py            # writes /tmp-style .sing files into ./_lean/
# Output .sing files go next to this script under _lean/ ; a bash driver runs
# them in WSL Singular.

import itertools, os, sys
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.getcwd(), '_lean')
os.makedirs(OUT, exist_ok=True)

y1, y2, y3, s = sp.symbols('y1 y2 y3 s')
Y = (y1, y2, y3)
PRIMES = sys.argv[1:] or ['32003']


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


def numden_clear(exprs):
    """clear only rational-number / c-denominators generator-wise, keep degree low."""
    out = []
    for e in exprs:
        n, _ = sp.fraction(sp.together(e))
        n = sp.expand(n)
        if n != 0:
            out.append(n)
    return out


def s_poly(p):
    return str(sp.expand(p)).replace('**', '^')


def emit(k, unknowns, Jpolys, Tpolys, prime):
    ring = ','.join(str(u) for u in unknowns)
    path = os.path.join(OUT, f'lean_b{k}_p{prime}.sing')
    with open(path, 'w', newline='\n') as f:
        f.write(f'// LEAN degree-5 rank-3 branch {k}, mod {prime}\n')
        f.write(f'ring R = {prime},({ring}),dp;\n')
        f.write('short=0; option(redSB);\n')
        f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jpolys) + ';\n')
        f.write('int t0 = timer;\n')
        f.write('ideal G = slimgb(J);\n')
        f.write('"branch ' + str(k) + ' mod ' + str(prime) + '";\n')
        f.write('"  slimgb(J) time 1/100s:", timer - t0, "   GB size:", size(G);\n')
        f.write('if (size(G) == 1 && G[1] == 1) {\n')
        f.write('  "  J = (1): residual INCONSISTENT -> branch EMPTY outright";\n')
        f.write('} else {\n')
        f.write('  "  dim V(J):", dim(std(G));\n')
        f.write('  ideal T = ' + ',\n    '.join(s_poly(p) for p in Tpolys) + ';\n')
        f.write('  int bad = 0; int i;\n')
        f.write('  for (i = 1; i <= size(T); i++) {\n')
        f.write(f'    ring Rw = {prime},({ring},rabw),dp;\n')
        f.write('    ideal G = imap(R, G);\n')
        f.write('    ideal T = imap(R, T);\n')
        f.write('    ideal Q = G, T[i]*rabw - 1;\n')
        f.write('    ideal g = slimgb(Q);\n')
        f.write('    if (size(g) != 1 || g[1] != 1) { bad = bad + 1; "    coeff", i, "NOT in sqrt(J)"; }\n')
        f.write('    setring R; kill Rw;\n')
        f.write('  }\n')
        f.write('  "  size(T)=", size(T), "  bad=", bad, "  EMPTY(mod p)=", (bad==0);\n')
        f.write('}\n')
        f.write('quit;\n')
    return path


b3g, cs = form(Y, 3, 'c')
eqs4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
sols_b = sp.solve(eqs4, cs, dict=True)
a5, as_ = form(Y, 5, 'a')

for k, sb in enumerate(sols_b):
    b3k = sp.expand(b3g.subs(sb))
    if b3k == 0:
        print(f'branch {k}: b3 == 0 -> EMPTY by domain argument (E0 = y1 det A); skip.')
        continue
    Es, A = E_components(a5, b3k)
    Jraw = []
    for e in Es:
        Jraw += coeffs_in_y(e)
    Jpolys = numden_clear(Jraw)
    Tpolys = numden_clear(coeffs_in_y(A.det(method='berkowitz')))
    unknowns = sorted({v for v in (set().union(*[e.free_symbols for e in Es]) | A.free_symbols)
                       if v not in set(Y)}, key=str)
    for prime in PRIMES:
        path = emit(k, unknowns, Jpolys, Tpolys, prime)
        print(f'branch {k}, mod {prime}: {len(unknowns)} vars, {len(Jpolys)} J-gens, '
              f'{len(Tpolys)} T-gens -> {path}')
