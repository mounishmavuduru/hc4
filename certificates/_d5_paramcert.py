# _d5_paramcert.py -- turn the parametric "J=(1) over Q(c)" result into a
# self-verifying Nullstellensatz certificate, and extract the exceptional
# locus Z.
#
# Step 1 (this script): build J in sympy, emit a Singular script that computes
#   matrix T = lift(J, ideal(1));   // over Q(c),  sum_i J[i]*T[i] = 1
# and prints the T[i] (rational functions in c) as strings.
#
# Step 2 (verify_paramcert.py): read the T[i], clear to a common denominator
# D(c), and check IN SYMPY, exactly, that  sum_i g_i(a,c) * (D*T_i)(c) == D(c)
# identically in Q[a,c].  Since D depends only on c, this proves:
#   for every c with D(c) != 0, the a-system J=0 is inconsistent (no a5).
# The exceptional locus is Z = V(D).
#
#   py _d5_paramcert.py        # writes _lean/param/cert_b{k}.sing and the
#                              # matching g-generators to cert_b{k}_gens.txt
import itertools, os
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.getcwd(), '_lean', 'param')
os.makedirs(OUT, exist_ok=True)

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
    Jp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jsym if p != 0]
    cvars = sorted({v for p in Jp for v in p.free_symbols if str(v).startswith('c')}, key=str)
    av = ','.join(str(a) for a in as_)
    par = ','.join(str(c) for c in cvars)
    # emit generator list for the sympy verifier (one per line, in J order)
    with open(os.path.join(OUT, f'cert_b{k}_gens.txt'), 'w', newline='\n') as f:
        f.write('# a-vars: ' + av + '\n')
        f.write('# c-params: ' + par + '\n')
        for p in Jp:
            f.write(str(sp.expand(p)) + '\n')
    # emit Singular lift script
    path = os.path.join(OUT, f'cert_b{k}.sing')
    with open(path, 'w', newline='\n') as f:
        f.write(f'// certificate: lift 1 into J over Q(c), branch {k}\n')
        f.write(f'ring R = (0,{par}),({av}),dp;\n')
        f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jp) + ';\n')
        f.write('matrix T = lift(J, ideal(1));\n')
        f.write('int n = size(J);\n')
        f.write('"NGEN", n;\n')
        f.write('int i;\n')
        f.write('for (i=1;i<=n;i++){ "T", i, "=", T[i,1]; }\n')
        f.write('// sanity: recompute the combination inside Singular too\n')
        f.write('poly chk = 0; for (i=1;i<=n;i++){ chk = chk + J[i]*T[i,1]; }\n')
        f.write('"SELFCHECK", chk;\n')
        f.write('quit;\n')
    print(f'branch {k}: {len(Jp)} gens, params {[str(c) for c in cvars]} -> {os.path.basename(path)}')
