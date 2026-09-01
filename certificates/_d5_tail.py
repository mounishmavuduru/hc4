# _d5_tail.py -- close the exceptional loci Z_k of the parametric decision.
# The generic-c proof (verify_paramcert.py) leaves the branches undecided only
# on Z_0 = V(c8)U V(c9)U V(3c0c9-c1c8)  and  Z_1 = V(c5)U V(2c0c5-c1c4).
# Here we take each irreducible component, substitute it into b3, rebuild the
# residual ideal J and det A in the a-variables, and emit a parametric Singular
# script over Q(remaining c) that decides the component:
#   J=(1)                        -> component EMPTY off a deeper locus
#   det A reduces to 0 mod std(J)-> a5 forced to a cone         -> EMPTY
#   else                         -> SURVIVOR / needs a deeper split
#
#   py _d5_tail.py     # writes _lean/tail/tail_<label>.sing
import itertools, os
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.getcwd(), '_lean', 'tail')
os.makedirs(OUT, exist_ok=True)

y1, y2, y3, s = sp.symbols('y1 y2 y3 s')
Y = (y1, y2, y3)
c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 = sp.symbols('c0:10')


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


def emit(label, avars, cvars, Jp, Tp):
    par = ','.join(str(c) for c in cvars) if cvars else ''
    av = ','.join(str(a) for a in avars)
    coeff = f'(0,{par})' if par else '0'
    path = os.path.join(OUT, f'tail_{label}.sing')
    Jp = [p for p in Jp if p != 0]
    Tp = [p for p in Tp if p != 0]
    with open(path, 'w', newline='\n') as f:
        f.write(f'// tail component {label}: k(c)={coeff}, vars=a\n')
        f.write(f'ring R = {coeff},({av}),dp;\n')
        f.write('ideal J = ' + (',\n  '.join(s_poly(p) for p in Jp) if Jp else '0') + ';\n')
        f.write('ideal G = std(J);\n')
        f.write(f'"COMPONENT {label}  GBsize", size(G);\n')
        f.write('if (size(G)==1 && G[1]==1) { "  RESULT EMPTY: J=(1)"; quit; }\n')
        f.write('"  dim", dim(G);\n')
        f.write('ideal T = ' + (',\n  '.join(s_poly(p) for p in Tp) if Tp else '0') + ';\n')
        f.write('ideal red = reduce(T, G);\n')
        f.write('int nz = size(simplify(red,2));\n')
        f.write('if (nz==0) { "  RESULT EMPTY: det A in J (a5 forced to a cone)"; }\n')
        f.write('else { "  RESULT UNDECIDED here: det A not in J, nz=", nz, " -- needs radical/deeper split"; }\n')
        f.write('quit;\n')
    return path, len(Jp), len(Tp)


# rebuild the two E4-branches
b3g, cs = form(Y, 3, 'c')
eqs4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
sols_b = sp.solve(eqs4, cs, dict=True)
a5, as_ = form(Y, 5, 'a')
branch_b3 = []
for sb in sols_b:
    b = sp.expand(b3g.subs(sb))
    branch_b3.append(b)

# component substitutions (generic point of each irreducible piece)
COMPONENTS = {
    'b0_c8eq0':      (0, {c8: 0}),
    'b0_c9eq0':      (0, {c9: 0}),
    'b0_3c0c9=c1c8': (0, {c1: 3 * c0 * c9 / c8}),   # c8 != 0 generic
    'b1_c5eq0':      (1, {c5: 0}),
    'b1_2c0c5=c1c4': (1, {c1: 2 * c0 * c5 / c4}),   # c4 != 0 generic
}

for label, (k, comp) in COMPONENTS.items():
    b3 = branch_b3[k]
    if b3 == 0:
        print(f'{label}: base b3==0 skip'); continue
    b3s = sp.expand(b3.subs(comp))
    Es, A = E_components(a5, b3s)
    Jsym = []
    for e in Es:
        Jsym += coeffs_in_y(e)
    Jp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jsym if p != 0]
    Tp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in coeffs_in_y(A.det(method='berkowitz')) if p != 0]
    cvars = sorted({v for p in (Jp + Tp) for v in p.free_symbols if str(v).startswith('c')}, key=str)
    path, nj, nt = emit(label, as_, cvars, Jp, Tp)
    print(f'{label}: params {[str(c) for c in cvars]}, {nj} J-gens, {nt} T-gens -> {os.path.basename(path)}')
