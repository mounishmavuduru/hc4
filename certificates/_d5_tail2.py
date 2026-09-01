# _d5_tail2.py -- finish the two components the first pass left open:
#  (A) b1_2c0c5=c1c4 : reduce(det A, std J) != 0, but det A may lie in sqrt(J).
#      Decide rigorously by Rabinowitsch over Q(c0,c2,c4,c5): det A coeff t is
#      in sqrt(J)  <=>  1 in <J, t*w-1>.  Component EMPTY (a5 forced cone)
#      <=> every det-A coefficient is in sqrt(J).
#  (B) monomial punctures c5=0 (branch1) and c9=0 (branch0): the solved branch
#      b3k has that parameter in a denominator, so direct substitution gives
#      NaN.  Re-derive from the GENERIC cubic: impose E4=0 AND the puncture,
#      re-solve for a parametrization, and run the parametric a-decision on
#      each resulting sub-branch.
#
#   py _d5_tail2.py   # writes _lean/tail/{radA_*.sing, punc_*.sing}
import itertools, os
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.getcwd(), '_lean', 'tail')
os.makedirs(OUT, exist_ok=True)

y1, y2, y3, s = sp.symbols('y1 y2 y3 s')
Y = (y1, y2, y3)
C = sp.symbols('c0:10')
c0, c1, c2, c3, c4, c5, c6, c7, c8, c9 = C


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


def build_JT(b3):
    Es, A = E_components(a5, b3)
    Jsym = []
    for e in Es:
        Jsym += coeffs_in_y(e)
    Jp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jsym if p != 0]
    Tp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in coeffs_in_y(A.det(method='berkowitz')) if p != 0]
    return Jp, Tp


a5, as_ = form(Y, 5, 'a')
b3g, cs = form(Y, 3, 'c')
eqs4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
sols_b = sp.solve(eqs4, cs, dict=True)
branch_b3 = [sp.expand(b3g.subs(sb)) for sb in sols_b]

AV = ','.join(str(a) for a in as_)


def emit_radical(label, b3, params):
    """Rabinowitsch: every det-A coeff in sqrt(J) over Q(params)?"""
    Jp, Tp = build_JT(b3)
    par = ','.join(str(p) for p in params)
    coeff = f'(0,{par})' if par else '0'
    path = os.path.join(OUT, f'radA_{label}.sing')
    with open(path, 'w', newline='\n') as f:
        f.write(f'// radical closure of component {label}; k(c)={coeff}\n')
        f.write(f'ring R = {coeff},({AV}),dp;\n')
        f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jp) + ';\n')
        f.write('ideal G = std(J);\n')
        f.write(f'"RADA {label} GBsize", size(G), " dim", dim(G);\n')
        f.write('ideal T = ' + ',\n  '.join(s_poly(p) for p in Tp) + ';\n')
        f.write('int bad=0; int i;\n')
        f.write('for(i=1;i<=size(T);i++){\n')
        f.write(f'  ring Rw={coeff},({AV},w),dp;\n')
        f.write('  ideal G=imap(R,G); ideal T=imap(R,T);\n')
        f.write('  ideal Q=G,T[i]*w-1; ideal g=std(Q);\n')
        f.write('  if(size(g)!=1||g[1]!=1){bad=bad+1;}\n')
        f.write('  setring R; kill Rw;\n')
        f.write('}\n')
        f.write('if(bad==0){"  RESULT EMPTY: every det A coeff in sqrt(J) -> a5 forced cone";}\n')
        f.write('else{"  RESULT SURVIVOR CANDIDATE: bad=",bad;}\n')
        f.write('quit;\n')
    return path, len(Jp), len(Tp)


# (A) branch1 component 2c0c5=c1c4  (c4 != 0 generic): c1 = 2c0c5/c4
b1 = branch_b3[1]
compA = {c1: 2 * c0 * c5 / c4}
b1A = sp.expand(b1.subs(compA))
pathA, njA, ntA = emit_radical('b1_2c0c5=c1c4', b1A, [c0, c2, c4, c5])
print(f'(A) radA_b1_2c0c5=c1c4: {njA} J-gens, {ntA} T-gens -> {os.path.basename(pathA)}')

# (B) punctures: re-derive from the generic cubic with the puncture imposed.
def emit_punctures(tag, base_k, extra_zero):
    """impose E4=0 AND extra_zero=0 on the generic cubic, re-solve, decide each sub-branch."""
    eqs = eqs4 + [extra_zero]
    sub_sols = sp.solve(eqs, cs, dict=True)
    outs = []
    for j, sb in enumerate(sub_sols):
        b3 = sp.expand(b3g.subs(sb))
        if b3 == 0:
            print(f'  {tag} sub {j}: b3==0 -> EMPTY (domain arg); skip'); continue
        Jp, Tp = build_JT(b3)
        params = sorted({v for p in (Jp + Tp) for v in p.free_symbols if str(v).startswith('c')}, key=str)
        par = ','.join(str(p) for p in params)
        coeff = f'(0,{par})' if par else '0'
        path = os.path.join(OUT, f'punc_{tag}_{j}.sing')
        with open(path, 'w', newline='\n') as f:
            f.write(f'// puncture {tag} sub-branch {j}; k(c)={coeff}\n')
            f.write(f'ring R = {coeff},({AV}),dp;\n')
            f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jp) + ';\n')
            f.write('ideal G = std(J);\n')
            f.write(f'"PUNC {tag} sub {j} GBsize", size(G);\n')
            f.write('if(size(G)==1&&G[1]==1){"  RESULT EMPTY: J=(1)";quit;}\n')
            f.write('"  dim", dim(G);\n')
            f.write('ideal T = ' + ',\n  '.join(s_poly(p) for p in Tp) + ';\n')
            f.write('ideal red = reduce(T,G); int nz=size(simplify(red,2));\n')
            f.write('if(nz==0){"  RESULT EMPTY: det A in J (cone)";}\n')
            f.write('else{"  RESULT UNDECIDED: nz=",nz;}\n')
            f.write('quit;\n')
        outs.append(path)
        print(f'  {tag} sub {j}: params [{par}], {len(Jp)} J-gens -> {os.path.basename(path)}')
    return outs


print('(B) branch1 puncture c5=0:')
emit_punctures('b1_c5', 1, c5)
print('(B) branch0 puncture c9=0:')
emit_punctures('b0_c9', 0, c9)
