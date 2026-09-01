# _d5_fiber.py -- fiberwise mod-p decision for the degree-5 rank-3 pivot-free
# branches.  We FIX the free b3 parameters (the c-unknowns) to random integers
# on the branch, collapsing J and T to ideals in the a-unknowns only (the 21
# coefficients of the quintic a5).  These are low-degree (<=3) with numeric
# coefficients, so std over F_p is fast.
#
# Per fiber we decide:  is  E0=E1=E2=E3=0  =>  det Hess_3(a5) == 0 ?
#   std(J) = (1)               -> residual inconsistent (no a5 at all) -> EMPTY
#   every T-gen in sqrt(J)     -> det A forced 0            -> fiber EMPTY
#   some  T-gen not in sqrt(J) -> a pivot-free a5 EXISTS at this b3 (SURVIVOR)
#
# Fiberwise EMPTY for many random b3 is strong structured evidence (it varies
# ALL of a exactly via a Groebner basis, per fixed b3 -- stronger than random
# full-point sampling).  A survivor is an exact counterexample candidate.
#
#   py _d5_fiber.py [nfibers] [prime]
# writes _lean/fiber/fib_b{k}_{j}_p{prime}.sing (fed to Singular via stdin).

import itertools, os, sys, random
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.getcwd(), '_lean', 'fiber')
os.makedirs(OUT, exist_ok=True)

NFIB = int(sys.argv[1]) if len(sys.argv) > 1 else 4
PRIME = sys.argv[2] if len(sys.argv) > 2 else '32003'
random.seed(20260830)

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


def emit(k, j, avars, Jpolys, Tpolys):
    ring = ','.join(str(u) for u in avars)
    path = os.path.join(OUT, f'fib_b{k}_{j}_p{PRIME}.sing')
    Jp = [p for p in Jpolys if p != 0]
    Tp = [p for p in Tpolys if p != 0]
    with open(path, 'w', newline='\n') as f:
        f.write(f'// fiber decision: branch {k}, fiber {j}, mod {PRIME}\n')
        f.write(f'ring R = {PRIME},({ring}),dp;\n')
        f.write('ideal J = ' + ',\n  '.join(s_poly(p) for p in Jp) + ';\n')
        f.write('ideal G = std(J);\n')
        f.write(f'"b{k} fib{j}: GBsize", size(G);\n')
        f.write('if (size(G)==1 && G[1]==1) { "  RESULT: J=(1) inconsistent -> fiber EMPTY"; quit; }\n')
        f.write('ideal T = ' + ',\n  '.join(s_poly(p) for p in Tp) + ';\n')
        f.write('int bad=0; int i;\n')
        f.write('for (i=1;i<=size(T);i++){\n')
        f.write(f'  ring Rw={PRIME},({ring},w),dp;\n')
        f.write('  ideal G=imap(R,G); ideal T=imap(R,T);\n')
        f.write('  ideal Q=G,T[i]*w-1; ideal g=std(Q);\n')
        f.write('  if(size(g)!=1||g[1]!=1){bad=bad+1;}\n')
        f.write('  setring R; kill Rw;\n')
        f.write('}\n')
        f.write('if(bad==0){"  RESULT: det A forced 0 -> fiber EMPTY";}\n')
        f.write('else{"  RESULT: SURVIVOR -- pivot-free a5 exists, bad=",bad;}\n')
        f.write('quit;\n')
    return path, len(Jp), len(Tp)


b3g, cs = form(Y, 3, 'c')
eqs4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
sols_b = sp.solve(eqs4, cs, dict=True)
a5, as_ = form(Y, 5, 'a')

manifest = []
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
    freec = sorted({v for e in Es for v in e.free_symbols if str(v).startswith('c')}, key=str)
    for j in range(NFIB):
        sub = {c: random.randint(-6, 6) for c in freec}
        # avoid the all-zero b3 fiber
        if all(val == 0 for val in sub.values()):
            sub[freec[0]] = 3
        Jp = [sp.expand(sp.nsimplify(p.subs(sub))) for p in Jsym]
        Tp = [sp.expand(p.subs(sub)) for p in Tsym]
        # clear rational denominators generator-wise
        Jp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jp]
        Tp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Tp]
        path, nj, nt = emit(k, j, as_, Jp, Tp)
        manifest.append(path)
        print(f'branch {k} fiber {j}: c={sub} | {len(as_)} a-vars, {nj} J-gens, {nt} T-gens -> {os.path.basename(path)}')

with open(os.path.join(OUT, 'manifest.txt'), 'w', newline='\n') as f:
    f.write('\n'.join(os.path.basename(p) for p in manifest) + '\n')
print(f'\n{len(manifest)} fiber scripts written to {OUT}')
