# _d5_close.py -- fully automated recursive closure of the degree-5 rank-3
# pivot-free case.  Drives WSL Singular from Windows sympy.
#
# Node = a set S of polynomial constraints on the coefficients c0..c9 of the
# generic cubic b3.  At each node we solve  E4=0 AND S  on the GENERIC cubic
# (so no branch-denominators -> no NaN), and for each resulting b3 we decide
# the a-system J = <coeffs_y E0,E1,E2,E3> over Q(free params):
#   * std(J) = (1)               -> EMPTY off a deeper locus V(D); D is read
#                                   from a lift certificate, factored, and we
#                                   recurse with S ∪ {f=0} for each irreducible
#                                   factor f (that is not a nonzero constant and
#                                   is not already implied).
#   * det A reduces to 0 mod G   -> a5 forced to a cone -> EMPTY (leaf).
#   * else Rabinowitsch: every det-A coeff in sqrt(J) -> cone -> EMPTY (leaf);
#                        otherwise -> SURVIVOR (a pivot-free potential!) -> stop.
#
# Termination: each recursion adds a constraint, dropping dimension; the tree
# is finite (Noetherian).  Depth/'node budget caps guard runaway.
#
#   py _d5_close.py            # prints the full decision tree
import itertools, os, re, subprocess, sys
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
NODEDIR = os.path.join(os.getcwd(), '_lean', 'close')
os.makedirs(NODEDIR, exist_ok=True)

WSL = ['wsl', '-d', 'Ubuntu-24.04', '--']
DEPTH_CAP = 8
NODE_CAP = 400

y1, y2, y3, s = sp.symbols('y1 y2 y3 s')
Y = (y1, y2, y3)
C = list(sp.symbols('c0:10'))


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


a5, as_ = form(Y, 5, 'a')
b3g, _cs = form(Y, 3, 'c')
E4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])   # a-free constraints (pure c)
AV = ','.join(str(a) for a in as_)


def run_singular(text):
    fn = os.path.join(NODEDIR, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    wpath = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    p = subprocess.run(WSL + ['bash', '-c', f'timeout 300 Singular -q < "{wpath}"'],
                       capture_output=True, text=True)
    return p.stdout + p.stderr


def build_JT(b3):
    Es, A = E_components(a5, b3)
    Jsym = []
    for e in Es:
        Jsym += coeffs_in_y(e)
    Jp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jsym if p != 0]
    Tp = [sp.expand(sp.fraction(sp.together(p))[0]) for p in coeffs_in_y(A.det(method='berkowitz')) if p != 0]
    return Jp, Tp


def params_of(polys):
    return sorted({v for p in polys for v in p.free_symbols if str(v).startswith('c')}, key=str)


def decide(Jp, Tp):
    """returns (verdict, data). verdict in EMPTY1/CONE/SURV/DEEPER; data=D for DEEPER."""
    pars = params_of(Jp + Tp)
    par = ','.join(str(p) for p in pars)
    coeff = f'(0,{par})' if par else '0'
    # phase 1: std, check (1); if (1), lift for D; else reduce; else rabinowitsch
    txt = [f'ring R = {coeff},({AV}),dp;']
    txt.append('ideal J = ' + ',\n'.join(s_poly(p) for p in Jp) + ';')
    txt.append('ideal G = std(J);')
    txt.append('if (size(G)==1 && G[1]==1){ "V1"; matrix Tc = lift(J, ideal(1)); int i;'
               ' for(i=1;i<=size(J);i++){ "L",i,"=",Tc[i,1]; } "ENDL"; quit; }')
    txt.append('"DIM", dim(G);')
    txt.append('ideal T = ' + ',\n'.join(s_poly(p) for p in Tp) + ';')
    txt.append('ideal red = reduce(T,G); "NZ", size(simplify(red,2));')
    out = run_singular('\n'.join(txt) + '\n')
    if 'V1' in out:
        # parse lift denominators -> D
        syms = {str(v): v for v in (as_ + C)}
        dens = []
        for m in re.finditer(r'^L\s+\d+\s*=\s*(.*)$', out, re.M):
            expr = m.group(1).strip()
            if expr in ('0', ''):
                continue
            e = sp.sympify(expr.replace('^', '**'), locals=syms)
            _, d = sp.fraction(sp.together(e))
            dens.append(sp.expand(d))
        D = sp.lcm([sp.Integer(1)] + dens) if dens else sp.Integer(1)
        return 'EMPTY1', D
    mnz = re.search(r'NZ\s+(\d+)', out)
    if mnz and int(mnz.group(1)) == 0:
        return 'CONE', None
    # radical test: every T coeff in sqrt(J)?
    txt2 = [f'ring R = {coeff},({AV}),dp;']
    txt2.append('ideal J = ' + ',\n'.join(s_poly(p) for p in Jp) + ';')
    txt2.append('ideal G = std(J);')
    txt2.append('ideal T = ' + ',\n'.join(s_poly(p) for p in Tp) + ';')
    txt2.append('int bad=0; int i;')
    txt2.append(f'for(i=1;i<=size(T);i++){{ ring Rw={coeff},({AV},w),dp;'
                ' ideal G=imap(R,G); ideal T=imap(R,T); ideal Q=G,T[i]*w-1; ideal g=std(Q);'
                ' if(size(g)!=1||g[1]!=1){bad=bad+1;} setring R; kill Rw; }')
    txt2.append('"BAD", bad;')
    out2 = run_singular('\n'.join(txt2) + '\n')
    mb = re.search(r'BAD\s+(\d+)', out2)
    if mb and int(mb.group(1)) == 0:
        return 'CONE', None
    return 'SURV', out2


def factors_of(D):
    fs = []
    for base, _ in sp.factor_list(D)[1]:
        if base.free_symbols:       # skip numeric constants
            fs.append(sp.expand(base))
    return fs


NODES = 0
LEAVES = []
SURV = []


def close(S, depth, tag):
    """S = list of constraint polynomials (=0). Recurse."""
    global NODES
    NODES += 1
    if NODES > NODE_CAP:
        raise RuntimeError('node budget exceeded')
    if depth > DEPTH_CAP:
        LEAVES.append((tag, 'DEPTHCAP', S)); print(f'{"  "*depth}[{tag}] DEPTH CAP'); return
    # solve E4=0 AND S on the generic cubic
    eqs = list(E4) + list(S)
    sols = sp.solve(eqs, C, dict=True)
    if not sols:
        LEAVES.append((tag, 'NOSOL', S)); print(f'{"  "*depth}[{tag}] no b3 (constraints inconsistent) -> EMPTY'); return
    for si, sb in enumerate(sols):
        b3 = sp.expand(b3g.subs(sb))
        subtag = f'{tag}.{si}'
        if b3 == 0:
            print(f'{"  "*depth}[{subtag}] b3==0 -> EMPTY (domain arg)'); LEAVES.append((subtag, 'B3ZERO', S)); continue
        Jp, Tp = build_JT(b3)
        if not Jp:
            print(f'{"  "*depth}[{subtag}] J empty -> check detA');
        verdict, data = decide(Jp, Tp)
        if verdict == 'CONE':
            print(f'{"  "*depth}[{subtag}] CONE (detA in sqrt J) -> EMPTY'); LEAVES.append((subtag, 'CONE', S))
        elif verdict == 'SURV':
            print(f'{"  "*depth}[{subtag}] *** SURVIVOR CANDIDATE ***'); SURV.append((subtag, S, data))
        elif verdict == 'EMPTY1':
            D = data
            fs = factors_of(D)
            print(f'{"  "*depth}[{subtag}] J=(1) EMPTY off V(D), D-factors: {[str(f) for f in fs]}')
            for f in fs:
                # skip if f already implied by S (numeric reduce)
                close(list(S) + [f], depth + 1, f'{subtag}|{f}')
        else:
            print(f'{"  "*depth}[{subtag}] {verdict}'); LEAVES.append((subtag, verdict, S))


if __name__ == '__main__':
    print('=== recursive closure: degree-5 rank-3 pivot-free ===')
    close([], 0, 'root')
    print()
    print(f'nodes visited: {NODES}')
    print(f'survivor candidates: {len(SURV)}')
    for t, S, _ in SURV:
        print('  SURVIVOR at', t, 'constraints', S)
    ok = len(SURV) == 0 and all(v in ('CONE', 'B3ZERO', 'NOSOL') for _, v, _ in LEAVES)
    print('ALL CHECKS PASSED' if ok else 'INCOMPLETE / see leaves')
    # report any non-clean leaves
    for t, v, S in LEAVES:
        if v not in ('CONE', 'B3ZERO', 'NOSOL'):
            print('  LEAF', t, v, S)
    sys.exit(0 if ok else 1)
