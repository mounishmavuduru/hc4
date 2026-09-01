# _d5_surv_target.py -- full-variety radical-membership test of a chosen TARGET
# polynomial g(a,y) against J, for the survivor family b3 = y1^2*(c2y1+c0y2+c1y3).
#
# Same unbiased single-combination Rabinowitsch as _d5_surv_rabin.py but for an
# arbitrary target: is g == 0 on V(J)?  <=> every y-coeff g_i of g in sqrt(J)
# <=> (generic lambda) t = sum lam_i g_i in sqrt(J) <=> 1 in <J, t*w-1>.
#
# TARGET selectable:
#   detA   : g = det Hess3(a5)           -- rank<3 test (expect SURV: proper)
#   iso    : g = v*^T Hess3(a5) v*        -- D^2_{v*} a5, isotropy of v*=(0,c1,-c0)
#   Av*    : g = (Hess3(a5) v*)_k         -- v* in kernel of A (stronger)
#
#   py -u _d5_surv_target.py iso
import os, re, sys, random
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, E_components, coeffs_in_y, hess, grad, AV, s_poly

NODEDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survtgt')
os.makedirs(NODEDIR, exist_ok=True)
y1, y2, y3 = Y
c0, c1, c2 = sp.symbols('c0 c1 c2')
PRIMES = [32003, 40009, 15013]
PTS = 2
LAMBDAS = 2


def run_here(text, tmo=270):
    fn = os.path.join(NODEDIR, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    wpath = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    import subprocess
    p = subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{wpath}"'],
                       capture_output=True, text=True)
    return p.stdout + p.stderr


def build(target, C0, C1, C2):
    b3 = C2 * y1**3 + C0 * y1**2 * y2 + C1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    vs = sp.Matrix([0, C1, -C0])
    if target == 'detA':
        gpoly = [sp.expand(A.det(method='berkowitz'))]
    elif target == 'iso':
        gpoly = [sp.expand((vs.T * A * vs)[0])]
    elif target == 'Av':
        gpoly = [sp.expand(x) for x in (A * vs)]
    else:
        raise SystemExit('target in {detA,iso,Av}')
    G = []
    for gp in gpoly:
        G += [q for q in coeffs_in_y(gp) if q != 0]
    return J, G


def rabin(p, J, G, lam):
    tcomb = ' + '.join(f'{l}*({s_poly(g)})' for l, g in zip(lam, G))
    txt = [f'ring R = {p},({AV},w),dp;',
           'ideal J = ' + ',\n'.join(s_poly(g) for g in J) + ';',
           f'poly t = {tcomb};',
           'ideal Q = J, t*w - 1; ideal g = std(Q);',
           'if (size(g)==1 && g[1]==1){ "INRAD"; } else { "NOTRAD"; }',
           'quit;']
    out = run_here('\n'.join(txt) + '\n')
    return ('IN' if 'INRAD' in out else 'NOT' if 'NOTRAD' in out else None), out


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else 'iso'
    rng = random.Random(23)
    print(f'=== survivor family: is target "{target}" == 0 on V(J)?  '
          '(full-variety radical) ===\n')
    allin, notin, n = True, [], 0
    for p in PRIMES:
        for _ in range(PTS):
            C0, C1, C2 = (rng.randrange(1, p) for _ in range(3))
            J, G = build(target, C0, C1, C2)
            for li in range(LAMBDAS):
                lam = [rng.randrange(1, p) for _ in G]
                code, out = rabin(p, J, G, lam)
                n += 1
                if code == 'IN':
                    print(f'  p={p} c=({C0},{C1},{C2}) lam#{li}: target in sqrt(J) '
                          '-> VANISHES on V(J)')
                elif code == 'NOT':
                    print(f'  p={p} c=({C0},{C1},{C2}) lam#{li}: target NOT in sqrt(J) '
                          '-> does NOT vanish on V(J)')
                    allin = False
                    notin.append((p, C0, C1, C2, li))
                else:
                    print(f'  p={p} c=({C0},{C1},{C2}) lam#{li}: PARSE/TIMEOUT '
                          f'{out[-120:].strip()!r}')
                    allin = False
    print()
    if allin and n:
        print(f'RESULT ({n} tests): "{target}" == 0 on V(J) for all primes/points.')
        if target == 'iso':
            print('=> every V(J)-solution a5 has isotropic direction v*=(0,c1,-c0):')
            print('   D^2_{v*} a5 == 0.  Such a5 has a LINEAR DIRECTION, so it is')
            print('   EXCLUDED from the "no-isotropic-direction" branch.')
            print('   => the survivor family carries NO in-branch point -> EMPTY.')
        return 0
    print(f'RESULT: "{target}" does NOT vanish on V(J) in {len(notin)} test(s): {notin}')
    return 2


if __name__ == '__main__':
    sys.exit(main())
