# _d5_surv_rabin.py -- unbiased modular CONE test for the survivor family
# b3 = y1^2*(c2 y1 + c0 y2 + c1 y3), over the FULL variety V(J) (no slicing).
#
# CONE  <=>  det Hess3(a5) == 0 on V(J)  <=>  every y-coeff T_i of det A in
# sqrt(J).  Test it with a SINGLE random combination t = sum lambda_i T_i and
# one Rabinowitsch:
#         1 in <J, t*w - 1>   <=>   t in sqrt(J).
# For generic lambda over F_p, t vanishes on a component W of V(J) iff every
# T_i does, so  t in sqrt(J)  <=>  all T_i in sqrt(J)  <=>  CONE.  This is the
# full-variety radical test (unlike the sliced/biased fallback), and it is one
# std instead of the 55 the full sat computes -- feasible where sat timed out.
#
# Rigor: several primes, several c, several independent lambda.  A false CONE
# would need t in sqrt(J) for generic lambda while some T_i is not -- impossible
# by the genericity argument; a real non-cone shows t not in sqrt(J) for almost
# every lambda.  Own node dir (_lean/survrabin).
#
#   py -u _d5_surv_rabin.py
import os, re, sys, random
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, as_, E_components, coeffs_in_y, AV, s_poly

NODEDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survrabin')
os.makedirs(NODEDIR, exist_ok=True)
y1, y2, y3 = Y
PRIMES = [32003, 40009, 15013]
PTS_PER_PRIME = 2
LAMBDAS = 2          # independent random combinations per (p,c)


def run_here(text, tmo=270):
    fn = os.path.join(NODEDIR, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    wpath = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    import subprocess
    p = subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{wpath}"'],
                       capture_output=True, text=True)
    return p.stdout + p.stderr


def build_JT(c0, c1, c2):
    b3 = c2 * y1**3 + c0 * y1**2 * y2 + c1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    T = [q for q in coeffs_in_y(sp.expand(A.det(method='berkowitz'))) if q != 0]
    return J, T


def rabin(p, J, T, lam):
    tcomb = ' + '.join(f'{l}*({s_poly(g)})' for l, g in zip(lam, T))
    txt = [f'ring R = {p},({AV},w),dp;',
           'ideal J = ' + ',\n'.join(s_poly(g) for g in J) + ';',
           f'poly t = {tcomb};',
           'ideal Q = J, t*w - 1;',
           'ideal g = std(Q);',
           'if (size(g)==1 && g[1]==1){ "TINRAD"; } else { "TNOTRAD"; }',
           'quit;']
    out = run_here('\n'.join(txt) + '\n')
    if 'TINRAD' in out:
        return 'CONE', out
    if 'TNOTRAD' in out:
        return 'SURV', out
    return None, out


def main():
    rng = random.Random(11)
    print('=== survivor family: full-variety single-combination Rabinowitsch ===')
    print('    b3 = y1^2*(c2 y1 + c0 y2 + c1 y3); t = sum lam_i * detA_coeff_i\n')
    allcone, survs, n = True, [], 0
    for p in PRIMES:
        for _ in range(PTS_PER_PRIME):
            c0, c1, c2 = (rng.randrange(1, p) for _ in range(3))
            J, T = build_JT(c0, c1, c2)
            for li in range(LAMBDAS):
                lam = [rng.randrange(1, p) for _ in T]
                code, out = rabin(p, J, T, lam)
                n += 1
                if code == 'CONE':
                    print(f'  p={p} c=({c0},{c1},{c2}) lam#{li}: t in sqrt(J) -> CONE')
                elif code == 'SURV':
                    print(f'  p={p} c=({c0},{c1},{c2}) lam#{li}: t NOT in sqrt(J) -> *** SURVIVOR ***')
                    allcone = False
                    survs.append((p, c0, c1, c2, li))
                else:
                    print(f'  p={p} c=({c0},{c1},{c2}) lam#{li}: PARSE/TIMEOUT; tail {out[-140:].strip()!r}')
                    allcone = False
    print()
    if allcone and n:
        print(f'RESULT ({n} tests): det Hess3(a5) in sqrt(J) on the FULL variety '
              'across all primes/points/lambdas.')
        print('=> survivor family b3=y1^2*(...) is a rank<3 CONE -> EMPTY (unbiased).')
        return 0
    if survs:
        print(f'RESULT: {len(survs)} GENUINE survivor test(s):', survs)
        return 2
    print('RESULT: some test timed out/inconclusive -- inspect above.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
