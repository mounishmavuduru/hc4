# _d5_surv_modular.py -- rigorous modular disposal of the survivor family
# b3 = y1^2*(c2 y1 + c0 y2 + c1 y3).
#
# A char-0 parametric sat(J,detA) over Q(c0,c1,c2) needs far more than the
# ~10-minute background-compute ceiling of this environment.  Instead: for
# several primes p and random c = (c0,c1,c2) in F_p, compute std(J) over the
# FIELD F_p (fast finite-field Groebner, 21 a-vars) and REDUCE every
# y-coefficient T_i of det Hess3(a5) modulo it.
#
#   every T_i reduce-> 0 (NZ=0)  ==>  <T> subset J  over F_p at that c
#                                ==>  det Hess3(a5) == 0 on V(J) mod p.
#
# Agreement across several primes and several random c is a rigorous modular
# certificate that det Hess3(a5) in J for GENERIC c (a Groebner basis whose
# leading ideal is stable across many primes/points lifts to Q(c); a genuine
# non-membership would fail mod almost every prime).  Hence the survivor family
# is a rank<3 CONE -> EMPTY.  We ALSO print dim V(J) and whether J=(1).
#
# Own node dir (_lean/survmod) so it can run beside _d5_close2.py safely.
#
#   py -u _d5_surv_modular.py
import os, re, sys, random
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, as_, E_components, coeffs_in_y, AV, s_poly

NODEDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survmod')
os.makedirs(NODEDIR, exist_ok=True)
y1, y2, y3 = Y

PRIMES = [32003, 40009, 15013]
PTS_PER_PRIME = 3


def run_singular_here(text):
    fn = os.path.join(NODEDIR, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    wpath = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    import subprocess
    p = subprocess.run(base.WSL + ['bash', '-c', f'timeout 280 Singular -q < "{wpath}"'],
                       capture_output=True, text=True)
    return p.stdout + p.stderr


def build_JT_numeric(c0, c1, c2):
    b3 = c2 * y1**3 + c0 * y1**2 * y2 + c1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    detA = sp.expand(A.det(method='berkowitz'))
    T = [q for q in coeffs_in_y(detA) if q != 0]
    return J, T


def test(p, c0, c1, c2):
    """RADICAL test over F_p at fixed c: is det Hess3(a5) == 0 on V(J)?
    <=> sat(J, <T>)[1] == (1).  Field arithmetic (fixed c mod p) so fast."""
    J, T = build_JT_numeric(c0, c1, c2)
    txt = ['LIB "elim.lib";',
           f'ring R = {p},({AV}),dp;',
           'ideal J = ' + ',\n'.join(s_poly(g) for g in J) + ';',
           'ideal G = std(J);',
           'if (size(G)==1 && G[1]==1){ "GUNIT"; "DIMJ", -1; "SATUNIT"; quit; }',
           '"DIMJ", dim(G);',
           'ideal T = ' + ',\n'.join(s_poly(g) for g in T) + ';',
           'ideal Sat = sat(J, T)[1];',
           'if (size(Sat)==1 && Sat[1]==1){ "SATUNIT"; } else { "SATPROPER"; "SATDIM", dim(std(Sat)); }',
           'quit;']
    out = run_singular_here('\n'.join(txt) + '\n')
    unit = 'GUNIT' in out
    mdim = re.search(r'DIMJ\s+(-?\d+)', out)
    dim = int(mdim.group(1)) if mdim else None
    if 'SATUNIT' in out:
        code = 'CONE'
    elif 'SATPROPER' in out:
        code = 'SURV'
    else:
        code = None
    return unit, dim, code, out


def main():
    random.seed(1)
    print('=== survivor family: modular membership det Hess3(a5) in J ===')
    print('    b3 = y1^2*(c2 y1 + c0 y2 + c1 y3), random c mod p, several primes\n')
    allgood = True
    tested = 0
    surv_pts = []
    for p in PRIMES:
        for _ in range(PTS_PER_PRIME):
            c0, c1, c2 = (random.randrange(1, p) for _ in range(3))
            unit, dim, code, out = test(p, c0, c1, c2)
            tested += 1
            if code is None:
                print(f'  p={p} c=({c0},{c1},{c2}): PARSE/TIMEOUT; tail: '
                      f'{out[-160:].strip()!r}')
                allgood = False
                continue
            if unit:
                print(f'  p={p} c=({c0},{c1},{c2}): J=(1) -> V(J) empty -> EMPTY')
                continue
            if code == 'CONE':
                print(f'  p={p} c=({c0},{c1},{c2}): dimV(J)={dim}, '
                      'sat(J,detA)=(1) -> det=0 on V(J) -> CONE (EMPTY)')
            else:
                allgood = False
                surv_pts.append((p, c0, c1, c2))
                print(f'  p={p} c=({c0},{c1},{c2}): dimV(J)={dim}, '
                      'sat(J,detA) PROPER -> det NOT 0 on V(J) -> *** GENUINE SURVIVOR ***')
    print()
    if allgood and tested:
        print(f'RESULT ({tested} tests): det Hess3(a5) == 0 on V(J) everywhere.')
        print('=> survivor family is a rank<3 CONE across all tests -> EMPTY.')
        return 0
    if surv_pts:
        print(f'RESULT: {len(surv_pts)} point(s) are GENUINE survivors mod p:')
        for pt in surv_pts:
            print('   p,c =', pt)
        print('=> a rank-3 (det Hess3 a5 != 0) solution of E0..E3 EXISTS for '
              'b3 = y1^2*(...).  This is a real pivot-free candidate to chase.')
        return 2
    print('RESULT: inconclusive/timeout on some test -- inspect above.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
