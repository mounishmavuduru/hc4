# ############################################################################
# # SUPERSEDED -- DO NOT USE.  This sliced test gives a FALSE cone on the      #
# # survivor family: the det A != 0 (rank-3) locus is a LOWER-dimensional      #
# # component of V(J), which a generic codim-6 slice misses, so the slice      #
# # wrongly reports CONE.  Use the full-variety Rabinowitsch instead           #
# # (_d5_surv_rabin.py / _d5_surv_target.py) or the exact char-0 sat test      #
# # (clean_cone_test in _d5_survivors.py).  Kept only as a documented pitfall. #
# ############################################################################
# _d5_surv_slice.py -- FALLBACK survivor test if sat(J,detA) over F_p is too
# heavy on the dim-6 variety.  Slice V(J) (dim 6 in 21 a-vars, fixed c mod p)
# by 6 GENERIC linear forms -> a 0-dimensional ideal J' with V(J') subset V(J).
# Then a single Rabinowitsch on a random combination t = sum lambda_i * T_i of
# the det Hess3(a5) y-coefficients:
#
#     1 in <J', t*w - 1>   <=>  t not vanishing nowhere... i.e. t in sqrt(J')
#                          <=>  t == 0 at every point of V(J')
#     (generic lambda)     <=>  every T_i == 0 at every point of V(J')
#                          <=>  det Hess3(a5) == 0 on V(J').
#
# If detA =/= 0 somewhere on the dim-6 V(J), a generic 6-plane meets that locus,
# so some slice point has detA =/= 0 -> t not in sqrt(J') -> SURVIVOR.  If detA
# == 0 on all of V(J), every slice point kills it -> CONE.  All computations are
# 0-dimensional (fast) over the FIELD F_p.
#
#   py -u _d5_surv_slice.py
import os, re, sys, random
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, as_, E_components, coeffs_in_y, AV, s_poly

NODEDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survslice')
os.makedirs(NODEDIR, exist_ok=True)
y1, y2, y3 = Y
PRIMES = [32003, 40009, 15013]
PTS_PER_PRIME = 3
SLICES = 6          # codim of V(J) in the 21-dim a-space (dim V(J)=6)


def run_here(text):
    fn = os.path.join(NODEDIR, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    wpath = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    import subprocess
    p = subprocess.run(base.WSL + ['bash', '-c', f'timeout 200 Singular -q < "{wpath}"'],
                       capture_output=True, text=True)
    return p.stdout + p.stderr


def build_JT_numeric(c0, c1, c2):
    b3 = c2 * y1**3 + c0 * y1**2 * y2 + c1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    T = [q for q in coeffs_in_y(sp.expand(A.det(method='berkowitz'))) if q != 0]
    return J, T


def test(p, c0, c1, c2, rng):
    J, T = build_JT_numeric(c0, c1, c2)
    lin = []
    for _ in range(SLICES):
        terms = [str(rng.randrange(1, p)) + '*' + str(v) for v in as_]
        lin.append(' + '.join(terms) + ' + ' + str(rng.randrange(1, p)))
    lam = [rng.randrange(1, p) for _ in T]
    tcomb = ' + '.join(f'{l}*({s_poly(g)})' for l, g in zip(lam, T))
    txt = [f'ring R = {p},({AV},w),dp;',
           'ideal J = ' + ',\n'.join(s_poly(g) for g in J) + ',\n' + ',\n'.join(lin) + ';',
           'ideal G = std(J);',
           'if (size(G)==1 && G[1]==1){ "SLICEEMPTY"; quit; }',
           '"VDIM", vdim(G);',
           f'poly t = {tcomb};',
           'ideal Q = J, t*w - 1;',
           'ideal g = std(Q);',
           'if (size(g)==1 && g[1]==1){ "CONE"; } else { "SURV"; }',
           'quit;']
    out = run_here('\n'.join(txt) + '\n')
    mvd = re.search(r'VDIM\s+(-?\d+)', out)
    vd = int(mvd.group(1)) if mvd else None
    if 'SLICEEMPTY' in out:
        return 'SLICEEMPTY', vd, out
    if 'CONE' in out:
        return 'CONE', vd, out
    if 'SURV' in out:
        return 'SURV', vd, out
    return None, vd, out


def main():
    rng = random.Random(7)
    print('=== survivor family: sliced radical test (fallback) ===')
    print('    b3 = y1^2*(c2 y1 + c0 y2 + c1 y3); 6-plane slice, then Rabinowitsch\n')
    allcone = True
    survs = []
    n = 0
    for p in PRIMES:
        for _ in range(PTS_PER_PRIME):
            c0, c1, c2 = (rng.randrange(1, p) for _ in range(3))
            code, vd, out = test(p, c0, c1, c2, rng)
            n += 1
            if code == 'CONE':
                print(f'  p={p} c=({c0},{c1},{c2}): vdim={vd}, detA=0 on slice -> CONE')
            elif code == 'SLICEEMPTY':
                print(f'  p={p} c=({c0},{c1},{c2}): slice empty (bad slice) -- retrying skipped')
                allcone = False
            elif code == 'SURV':
                print(f'  p={p} c=({c0},{c1},{c2}): vdim={vd}, detA != 0 on a slice point -> *** SURVIVOR ***')
                allcone = False
                survs.append((p, c0, c1, c2))
            else:
                print(f'  p={p} c=({c0},{c1},{c2}): PARSE/TIMEOUT; tail {out[-140:].strip()!r}')
                allcone = False
    print()
    if allcone and n:
        print(f'RESULT ({n} tests): detA=0 on every slice -> survivor family is a CONE -> EMPTY.')
        return 0
    if survs:
        print(f'RESULT: {len(survs)} GENUINE survivor point(s):', survs)
        return 2
    print('RESULT: inconclusive -- inspect.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
