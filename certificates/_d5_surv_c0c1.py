# _d5_surv_c0c1.py -- decide the c0=c1=0 sub-locus of the R-D5-TAIL survivor
# family, where b3 = c2*y1^3 and the general excluding direction v*=(0,c1,-c0)
# VANISHES (so M2's argument is vacuous there).  Hess(b3) has rank 1 with a
# 2-dimensional kernel span{e2,e3}.  Decide, over F_p at random c2 (multi-prime):
#   (0) is J = (1)?  -> no a5 solution at all -> EMPTY.
#   (1) det Hess3(a5) in sqrt(J)? -> every solution rank<3 (a cone) -> EMPTY.
#   (2) does every solution carry an isotropic/linear direction among the kernel
#       e2,e3?  D^2_{e2}a5 = e2^T A e2 in sqrt(J)  or  e3^T A e3 in sqrt(J)?
#       (either one being forced would exclude from the no-linear-direction
#        branch, closing this sub-locus like the general family.)
# Full-variety single-combination Rabinowitsch, 3 primes.
#
#   py -u _d5_surv_c0c1.py
import os, re, random, subprocess, sys
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, E_components, coeffs_in_y, hess, AV, s_poly

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'c0c1')
os.makedirs(ND, exist_ok=True)
PRIMES = [32003, 40009, 15013]


def run(text, tmo=260):
    fn = os.path.join(ND, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


def build(c2v):
    b3 = c2v * y1**3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    return J, A


def radical_vanishes(p, J, G, rng):
    """is every y-coeff of the polynomial(s) G in sqrt(J)? single-combo Rabinowitsch."""
    Gc = []
    for g in G:
        Gc += [q for q in coeffs_in_y(sp.expand(g)) if q != 0]
    if not Gc:
        return True, ''
    lam = [rng.randrange(1, p) for _ in Gc]
    t = ' + '.join(f'{l}*({s_poly(g)})' for l, g in zip(lam, Gc))
    txt = [f'ring R={p},({AV},w),dp;',
           'ideal J=' + ',\n'.join(s_poly(g) for g in J) + ';',
           'ideal S=std(J); if(size(S)==1 && S[1]==1){"JUNIT";}',
           f'poly t={t};', 'ideal Q=J,t*w-1; ideal g=std(Q);',
           'if(size(g)==1 && g[1]==1){"INRAD";}else{"NOTRAD";}', 'quit;']
    out = run('\n'.join(txt) + '\n')
    if 'JUNIT' in out:
        return 'JUNIT', out
    return ('INRAD' in out), out


def main():
    rng = random.Random(31)
    print('=== c0=c1=0 sub-locus  b3 = c2*y1^3  (v* vanishes) ===\n')
    verdict = 'EMPTY'
    for p in PRIMES:
        c2v = rng.randrange(1, p)
        J, A = build(c2v)
        detA = A.det(method='berkowitz')
        e2, e3 = sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])
        d22 = (e2.T * A * e2)[0]      # D^2_{e2} a5
        d33 = (e3.T * A * e3)[0]      # D^2_{e3} a5
        rd, _ = radical_vanishes(p, J, [detA], rng)
        r22, _ = radical_vanishes(p, J, [d22], rng)
        r33, _ = radical_vanishes(p, J, [d33], rng)
        if rd == 'JUNIT' or rd is True:
            tag = 'J=(1) -> EMPTY' if rd == 'JUNIT' else 'detA in sqrt(J) -> rank<3 cone -> EMPTY'
        elif r22 is True or r33 is True:
            tag = f'rank-3 exist but isotropic: D2_e2 in radJ={r22}, D2_e3 in radJ={r33} -> EXCLUDED -> EMPTY'
        else:
            tag = (f'detA in radJ={rd}, D2_e2 in radJ={r22}, D2_e3 in radJ={r33} '
                   '-> NEITHER cone nor forced isotropic -> SURVIVOR (needs analysis)')
            verdict = 'SURVIVOR'
        print(f'  p={p} c2={c2v}: {tag}')
    print()
    print('SUB-LOCUS VERDICT:', verdict)
    return 0 if verdict == 'EMPTY' else 2


if __name__ == '__main__':
    sys.exit(main())
