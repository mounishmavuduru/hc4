# _iso_power.py -- the F_p POWER CERTIFICATE behind R-D5-TAIL fact M2:
#     (v*^T A v*)^2  in  <J>      while      v*^T A v*  is NOT in <J>,
# for the survivor family b3 = y1^2*(c2 y1 + c0 y2 + c1 y3), v* = (0,c1,-c0),
# A = Hess3 a5, J = <coeff_y E0..E3>.
#
# (v*^T A v*)^2 in <J>  =>  every y-coefficient of (v*^T A v*)^2 lies in J  =>
# D^2_{v*} a5 == 0 on V(J): a DEFINITE (non-probabilistic) witness that
# iso := v*^T A v* is in sqrt(J) at that specialization, i.e. every solution's
# a5 carries the isotropic direction v*.
#
# Upgraded (optional-cert item from the endgame audit): it was a single-prime,
# single-point print-only probe. Now it loops over several primes AND several
# random c-points and ASSERTS (fail-closed) that iso^2 reduces to 0 mod std(J)
# at every one, while also recording that iso itself does not (so the power is
# genuinely needed). Reduction is done in the ring F_p[a-vars, y1,y2,y3], where
# J (in a-vars only) generates J*F_p[a,y], so reduce(iso^k, std(J)) == 0 iff
# every y-coefficient of iso^k is in J. Needs WSL Singular.
#
#   py -u _iso_power.py
import os, re, random, subprocess, sys
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, E_components, coeffs_in_y, AV, s_poly

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survtgt')
os.makedirs(ND, exist_ok=True)
PRIMES = [32003, 40009, 15013]
PTS_PER_PRIME = 2


def run(t, tmo=250):
    fn = os.path.join(ND, 'pw.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(t)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


def test(p, C0, C1, C2):
    b3 = C2 * y1**3 + C0 * y1**2 * y2 + C1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    vs = sp.Matrix([0, C1, -C0])
    iso = sp.expand((vs.T * A * vs)[0])
    txt = [f'ring R={p},({AV},y1,y2,y3),dp;',
           'ideal J=' + ',\n'.join(s_poly(g) for g in J) + ';',
           'ideal S=std(J);',
           f'poly iso={s_poly(iso)};',
           'poly r1=reduce(iso,S);   "NZ1", (r1!=0);',
           'poly r2=reduce(iso^2,S); "NZ2", (r2!=0);',
           'quit;']
    out = run('\n'.join(txt) + '\n')
    m1 = re.search(r'NZ1\s+(\d+)', out)
    m2 = re.search(r'NZ2\s+(\d+)', out)
    if not (m1 and m2):
        return None, None, out
    return int(m1.group(1)), int(m2.group(1)), out


def main():
    rng = random.Random(9)
    print('=== M2 power certificate: (v*^T A v*)^2 in <J> at random c, several primes ===')
    allok, n = True, 0
    for p in PRIMES:
        for _ in range(PTS_PER_PRIME):
            C0, C1, C2 = (rng.randrange(1, p) for _ in range(3))
            nz1, nz2, out = test(p, C0, C1, C2)
            n += 1
            if nz1 is None:
                print(f'  p={p} c=({C0},{C1},{C2}): PARSE/TIMEOUT {out[-120:].strip()!r}')
                allok = False
                continue
            ok = (nz2 == 0)
            allok = allok and ok
            print(f'  p={p} c=({C0},{C1},{C2}): iso in J? {"no" if nz1 else "yes"};  '
                  f'iso^2 in J? {"YES" if ok else "NO"}')
            assert ok, f'p={p} c=({C0},{C1},{C2}): iso^2 NOT in <J> -- M2 power cert FAILS'
    print()
    assert allok and n == len(PRIMES) * PTS_PER_PRIME
    print(f'ALL CHECKS PASSED: (v*^T A v*)^2 in <J> at all {n} (prime, point) tests '
          '(iso itself not in J, so the square is needed) => D^2_{v*} a5 == 0 on V(J).')
    sys.exit(0)


if __name__ == '__main__':
    main()
