# d5_lindir_decide.py
#
# The degree-5 residual branch: rank-3 leading form f5 WITH a linear direction.
# Rank-3 reduction gives f = a(x') + x4 b(x') + (1/2) x4^2 tau(x'); the weighted
# leading form is F = a5(y) + x4 b3(y) + (1/2) x4^2 y1 (y1 = tau-direction), and
# det Hess_4 f in C^x forces det Hess_4 F = 0, i.e. the graded tower
# E4 = E3 = E2 = E1 = E0 = 0.  E4 = 0 (pure b3) splits the cubic into two
# rational branches, exactly as in the pivot-free case.
#
# Here a5 is NOT generic: it has a linear direction.  Two sub-cases by whether
# that direction is the tau-direction y1:
#   A) a5 = P(y2,y3) + y1 Q(y2,y3)   (linear direction e1 == y1)
#   B) a5 = R(y1,y2) + y3 S(y1,y2)   (linear direction e3 =/= y1)
# each with 11 a-parameters (vs 21 generic).
#
# QUESTION decided here: does a rank-3 member exist?  With a5 constrained to the
# linear-direction form, is det Hess3(a5) == 0 forced on the solution variety
# V(J), J = <coeff_y E0,E1,E2,E3>?
#   * det A in sqrt(J)  => every solution is rank<3 => NO rank-3 member with a
#     linear direction of THIS shape => that (sub-case, b3-branch) is EMPTY.
#   * det A not in sqrt(J) => rank-3 members exist => the branch is NON-empty
#     (then pivot-freeness / further gates must be examined).
#
# Test: unbiased full-variety single-combination Rabinowitsch over F_p at random
# c on the E4=0 branch (fixed a-ring over F_p is fast).  Multi-prime.
#
#   py -u d5_lindir_decide.py
import os, re, random, subprocess, sys
import sympy as sp
import _d5_close as base
from _d5_close import Y, E_components, coeffs_in_y, form, hess, s_poly

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'lindir')
os.makedirs(ND, exist_ok=True)
PRIMES = [32003, 40009, 15013]


def homog2(vs, d, tag):
    cs = sp.symbols(f'{tag}0:{d + 1}')
    return sum(c * vs[0]**(d - i) * vs[1]**i for i, c in enumerate(cs)), list(cs)


b3g, C = form(Y, 3, 'c')
E4 = coeffs_in_y(hess(b3g, Y).adjugate()[0, 0])
B3_BRANCHES = sp.solve(E4, C, dict=True)     # the two E4=0 families


def a5_form(sub):
    if sub == 'A':
        P, Pc = homog2((y2, y3), 5, 'p'); Q, Qc = homog2((y2, y3), 4, 'q')
        return sp.expand(P + y1 * Q), Pc + Qc
    R, Rc = homog2((y1, y2), 5, 'r'); S, Sc = homog2((y1, y2), 4, 's')
    return sp.expand(R + y3 * S), Rc + Sc


def s_poly_modp(expr, p, avars):
    """Emit expr as a Singular poly string with every integer coefficient
    reduced mod p (so no literal exceeds p -- large literals crash the parser)."""
    P = sp.Poly(sp.expand(expr), *avars)
    terms = []
    for mon, co in P.terms():
        c = int(co) % p
        if c == 0:
            continue
        factors = [str(c)] + [f'{v}^{e}' if e > 1 else str(v)
                              for v, e in zip(avars, mon) if e > 0]
        terms.append('*'.join(factors))
    return ' + '.join(terms) if terms else '0'


def run(text, tmo=270):
    fn = os.path.join(ND, 'node.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


def decide(sub, bi, sb, rng):
    a5, av = a5_form(sub)
    b3 = sp.expand(b3g.subs(sb))
    # remaining free c-params after the E4=0 substitution
    cfree = sorted({v for v in b3.free_symbols if str(v).startswith('c')}, key=str)
    # random values for the free c's (mod p handled in-Singular via integer coeffs)
    results = []
    AV = ','.join(str(a) for a in av)
    for p in PRIMES:
        cs = {c: rng.randrange(1, p) for c in cfree}
        b3n = sp.expand(b3.subs(cs))
        (E3, E2, E1, E0), A = E_components(a5, b3n)
        J = []
        for e in (E3, E2, E1, E0):
            J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
        T = [q for q in coeffs_in_y(sp.expand(A.det(method='berkowitz'))) if q != 0]
        if not T:
            results.append((p, 'DETA0')); continue
        lam = [rng.randrange(1, p) for _ in T]
        tcomb = ' + '.join(f'{l % p}*({s_poly_modp(g, p, av)})' for l, g in zip(lam, T))
        txt = [f'ring R={p},({AV},w),dp;',
               'ideal J=' + ',\n'.join(s_poly_modp(g, p, av) for g in J) + ';',
               'ideal S=std(J); "DIMJ",dim(S);',
               f'poly t={tcomb};',
               'ideal Q=J,t*w-1; ideal g=std(Q);',
               'if(size(g)==1 && g[1]==1){"CONE";}else{"SURV";}',
               'quit;']
        out = run('\n'.join(txt) + '\n')
        md = re.search(r'DIMJ\s+(-?\d+)', out)
        code = 'CONE' if 'CONE' in out else 'SURV' if 'SURV' in out else 'ERR'
        results.append((p, code, md.group(1) if md else '?'))
    return results


def main():
    rng = random.Random(3)
    print('=== degree-5 rank-3 WITH linear direction: do rank-3 members exist? ===')
    print(f'    E4=0 has {len(B3_BRANCHES)} b3-branches; sub-cases A,B; primes {PRIMES}\n')
    verdict = {}
    for sub in ('A', 'B'):
        for bi, sb in enumerate(B3_BRANCHES):
            res = decide(sub, bi, sb, rng)
            codes = [r[1] for r in res]
            tag = f'sub {sub}, b3-branch {bi}'
            if all(c == 'CONE' for c in codes):
                print(f'[{tag}] det A in sqrt(J) all primes -> NO rank-3 member -> EMPTY   {res}')
                verdict[(sub, bi)] = 'EMPTY'
            elif any(c == 'SURV' for c in codes):
                print(f'[{tag}] det A NOT in sqrt(J) -> rank-3 members EXIST -> NON-EMPTY   {res}')
                verdict[(sub, bi)] = 'NONEMPTY'
            else:
                print(f'[{tag}] inconclusive {res}')
                verdict[(sub, bi)] = 'INCONCLUSIVE'
    print()
    print('SUMMARY:', verdict)
    return 0


if __name__ == '__main__':
    sys.exit(main())
