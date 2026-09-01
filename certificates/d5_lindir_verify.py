# d5_lindir_verify.py -- confirm and (char-0) strengthen the finding that the
# rank-3-WITH-linear-direction weighted tower is EMPTY at generic b3.
#
# For each sub-case (A: a5=P(y2,y3)+y1 Q(y2,y3); B: a5=R(y1,y2)+y3 S(y1,y2)) and
# each E4=0 b3-branch, fix the free b3-coefficients to generic integers and test
# whether the a-system J = <coeff_y E0,E1,E2,E3> (11 a-unknowns) is the UNIT
# ideal, over F_p AND over Q.  std(J) = (1)  <=>  no a5 solves the tower for
# that b3  <=>  that (sub-case, branch) is EMPTY at generic b3.
#
#   py -u d5_lindir_verify.py
import os, re, subprocess, sys, random
import sympy as sp
import _d5_close as base
from _d5_close import Y, E_components, coeffs_in_y, form, hess, s_poly
from d5_lindir_decide import a5_form, b3g, B3_BRANCHES, s_poly_modp

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'lindir')
os.makedirs(ND, exist_ok=True)


def run(text, tmo=400):
    fn = os.path.join(ND, 'verify.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


def build_J(sub, sb, cvals):
    a5, av = a5_form(sub)
    b3 = sp.expand(b3g.subs(sb).subs(cvals))
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    return J, av


def is_unit(sub, sb, cvals, field):
    J, av = build_J(sub, sb, cvals)
    AV = ','.join(str(a) for a in av)
    if field == 'Q':
        gens = ',\n'.join(s_poly(g) for g in J)
        ring = f'ring R=0,({AV}),dp;'
    else:
        p = field
        gens = ',\n'.join(s_poly_modp(g, p, av) for g in J)
        ring = f'ring R={p},({AV}),dp;'
    txt = [ring, f'ideal J={gens};', 'ideal S=std(J);',
           '"SIZE",size(S); if(size(S)==1 && S[1]==1){"UNIT";}else{"PROPER","dim",dim(S);}',
           'quit;']
    out = run('\n'.join(txt) + '\n')
    return ('UNIT' in out), out


def main():
    rng = random.Random(101)
    print('=== rank-3 WITH linear direction: is the tower EMPTY at generic b3? ===\n')
    allunit = True
    for sub in ('A', 'B'):
        for bi, sb in enumerate(B3_BRANCHES):
            cfree = sorted({v for v in sp.expand(b3g.subs(sb)).free_symbols
                            if str(v).startswith('c')}, key=str)
            cvals = {c: rng.randrange(2, 40) for c in cfree}
            # F_p first (fast), then Q
            up, _ = is_unit(sub, sb, cvals, 32003)
            uq, outq = is_unit(sub, sb, cvals, 'Q')
            allunit &= (up and uq)
            print(f'[sub {sub}, b3-branch {bi}] c={cvals}')
            print(f'    F_p: {"J=(1) UNIT" if up else "PROPER"} | '
                  f'Q(char 0): {"J=(1) UNIT" if uq else "PROPER -- "+outq.strip()[-80:]}')
    print()
    if allunit:
        print('RESULT: J=(1) over F_p AND over Q for all (sub-case, b3-branch) at')
        print('generic b3 => the rank-3-with-linear-direction weighted tower has NO')
        print('a5 solution at generic b3 (char 0).  Combined with R-D5-GEN, the ONLY')
        print('degree-5 rank-3 leading forms live at b3 = y1^2*(c2y1+c0y2+c1y3)')
        print('(the R-D5-TAIL survivors, a5 rank-3 with isotropic v*).')
        return 0
    print('RESULT: some case not the unit ideal -- inspect (a genuine solution locus).')
    return 1


if __name__ == '__main__':
    sys.exit(main())
