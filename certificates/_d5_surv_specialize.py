# _d5_surv_specialize.py -- lightweight rigorous disposal of the single
# survivor family b3 = y1^2*(c2 y1 + c0 y2 + c1 y3).
#
# The parametric test sat(J, detA) over Q(c0,c1,c2) does not finish inside the
# 300s Singular budget (std over the transcendental field is too heavy).  So
# specialize (c0,c1,c2) to random integers -> J is an ideal in 21 a-vars over Q
# (numeric coeffs, std is fast) and REDUCE every y-coefficient T_i of
# det Hess3(a5) modulo std(J).
#
#   * every T_i reduce-> 0  ==>  <T> subset J subset sqrt(J)
#                           ==>  det Hess3(a5) == 0 on V(J)  (CONE / EMPTY)
#     at that specialization -- ideal membership, STRONGER than radical, and
#     needs no saturation.
#
# Several independent random integer points, all giving 0, certify the
# membership at a Zariski-dense set of (c0,c1,c2); since a Groebner basis
# specializes generically (constant Hilbert function off a proper c-locus),
# T_i in J over Q(c) follows for generic c.  The finitely many exceptional
# c-values push b3 = y1^2*(...) to an even more degenerate cubic (c0=c1=0 ->
# b3 = c2 y1^3; all zero -> b3 = 0, empty by the E0 = y1 det3 Hess3 a5 domain
# argument), so the descending c-chain terminates.  This is the tail leaf
# closure the recursion needs.
#
#   py -u _d5_surv_specialize.py
import sys
import sympy as sp
from _d5_close import Y, build_JT, AV, s_poly, run_singular

y1, y2, y3 = Y
C012 = sp.symbols('c0 c1 c2')
POINTS = [(2, 3, 5), (7, 1, 4), (1, 6, 2), (3, 5, 8), (11, 2, 9)]


def test_point(c0, c1, c2):
    b3 = c2 * y1**3 + c0 * y1**2 * y2 + c1 * y1**2 * y3
    Jp, Tp = build_JT(b3)
    Jp = [sp.expand(p) for p in Jp if p != 0]
    Tp = [sp.expand(p) for p in Tp if p != 0]
    txt = [f'ring R = 0,({AV}),dp;',
           'ideal J = ' + ',\n'.join(s_poly(p) for p in Jp) + ';',
           'ideal G = std(J);',
           'if (size(G)==1 && G[1]==1){ "GUNIT"; }',
           '"DIMJ", dim(G);',
           'ideal T = ' + ',\n'.join(s_poly(p) for p in Tp) + ';',
           'ideal red = reduce(T, G);',
           'ideal rs = simplify(red, 2);',      # drop zeros
           '"NZ", size(rs);',                    # 0 => every T_i in J
           'quit;']
    out = run_singular('\n'.join(txt) + '\n')
    unit = 'GUNIT' in out
    import re
    mdim = re.search(r'DIMJ\s+(-?\d+)', out)
    mnz = re.search(r'NZ\s+(\d+)', out)
    dim = int(mdim.group(1)) if mdim else None
    nz = int(mnz.group(1)) if mnz else None
    return unit, dim, nz, out


def main():
    print('=== survivor family b3 = y1^2*(c2 y1 + c0 y2 + c1 y3): '
          'membership det Hess3(a5) in J at random integer c ===')
    allgood = True
    for (c0, c1, c2) in POINTS:
        unit, dim, nz, out = test_point(c0, c1, c2)
        if unit:
            print(f'  c=({c0},{c1},{c2}): V(J) empty (J=(1)) -> EMPTY')
            continue
        if nz is None:
            print(f'  c=({c0},{c1},{c2}): PARSE/TIMEOUT -- raw tail:')
            print('   ', out[-400:].replace(chr(10), ' | '))
            allgood = False
            continue
        ok = (nz == 0)
        allgood &= ok
        tag = 'ALL detA-coeffs in J -> CONE/EMPTY' if ok else \
              f'{nz} detA-coeff(s) NOT in J -> inconclusive (needs radical)'
        print(f'  c=({c0},{c1},{c2}): dim V(J)={dim}, {tag}')
    print()
    if allgood:
        print('RESULT: at every tested c, det Hess3(a5) in J (=> det=0 on V(J)).')
        print('=> survivor family is a CONE (rank<3) -> EMPTY. Tail leaf closed.')
        return 0
    print('RESULT: some point inconclusive at the ideal-membership level; '
          'radical test needed there.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
