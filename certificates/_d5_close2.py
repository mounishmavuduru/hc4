# _d5_close2.py -- authoritative re-run of the degree-5 rank-3 pivot-free
# closure with a FAST, RIGOROUS cone test.
#
# _d5_close.py's decide() certifies "det Hess3(a5) == 0 on V(J)" with a
# Rabinowitsch loop that runs one std() per y-coefficient of det A (55 of them)
# under a 300s ceiling.  That loop is slow and, worse, a timeout on any single
# generator is parsed as "not a cone" -> a FALSE survivor.  Every survivor the
# first run reported is the single degenerate family b3 = y1^2*(c2 y1+c0 y2+c1 y3)
# (verified in _d5_survivors.py), which is exactly where that loop times out.
#
# Here decide() instead does ONE saturation:
#     sat(J, <det A coeffs>)[1] == (1)  <=>  V(J) subset V(det A)
#                                       <=>  det Hess3(a5) == 0 on V(J)  (CONE).
# One Groebner/saturation, no per-generator loop, no timeout blow-up.
#
# Writes to _lean/close2 (its own node file + log) so it can run CONCURRENTLY
# with a still-running _d5_close.py without clobbering its node.sing.
#
#   py -u _d5_close2.py
import os, re, sys
import sympy as sp
import _d5_close as base
from _d5_close import AV, s_poly, as_, C, run_singular

base.NODEDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'close2')
os.makedirs(base.NODEDIR, exist_ok=True)


def decide(Jp, Tp):
    """Fast rigorous replacement. Returns (verdict, data):
       EMPTY1 (J=(1), off V(D))  |  CONE (detA==0 on V(J))  |  SURV."""
    pars = sorted({v for p in Jp + Tp for v in p.free_symbols
                   if str(v).startswith('c')}, key=str)
    par = ','.join(str(p) for p in pars)
    coeff = f'(0,{par})' if par else '0'
    txt = [f'ring R = {coeff},({AV}),dp;',
           'ideal J = ' + ',\n'.join(s_poly(p) for p in Jp) + ';',
           'ideal G = std(J);',
           'if (size(G)==1 && G[1]==1){ "V1"; matrix Tc = lift(J, ideal(1)); int i;'
           ' for(i=1;i<=size(J);i++){ "L",i,"=",Tc[i,1]; } "ENDL"; quit; }',
           'ideal T = ' + ',\n'.join(s_poly(p) for p in Tp) + ';',
           'ideal Sat = sat(J, T)[1];',
           'if (size(Sat)==1 && Sat[1]==1){ "VCONE"; } else { "VSURV"; }',
           'quit;']
    out = run_singular('LIB "elim.lib";\n' + '\n'.join(txt) + '\n')
    if 'V1' in out:
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
    if 'VCONE' in out:
        return 'CONE', None
    if 'VSURV' in out:
        return 'SURV', out
    return 'SURV', out          # fail-closed: unparsed -> treat as survivor


base.decide = decide            # close() reads decide as a module global

if __name__ == '__main__':
    print('=== authoritative closure (fast sat cone test): deg-5 rank-3 pivot-free ===')
    base.close([], 0, 'root')
    print()
    print(f'nodes visited: {base.NODES}')
    print(f'survivor candidates: {len(base.SURV)}')
    for t, S, _ in base.SURV:
        print('  SURVIVOR at', t, 'constraints', S)
    ok = len(base.SURV) == 0 and all(v in ('CONE', 'B3ZERO', 'NOSOL')
                                     for _, v, _ in base.LEAVES)
    print('ALL CHECKS PASSED' if ok else 'INCOMPLETE / see leaves')
    for t, v, S in base.LEAVES:
        if v not in ('CONE', 'B3ZERO', 'NOSOL'):
            print('  LEAF', t, v, S)
    sys.exit(0 if ok else 1)
