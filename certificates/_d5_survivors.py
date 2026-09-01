# _d5_survivors.py -- triage the SURVIVOR CANDIDATE leaves emitted by
# _d5_close.py.  A "survivor" is a node where the a-system J did NOT reduce to
# (1) and the Rabinowitsch cone test did NOT certify det A in sqrt(J).  Before
# treating any such node as a genuine pivot-free candidate we must rule out the
# two ways it can be spurious:
#
#   (R) RANK COLLAPSE.  Each survivor sits under a long chain of constraints
#       c_i = 0, i.e. a very degenerate cubic b3.  If Hess(b3) has rank < 3 at
#       that b3, the node lies OUTSIDE the degree-5 rank-3 stratum that this
#       whole branch is defined by -- it belongs to the rank<=2 strata handled
#       separately (ledger AP0/AP1/rank-reduction).  Not a counterexample here.
#
#   (T) COMPUTE TIMEOUT.  decide()'s Rabinowitsch loop runs under `timeout 300`.
#       A timeout parses as "not the unit ideal" -> BAD>0 -> false SURV.  Re-run
#       the cone test at a longer wall clock to see if it actually certifies.
#
# Phase 1 (rank-only, PURE sympy, no WSL -- safe to run while _d5_close is still
# going): reconstruct each survivor's b3 from its tag path and print rank Hess.
# Phase 2 (--full, uses WSL Singular -- run only AFTER the closer finishes so it
# does not contend): for any rank-3 survivor, re-run the cone test cleanly.
#
#   py _d5_survivors.py                 # rank triage from the closelog
#   py _d5_survivors.py --full          # + clean cone re-test on rank-3 nodes
import os, re, sys
import sympy as sp
from _d5_close import (Y, C, E4, b3g, hess, build_JT, AV, s_poly, run_singular)


def clean_cone_test(Jp, Tp):
    """Rigorous single-shot: is det Hess3(a5) == 0 on V(J)?  i.e. does every
    y-coefficient of det A vanish on V(J)?  <=> sat(J, <T>)[1] == (1)
    (saturation J:<T>^inf = (1) iff V(J) subset of V(<T>)).  No per-generator
    Rabinowitsch loop, so no timeout blow-up.  Returns 'CONE' or 'SURV' (or
    'INCONS' if J itself is already the unit ideal -> V(J) empty -> EMPTY)."""
    pars = sorted({v for p in Jp + Tp for v in p.free_symbols
                   if str(v).startswith('c')}, key=str)
    par = ','.join(str(p) for p in pars)
    coeff = f'(0,{par})' if par else '0'
    txt = [f'ring R = {coeff},({AV}),dp;',
           'ideal J = ' + ',\n'.join(s_poly(p) for p in Jp) + ';',
           'ideal G = std(J);',
           'if (size(G)==1 && G[1]==1){ "VERDICT INCONS"; quit; }',
           'ideal T = ' + ',\n'.join(s_poly(p) for p in Tp) + ';',
           'ideal Sat = sat(J, T)[1];',
           'if (size(Sat)==1 && Sat[1]==1){ "VERDICT CONE"; } else { "VERDICT SURV"; "SATDIM", dim(std(Sat)); }',
           'quit;']
    out = run_singular('LIB "elim.lib";\n' + '\n'.join(txt) + '\n')
    if 'VERDICT INCONS' in out:
        return 'INCONS', out
    if 'VERDICT CONE' in out:
        return 'CONE', out
    if 'VERDICT SURV' in out:
        return 'SURV', out
    return 'ERROR', out

LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '_lean', 'close', 'closelog.txt')


def survivor_tags(logtext):
    tags = []
    for m in re.finditer(r'^\s*\[([^\]]+)\]\s*\*\*\* SURVIVOR CANDIDATE', logtext, re.M):
        tags.append(m.group(1))
    return tags


def constraints_from_tag(tag):
    """tag = 'root.k|f1.i1|f2.i2|...'; each fj (a factor string) is a =0
    constraint.  Return the list of sympy constraint polynomials."""
    segs = tag.split('|')[1:]                      # drop 'root.k'
    S = []
    for seg in segs:
        fac = seg.rsplit('.', 1)[0]                # strip the '.idx' suffix
        S.append(sp.sympify(fac.replace('^', '**')))
    return S


def rank_hess(b3):
    return hess(b3, Y).rank()


def main():
    full = '--full' in sys.argv
    with open(LOG, encoding='utf-8', errors='replace') as f:
        txt = f.read()
    tags = survivor_tags(txt)
    if not tags:
        print('no SURVIVOR CANDIDATE lines in closelog yet.')
        return 0
    # dedupe on the reconstructed b3 (the real object), not the tag path:
    # different constraint chains can land on the same degenerate b3.
    uniq = {}   # b3-string -> (tag, S, b3)
    for t in tags:
        S = constraints_from_tag(t)
        sols = sp.solve(list(E4) + list(S), C, dict=True)
        for sb in sols:
            b3 = sp.expand(b3g.subs(sb))
            if b3 == 0:
                continue
            uniq.setdefault(str(b3), (t, S, b3))
    print(f'{len(tags)} survivor leaves -> {len(uniq)} distinct non-zero b3.\n')

    for _, (t, S, b3) in uniq.items():
        print(f'[{t}]')
        print(f'   b3 = {b3}')
        print(f'   b3 / y1^2 = {sp.simplify(b3 / Y[0]**2)}   (rank Hess b3 = {rank_hess(b3)})')

    if not full:
        print('\nphase 1 (rank/shape info) only.  Re-run with --full AFTER '
              '_d5_close finishes for the rigorous single-shot cone test.')
        return 1

    # Phase 2: rigorous sat-based cone test (one saturation per distinct b3).
    print('\n=== phase 2: rigorous cone test  sat(J,detA)==(1) ? ===')
    genuine = []
    for _, (t, S, b3) in uniq.items():
        Jp, Tp = build_JT(b3)
        verdict, out = clean_cone_test(Jp, Tp)
        note = {'CONE': 'det Hess3(a5)=0 on V(J) -> EMPTY',
                'INCONS': 'V(J) empty -> EMPTY',
                'SURV': 'det Hess3(a5) NOT forced 0 -> GENUINE candidate',
                'ERROR': 'parse/compute error -- inspect'}.get(verdict, verdict)
        print(f'[{t}] b3=y1^2*({sp.simplify(b3/Y[0]**2)}) -> {verdict}: {note}')
        if verdict not in ('CONE', 'INCONS'):
            genuine.append((t, S, b3, verdict, out))
    print()
    if not genuine:
        print('ALL survivor families certify EMPTY (CONE/INCONS) on the '
              'rigorous test.\n=> degree-5 rank-3 pivot-free branch is EMPTY. '
              'Closure complete.')
        return 0
    print(f'*** {len(genuine)} GENUINE survivor(s) -- potential pivot-free point ***')
    for t, S, b3, v, _ in genuine:
        print('  ', t, '| b3 =', b3, '|', v)
    return 2


if __name__ == '__main__':
    sys.exit(main())
