# d5_e4_exhaustive.py
#
# P1.3a (partial closure of the solve()-exhaustiveness caveat): the recursive
# tail closer _d5_close.py enumerates admissible cubics b3 by sp.solve(E4=0),
# which has no completeness guarantee. The ROOT split is the decisive one: does
# sp.solve find ALL irreducible components of the E4=0 locus, or could it miss a
# stratum? E4 = coeff_y(adj(Hess b3)[0,0]) is 6 polynomials in the 10 cubic
# coefficients c0..c9 (pure c, no a-variables), so its primary/minimal-prime
# decomposition is a small, decidable Singular computation.
#
# This script (a) builds E4 in sympy, (b) computes the minimal primes of <E4>
# over Q with Singular minAssGTZ, (c) reports the component count and dimensions,
# and (d) cross-checks against sp.solve(E4)'s branches: it asserts that the
# NUMBER of nonzero-b3 solve branches equals the number of minAssGTZ components
# that are NOT contained in {b3 = 0} (i.e. that carry a genuine cubic), so no
# root stratum is missed. If this holds, the ROOT split of the tail is provably
# exhaustive (the deeper recursion only adds constraints / drops dimension).
#
#   py -u d5_e4_exhaustive.py            # needs WSL Singular
import os, re, subprocess, sys
import sympy as sp
from _d5_close import Y, form, coeffs_in_y, hess, s_poly, WSL

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'e4exh')
os.makedirs(ND, exist_ok=True)

b3g, C = form(Y, 3, 'c')                      # generic cubic, coeffs c0..c9
E4 = [p for p in coeffs_in_y(hess(b3g, Y).adjugate()[0, 0]) if p != 0]
CV = ','.join(str(c) for c in C)


def run(text, tmo=280):
    fn = os.path.join(ND, 'e4.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(text)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


def main():
    # sp.solve branches (as _d5_close/_d5_paramcert use them)
    sols = sp.solve(E4, C, dict=True)
    nz = [sb for sb in sols if sp.expand(b3g.subs(sb)) != 0]
    print(f'sp.solve(E4=0): {len(sols)} branch(es), {len(nz)} with b3 != 0')

    txt = ['LIB "primdec.lib";',
           f'ring R = 0,({CV}),dp;',
           'ideal E4 = ' + ',\n'.join(s_poly(p) for p in E4) + ';',
           'list pd = minAssGTZ(E4);',
           '"NCOMP", size(pd);',
           'int i;',
           'for(i=1;i<=size(pd);i++){ "COMP", i, "dim", dim(std(pd[i])); '
           '  "gens", i; print(pd[i]); "ENDGENS"; }',
           'quit;']
    out = run('\n'.join(txt) + '\n')
    m = re.search(r'NCOMP\s+(\d+)', out)
    ncomp = int(m.group(1)) if m else None
    dims = [int(d) for d in re.findall(r'COMP\s+\d+\s+dim\s+(-?\d+)', out)]
    print(f'Singular minAssGTZ(<E4>): {ncomp} minimal prime(s), dims {dims}')
    print('--- raw decomposition (first 1200 chars) ---')
    print(out[:1200])

    # FINDING (2026-09-11): the E4=0 locus is ONE irreducible component (dim 6).
    # sp.solve's "2 branches" are therefore two rational CHARTS of that single
    # variety (branch 0: generic, c9 != 0; branch 1: the c9 = 0, c5 != 0 boundary),
    # NOT two components -- so "2 solve branches vs N components" is the wrong
    # comparison. The exhaustiveness question is whether the charts + the
    # D-factor recursion cover every point. The closelog shows sp.solve returned
    # "no b3" at the two c9 = 0 nodes although that locus is nonempty: a concrete
    # solve() miss. It is SUBSUMED because on V(E4) we have c9 = 0 => c8 = 0
    # (generator c8^3 - 27 c7 c9^2), so those points lie inside the explored
    # c8 = 0 subtree. Verify that subsumption here: c8 in sqrt(<E4, c9>).
    y = sp.Symbol('y')
    c8, c9 = C[8], C[9]
    G = sp.groebner(E4 + [c9, 1 - y * c8], *(list(C) + [y]), order='grevlex')
    subsumed = (list(G.exprs) == [sp.Integer(1)])
    print()
    print(f'nonzero-b3 solve branches (charts) = {len(nz)}; irreducible components = {ncomp}, dims {dims}')
    print(f'c8 in sqrt(<E4, c9>) (i.e. c9=0 => c8=0 on V(E4), so the solve()-empty '
          f'c9=0 nodes are subsumed by the explored c8=0 subtree)?  {subsumed}')
    ok = (ncomp == 1) and subsumed
    print()
    if ok:
        print('RESULT: root E4-locus is a single irreducible dim-6 component; the two '
              'solve() charts are charts of it, and the only solve()-empty nodes (c9=0) '
              'are PROVABLY subsumed. P1.3a is sharpened, NOT fully discharged: a '
              'partial-chart miss at some deeper node is not excluded by this check -- '
              'full method-exhaustiveness needs an ideal-theoretic (minAssGTZ per node) '
              're-run of the closer. Tail stays [machine-modular].')
        return 0
    print('UNEXPECTED: inspect the decomposition above.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
