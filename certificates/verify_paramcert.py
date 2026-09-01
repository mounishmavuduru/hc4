#!/usr/bin/env python
# verify_paramcert.py -- INDEPENDENT (Singular-free) verification of the
# parametric Nullstellensatz certificate for the degree-5 rank-3 pivot-free
# branches, plus extraction of the exceptional locus Z = V(D).
#
# Reads, per branch k, the files produced by _d5_paramcert.py + _run_cert.sh:
#   _lean/param/cert_b{k}_gens.txt   -- the J generators g_i (sympy syntax)
#   _lean/param/cert_b{k}.out        -- Singular's lift coefficients T_i in Q(c)
#
# Verifies EXACTLY in sympy that  sum_i g_i * T_i == 1  identically in Q(c)[a].
# This is a Positivstellensatz-style membership witness: it proves 1 in J over
# Q(c), hence for every c with D(c) != 0 (D = common denominator) the a-system
# E0=E1=E2=E3=0 has NO solution -> the branch has no rank-3 pivot-free
# potential with that b3.  The tail is Z = V(D).
#
# Fail-closed: every check is an assert; exit 0 <=> all held.
import os, re, sys
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(os.getcwd(), '_lean', 'param')


def load_gens(k):
    lines = open(os.path.join(BASE, f'cert_b{k}_gens.txt')).read().splitlines()
    avars = lines[0].split(':', 1)[1].strip()
    cpars = lines[1].split(':', 1)[1].strip()
    gens = [ln for ln in lines[2:] if ln.strip() and not ln.startswith('#')]
    syms = {s.name: s for s in sp.symbols(avars.replace(' ', '') + ',' + cpars.replace(' ', ''))}
    G = [sp.sympify(g, locals=syms) for g in gens]
    return G, syms


def load_T(k, syms):
    txt = open(os.path.join(BASE, f'cert_b{k}.out')).read()
    n = int(re.search(r'NGEN\s+(\d+)', txt).group(1))
    T = [sp.Integer(0)] * (n + 1)
    for m in re.finditer(r'^T\s+(\d+)\s*=\s*(.*)$', txt, re.M):
        i = int(m.group(1))
        expr = m.group(2).strip().replace('^', '**')
        T[i] = sp.sympify(expr, locals=syms)
    return n, T[1:]


def main():
    allok = True
    for k in (0, 1):
        G, syms = load_gens(k)
        n, T = load_T(k, syms)
        assert len(G) == n, f'branch {k}: gens {len(G)} != NGEN {n}'
        # combination sum g_i * T_i, as a single rational function
        S = sp.together(sum(g * t for g, t in zip(G, T)))
        num, den = sp.fraction(S)
        num = sp.expand(num)
        den = sp.expand(den)
        # S == 1  <=>  num - den == 0
        diff = sp.expand(num - den)
        ok = diff == 0
        print(f'branch {k}: sum g_i T_i == 1 ?  {ok}')
        assert ok, f'branch {k}: certificate FAILED, num-den = {diff}'
        # exceptional locus: common denominator D(c) (depends only on c)
        cset = {s for name, s in syms.items() if name.startswith('c')}
        assert den.free_symbols <= cset, f'branch {k}: denominator involves a-vars: {den.free_symbols - cset}'
        Dfac = sp.factor(den)
        print(f'  D(c) = {Dfac}')
        print(f'  => branch {k} EMPTY for all b3 with D(c) != 0; exceptional locus Z_{k} = V(D).')
        allok = allok and ok
    print()
    print('ALL CHECKS PASSED' if allok else 'FAILED')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
