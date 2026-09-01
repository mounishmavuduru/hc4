#!/usr/bin/env python
# verify_paramcert.py -- INDEPENDENT (Singular-free) verification of the
# parametric Nullstellensatz certificate for the degree-5 rank-3 pivot-free
# branches, plus extraction of the exceptional locus Z = V(D).
#
# Reads, per branch k, the certificate files (tracked under certs/, or the
# scratch _lean/param/ if certs/ is absent):
#   cert_b{k}_gens.txt   -- the J generators g_i (sympy syntax, one per line)
#   cert_b{k}.out        -- Singular's lift cofactors T_i in Q(c)
#
# Two independent checks, both fail-closed:
#   (1) GENERATOR RE-DERIVATION: rebuild J = <coeff_y(E0..E3)> for branch k
#       from the model here in sympy (NOT trusting the gens file) and assert the
#       loaded g_i equal the re-derived generators as a multiset up to nonzero
#       rational scalar. This proves the file's generators really are the
#       E-system's generators (closes the "trust the gens file" gap).
#   (2) MEMBERSHIP: verify EXACTLY in sympy that  sum_i g_i * T_i == 1
#       identically in Q(c)[a]. Together with (1) this proves 1 in J over Q(c),
#       hence for every c with D(c) != 0 (D = common denominator) the a-system
#       E0=E1=E2=E3=0 has NO solution -> the branch has no rank-3 pivot-free
#       potential with that b3. The tail is Z = V(D).
#
# Fail-closed: every check is an assert; exit 0 <=> all held.
import itertools, os, re, sys
from collections import Counter
import sympy as sp

os.chdir(os.path.dirname(os.path.abspath(__file__)))
# prefer the tracked witnesses; fall back to the scratch dir if regenerated
BASE = 'certs' if os.path.exists(os.path.join('certs', 'cert_b0_gens.txt')) \
       else os.path.join('_lean', 'param')

y1, y2, y3, s = sp.symbols('y1 y2 y3 s')
Y = (y1, y2, y3)


# --- model (replicated inline so the verifier depends on neither Singular nor
#     the other project scripts) ---------------------------------------------
def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def grad(g, vs):
    return sp.Matrix([sp.diff(g, v) for v in vs])


def form(vs, deg, tag):
    mons = sorted({sp.prod(t) for t in itertools.combinations_with_replacement(vs, deg)}, key=str)
    cs = sp.symbols(f'{tag}0:{len(mons)}')
    return sp.expand(sum(c * m for c, m in zip(cs, mons))), list(cs)


def coeffs_in_y(expr):
    e = sp.expand(expr)
    return [] if e == 0 else [sp.expand(c) for c in sp.Poly(e, *Y).coeffs()]


def E_components(a5, b3):
    A, B, q = hess(a5, Y), hess(b3, Y), grad(b3, Y)
    adjA, adjB = A.adjugate(), B.adjugate()
    L = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff((A + s * B).adjugate()[i, j], s).subs(s, 0)))
    E2 = y1 * sp.trace(A * adjB) - (q.T * adjB * q)[0, 0] - 2 * (L * q)[0] - adjA[0, 0]
    E1 = y1 * sp.trace(adjA * B) - (q.T * L * q)[0, 0] - 2 * (adjA * q)[0]
    E0 = y1 * A.det(method='berkowitz') - (q.T * adjA * q)[0, 0]
    E3 = y1 * B.det(method='berkowitz') - 2 * (adjB * q)[0] - L[0, 0]
    return [sp.expand(e) for e in (E3, E2, E1, E0)], A


# reproduce the two E4=0 branches exactly as _d5_paramcert.py does
_b3g, _cs = form(Y, 3, 'c')
_sols_b = sp.solve(coeffs_in_y(hess(_b3g, Y).adjugate()[0, 0]), _cs, dict=True)
_a5, _as = form(Y, 5, 'a')


def rederive_gens(k):
    """rebuild J = <coeff_y(E0..E3)> for branch k, numerators, in J order."""
    b3k = sp.expand(_b3g.subs(_sols_b[k]))
    Es, _A = E_components(_a5, b3k)
    Jsym = []
    for e in Es:
        Jsym += coeffs_in_y(e)
    return [sp.expand(sp.fraction(sp.together(p))[0]) for p in Jsym if p != 0]


def canon(expr, allsyms):
    """canonical form up to nonzero rational scalar: primitive part, leading
    coefficient made positive; returned as a frozenset of (monomial, coeff)."""
    P = sp.Poly(sp.expand(expr), *allsyms)
    if P.is_zero:
        return frozenset()
    P = P.primitive()[1]                 # strip integer content
    if P.LC() < 0:
        P = -P
    return frozenset(P.as_dict().items())


def multiset_up_to_scalar(polys, allsyms):
    return Counter(canon(p, allsyms) for p in polys)


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
    print(f'reading certificates from {BASE}/')
    allok = True
    for k in (0, 1):
        G, syms = load_gens(k)
        n, T = load_T(k, syms)
        assert len(G) == n, f'branch {k}: gens {len(G)} != NGEN {n}'

        # (1) generator re-derivation: loaded gens == E-system gens (up to scalar)
        Gref = rederive_gens(k)
        allsyms = sorted(set().union(*[g.free_symbols for g in G + Gref]),
                         key=lambda t: (t.name[0], int(t.name[1:])))
        got = multiset_up_to_scalar(G, allsyms)
        want = multiset_up_to_scalar(Gref, allsyms)
        gens_ok = (got == want)
        print(f'branch {k}: loaded {len(G)} gens == re-derived coeff_y(E0..E3) '
              f'(up to scalar)?  {gens_ok}')
        assert gens_ok, (f'branch {k}: loaded generators are NOT the E-system '
                         f'generators (missing {want - got}, extra {got - want})')

        # (2) membership: sum g_i T_i == 1
        S = sp.together(sum(g * t for g, t in zip(G, T)))
        num, den = sp.fraction(S)
        diff = sp.expand(sp.expand(num) - sp.expand(den))
        mem_ok = diff == 0
        print(f'branch {k}: sum g_i T_i == 1 ?  {mem_ok}')
        assert mem_ok, f'branch {k}: certificate FAILED, num-den = {diff}'

        # exceptional locus: common denominator D(c) (depends only on c)
        den = sp.expand(den)
        cset = {sym for name, sym in syms.items() if name.startswith('c')}
        assert den.free_symbols <= cset, \
            f'branch {k}: denominator involves a-vars: {den.free_symbols - cset}'
        print(f'  D(c) = {sp.factor(den)}')
        print(f'  => branch {k} EMPTY for all b3 with D(c) != 0; exceptional locus Z_{k} = V(D).')
        allok = allok and gens_ok and mem_ok
    print()
    print('ALL CHECKS PASSED' if allok else 'FAILED')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
