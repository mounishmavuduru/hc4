# PE  --  PIVOT EXISTENCE, part 1: the tool and the baseline.
#
# Definition.  f in C[x_1..x_n] admits a PIVOT if there is v != 0 with
# D_v^2 f  identically constant.   Since D_v^2 f = v^T (Hess f) v and
# f = sum_k f_k (homogeneous parts), D_v^2 f is constant  <=>  D_v^2 f_k = 0
# for every k >= 3.  Writing N(x) := Hess f - Hess f(0) = sum_{alpha != 0}
# x^alpha N_alpha  (N_alpha constant symmetric matrices), this reads
#
#     PIVOT CONE  P(f) := { v : v^T N_alpha v = 0 for all alpha != 0 }
#                       = common zero cone of the quadratic forms
#                         q_alpha(v) = v^T N_alpha v.
#
# So "f admits a pivot" <=> P(f) != {0} <=> the linear system of quadrics
# {q_alpha} has a base point in P^{n-1}.
#
# This file provides
#   pivot_quadrics(f, xs, vs)   -> the list of quadrics q_alpha
#   variation_space_dim(f, xs)  -> dim span{N_alpha} in Sym_n(C)
#   cone_is_trivial(quads, vs)  -> exact Groebner decision  P(f) = {0}?
#   witness(quads, vs)          -> an explicit nonzero v when one is found
# and runs the baseline: every 4-variable constant-Hessian potential known to
# this project has a nonempty pivot cone, while the 5-variable Meng-Yang
# potential has an empty one (independent re-derivation of D1 with this tool).
#
# Exact arithmetic only.  Fail-closed: any assertion failure exits nonzero.

import sympy as sp
from itertools import combinations

# ----------------------------------------------------------------------
# core


def hess(f, xs):
    return sp.Matrix(len(xs), len(xs), lambda i, j: sp.diff(f, xs[i], xs[j]))


def pivot_quadrics(f, xs, vs):
    """Quadratic forms in vs whose common zero cone is the pivot cone of f."""
    D2 = sp.expand(sum(vs[i] * vs[j] * sp.diff(f, xs[i], xs[j])
                       for i in range(len(xs)) for j in range(len(xs))))
    P = sp.Poly(D2, *xs)
    quads = []
    for monom, coeff in P.terms():
        if sum(monom) > 0:
            c = sp.expand(coeff)
            if c != 0:
                quads.append(c)
    # sanity: every generator is a homogeneous quadratic in v
    for q in quads:
        pq = sp.Poly(q, *vs)
        assert pq.is_homogeneous and pq.total_degree() == 2, q
    return quads


def variation_space_dim(f, xs):
    """dim of W = span{ N_alpha } inside Sym_n(C), N = Hess f - Hess f(0)."""
    n = len(xs)
    H = hess(f, xs)
    N = sp.Matrix(n, n, lambda i, j: sp.expand(H[i, j] - H[i, j].subs(
        {x: 0 for x in xs})))
    mons = set()
    for i in range(n):
        for j in range(n):
            if N[i, j] != 0:
                for m in sp.Poly(N[i, j], *xs).monoms():
                    if sum(m) > 0:
                        mons.add(m)
    mons = sorted(mons)
    rows = []
    for m in mons:
        row = []
        for i in range(n):
            for j in range(i, n):
                row.append(sp.Poly(N[i, j], *xs).coeff_monomial(m)
                           if N[i, j] != 0 else 0)
        rows.append(row)
    if not rows:
        return 0
    return sp.Matrix(rows).rank()


def cone_is_trivial(quads, vs, order='grevlex'):
    """True iff the common zero cone of `quads` is {0} (exact Groebner test).

    All generators are homogeneous, so the zero set is a cone; a cone is {0}
    iff it is a finite set iff the leading-term ideal contains a pure power of
    each variable."""
    if not quads:
        return False
    G = sp.groebner(quads, *vs, order=order)
    if list(G.exprs) == [sp.Integer(1)]:
        return True          # (cannot happen for homogeneous quadrics)
    pure = {i: False for i in range(len(vs))}
    for g in G.exprs:
        m = sp.Poly(g, *vs).monoms(order=order)[0]
        nz = [i for i, e in enumerate(m) if e > 0]
        if len(nz) == 1:
            pure[nz[0]] = True
    return all(pure.values())


def witness(quads, vs, extra_tries=()):
    """Try to exhibit an explicit nonzero common zero (a pivot direction)."""
    n = len(vs)
    cands = [sp.Matrix([1 if i == k else 0 for i in range(n)])
             for k in range(n)]
    cands += [sp.Matrix(list(t)) for t in extra_tries]
    # also all 0/+-1/+-I vectors with small support
    from itertools import product
    vals = [0, 1, -1, sp.I, -sp.I]
    for supp in range(1, min(n, 3) + 1):
        for pos in combinations(range(n), supp):
            for combo in product(vals[1:], repeat=supp):
                vec = [0] * n
                for p, c in zip(pos, combo):
                    vec[p] = c
                cands.append(sp.Matrix(vec))
    for v in cands:
        if all(x == 0 for x in v):
            continue
        sub = {vs[i]: v[i] for i in range(n)}
        if all(sp.simplify(q.subs(sub)) == 0 for q in quads):
            return v
    return None


def report(name, f, xs, vs, expect_trivial=None, quiet=False):
    H = hess(f, xs)
    d = sp.simplify(H.det(method='berkowitz'))
    quads = pivot_quadrics(f, xs, vs)
    triv = cone_is_trivial(quads, vs) if quads else False
    w = witness(quads, vs) if not triv else None
    if not quiet:
        print(f"  {name}")
        print(f"    det Hess = {d}   (#quadrics = {len(quads)}, "
              f"dim W = {variation_space_dim(f, xs)})")
        print(f"    pivot cone trivial? {triv}"
              + (f"   pivot witness v = {list(w)}" if w is not None else ""))
    if expect_trivial is not None:
        assert triv == expect_trivial, (name, triv, expect_trivial)
    return d, quads, triv, w


# ----------------------------------------------------------------------
if __name__ == "__main__":
    x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
    xs4 = [x1, x2, x3, x4]
    vs4 = list(sp.symbols('v1:5'))

    print("=" * 70)
    print("BASELINE A: four-variable constant-Hessian potentials")
    print("=" * 70)

    # --- de Bondt's dillen4 example (a Meng doubling of a planar Keller map)
    b = x1 + x2**2
    e = x2 + (x1 + x2**2)**2
    f_dillen4 = b * x3 + e * x4
    d, _, triv, _ = report("dillen4 = (x1+x2^2)x3 + (x2+(x1+x2^2)^2)x4",
                           f_dillen4, xs4, vs4, expect_trivial=False)
    assert d.is_number and d != 0

    # --- general Meng doubling with an inert planar potential
    a = x1**5 * x2 - 3 * x1 * x2**3 + x1**2
    f_doub = a + b * x3 + e * x4
    d, _, _, _ = report("doubling + inert a(x1,x2)", f_doub, xs4, vs4,
                        expect_trivial=False)
    assert d.is_number and d != 0

    # --- the project's AP1 family  a + x3(-k x1 + beta(xi)) + x4 xi
    xi = x2 + x1**2 / 2
    kap = sp.Integer(2)
    beta = xi**3 - 5 * xi
    f_ap1 = (x1**4 * x2 + x2**3) + x3 * (-kap * x1 + beta) + x4 * xi
    d, _, _, _ = report("AP1 classified family", f_ap1, xs4, vs4,
                        expect_trivial=False)
    assert sp.simplify(d - kap**2) == 0

    # --- a quadratic-pivot potential (Theorem A(i) class):
    #     f = e0(x') + x4 e1(x') + x4^2/2  with the forced pencil identity.
    #     take e1 = x1^2/2 + x3 and e0 = e1^2/2 + (planar const-Hess in x1,x2)
    #     ~e0 = x1 x3 + x2^2/2 + x1^5 has Hess_3 with cofactor C_11 = 0 and
    #     det_3 = -1, so det_3(Hess_3 ~e0 + s K) = -1 for every s (K = E_11).
    e1 = x1**2 / 2 + x3
    e0 = e1**2 / 2 + (x1 * x3 + x2**2 / 2 + x1**5)
    f_q2 = e0 + x4 * e1 + x4**2 / 2
    d, _, _, _ = report("quadratic-pivot instance", f_q2, xs4, vs4,
                        expect_trivial=False)
    assert d.is_number and d != 0

    # --- pure top form:  f = q + f_5  with Hess f_5 pointwise nilpotent
    #     (cone in the x4 direction => automatically a pivot)
    q = (x1**2 + x2**2 + x3**2 + x4**2) / 2
    f_cone5 = q + (x1 + sp.I * x2)**5          # isotropic linear form
    d, _, _, _ = report("q + (x1+i x2)^5", f_cone5, xs4, vs4,
                        expect_trivial=False)
    assert sp.simplify(d - 1) == 0

    print()
    print("=" * 70)
    print("BASELINE B: negative controls")
    print("=" * 70)

    # a generic-looking quartic (NOT constant Hessian): cone must be empty
    f_gen = q + x1**3 + x2**3 * x3 + x1 * x2 * x4 + x3**4 + x4**3 * x1
    quads = pivot_quadrics(f_gen, xs4, vs4)
    triv = cone_is_trivial(quads, vs4)
    print(f"  random non-Keller potential: pivot cone trivial? {triv} "
          f"(#quadrics {len(quads)})")
    assert triv, "the test must be able to certify an EMPTY pivot cone"

    # 5-variable Meng-Yang counterexample: empty pivot cone (re-derivation)
    y1, y2, y3 = sp.symbols('y1 y2 y3')
    w5 = [x1, x2, y1, y2, y3]
    vs5 = list(sp.symbols('u1:6'))
    u = 1 + x1 * x2
    A = y1 * u**3 + 3 * x1 * y2 * u**2 - x1**3 * y3
    B = (y1 * x2**2 * u * (4 + 3 * x1 * x2)
         + y2 * (x2 + 3 * x1 * x2**2 * (4 + 3 * x1 * x2))
         + y3 * (2 * x1 - 3 * x1**2 * x2))
    lam, mu = sp.Integer(2), sp.Integer(13)
    Psi = B + (lam / 2) * A**2 + mu * A
    quads5 = pivot_quadrics(Psi, w5, vs5)
    triv5 = cone_is_trivial(quads5, vs5)
    print(f"  Meng-Yang Psi (5 vars): pivot cone trivial? {triv5} "
          f"(#quadrics {len(quads5)}, dim W = {variation_space_dim(Psi, w5)})")
    assert triv5, "Meng-Yang must have an EMPTY pivot cone (lemma D1)"

    print()
    print("ALL BASELINE CHECKS PASS.")
    print("Every 4-variable constant-Hessian potential in the project's")
    print("catalogue admits a pivot; the 5-variable Meng-Yang potential does")
    print("not.  The tool certifies both directions.")
