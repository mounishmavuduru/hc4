# r5_rank3_cascade.py  --  key: r5  (HC_4 in DEGREE 5, rank-3 branch: the
#                                     weighted cascade below the top form)
#
# INPUT (r5_rank3_topform.py, Theorem R3-TOP).  In the rank-3 branch of
# degree 5 with NON-constant x4^2-coefficient, after affine normalisation
#      f = G_I + a4(x') + x4 b2(x') + a3(x') + x4 b1(x') + a2(x'),
#      G_I = x1^3 x2^2/2 + x1^4 x3 + x1^2 x2 x4 + x1 x4^2/2,
# with a4, b2, a3, b1, a2 ARBITRARY forms of the indicated degrees in
# x' = (x1,x2,x3)  (deg a_k = k, deg b_k = k).  [No x4^3 term: forced by
# [det]_10; no x4^2 term below G_I: the affine normalisation tau = x1 kills
# the constant part of tau; f1, f0 dropped (irrelevant to Hess).]
#
# THEOREM R3-CASCADE (proved here).  For every such f, det Hess f is NOT a
# nonzero constant.  Consequently, in the rank-3 branch of degree 5 the
# x4^2-coefficient tau is constant, i.e. e4 IS A PIVOT of f.
#
# PROOF.  Weights wt(x') = 1, wt(x4) = 2; [det Hess f]_{wt k} =: D_k, with
# D_10 = det Hess G_I = 0.  If det Hess f = c in C^x then D_k = 0 for
# 1 <= k <= 9 and D_0 = det Hess f2 = c.  We show, level by level, that the
# coefficient equations of D_9, D_8, D_7, D_6 force
#      a4 = x1^2 ( a1111 x1^2 + a1112 x1 x2 + a1113 x1 x3 + b12 x2^2 + b13 x2 x3 ),
#      b2 = x1 ( b11 x1 + b12 x2 + b13 x3 ),
#      a3 = x1 ( c111 x1^2 + c112 x1 x2 + c113 x1 x3 + (b12 x2 + b13 x3)^2/2 ),
#      b1 = d1 x1,        a2 = x1 ( e11 x1 + e12 x2 + e13 x3 ),
# and then f2 = a2 + x4 b1 = x1 * (linear) has det Hess f2 = 0 = D_0,
# contradicting c != 0.  Rigour of each level:
#   D_9: LINEAR in (a4,b2) (weight count: 3 pieces of weight 5, one of
#        weight 4); the coefficient matrix has rank 13 on 21 unknowns and the
#        displayed 8-parameter family satisfies it, so it IS the solution set.
#   D_8: after the D_9 substitution, the ideal I8 of its coefficients has
#        radical equal to the ideal J8 = (b22, b23, c222, c223, c233, c333,
#        c122 - b12^2/2 - d2, c123 - b12 b13 - d3, c133 - b13^2/2):
#        I8 c J8 (reduction to zero) and b22^2, b23^2, and the linear
#        generators lie in I8 (so J8 c rad I8; J8 is prime).
#   D_7: LINEAR in (e22,e23,e33), rank 3, solution e22 = b12 d2,
#        e23 = b12 d3 + b13 d2, e33 = b13 d3.
#   D_6: coefficients generate an ideal I6 with d2^2, d3^2 in I6 and
#        I6 c (d2,d3): rad I6 = (d2,d3).
#   Then D_0 = det Hess f2 reduces to 0.
# (D_5..D_1 impose no further condition; they are not needed.)
#
# All arithmetic exact.  Fail-closed.

import itertools
import sys
import time
import sympy as sp

x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
Y = (x1, x2, x3)
V = (x1, x2, x3, x4)
W = (1, 1, 1, 2)
T0 = time.time()


def say(msg):
    print('[%6.1fs] %s' % (time.time() - T0, msg))
    sys.stdout.flush()


def hess(g, vs):
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.diff(g, vs[i], vs[j]))


def gen(vs, d, tag):
    e = sp.Integer(0)
    ks = []
    for m in itertools.combinations_with_replacement(range(len(vs)), d):
        s = sp.Symbol(tag + ''.join(str(i + 1) for i in m))
        ks.append(s)
        t = s
        for i in m:
            t *= vs[i]
        e += t
    return sp.expand(e), ks


def wparts(pol):
    P = sp.Poly(sp.expand(pol), *V)
    out = {}
    for mon, co in P.terms():
        d = sum(w * e for w, e in zip(W, mon))
        t = co
        for v, e in zip(V, mon):
            t *= v**e
        out[d] = out.get(d, 0) + t
    return {k: sp.expand(v) for k, v in out.items()}


def coeffs(e):
    e = sp.expand(e)
    if e == 0:
        return []
    return [c for c in sp.Poly(e, *V).coeffs() if c != 0]


def S(name):
    return sp.Symbol(name)


GI = x1**3 * x2**2 / 2 + x1**4 * x3 + x1**2 * x2 * x4 + x1 * x4**2 / 2
assert sp.expand(GI - x1 * ((x4 + x1 * x2)**2 / 2 + x1**3 * x3)) == 0
a4, ak = gen(Y, 4, 'a')
b2, bk = gen(Y, 2, 'b')
a3, ck = gen(Y, 3, 'c')
b1, dk = gen(Y, 1, 'd')
a2, ek = gen(Y, 2, 'e')
allsyms = ak + bk + ck + dk + ek
f = GI + a4 + x4 * b2 + a3 + x4 * b1 + a2

say('computing det Hess f with %d generic coefficients ...' % len(allsyms))
D = sp.expand(hess(f, V).det(method='berkowitz'))
parts = wparts(D)
say('weights present: %s' % sorted(parts))
assert max(parts) <= 10 and parts.get(10, 0) == 0
assert sp.expand(parts[0] - hess(a2 + x4 * b1, V).det(method='berkowitz')) == 0
say('D_10 = 0 (= det Hess G_I) and D_0 = det Hess f2                    OK')

# --------------------------------------------------------------- level 9
say('LEVEL 9')
eq9 = coeffs(parts[9])
unk9 = ak + bk
for e in eq9:
    assert sp.Poly(e, *unk9).total_degree() == 1 and e.subs({s: 0 for s in unk9}) == 0
Mat9, _ = sp.linear_eq_to_matrix(eq9, unk9)
r9 = Mat9.rank()
sub9 = {S('a1122'): S('b12'), S('a1123'): S('b13') + 4 * S('b22'), S('a1133'): 4 * S('b23'),
        S('a1222'): S('b22'), S('a1223'): S('b23'), S('a1233'): 0, S('a1333'): 0,
        S('a2222'): 0, S('a2223'): 0, S('a2233'): 0, S('a2333'): 0, S('a3333'): 0,
        S('b33'): 0}
assert r9 == len(sub9) == 13
assert all(sp.expand(e.subs(sub9)) == 0 for e in eq9)
say('  D_9: homogeneous linear system of rank 13 in 21 unknowns; the 8-parameter')
say('       family sub9 satisfies it  =>  sub9 is the full solution set     PROVED')

# --------------------------------------------------------------- level 8
say('LEVEL 8')
w8 = sp.expand(parts[8].subs(sub9))
eq8 = coeffs(w8)
syms8 = sorted(set().union(*[e.free_symbols for e in eq8]), key=str)
say('  D_8: %d coefficient equations in %s' % (len(eq8), syms8))
J8gens = [S('b22'), S('b23'), S('c222'), S('c223'), S('c233'), S('c333'),
          S('c122') - S('b12')**2 / 2 - S('d2'),
          S('c123') - S('b12') * S('b13') - S('d3'),
          S('c133') - S('b13')**2 / 2]
ring8 = syms8
GJ8 = sp.groebner(J8gens, *ring8, order='grevlex', domain='QQ')
# I8 subset J8
for e in eq8:
    assert GJ8.reduce(e)[1] == 0
GI8 = sp.groebner(eq8, *ring8, order='grevlex', domain='QQ')
# J8 subset rad(I8): powers of the generators reduce to 0 mod I8
for g in J8gens:
    ok = False
    for k in (1, 2, 3):
        if GI8.reduce(sp.expand(g**k))[1] == 0:
            ok = True
            break
    assert ok, g
say('  D_8: I8 c J8 and every generator of J8 has a power in I8 => rad I8 = J8')
sub8 = {S('b22'): 0, S('b23'): 0, S('c222'): 0, S('c223'): 0, S('c233'): 0, S('c333'): 0,
        S('c122'): S('b12')**2 / 2 + S('d2'), S('c123'): S('b12') * S('b13') + S('d3'),
        S('c133'): S('b13')**2 / 2}
assert all(sp.expand(e.subs(sub8)) == 0 for e in eq8)
sub98 = {k: sp.expand(v.subs(sub8)) if hasattr(v, 'subs') else v for k, v in sub9.items()}
sub98.update(sub8)
say('       => b22 = b23 = 0, c2xx = c3xx = 0, c122 = b12^2/2 + d2, c123 = b12 b13 + d3,')
say('          c133 = b13^2/2                                                 PROVED')

# --------------------------------------------------------------- level 7
say('LEVEL 7')
w7 = sp.expand(parts[7].subs(sub98))
eq7 = coeffs(w7)
unk7 = [S('e22'), S('e23'), S('e33')]
syms7 = sorted(set().union(*[e.free_symbols for e in eq7]), key=str)
say('  D_7: %d coefficient equations in %s' % (len(eq7), syms7))
for e in eq7:
    assert sp.Poly(e, *unk7).total_degree() == 1
Mat7, rhs7 = sp.linear_eq_to_matrix(eq7, unk7)
assert Mat7.rank() == 3
sub7 = {S('e22'): S('b12') * S('d2'), S('e23'): S('b12') * S('d3') + S('b13') * S('d2'),
        S('e33'): S('b13') * S('d3')}
assert all(sp.expand(e.subs(sub7)) == 0 for e in eq7)
sub987 = dict(sub98)
sub987.update(sub7)
say('  D_7: linear in (e22,e23,e33) of rank 3; unique solution e22 = b12 d2,')
say('       e23 = b12 d3 + b13 d2, e33 = b13 d3                              PROVED')

# --------------------------------------------------------------- level 6
say('LEVEL 6')
w6 = sp.expand(parts[6].subs(sub987))
eq6 = coeffs(w6)
syms6 = sorted(set().union(*[e.free_symbols for e in eq6]), key=str)
say('  D_6: %d coefficient equations in %s' % (len(eq6), syms6))
assert set(syms6) <= {S('d2'), S('d3')}
GI6 = sp.groebner(eq6, S('d2'), S('d3'), order='grevlex', domain='QQ')
assert GI6.reduce(S('d2')**2)[1] == 0 and GI6.reduce(S('d3')**2)[1] == 0
for e in eq6:
    assert sp.expand(e.subs({S('d2'): 0, S('d3'): 0})) == 0
say('  D_6: d2^2, d3^2 in I6 and I6 c (d2,d3)  =>  rad I6 = (d2,d3): d2 = d3 = 0   PROVED')
final = {k: sp.expand(v.subs({S('d2'): 0, S('d3'): 0})) if hasattr(v, 'subs') else v
         for k, v in sub987.items()}
final.update({S('d2'): 0, S('d3'): 0})

# --------------------------------------------------------------- level 0
say('LEVEL 0')
w0 = sp.expand(parts[0].subs(final))
assert w0 == 0
f2fin = sp.expand((a2 + x4 * b1).subs(final))
assert sp.expand(f2fin - x1 * (S('e11') * x1 + S('e12') * x2 + S('e13') * x3 + S('d1') * x4)) == 0
say('  D_0 = det Hess f2 = 0 identically on the solution set: f2 = x1 * (linear)')
# and the surviving shape of the whole f, for the record
ffin = sp.expand(f.subs(final))
assert sp.expand(ffin - x1 * sp.cancel(ffin / x1)) == 0
say('  (indeed the entire surviving f is divisible by x1)')
# (levels 5..1 are automatically satisfied on the solution set; recorded, not needed)
for k in range(5, 0, -1):
    assert sp.expand(parts[k].subs(final)) == 0
say('  D_5..D_1 vanish identically on the solution set (not needed)')

print()
print('THEOREM R3-CASCADE: no f = G_I + (lower weight) has det Hess f in C^x.')
print('Hence in the rank-3 branch of degree 5 the x4^2-coefficient is constant:')
print('e4 is a pivot, and Theorem A / the pivot trichotomy + Moh settle HC_4 there.')
print('r5_rank3_cascade: ALL PASS')
