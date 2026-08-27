# RED TEAM: independent re-parse of Gao's F4/F5 from the LaTeX source, plus a
# hostile audit of the parser used by gradient_equivalence_test.py.
#
# Method: a STRICT tokenizing scanner that must consume the whole term string.
# Grammar per term:  [+-]?  (\tfrac{a}{b} | integer)?  (var (^d | ^{dd})?)*
# Anything it cannot consume is a hard error (no silent fallback).
# Then compare, monomial by monomial, against the production parser's output.

import re
import sys
import sympy as sp

x, y, z, t = sp.symbols('x y z t')
VARS = {'x': x, 'y': y, 'z': z, 't': t}
TEXPATH = r'C:\Users\mouni\hc4\lit\gao\Jacobian_CE.tex'
TEX = open(TEXPATH, encoding='utf-8').read()

FAIL = []


def check(cond, msg):
    print(('  OK   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAIL.append(msg)


# ---------------------------------------------------------------- strict parser
TERM_RE = re.compile(r'''
    (?P<sign>[+-])?
    (?:\\tfrac\{(?P<num>-?\d+)\}\{(?P<den>\d+)\}|(?P<int>\d+))?
    (?P<vars>(?:[xyzt](?:\^(?:\d|\{\d+\}))?)*)
    $''', re.X)

VAR_RE = re.compile(r'([xyzt])(?:\^(?:(\d)|\{(\d+)\}))?')


def strict_parse_component(src):
    """src: raw LaTeX of one component RHS. Returns sympy expr."""
    s = src
    s = re.sub(r'\\\\(\[\d+pt\])?', '', s)      # line breaks
    s = s.replace('&', '').replace('\n', '')
    s = re.sub(r'\s+', '', s)
    if not s:
        raise ValueError('empty component')
    # split into terms keeping the sign
    pieces = re.findall(r'[+-]?[^+-]+', s)
    # careful: '-' inside \tfrac{-7}{9} would split wrongly, so re-glue
    terms = []
    buf = ''
    for p in pieces:
        buf += p
        # balanced braces?
        if buf.count('{') == buf.count('}'):
            terms.append(buf)
            buf = ''
    if buf:
        raise ValueError('unbalanced braces near %r' % buf[:60])
    total = sp.Integer(0)
    for term in terms:
        m = TERM_RE.match(term)
        if m is None:
            raise ValueError('UNPARSED TERM: %r' % term)
        sign = -1 if m.group('sign') == '-' else 1
        if m.group('num') is not None:
            coeff = sp.Rational(int(m.group('num')), int(m.group('den')))
        elif m.group('int') is not None:
            coeff = sp.Integer(int(m.group('int')))
        else:
            coeff = sp.Integer(1)
        mono = sp.Integer(1)
        vs = m.group('vars')
        # verify VAR_RE tiles vs exactly
        consumed = 0
        for vm in VAR_RE.finditer(vs):
            if vm.start() != consumed:
                raise ValueError('gap in var block %r' % vs)
            consumed = vm.end()
            e = 1
            if vm.group(2) is not None:
                e = int(vm.group(2))
            elif vm.group(3) is not None:
                e = int(vm.group(3))
            mono *= VARS[vm.group(1)] ** e
        if consumed != len(vs):
            raise ValueError('trailing junk in var block %r' % vs)
        if coeff == 1 and mono == 1 and m.group('int') is None and m.group('num') is None and vs == '':
            raise ValueError('empty term %r' % term)
        total += sign * coeff * mono
    return sp.expand(total)


def blocks(marker):
    i = TEX.index(marker)
    j = TEX.index(r'\begin{align*}', i)
    k = TEX.index(r'\end{align*}', j)
    return TEX[j + len(r'\begin{align*}'):k]


def strict_components(block, d):
    parts = re.split(r'F_\{%d,\d\}\s*=\{\}&' % d, block)
    assert parts[0].strip() == '', parts[0][:80]
    return [strict_parse_component(p) for p in parts[1:]]


print('=' * 74)
print('A. Hostile unit tests of the PRODUCTION parser (gradient_equivalence_test)')
print('=' * 74)
sys.path.insert(0, r'C:\Users\mouni\hc4\certificates')
prod_src = open(r'C:\Users\mouni\hc4\certificates\gradient_equivalence_test.py',
                encoding='utf-8').read()
# lift the exact transformation chain out of the production parser
def prod_transform(s):
    s = re.sub(r'\\\\(\[\d+pt\])?', '', s)
    s = s.replace('&', '').replace('\n', '').replace(' ', '')
    s = re.sub(r'\\tfrac\{(-?\d+)\}\{(\d+)\}', r'(\1/\2)*', s)
    s = re.sub(r'\^\{(\d+)\}', r'**\1', s)
    s = re.sub(r'\^(\d)', r'**\1', s)
    s = re.sub(r'(?<![\*/\d])(\d+)(?=[xyzt])', r'\1*', s)
    s = re.sub(r'(?<=[xyzt\)])(?=[xyzt])', '*', s)
    s = re.sub(r'(\*\*\d+)(?=[xyzt])', r'\1*', s)
    return s


CRAFTED = [
    # (latex, correct sympy expr or None meaning "must raise")
    (r'x^{10}',      x**10),
    (r'x^{12}yz^4t^4', x**12*y*z**4*t**4),
    (r'x^1 0',       None),      # LaTeX: x^1 * 0 = 0  -> parser gives x**10 (WRONG)
    (r'x^2 3',       None),      # LaTeX: x^2 * 3      -> ?
    (r'\tfrac{22}{9}x^6yz^2t^2', sp.Rational(22, 9)*x**6*y*z**2*t**2),
    (r'-\tfrac{7}{9}x^3t^2', -sp.Rational(7, 9)*x**3*t**2),
    (r'14x^4yz^2t',  14*x**4*y*z**2*t),
    (r'2xy',         2*x*y),
    (r'x^{2}y^{3}',  x**2*y**3),
]
for tex_s, want in CRAFTED:
    got = None
    err = None
    try:
        got = sp.expand(sp.sympify(prod_transform(tex_s), locals=VARS))
    except Exception as e:                                    # noqa
        err = e
    print(f'  {tex_s!r:32s} -> {got if err is None else "ERROR:"+type(err).__name__}'
          + ('' if want is None else f'   (want {want})'))
    if want is not None:
        check(err is None and sp.simplify(got - want) == 0, f'crafted {tex_s!r}')

print()
print('  KNOWN HAZARD: "x^1 0" (LaTeX for x^1 followed by literal 0) is turned')
print('  into x**10 because spaces are deleted BEFORE the ^(\\d) rewrite.')
print('  Impact depends on whether such a pattern occurs in the real source:')
for nm, d in (('F4', 4), ('F5', 5)):
    b = blocks('$F_%d=(F_{%d,1}' % (d, d))
    hz = re.findall(r'\^\s*\d\s*\d', b)
    check(not hz, f'{nm}: no "^ digit digit" (unbraced multi-digit exponent) in source: {hz}')
    cmds = sorted(set(re.findall(r'\\[a-zA-Z]+', b)))
    print(f'  {nm}: LaTeX commands present in block = {cmds}')
    check(set(cmds) <= {'tfrac'}, f'{nm}: only \\tfrac appears (no \\frac/\\cdot/\\left/...)')
    check('(' not in b and ')' not in b, f'{nm}: no parentheses in source (no implicit-mult hazard)')

print()
print('=' * 74)
print('B. Independent strict re-parse and monomial-by-monomial comparison')
print('=' * 74)

# production parse (verbatim copy of production function)
def prod_components(block):
    parts = re.split(r'F_\{\d,\d\}\s*=\{\}&', block)
    comps = []
    for chunk in parts[1:]:
        comps.append(sp.expand(sp.sympify(prod_transform(chunk), locals=VARS)))
    return comps


for nm, d, expdet, expdeg in (('F4', 4, sp.Rational(-44, 9), [4, 11, 12, 21]),
                              ('F5', 5, sp.Rational(160, 29), [3, 12, 14, 16])):
    b = blocks('$F_%d=(F_{%d,1}' % (d, d))
    P = prod_components(b)
    S = strict_components(b, d)
    check(len(P) == 4 and len(S) == 4, f'{nm}: 4 components from both parsers')
    same = all(sp.expand(P[i] - S[i]) == 0 for i in range(4))
    check(same, f'{nm}: strict parser agrees with production parser on all 4 components')
    if not same:
        for i in range(4):
            dd = sp.expand(P[i] - S[i])
            if dd != 0:
                print('    component', i + 1, 'differs by', dd)
    nt = [len(sp.Poly(c, x, y, z, t).terms()) for c in S]
    dg = [sp.Poly(c, x, y, z, t).total_degree() for c in S]
    print(f'  {nm}: #terms per component (strict parse) = {nt}; degrees = {dg}')
    check(dg == expdeg, f'{nm}: component degrees {dg} == claimed {expdeg}')
    # term-count cross-check against the raw source: count +/- separators
    raw = re.sub(r'\\\\(\[\d+pt\])?', '', b).replace('&', '').replace('\n', '')
    raw = re.sub(r'\s+', '', raw)
    chunks = re.split(r'F_\{%d,\d\}=\{\}' % d, raw)[1:]
    for i, ch in enumerate(chunks):
        # count top-level +/- outside braces  => number of source terms
        cnt, depth = 1, 0
        for ci, c in enumerate(ch):
            if c == '{':
                depth += 1
            elif c == '}':
                depth -= 1
            elif c in '+-' and depth == 0 and ci > 0:
                cnt += 1
        check(cnt == nt[i], f'{nm},{i+1}: source term count {cnt} == parsed monomial count {nt[i]}'
                            ' (no cancellation / no dropped or merged terms)')

    J = sp.Matrix([[sp.diff(S[i], v) for v in (x, y, z, t)] for i in range(4)])
    dsym = sp.together(J.det(method='berkowitz'))
    check(sp.simplify(dsym - expdet) == 0,
          f'{nm}: det J (independent parse, Berkowitz algorithm) == {expdet}')
    # numeric spot check at exact rational points, independent of symbolic det
    import random
    random.seed(7)
    ok = True
    for _ in range(4):
        pt = {x: sp.Rational(random.randint(-9, 9), random.randint(1, 7)),
              y: sp.Rational(random.randint(-9, 9), random.randint(1, 7)),
              z: sp.Rational(random.randint(-9, 9), random.randint(1, 7)),
              t: sp.Rational(random.randint(-9, 9), random.randint(1, 7))}
        ok &= (J.subs(pt).det() == expdet)
    check(ok, f'{nm}: det J == {expdet} at 4 exact random rational points')

print()
print('FAILURES:', len(FAIL))
for f_ in FAIL:
    print('  -', f_)
print('PARSER AUDIT', 'PASS' if not FAIL else 'FAIL')
