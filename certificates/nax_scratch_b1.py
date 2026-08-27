import re
import sympy as sp

x, y, z, t = sp.symbols('x y z t')
LOCALS = {'x': x, 'y': y, 'z': z, 't': t}
TEX = open(r'C:\Users\mouni\hc4\lit\gao\Jacobian_CE.tex', encoding='utf-8').read()


def align_block_after(marker):
    i = TEX.index(marker)
    j = TEX.index(r'\begin{align*}', i)
    k = TEX.index(r'\end{align*}', j)
    return TEX[j + len(r'\begin{align*}'):k]


def parse_components(block):
    parts = re.split(r'F_\{\d,\d\}\s*=\{\}&', block)
    comps = []
    for chunk in parts[1:]:
        s = chunk
        s = re.sub(r'\\\\(\[\d+pt\])?', '', s)
        s = s.replace('&', '').replace('\n', '').replace(' ', '')
        s = re.sub(r'\\tfrac\{(-?\d+)\}\{(\d+)\}', r'(\1/\2)*', s)
        s = re.sub(r'\^\{(\d+)\}', r'**\1', s)
        s = re.sub(r'\^(\d)', r'**\1', s)
        s = re.sub(r'(?<![\*/\d])(\d+)(?=[xyzt])', r'\1*', s)
        s = re.sub(r'(?<=[xyzt\)])(?=[xyzt])', '*', s)
        s = re.sub(r'(\*\*\d+)(?=[xyzt])', r'\1*', s)
        comps.append(sp.expand(sp.sympify(s, locals=LOCALS)))
    return comps


F4 = parse_components(align_block_after(r'$F_4=(F_{4,1}'))
F5 = parse_components(align_block_after(r'$F_5=(F_{5,1}'))
for nm, FF in (('F4', F4), ('F5', F5)):
    print(nm, 'degrees', [sp.Poly(f, x, y, z, t).total_degree() for f in FF],
          'terms', [len(sp.Poly(f, x, y, z, t).terms()) for f in FF])
J4 = sp.Matrix(4, 4, lambda i, j: sp.diff(F4[i], [x, y, z, t][j]))
print('det J4 =', sp.factor(J4.det(method='domain-ge')))
J5 = sp.Matrix(4, 4, lambda i, j: sp.diff(F5[i], [x, y, z, t][j]))
print('det J5 =', sp.factor(J5.det(method='domain-ge')))
