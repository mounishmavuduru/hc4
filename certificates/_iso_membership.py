import os, re, random, subprocess
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, E_components, coeffs_in_y, AV, s_poly

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survtgt')
os.makedirs(ND, exist_ok=True)


def run(t):
    fn = os.path.join(ND, 'm.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(t)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout 200 Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


random.seed(5)
for p in (32003, 40009, 15013):
    C0, C1, C2 = (random.randrange(1, p) for _ in range(3))
    b3 = C2 * y1**3 + C0 * y1**2 * y2 + C1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    vs = sp.Matrix([0, C1, -C0])
    iso = sp.expand((vs.T * A * vs)[0])
    G = [q for q in coeffs_in_y(iso) if q != 0]
    txt = [f'ring R={p},({AV}),dp;',
           'ideal J=' + ',\n'.join(s_poly(g) for g in J) + ';',
           'ideal S=std(J);',
           'ideal T=' + ',\n'.join(s_poly(g) for g in G) + ';',
           '"NZ",size(simplify(reduce(T,S),2));',
           'quit;']
    out = run('\n'.join(txt) + '\n')
    m = re.search(r'NZ\s+(\d+)', out)
    print(f'p={p} c=({C0},{C1},{C2}): iso has {len(G)} coeffs; '
          f'NZ(not in J)= {m.group(1) if m else "??" }')
