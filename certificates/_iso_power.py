import os, re, random, subprocess
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, E_components, coeffs_in_y, AV, s_poly

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survtgt')
os.makedirs(ND, exist_ok=True)


def run(t):
    fn = os.path.join(ND, 'pw.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(t)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout 250 Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


random.seed(9)
p = 32003
C0, C1, C2 = (random.randrange(1, p) for _ in range(3))
b3 = C2 * y1**3 + C0 * y1**2 * y2 + C1 * y1**2 * y3
(E3, E2, E1, E0), A = E_components(a5, b3)
J = []
for e in (E3, E2, E1, E0):
    J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
vs = sp.Matrix([0, C1, -C0])
iso = sp.expand((vs.T * A * vs)[0])
# reduce iso^k in the ring where iso is a single poly in a,y -> need y as ring vars too.
# Work in ring over F_p with BOTH a-vars and y-vars; J is in a-vars (still valid).
YV = 'y1,y2,y3'
txt = [f'ring R={p},({AV},{YV}),dp;',
       'ideal J=' + ',\n'.join(s_poly(g) for g in J) + ';',
       'ideal S=std(J);',
       f'poly iso={s_poly(iso)};',
       'poly r1=reduce(iso,S);   "NZ1", (r1!=0);',
       'poly r2=reduce(iso^2,S); "NZ2", (r2!=0);',
       'poly r3=reduce(iso^3,S); "NZ3", (r3!=0);',
       'quit;']
out = run('\n'.join(txt) + '\n')
print(f'p={p} c=({C0},{C1},{C2})')
for k in (1, 2, 3):
    m = re.search(rf'NZ{k}\s+(\d+)', out)
    print(f'  iso^{k} reduces to nonzero (1) / zero(0): {m.group(1) if m else "??"}')
print('--- raw tail ---')
print(out[-300:])
