# char-0 point check: iso^2 in J over Q at a fixed generic integer c.
import os, re, subprocess
import sympy as sp
import _d5_close as base
from _d5_close import Y, a5, E_components, coeffs_in_y, AV, s_poly

y1, y2, y3 = Y
ND = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_lean', 'survtgt')
os.makedirs(ND, exist_ok=True)


def run(t, tmo=550):
    fn = os.path.join(ND, 'q.sing')
    with open(fn, 'w', newline='\n') as f:
        f.write(t)
    w = '/mnt/c/' + os.path.abspath(fn)[3:].replace('\\', '/')
    return subprocess.run(base.WSL + ['bash', '-c', f'timeout {tmo} Singular -q < "{w}"'],
                          capture_output=True, text=True).stdout


for (C0, C1, C2) in [(2, 3, 5), (7, 4, 9)]:
    b3 = C2 * y1**3 + C0 * y1**2 * y2 + C1 * y1**2 * y3
    (E3, E2, E1, E0), A = E_components(a5, b3)
    J = []
    for e in (E3, E2, E1, E0):
        J += [q for q in coeffs_in_y(sp.expand(e)) if q != 0]
    vs = sp.Matrix([0, C1, -C0])
    iso = sp.expand((vs.T * A * vs)[0])
    txt = [f'ring R=0,({AV},y1,y2,y3),dp;',
           'ideal J=' + ',\n'.join(s_poly(g) for g in J) + ';',
           'ideal S=std(J);',
           f'poly iso={s_poly(iso)};',
           'poly r2=reduce(iso^2,S); "ISO2_NONZERO", (r2!=0);',
           'quit;']
    out = run('\n'.join(txt) + '\n')
    m = re.search(r'ISO2_NONZERO\s+(\d+)', out)
    print(f'Q c=({C0},{C1},{C2}): iso^2 in J ? -> {"NO" if (m and m.group(1)==chr(49)) else ("YES" if m else "?? "+out[-120:])}')
