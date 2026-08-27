# scratch: master formula at the level of 2-JETS (independent symbols) => complete proof.
import sympy as sp

z = sp.symbols('z1 z2 z3')
lam, mu = sp.symbols('lam mu')

def jet(nm):
    return {k: sp.Symbol(f'{nm}_{k}') for k in ('', '1', '2', '11', '12', '22')}

a = [jet(f'a{i}') for i in range(3)]
b = [jet(f'b{i}') for i in range(3)]
d1 = {'': '1', '1': '11', '2': '12'}
d2 = {'': '2', '1': '12', '2': '22'}

S = sum(a[i][''] * z[i] for i in range(3))
Sp = [sum(a[i][d[''][0] if False else ('1' if d is d1 else '2')] * z[i] for i in range(3)) for d in (d1, d2)]
Sp = [sum(a[i]['1']*z[i] for i in range(3)), sum(a[i]['2']*z[i] for i in range(3))]
Spq = {(0, 0): sum(a[i]['11']*z[i] for i in range(3)),
       (0, 1): sum(a[i]['12']*z[i] for i in range(3)),
       (1, 1): sum(a[i]['22']*z[i] for i in range(3))}
Bpq = {(0, 0): sum(b[i]['11']*z[i] for i in range(3)),
       (0, 1): sum(b[i]['12']*z[i] for i in range(3)),
       (1, 1): sum(b[i]['22']*z[i] for i in range(3))}
pk = ['1', '2']
H = sp.zeros(5, 5)
for p in range(2):
    for q in range(2):
        key = (min(p, q), max(p, q))
        H[p, q] = sp.expand(Bpq[key] + lam*(Sp[p]*Sp[q] + S*Spq[key]) + mu*Spq[key])
for p in range(2):
    for j in range(3):
        H[p, 2+j] = H[2+j, p] = sp.expand(b[j][pk[p]] + lam*(a[j][pk[p]]*S + a[j]['']*Sp[p]) + mu*a[j][pk[p]])
for j in range(3):
    for k in range(3):
        H[2+j, 2+k] = sp.expand(lam*a[j]['']*a[k][''])

lhs = sp.expand(H.det(method='berkowitz'))
tau = lam*S + mu
JF = sp.Matrix(3, 3, lambda i, j: [sp.expand(b[i]['1'] + tau*a[i]['1']),
                                   sp.expand(b[i]['2'] + tau*a[i]['2']),
                                   a[i]['']][j])
rhs = sp.expand(lam*sp.expand(JF.det())**2)
print('MASTER FORMULA (2-jet, fully generic) holds:', sp.expand(lhs - rhs) == 0)
