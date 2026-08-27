with(LinearAlgebra):
Mx:=convert(
[
[0,f(x+a),f(x+b),f(x+c)],
[f(x+a),f(x+a+a),f(x+a+b),f(x+a+c)],
[f(x+b),f(x+b+a),f(x+b+b),f(x+b+c)],
[f(x+c),f(x+c+a),f(x+c+b),f(x+c+c)]
]
,Matrix);

Ma:=subs(x=a,Mx);
Mb:=subs(x=b,Mx);
Mc:=subs(x=c,Mx);

DF3:=Determinant(da*Ma+db*Mb+dc*Mc):

p:=[a,b,c];
dp:=[da,db,dc];
d4f:=add(add(add(add(f(p[i]+p[j]+p[k]+p[l])*
   dp[i]*dp[j]*dp[k]*dp[l],i=1..3),j=1..3),k=1..3),l=1..3):

d3f:=add(add(add(f(p[i]+p[j]+p[k])*
   dp[i]*dp[j]*dp[k],i=1..3),j=1..3),k=1..3):

dxHess :=
   - f(x + 2*a)*f(b + c)^2 - f(x + 2*b)*f(a + c)^2 - f(x + 2*c)*f(a + b)^2
   + 2*f(a + b)*f(a + c)*f(x + b + c) + 2*f(a + b)*f(b + c)*f(x + a + c)
   + 2*f(a + c)*f(b + c)*f(x + a + b) - 2*f(2*a)*f(b + c)*f(x + b + c) + 
   f(2*a)*f(2*b)*f(x + 2*c) + f(2*a)*f(2*c)*f(x + 2*b) - 2*f(2*b)*f(a + c)*
  f(x + a + c) + f(2*b)*f(2*c)*f(x + 2*a) - 2*f(2*c)*f(a + b)*f(x + a + b):

dH:=subs(x=a, dxHess)*da+subs(x=b, dxHess)*db+subs(x=c, dxHess)*dc:

read "f4ord.mpl":

test:=expand(d4f*Hess- d3f*dH - 3* DF3);
