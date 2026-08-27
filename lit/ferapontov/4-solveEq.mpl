####   This programs reads in the output of the previous
####   FORM program (slightly edited) and solves the resulting system
####   in the derivatives of the Lagrangian of order 4

read "3-Lagr3dim4maple":

####  Collect the terms with like powers of \lambda_s, \mu_s,
####  order them lexicographically and transform into a set of 
####  coefficients - these are the equations for 
####   the derivatives of the Lagrangian of order 4

sys:=collect(expand(all),[la(i),la(j),la(k),mu(i),mu(j),mu(k)],ff):
sys:=sort(expand(sys),[la(i),la(j),la(k),mu(i),mu(j),mu(k)],plex):

sys:=subs(la(i)=1,la(j)=1,la(k)=1,mu(i)=1,mu(j)=1,mu(k)=1,
            convert(sys,set)):
nops(sys);
ff:=x->x:
sys:=sys:

Hess:=(- (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) ):

sys:=map(expand,sys):

####  Solve the system

sol:=solve(sys,
{f(4*a),f(4*b),f(4*c),f(3*a+b),f(3*a+c),f(3*b+a),f(3*b+c),
f(3*c+a),f(3*c+b),f(2*a+2*b),f(2*c+2*b),f(2*a+2*c),f(2*a+b+c)
,f(2*b+a+c),f(2*c+b+a)}):

sol[1];

######  Here we test that the solution is correct - just in case...

sys1:=subs(sol,sys):
map(normal,sys1);

####  We introduce a multiplier Hess1 as a substitute  of 1/Hess
####    and save the resulting solution set in file "f4ord.out"
####   - this is to facilitate subsequent FORM processing
sol2:=NULL:
for i in sol  do 
sol2:=sol2,lhs(i)=normal(Hess*Hess1*rhs(i)): end do:
save sol2,"f4ord.out";

quit;

