#: SmallSize 30000000 
#: LargeSize 100000000 
#: TermsInSmall 1000000 
* off statistics;
CFunctions F,G,f,D,diff,mu,la,B,nB,dB,U;
Symbols xx,a,b,c,d,i,j,k,l,[1/Hess],[Hess];
S [1/det20],x,[1/x],y,[1/y],z,[1/z ],[1/f(2*b)],[1/F(b)];
L Hess= - (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2);
Local U1=F(d)*G(a)-F(a)*G(d)-f(a+a)*[1/f(2*b)]*F(b);
Local U2=F(d)*G(b)-F(b)*G(d)-2*f(a+b)*[1/f(2*b)]*F(b)+F(a);
Local U3=F(d)*G(c)-F(c)*G(d)-2*f(a+c)*[1/f(2*b)]*F(b)-G(a);
Local U4=G(c)+f(c+c)*[1/f(2*b)]*F(b);
Local U5=G(b)+2*f(b+c)*[1/f(2*b)]*F(b)-F(c);
L [det20] = [Hess] * (  - F(b)^2 );
.sort

**** We differentiate the equations (20) for
**** the first order derivatives of F(a,b,c,d) and G(a,b,c,d)
#do i={1,2,3,4,5}
#do j={a,b,c,d}
L U`i'`j'= diff(`j')*U`i';
#enddo
#enddo
repeat;
id diff(i?)*F(j?)=F(j+i)+  D(F(j))*diff(i);
id diff(i?)*G(j?)=G(j+i)+  D(G(j))*diff(i);
id diff(i?{a,b,c})*f(j?)=f(j+i)+  D(f(j))*diff(i);
** !!!!! ****
id diff(i?{a,b,c})*[1/f(2*b)] = - [1/f(2*b)]^2*f(2*b+i)
    +  D([1/f(2*b)])*diff(i);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
id D(a?,a)=1;
.sort

**** We have to eliminate ALL derivatives of G from
**** the obtained expressions, for this we first
**** define a preprocessor macroprocedure "sol"
**** which solves a linear equation u in one unknown p
**** and substitutes the obtained expression into
**** all equations
#procedure sol(u,p)
skip;
L tsol=`u';
id `p'=0;
.sort
skip;
nskip tsol,co;
L co=`u'-tsol;
id `p'=1;
.sort;
G [`p']=-(tsol)/co;
id `p'=[`p'];
.sort
repeat; 
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort
#endprocedure


**** now we solve first-order equations for G
**** so in fact we get the expressions (22)
**** (and substitute the results into all equations)
**** 
#call sol(U3,G(a))
#call sol(U4,G(c))
#call sol(U5,G(b))

**** it's better to solve the next equation "by hand"
*#call sol(U2,G(d))
G [G(d)] = (+ F(a) - 2*F(b)*F(d)*f(b + c)*[1/f(2*b)]
    - 2*F(b)*f(a + b)*[1/f(2*b)] + F(c)*F(d)) * [1/F(b)];
**** This expression is in fact the relation (23)
L U1F= -U1*f(2*b)*F(b);
id G(d) = [G(d)];
repeat; 
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort


**** now we solve second-order equations for G
#call sol(U3a,G(2*a))
#call sol(U3b,G(a+b))
#call sol(U3c,G(a+c))
#call sol(U3d,G(a+d))

#call sol(U4b,G(b+c))
#call sol(U4c,G(2*c))
#call sol(U4d,G(c+d))

#call sol(U5b,G(2*b))
#call sol(U5d,G(b+d))

**** again it's better to do "by hand"
L [G(2*d)] =    (F(a + d) - F(a)*F(b + d)*[1/F(b)]
       - 2*F(b)*F(2*d)*f(b + c)*[1/f(2*b)]  + F(c)*F(2*d)
       - F(c)*F(d)*F(b + d)*[1/F(b)] + F(d)*F(c + d))*[1/F(b)];
id G(2*d) = [G(2*d)];
repeat; 
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

**** Now we input the simplified expressions for
**** the second-order derivatives of F(a,b,c,d)
**** (in terms of its first-order derivatives and the derivatives
**** of the Lagrangian)
skip;
#include- F2ordReduced;
.sort

**** Here we assign the second-order derivatives of F(a,b,c,d)
**** to their values given in "F2ordReduced"
**** and substitute them into the equations Uij
skip;
nskip
#do i={1,2,3,4,5}
#do j={a,b,c,d}
U`i'`j'
#enddo
#enddo
;
#do j={F(2*a),F(2*b),F(2*c),F(2*d),F(a+b),F(a+c),F(a+d),F(b+c),F(b+d),F(c+d)}
id `j'=[`j'];
#enddo
.sort

**** We get rid of the common denominator to make them
**** easier to inspect
skip;
nskip
#do i={1,2,3,4,5}
#do j={a,b,c,d}
U`i'`j'
#enddo
#enddo
;
Multiply [det20]^3;
repeat; 
id F(b)*[1/F(b)]=1;
id [Hess]*[1/Hess]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

**** Now comes the first simplification:
**** we cancel the Hessian (of the Lagrangian)
**** if it is present 
skip;
nskip
#do i={1,2,3,4,5}
#do j={a,b,c,d}
U`i'`j'
#enddo
#enddo
;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

****  the main simplification: we take (23)
**** into account and print the resulting Uij
**** - they are zeros as expected
skip;
nskip
#do i={1,2,3,4,5}
#do j={a,b,c,d}
U`i'`j'
#enddo
#enddo
;
repeat; 
id [Hess]*[1/Hess]=1;
id F(a)^2= - [1/f(2*b)]*
( - 2*F(a)*F(b)*F(d)*f(b + c) - 2*F(a)*F(b)*f(a + b) + 2*
 F(a)*F(c)*F(d)*f(2*b) + F(b)^2*F(d)^2*f(2*c) + 2*F(b)^2*F(d)*f(a + c)
 + F(b)^2*f(2*a) - 2*F(b)*F(c)*F(d)^2*f(b + c) - 2*F(b)*F(c)*F(d)*f(a
          + b) + F(c)^2*F(d)^2*f(2*b));
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
print;
.end

