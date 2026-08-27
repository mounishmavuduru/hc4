*** This FORM program simplifies the output of a MAPLE program
*** which produced the expressions of the second order derivatives 
*** of F(a,b,c,d)
***   The idea was very simple: take into consideration
***   the extra relation (23) for the derivatives
***     F_a, F_b, F_c, F_\xi 
***   and it turns out that this DRAMATICALLY reduces the size
***    of the expressions, and it is possible to cancel one of
***    the Hessians in the denominators.
***  The resulting simplified expressions are in the file
***     "F2ordReduced"

***  Notations:
***     f - the Lagrangian, it depends on a,b,c, which are:
***     a = u_x   (as in the paper)
***     b = u_y
***     c = u_t
***     d = S_x = \xi  (as in the paper)
***     F - The r.h.s. of the first equation for the pseudopotential
***            it depends on a,b,c,d (eq. (19) in the paper)


CFunctions F,f;
Symbols a,b,c,d,i,j,k,l,[1/Hess],[Hess];
Symbols [1/det20],x,[1/x],y,[1/y],z,[1/z ],[1/f(2*b)],[1/F(b)];

L Hess= - (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2);

repeat; 
id F(a)^2= - [1/f(2*b)]*
( - 2*F(a)*F(b)*F(d)*f(b + c) - 2*F(a)*F(b)*f(a + b) + 2*
 F(a)*F(c)*F(d)*f(2*b) + F(b)^2*F(d)^2*f(2*c) + 2*F(b)^2*F(d)*f(a + c)
 + F(b)^2*f(2*a) - 2*F(b)*F(c)*F(d)^2*f(b + c) - 2*F(b)*F(c)*F(d)*f(a
          + b) + F(c)^2*F(d)^2*f(2*b));
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
print;
.sort
   
skip;
L H=Hess;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
print H;
.sort


**** Here we load HUGE expressions profuced by Maple
****  and saved later in FORM internal format.

load F2ord1.sav;
.sort

load F2ord2.sav;
.sort

**** unfortunately the available version of FORM can not load more
****  than 8 expressions simultaneously so we split processing 
**** into several runs commenting the unnecessary load instructions
*load F2ord3.sav;
*.sort

**** I do not know how to make the loaded global expression
****  to work except to make them appear in r.h.s of other expressions
skip;
#do j={F(2*a),F(2*b),F(2*c),F(2*d),F(a+b),F(a+c),F(a+d)}
* ,F(b+c),F(b+d),F(c+d)}
L [N`j']=[`j'];
#enddo
id [1/det20] = [1/Hess] * (  - [1/F(b)]^2 );
.sort

****  Now - simplification! First the extra relation (23)
****    of the paper
repeat; 
id F(a)^2= - [1/f(2*b)]*
( - 2*F(a)*F(b)*F(d)*f(b + c) - 2*F(a)*F(b)*f(a + b) + 2*
 F(a)*F(c)*F(d)*f(2*b) + F(b)^2*F(d)^2*f(2*c) + 2*F(b)^2*F(d)*f(a + c)
 + F(b)^2*f(2*a) - 2*F(b)*F(c)*F(d)^2*f(b + c) - 2*F(b)*F(c)*F(d)*f(a
          + b) + F(c)^2*F(d)^2*f(2*b));
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Bracket [1/f(2*b)],[1/F(b)]; 
*print;
.sort

**** Now we cancel (one of two) Hessians in numerators/denominators
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Bracket [1/f(2*b)],[1/F(b)],[Hess],[1/Hess]; 
print;
.end


