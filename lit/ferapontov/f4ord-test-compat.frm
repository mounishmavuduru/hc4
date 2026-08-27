****   This program tests the compatibility conditions for the derivatives
****   of the Lagrangian of order 4 found by Maple.
****   The derivatives are read from the file f4ord.frm


CFunctions f,D,diff;
Symbols xx,a,b,c,i,j,k,[1/Hess],[Hess],[1/f(2*b)];
L Hess:=(- (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) );
L dxxHess=diff(xx)*Hess;
repeat; 
id diff(i?)*f(j?)=f(j+i)+  D(f(j))*diff(i);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
.sort

#include- f4ord.frm;
.sort


****   This module forms the compatibility conditions to be checked.
****   They are printed in the end of the module for inspection.
#define p1 "a"
#define p2 "b"
#define p3 "c"
#$nmr=0;
skip;
#do k ={3*a,3*b,3*c,2*a+b,2*a+c,2*b+a,2*b+c,2*c+a,2*c+b,a+b+c}
#do i=1,2
#do j={`i'+1},3
#$nmr=`$nmr'+1;
L compat`$nmr' = diff(`p`j'')*f(`k'+`p`i'')-diff(`p`i'')*f(`k'+`p`j'');
#enddo
#enddo
#enddo
print;
.sort

skip;
nskip 
#do i=1,`$nmr'
compat`i'
#enddo
;
#do i ={f(4*a),f(4*b),f(4*c),f(3*a+b),f(3*a+c),f(3*b+a),f(3*b+c),f(3*c+a),f(3*c+b),f(2*a+2*b),f(2*b+2*c),f(2*a+2*c),f(2*a+b+c),f(2*b+a+c),f(2*c+b+a)}
 id `i'=[`i'];
#enddo
.sort

skip;
nskip 
#do i=1,`$nmr'
compat`i'
#enddo
;
repeat; 
id diff(i?)*f(j?)=f(j+i)+  D(f(j))*diff(i);
id diff(i?)*[1/f(2*b)] = - [1/f(2*b)]^2*f(2*b+i) +  D([1/f(2*b)])*diff(i);
id diff(xx?)*[1/Hess] = - [1/Hess]^2*dxxHess +  D([1/Hess])*diff(xx);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
.sort

skip;
nskip 
#do i=1,`$nmr'
compat`i'
#enddo
;
repeat; 
#do i ={f(4*a),f(4*b),f(4*c),f(3*a+b),f(3*a+c),f(3*b+a),f(3*b+c),f(3*c+a),f(3*c+b),f(2*a+2*b),f(2*b+2*c),f(2*a+2*c),f(2*a+b+c),f(2*b+a+c),f(2*c+b+a)}
 id `i'=[`i'];
#enddo
endrepeat;
.sort

skip;
nskip 
#do i=1,`$nmr'
compat`i'
#enddo
;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
print;
.end

