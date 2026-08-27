****  This program tests that the expressions for the derivatives
****   of the Lagrangian of order 4 found from the requirement
****   of existence of hydrodynamic reductions
****   imply existence of a pseudopotential.
****   These 4th-order eqns for f(a,b,c) are read from f4ord.frm,
****   the 2nd-order eqns for F(a,b,c,d) are read from F2ordReduced.
****   We form the 3rd-order compatibility conditions for F
****   and check that they are satisfied modulo the 4th-order
****    relations for f.

#: SmallSize 10000000 
#: LargeSize 30000000 
#: TermsInSmall 1000000 
* off statistics;
CFunctions F,G,f,D,diff;
Symbols xx,a,b,c,d,i,j,k,l,[1/Hess],[Hess],[1/det20],[1/f(2*b)],[1/F(b)];

#define p1 "a"
#define p2 "b"
#define p3 "c"
#define p4 "d"

#$nmr=1;

#do i=1,3
#do j={`i'+1},4
L rez`$nmr' = diff(`p`j'')*F(2*`p`i'')-diff(`p`i'')*F(`p`i''+`p`j'');
*#redefine nmr {`nmr'+1}
#$nmr=`$nmr'+1;
#enddo
#enddo

#do i=1,3
#do j={`i'+1},4
L rez`$nmr' = diff(`p`i'')*F(2*`p`j'')-diff(`p`j'')*F(`p`i''+`p`j'');
#$nmr=`$nmr'+1;
#enddo
#enddo

#do i=1,2
#do j={`i'+1},3
#do k={`j'+1},4
L rez`$nmr' = diff(`p`j'')*F(`p`i''+`p`k'')-diff(`p`k'')*F(`p`i''+`p`j'');
#$nmr=`$nmr'+1;
L rez`$nmr' = diff(`p`i'')*F(`p`j''+`p`k'')-diff(`p`k'')*F(`p`i''+`p`j'');
#$nmr=`$nmr'+1;
L rez`$nmr' = diff(`p`j'')*F(`p`i''+`p`k'')-diff(`p`i'')*F(`p`k''+`p`j'');
#$nmr=`$nmr'+1;
#enddo
#enddo
#enddo

print;
.sort

skip;
L Hess= - (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2);
L [det20] = [Hess] * (  - F(b)^2 );
.sort

skip;
L d1Hess=diff(i)*Hess;
repeat;
id diff(i?)*f(j?)=f(j+i)+  D(f(j))*diff(i);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
id D(a?,a)=1;
*print d1Hess;
.sort

skip;
nskip d1Hess;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
*print d1Hess;
.sort

skip;
#include- F2ordReduced;
*Bracket [1/f(2*b)],[1/F(b)],[1/det20],[Hess],[1/Hess]; 
.sort

skip;
nskip 
#do i=1,`$nmr'
rez`i'
#enddo
;
repeat; 
#do j={F(2*a),F(2*b),F(2*c),F(2*d),F(a+b),F(a+c),F(a+d),F(b+c),F(b+d),F(c+d)}
id `j'=[`j'];
#enddo
id F(b)*[1/F(b)]=1;
id [Hess]*[1/Hess]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

skip;
nskip 
#do i=1,`$nmr'
rez`i'
#enddo
;
repeat;
id diff(i?)*F(j?)=F(j+i)+  D(F(j))*diff(i);
id diff(i?{a,b,c})*f(j?)=f(j+i)+  D(f(j))*diff(i);
** !!!!! ****
id diff(i?{a,b,c})*[1/f(2*b)] = - [1/f(2*b)]^2*f(2*b+i)
    +  D([1/f(2*b)])*diff(i);
id diff(i?)*[1/F(b)] = - [1/F(b)]^2*F(b+i)
    +  D([1/F(b)])*diff(i);
id diff(i?{a,b,c})*[1/Hess] = - [1/Hess]^2*d1Hess
    +  D([1/Hess])*diff(i);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
id D(a?,a)=1;
*print;
.sort


skip;
nskip 
#do i=1,`$nmr'
rez`i'
#enddo
;
repeat; 
#do j={F(2*a),F(2*b),F(2*c),F(2*d),F(a+b),F(a+c),F(a+d),F(b+c),F(b+d),F(c+d)}
id `j'=[`j'];
#enddo
id F(b)*[1/F(b)]=1;
id [Hess]*[1/Hess]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

skip;
nskip 
#do i=1,`$nmr'
rez`i'
#enddo
;
repeat; 
id [Hess]*[1/Hess]=1;
id F(a)^2= - [1/f(2*b)]*
( - 2*F(a)*F(b)*F(d)*f(b + c) - 2*F(a)*F(b)*f(a + b) + 2*
 F(a)*F(c)*F(d)*f(2*b) + F(b)^2*F(d)^2*f(2*c) + 2*F(b)^2*F(d)*f(a + c)
 + F(b)^2*f(2*a) - 2*F(b)*F(c)*F(d)^2*f(b + c) - 2*F(b)*F(c)*F(d)*f(a
          + b) + F(c)^2*F(d)^2*f(2*b));
id F(b)*[1/F(b)]=1;
id f(2*b)*[1/f(2*b)]=1;
id [Hess]*[1/Hess]=1;
endrepeat;
.sort

skip;
nskip 
#do i=1,`$nmr'
rez`i'
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

skip;
#include- f4ord.frm
.sort

#do i ={f(a+3*b),f(3*a+c),f(a+b+2*c),f(a+3*c),f(4*c),f(a+2*b+c),f(2*a+b+c),f(2*a+2*b),f(3*a+b),f(b+3*c),f(4*a),f(2*b+2*c),f(4*b),f(2*a+2*c),f(3*b+c)}
 id `i'=[`i'];
#enddo
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

skip;
nskip 
#do i=1,`$nmr'
rez`i'
#enddo
;
print;
.end

