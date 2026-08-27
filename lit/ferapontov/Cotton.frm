CFunction Gamma,RR,f,fo,D,diff,Haan,Haan1,Nijen,
  u,v,dv,f2,fo2,[1/sqrtHess];
Symbols xx,a,b,c,[1/Hess],[Hess],[1/f(2*b)];
Symbols i,j,k,l1,i1,j1,k1,l2,i2,j2,k2,m,r,p,q,l,s,s1,s2,s3;

#include- f4ord.frm;
.sort

skip;
L Hess:=(- (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) );
L dxxHess=diff(xx)*Hess;
repeat; 
id diff(i?)*f(j?)=f(j+i)+  D(f(j))*diff(i);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
.sort

skip;
**** hereafter we set Gkij = Gamma^k1_{i1 j1}  !!!!!!!
L Gkij=sum_(s1,1,3,fo(k1,s1)*f(s1,i1,j1)/2);
.sort

skip;
**** hereafter we set Rijkl = R_{ijkl}  !!!!!!!
L Rijkl=sum_(s1,1,3,sum_(s2,1,3,
   fo(s1,s2)*( f(i,s1,k)*f(s2,j,l)-f(i,s1,l)*f(s2,j,k) )/4));
.sort

skip;
L Rjl=sum_(s1,1,3,sum_(s2,1,3,fo(s1,s2)*RR(s1,j1,s2,l1)));
id RR(i?,j?,k?,l?)=Rijkl;
.sort

skip;
L R=sum_(s1,1,3,sum_(s2,1,3,fo(s1,s2)*RR(s1,s2)));
id RR(j1?,l1?)=Rjl;
repeat;
****  this substitutes the expressions for the elements of the inverse 
****  matrix  (f_{ij})^(-1)
#include- fo.frm
id f(i1?,j1?,k1?)=f(u(i1)+u(j1)+u(k1));
id f(i1?,j1?)=f(u(i1)+u(j1));
argument;
id u(1)=a;
id u(2)=b;
id u(3)=c;
endargument;
endrepeat;
*print;
.sort

skip;
L Pij=RR(i,j) - R*f(i,j)/4;
id RR(j1?,l1?)=Rjl;
*print;
.sort

skip;
L dkPij=diff(u(k))*Pij - sum_(s1,1,3,Gamma(s1,k,i)*RR(s1,j))
    - sum_(s1,1,3,Gamma(s1,k,j)*RR(s1,i));
repeat; 
id Gamma(k1?,i1?,j1?)=Gkij;
id RR(i?,j?)=Pij;
endrepeat;
.sort


skip;
L rez=RR(1,2,3) - RR(2,1,3);
id RR(k?,i?,j?)=dkPij;
repeat;
#include- fo.frm
id f(i1?,j1?,k1?,l1?)=f(u(i1)+u(j1)+u(k1)+u(l1));
id f(i1?,j1?,k1?)=f(u(i1)+u(j1)+u(k1));
id f(i1?,j1?)=f(u(i1)+u(j1));
argument;
id u(1)=a;
id u(2)=b;
id u(3)=c;
endargument;
endrepeat;
.sort

skip;
nskip rez;
repeat; 
id diff(i?)*f(j?)=f(j+i)+  D(f(j))*diff(i);
id diff(xx?)*[1/Hess] = - [1/Hess]^2*dxxHess +  D([1/Hess])*diff(xx);
endrepeat;
id diff(?a)=0;
id D(a?)=a;
.sort

skip;
nskip rez;
repeat; 
#do i ={f(4*a),f(4*b),f(4*c),f(3*a+b),f(3*a+c),f(3*b+a),f(3*b+c),f(3*c+a),f(3*c+b),f(2*a+2*b),f(2*b+2*c),f(2*a+2*c),f(2*a+b+c),f(2*b+a+c),f(2*c+b+a)}
 id `i'=[`i'];
#enddo
endrepeat;
.sort

skip;
nskip rez;
repeat; 
id f(a + b)^2*f(2*c)*[1/Hess] = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2)*[1/Hess] - 1;
endrepeat;
Bracket [1/Hess];
print;
.end



**************  the rest is just a collection of useful ticks I did
**************   on some stages of computations


*** THIS DOES NOT WORK !!!!
*** because one needs "argument" environment, see below
*skip;
*nskip rez;
*id i=1;
*id j=2;
*id k=3;
*print;
*.sort

skip;
nskip rez;
repeat;
id f(i1?,j1?,k1?)=f(u(i1)+u(j1)+u(k1));
id f(i1?,j1?)=f(u(i1)+u(j1));
argument;
id u(1)=a;
id u(2)=b;
id u(3)=c;
endargument;
endrepeat;
Bracket [1/Hess];
print;
.end




skip;
L fijkl1 = (  sum_(p,1,3, sum_(q,1,3,
             f(i1,j1,k1)*fo(p,q)*f(p,q,l1) +
             f(i1,j1,l1)*fo(p,q)*f(p,q,k1) +
             f(i1,l1,k1)*fo(p,q)*f(p,q,j1) +
             f(j1,l1,k1)*fo(p,q)*f(p,q,i1) )))/4;
L fijkl2 =    3/2*
 ( sum_(s1,1,3, 
  sum_(s2,1,3,
  sum_(p,1,3,
  sum_(q,1,3,
  sum_(r,1,3,
  sum_(s,1,3,
     f(p,q,i1)*f(r,s,j1)*f(s1,k1)*f(s2,l1)*e_(p,r,s1)*e_(q,s,s2) +
     f(p,q,i1)*f(r,s,k1)*f(s1,j1)*f(s2,l1)*e_(p,r,s1)*e_(q,s,s2) +
     f(p,q,i1)*f(r,s,l1)*f(s1,j1)*f(s2,k1)*e_(p,r,s1)*e_(q,s,s2) +
     f(p,q,j1)*f(r,s,k1)*f(s1,i1)*f(s2,l1)*e_(p,r,s1)*e_(q,s,s2) +
     f(p,q,j1)*f(r,s,l1)*f(s1,i1)*f(s2,k1)*e_(p,r,s1)*e_(q,s,s2) +
     f(p,q,k1)*f(r,s,l1)*f(s1,i1)*f(s2,j1)*e_(p,r,s1)*e_(q,s,s2) 
     )))))))*[1/Hess]/6;
L fijkl=fijkl1-fijkl2;
repeat;
id e_(1,2,3)=1;
endrepeat;
.sort

skip;
nskip rez;
repeat;
id [1/Hess]*f(1,1)*f(2,2)*f(3,3)= 1-
   [1/Hess]*(- ( f(1,2)^2*f(3,3)
 + f(2,2)*f(1,3)^2 - 2*f(1,2)*f(1,3)*f(2,3) + f(1,1)*f(2,3)^2) );
id f(i1?,j1?,k1?,l1?)=fijkl;
endrepeat;
.sort

skip;
nskip rez;
repeat;
#include- fo.frm
endrepeat;
.sort

skip;
nskip rez;
repeat;
id [1/Hess]*f(1,1)*f(2,2)*f(3,3)= 1-
   [1/Hess]*(- ( f(1,2)^2*f(3,3)
 + f(2,2)*f(1,3)^2 - 2*f(1,2)*f(1,3)*f(2,3) + f(1,1)*f(2,3)^2) );
endrepeat;
.sort

