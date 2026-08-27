***** a FORM program for writing down the conditions
*****    of a Lagrangian f(u_x,u_y,u_t) to produce
*****     an integrable (in the sense described in the paper) equation.

*****   Notations: f - the Lagrangian
*****                  hereafter we denote as
*****                  f(a+b) = \partial^2 f/\partial a \partial b
*****                  f(2*a) = \partial^2 f/\partial a^2
*****                  f(2*a+c) = \partial^3 f/\partial a^2 \partial c
*****                     etc..
*****         DispEq - l.h.s of the dispersion equation (8)
*****              a - u_x
*****              b - u_y
*****              c - u_t
*****          mu(i) - \mu^i
*****          la(i) - \la^i
*****     B(i,j)=Bij - B_{ij}
*****            Dij - the denominator of B_{ij}
*****            Nij - the numerator of B_{ij}
*****            Uij - the "inverted" (modulo DispEq) Dij:
*****                 1/Dij=Uij=(...)/((\la^i-\la^j)^{2}*Hess) (mod DispEq)
*****         [Hess] - Hessian of the Lagrangian f (the determinant)
*****       [1/Hess] - (Hessian)^{-1} of the Lagrangian f
*****          [1/x] - (\la^i-\la^j)^{-1} = [1/lji](i+j)
*****          [1/y] - (\la^j-\la^k)^{-1} = [1/lji](k+j)
*****          [1/z] - (\la^k-\la^i)^{-1} = [1/lji](i+k)

*off statistics;
CFunctions f,D,diff,mu,la,B,U;
Symbols a,b,c,i,j,k,[1/Hess],[Hess];
Symbols x,[1/x],y,[1/y],z,[1/z],[1/f(2*b)];
CFunctions [lji],[1/lji];
Local DispEq=f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+ 
  f(b+b)*mu(i)^2+f(c+c)*la(i)^2+2*f(c+b)*mu(i)*la(i);

**** Here we form the numerator Uij of the "inverted" (modulo DispEq) Dij ****
**** the denominator of Uij iz dUij=1/((\la^i-\la^j)^2*Hess)=
****    =[1/lji](i+j)^2*[1/Hess] to avoid (formally) denominators
****        - FORM is NOT strong in working with denominators...
L dUij=([1/lji](j+i)^2)*[1/Hess];
L xU =   (-f(a + b) - f(b + c)*la(i))*
     (f(a + b) + f(b + c)*la(j)) *dUij+
   (f(a + b) + f(b + c)*la(i))*(-f(a + b) -
      f(b + c)*la(j))*dUij +
   (f(2*a)*f(2*b) + f(2*b)*f(a + c)*la(i) +
     f(2*b)*f(a + c)*la(j) + f(2*b)*f(2*c)*la(i)*
      la(j))*dUij;
L yU = -(f(2*b)*(f(a + b) + f(b + c)*la(j)))*dUij;
L zU = -(f(2*b)*(f(a + b) + f(b + c)*la(i)))*dUij;
L tU = -f(2*b)^2*dUij;
L Uij=(xU+mu(i)*yU+mu(j)*zU+mu(i)*mu(j)*tU)/2;

L Hess= - (-(f(2*a)*f(2*b)*f(2*c)) + f(a + b)^2*f(2*c)
 + f(2*b)*f(a + c)^2 - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2);

**** The numerator and the denominator of B_{ij},
****   in fact they were found differentiating the DispEq,
****  in the next FORM module we test this.
L   Dij =
       - 2*f(a + b)*mu(i) - 2*f(a + b)*mu(j) - 2*f(a + c)*la(i)
      - 2*f(a + c)*la(j)  - 2*f(2*a) - 2*f(b + c)*mu(i)*la(j)
   - 2*f(b + c)*mu(j)*la(i) - 2*f(2*b)*mu(i)*mu(j) - 2*f(2*c)*la(i)*la(j);
L   Nij =
      2*f(a + b)*f(a + 2*b)*mu(i)*[1/f(2*b)] + 2*f(a + b)*f(2*b + c)*mu(i)*la(
      j)*[1/f(2*b)] - 2*f(a + b + c)*mu(i)*la(i) - 2*f(a + b + c)*mu(i)*la(j)
       - 2*f(a + b + c)*mu(j)*la(i) + 2*f(a + 2*b)*f(a + c)*la(i)*[1/f(2*b)]
       + 2*f(a + 2*b)*f(b + c)*mu(i)*la(i)*[1/f(2*b)] - 2*f(a + 2*b)*mu(i)*mu(
      j) + 2*f(a + c)*f(2*b + c)*la(i)*la(j)*[1/f(2*b)] - f(a + 2*c)*la(i)^2
       - 2*f(a + 2*c)*la(i)*la(j) + f(2*a)*f(a + 2*b)*[1/f(2*b)] + f(2*a)*f(2*
      b + c)*la(j)*[1/f(2*b)] + f(2*a)*f(3*b)*mu(j)*[1/f(2*b)] - 2*f(2*a + b)*
      mu(i) - f(2*a + b)*mu(j) - 2*f(2*a + c)*la(i) - f(2*a + c)*la(j) - f(3*a
      ) + 2*f(b + c)*f(2*b + c)*mu(i)*la(i)*la(j)*[1/f(2*b)] - 2*f(b + 2*c)*
      mu(i)*la(i)*la(j) - f(b + 2*c)*mu(j)*la(i)^2 - 2*f(2*b + c)*mu(i)*mu(j)*
      la(i) + 2*f(3*b)*f(a + b)*mu(i)*mu(j)*[1/f(2*b)] + 2*f(3*b)*f(a + c)*mu(
      j)*la(i)*[1/f(2*b)] + 2*f(3*b)*f(b + c)*mu(i)*mu(j)*la(i)*[1/f(2*b)] + 
      f(3*b)*f(2*c)*mu(j)*la(i)^2*[1/f(2*b)] + f(2*c)*f(a + 2*b)*la(i)^2*
      [1/f(2*b)] + f(2*c)*f(2*b + c)*la(i)^2*la(j)*[1/f(2*b)] - f(3*c)*la(i)^2
      *la(j);
.sort

**** From now on we set B_{ij}=Uij*Nij , instead of B_{ij}=Nij/Dij
skip;
L Bij=Nij*Uij;
*** this forces FORM calculate modulo DispEq
repeat; 
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
***  and cancel in terms where it is in both the numerator and denominator
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

****  in the next three modules we test that Uij,Nij,Dij are correct
skip;
L pr=Uij*Dij;
L dkDispEq=diff(j)*DispEq;
**** we use the standard trick to define differentiation
**** which works fine for sums/products/powers using the id statements
**** see e.g. FORM tutorials
print;
.sort
skip;
nskip pr, dkDispEq;
repeat; 
id diff(j?)*mu(i?)=D(j,a)*B(i,j)*(mu(i)-mu(j))+ D(mu(i))*diff(j);
id diff(j?)*la(i?)=D(j,a)*B(i,j)*(la(i)-la(j))+ D(la(i))*diff(j);
id diff(i?)*f(j?)=f(j+a)*D(i,a) +f(j+b)*D(i,a)*mu(i) +f(j+c)*D(i,a)*la(i)
   +  D(f(j))*diff(i);
id diff(i?)*[1/f(2*b)] = - [1/f(2*b)]^2*(f(2*b+a)*D(i,a) 
   +f(2*b+b)*D(i,a)*mu(i) +f(2*b+c)*D(i,a)*la(i)) +  D([1/f(2*b)])*diff(i);
endrepeat;
id B(i?,j?)=Bij;
id diff(?a)=0;
id D(a?)=a;
id D(a?,a)=1;
.sort
skip;
nskip pr, dkDispEq;
repeat; 
*** The following id just wraps the Hessian back to [Hess] (if it is present...)
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
*** this will force cancellation of \lambda's in numerator/denominator
id la(i)^2*[1/lji](i+j)^2 = 1 - [1/lji](i+j)^2*(-2*la(i)*la(j) + la(j)^2);
id [Hess]*[1/Hess]=1;
endrepeat;
Bracket [1/Hess],[1/lji];
print pr, dkDispEq;
.sort

*** for convenience we introduce the following short names for
*** all other B_{??} we need
L Bik=B(i,k);
L Bkj=B(k,j);
L Bji=B(j,i);
L Bki=B(k,i);
L Bjk=B(j,k);
repeat;
**** Since in Bij,Nij,... (r.h.s. of the folowing id statements) i and j are
****  EXPLICITELY present inside mu(.), la(.) etc. - it is possible
****  to use the r.h.s. as a universal definition for arbitrary indices.
****  Please note that i,j,k are NOT declared as FORM indices - they
****  are LOGICAL indices for us.
id B(i?,j?)=Bij;
endrepeat;
.sort

************  below we form the terms of the sum
** - d_k N_{ij} + N_{ij}( 1/D_{ij} d_k D_{ij} + B_{kj} + B_{ik} ) - D_{ij}B_{kj}B_{ik}
*** one by one in the course of many FORM modules...
*** we will denote this complete sum later as "all"


**** first comes the last term :-) 
**** we split forming of this term D_{ij}B_{kj}B_{ik}
**** into two steps for efficiency reasons
skip;
L bb1=Dij*Bkj;
.sort
skip;
drop bb1;
L nbb=bb1*Bik;
.sort

**** here we check that nbb has the expected denominator,
**** later in will partially cancel out
skip;
nskip nbb;
Brackets [1/lji],[1/Hess],[Hess];
print[] nbb;
.sort

*** simplify module DispEq
skip;
nskip nbb;
repeat;
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

***  The following change is necessary: until this point we enjoyed the fact
***  that the denominators had (\lambda_{?}-\lambda_{?}) SQUARED
***  so we could easily write this as [1/lji](?+?)^2,
***  but now we will try to cancel gradually such denominators with the numerators
***  step by step, so we FIX:
*****          [1/x] = (\la^i-\la^j)^{-1} = [1/lji](i+j)
*****          [1/y] = (\la^j-\la^k)^{-1} = [1/lji](k+j)
*****          [1/z] = (\la^k-\la^i)^{-1} = [1/lji](i+k)
skip;
nskip nbb;
repeat;
id [1/lji](j+i)=[1/x];
id [1/lji](k+i)=[1/z];
id [1/lji](j+k)=[1/y];
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] nbb;
.sort

*** now one of the tricks: we cancel some of the 
*** (\lambda_{?}-\lambda_{?}) in the denominators with that in the numerators
skip;
nskip nbb;
repeat;
id la(j)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]+[1/z];
id la(i)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]-[1/y];

id la(k)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]-[1/x];
id la(i)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]+[1/y];

id la(j)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]-[1/z];
id la(k)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]+[1/x];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
*** now you can see how this works!
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] nbb;
.sort

*** a similar trick with the Hessian
skip;
nskip nbb;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
*** now you can see how this works!
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] nbb;
.sort

*** form the derivatives d_k N_{ij} and d_k D_{ij}
skip;
L dd=diff(k)*Dij;
L dn=diff(k)*Nij;
repeat; 
id diff(j?)*mu(i?)=D(j,a)*B(i,j)*(mu(i)-mu(j))+ D(mu(i))*diff(j);
id diff(j?)*la(i?)=D(j,a)*B(i,j)*(la(i)-la(j))+ D(la(i))*diff(j);
id diff(i?)*f(j?)=f(j+a)*D(i,a) +f(j+b)*D(i,a)*mu(i) +f(j+c)*D(i,a)*la(i)
   +  D(f(j))*diff(i);
id diff(i?)*[1/f(2*b)] = - [1/f(2*b)]^2*(f(2*b+a)*D(i,a) 
   +f(2*b+b)*D(i,a)*mu(i) +f(2*b+c)*D(i,a)*la(i)) +  D([1/f(2*b)])*diff(i);
id B(i?,j?)=Bij;
endrepeat;
id diff(?a)=0;
id D(a?)=a;
id D(a?,a)=1;
.sort

repeat;
id [1/lji](j+i)=[1/x];
id [1/lji](k+i)=[1/z];
id [1/lji](j+k)=[1/y];
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] dd,dn;
.sort

skip;
nskip dd,dn;
repeat;
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
*** just to see how many terms are left and how big are the denomiantors
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] dd,dn;
.sort

**** again we form the term N_{ij} 1/D_{ij} d_k D_{ij}
**** in two steps for efficiency 
skip;
L NUdD1=Nij*dd;
repeat;
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
.sort

drop NUdD1;
skip;
L NUdD=NUdD1*Uij;
repeat;
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
id [1/lji](j+i)=[1/x];
id [1/lji](k+i)=[1/z];
id [1/lji](j+k)=[1/y];
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] NUdD;
.sort

****  again simplify numerators/denominators
skip;
nskip NUdD;
repeat;
id la(j)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]+[1/z];
id la(i)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]-[1/y];

id la(k)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]-[1/x];
id la(i)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]+[1/y];

id la(j)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]-[1/z];
id la(k)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]+[1/x];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] NUdD;
.sort

skip;
nskip NUdD;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
*id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
*   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] NUdD;
.sort

**** finally in four modules we combine all terms to form the
**** complete sum "all"
skip;
L all1=Nij*(Bkj + Bik);
repeat; 
id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all1;
.sort

skip;
nskip all1;
repeat;
id la(j)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]+[1/z];
id la(i)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]-[1/y];

id la(k)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]-[1/x];
id la(i)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]+[1/y];

id la(j)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]-[1/z];
id la(k)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]+[1/x];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all1;
.sort


skip;
nskip all1;
repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
*id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
*   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all1;
.sort

*** here comes the complete sum, we drop all others as unnecessary
drop;
Global all= - dn + all1 - nbb + NUdD;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all;
.sort

**** again cancellations in numerator/denominator
repeat;
id la(j)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]+[1/z];
id la(i)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]-[1/y];

id la(k)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]-[1/x];
id la(i)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]+[1/y];

id la(j)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]-[1/z];
id la(k)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]+[1/x];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all;
.sort

repeat; 
id [Hess]*[1/Hess]=1;
id f(a + b)^2*f(2*c) = -( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
*id mu(i?)^2=-(f(a+a)+2*f(a+b)*mu(i)+2*f(a+c)*la(i)+f(c+c)*la(i)^2
*   +2*f(c+b)*mu(i)*la(i))*[1/f(2*b)];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all;
.sort

*** we shall use this identity to get rid of the products
*** of different (\lambda_{?}-\lambda_{?}) in the denominators,
*** the two following modules result in the sum "all"
*** which has NOT MORE THAN ONE of [1/x],[1/y],[1/z] for each term
repeat;
id [1/x]*[1/y]=-[1/x]*[1/z]-[1/z]*[1/y];
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess], [Hess];
print[] all;
.sort
repeat;
id la(j)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]+[1/z];
id la(i)*[1/z]*[1/y]=la(k)*[1/z]*[1/y]-[1/y];
id la(k)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]-[1/x];
id la(i)*[1/x]*[1/y]=la(j)*[1/x]*[1/y]+[1/y];
id la(j)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]-[1/z];
id la(k)*[1/z]*[1/x]=la(i)*[1/z]*[1/x]+[1/x];
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all;
.sort

**** finally these two modules COMPLETELY REMOVE [1/x], [1/y], [1/z]
**** from the denominators and only [1/Hess] (not [1/Hess]^2 as in the beginning)
**** is left in the denominators
repeat;
id f(a + b)^2*f(2*c)= - ( - f(2*a)*f(2*b)*f(2*c) + f(2*b)*f(a + c)^2
     - 2*f(a + b)*f(a + c)*f(b + c) + f(2*a)*f(b + c)^2) - [Hess];
id [Hess]*[1/Hess]=1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess], [Hess];
print[] all;
.sort
repeat;
id la(j)*[1/y]=la(k)*[1/y]+1;
id la(i)*[1/x]=la(j)*[1/x]+1;
id la(i)*[1/z]=la(k)*[1/z]-1;
id f(2*b)*[1/f(2*b)]=1;
endrepeat;
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all;
.sort

**** We print the resulting sum to feed it afterwards
**** to Maple and save the sum (in some internal FORM format)
**** in the file all.sav
Brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print all;
.sort

brackets [1/x],[1/y],[1/z],[1/Hess],[Hess];
print[] all;
.store

save all.sav all;
.end

*******************   THE END  ******************

