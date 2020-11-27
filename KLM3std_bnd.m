function [K,A,B,C,iter]=KLM3std_bnd(T,L,bsen,m,numit,K,A,B,C)
%
% Krylov-Levenberg-Marquardt method for Tucker decomposition of order-3
% tensors. Input:
%
%    T... tensor to be decomposed
%    L .. binary tensor indicating zeros and nonzeros of the desired
%          kernel tensor K
%    bsen .. bound on sensitivity
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    K .... core tensor
%    A,B,C .... initial estimate of factor matrices (if available)
%
% Programmed by Petr Tichavsky, November 2020
%
[Ia,Ib,Ic]=size(T); 
Iabc=Ia+Ib+Ic;
RR=size(L);
fx=find(L(:)==1);
lx=length(fx); rx=prod(RR);
R1=RR(1); R2=RR(2); R3=RR(3);
ia1=lx; ib1=ia1+R1*Ia; ic1=ib1+R2*Ib; id1=ic1+R3*Ic;
if nargin<9
    A=randn(Ia,R1);
    B=randn(Ib,R2);
    C=randn(Ic,R3);
    K=zeros(R1,R2,R3);
    K(fx)=randn(lx,1);
    sen=sensitSTD(L,K,A,B,C);
    fac=(bsen/sen)^(1/6);
    K=K*fac; A=A*fac; B=B*fac; C=C*fac;
  %  [bsen sensitSTD(L,K,A,B,C)]
end    
%Y1=reshape(T,Ia,Ib*Ic);
%Y2=reshape(permute(T,[2,1,3]),Ib,Ia*Ic);
%Y3=reshape(permute(T,[3,1,2]),Ic,Ia*Ib);
iter=zeros(1,numit);
% initial parameter mu
nu=2;
err=chybaT(T,K,A,B,C);              %%% computes the error of the current approximation
for it=1:numit
%   AA=A'*A;
%   BB=B'*B;
%   CC=C'*C;
  E=T-multi(K,A,B,C);
  gK=multi(E,A',B',C');
  gA=reshape(E,Ia,Ib*Ic)*reshape(multi(K,eye(R1),B,C),R1,Ib*Ic)';
  gB=reshape(permute(E,[2,3,1]),Ib,Ia*Ic)*reshape(permute(multi(K,A,eye(R2),C),[2,3,1]),R2,Ia*Ic)';
  gC=reshape(permute(E,[3,1,2]),Ic,Ia*Ib)*reshape(permute(multi(K,A,B,eye(R3)),[3,1,2]),R3,Ia*Ib)';
%   J=dejtJ(K,A,B,C);
%   J'*E(:);
  g=[gK(fx); gA(:); gB(:); gC(:)];
  Y=zeros(id1,m+2); %  Krylov subspace
  ng=norm(g);
  Y(:,1)=g/ng;
  W=zeros(m+1,m+1); 
  zz=g/ng;
  for i=1:m+1  
    zz1=[zeros(rx,1); zz(lx+1:id1)];
    zz1(fx)=zz(1:lx);
    aux=Hx(K,A,B,C,zz1);
    zz=[aux(fx); aux(rx+1:end)];
    W(1:i,i)=Y(:,1:i)'*zz;
    Y(:,i+1)=zz-Y(:,1:i)*W(1:i,i);
    Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
    zz=Y(:,i+1);
end
W=W+W'-diag(diag(W));
  Z=zeros(id1,m+2); %  The second Krylov subspace
  u=GRADsenSTD(L,K,A,B,C); 
    X=zeros(m+1,m+1);
    zz=u/nu;
    for i=1:m+1
        zz1=[zeros(rx,1); zz(lx+1:id1)];
        zz1(fx)=zz(1:lx);
        aux=Hx(K,A,B,C,zz1);
        zz=[aux(fx); aux(rx+1:end)];
        X(1:i,i)=Z(:,1:i)'*zz;
        Z(:,i+1)=zz-Z(:,1:i)*X(1:i,i);
        Z(:,i+1)=Z(:,i+1)/norm(Z(:,i+1));
        zz=Z(:,i+1);
    end
    X=X+X'-diag(diag(X));
if it==1
   mu=max(diag(W))/2;
end   
    d1=g/mu-Y(:,1:m+1)*(W*((eye(m+1)+W/mu)\(Y(:,1:m+1)'*g)))/mu^2;
    d2=u/mu-Z(:,1:m+1)*(X*((eye(m+1)+X/mu)\(Z(:,1:m+1)'*u)))/mu^2;
    lam=u'*d1/(u'*d2);
    d=d1-lam*d2;
A1=A; B1=B; C1=C; K1=K;
K(fx)=K(fx)+d(1:ia1); 
A=A+reshape(d(ia1+1:ib1),Ia,R1); 
B=B+reshape(d(1+ib1:ic1),Ib,R2); 
C=C+reshape(d(1+ic1:id1),Ic,R3); 
sen=sensitSTD(L,K,A,B,C);
fac=(bsen/sen)^(1/6);
K=K*fac; A=A*fac; B=B*fac; C=C*fac;
err2=chybaT(T,K,A,B,C);
if err2>err                                     %%% step is not accepted
      mu=mu*nu; nu=2*nu; err2=err;    
      A=A1; B=B1; C=C1; K=K1;                   %%% return to the previous solution
   else
      rho=(err-err2)/(d'*(g+mu*d));
      err = err2;
      nu=2;
      mu=mu*max([1/3 1-(2*rho-1)^3]);
end   
iter(it)=err2;
if rem(it,10)==0
       [it -log10(err2)]  % to monitor decrease of the cost function each 10 iterations
end    
end
semilogy(iter)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chybaT(X,K,A,B,C)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
Y=multi(K,A,B,C);
err=sum((Y(:)-X(:)).^2);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function AB = krb(A,B)
%KRB Khatri-Rao product
%
[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in krb.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end
end
%
function y=Hx(K,A,B,C,x)
%
% to compute the matrix-vector product y=Hx for Tucker tensor decomposition
%
%[Ra,Rb,Rc]=size(K);
[Ia,Ra]=size(A);
[Ib,Rb]=size(B);
[Ic,Rc]=size(C);
%R=Ra;
ia1=Ra*Rb*Rc;
ib1=ia1+Ia*Ra;
ic1=ib1+Ib*Rb;
XK=reshape(x(1:Ra*Rb*Rc),Ra,Rb,Rc);
XA=reshape(x(ia1+1:ia1+Ia*Ra),Ia,Ra);
XB=reshape(x(ib1+1:ib1+Ib*Rb),Ib,Rb);
XC=reshape(x(ic1+1:ic1+Ic*Rc),Ic,Rc);
AA=A'*A; BB=B'*B; CC=C'*C;
YK=multi(XK,AA,BB,CC)+multi(K,A'*XA,BB,CC)+multi(K,AA,B'*XB,CC)+multi(K,AA,BB,C'*XC);
K2A=multi(K,eye(Ra),BB,CC);
K2B=multi(K,AA,eye(Rb),CC);
K2C=multi(K,AA,BB,eye(Rc));
K3A=multi(K,eye(Ra),XB'*B,CC);
K4A=multi(K,eye(Ra),BB,XC'*C);
%
K3B=multi(K,XA'*A,eye(Rb),CC);
K4B=multi(K,AA,eye(Rb),XC'*C);
%
K3C=multi(K,XA'*A,BB,eye(Rc));
K4C=multi(K,AA,XB'*B,eye(Rc));
%
YA=A*(reshape(XK,Ra,Rb*Rc)*reshape(K2A,Ra,Rb*Rc)')+XA*(reshape(K2A,Ra,Rb*Rc)*reshape(K,Ra,Rb*Rc)')...
+A*(reshape(K,Ra,Rb*Rc)*reshape(K3A+K4A,Ra,Rb*Rc)'); %+A*(reshape(K,R,R^2)*reshape(K4A,R,R^2)');
%
YB=B*(reshape(permute(XK,[2,3,1]),Rb,Ra*Rc)*reshape(permute(K2B,[2,3,1]),Rb,Ra*Rc)')+XB*(reshape(permute(K2B,[2,3,1]),Rb,Ra*Rc)*reshape(permute(K,[2,3,1]),Rb,Ra*Rc)')...
+B*(reshape(permute(K,[2,3,1]),Rb,Ra*Rc)*reshape(permute(K3B+K4B,[2,3,1]),Rb,Ra*Rc)');
%
YC=C*(reshape(permute(XK,[3,1,2]),Rc,Ra*Rb)*reshape(permute(K2C,[3,1,2]),Rc,Ra*Rb)')+XC*(reshape(permute(K2C,[3,1,2]),Rc,Ra*Rb)*reshape(permute(K,[3,1,2]),Rc,Ra*Rb)')...
+C*(reshape(permute(K,[3,1,2]),Rc,Ra*Rb)*reshape(permute(K3C+K4C,[3,1,2]),Rc,Ra*Rb)');
y=[YK(:); YA(:); YB(:); YC(:)];
end
%
function T2=multi(T,A,B,C)
[n1 n2 n3]=size(T);
r1=size(A,1); r2=size(B,1); r3=size(C,1);
T1=A*reshape(T,n1,n2*n3);
T2=reshape(T1,r1*n2,n3)*C.';
T2=permute(reshape(B*reshape(permute(reshape(T2,r1,n2,r3),[2 1 3]),n2,r1*r3),r2,r1,r3),[2 1 3]);
end