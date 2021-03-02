function [A,B,C,iter]=KLM3cp_bnd(T,R,bsen,m,numit,A,B,C)
%
% Krylov-Levenberg-Marquardt method for Tucker decomposition of order-3
% tensors (complex version). Input:
%
%    T... tensor to be decomposed
%    L .. binary tensor indicating zeros and nonzeros of the desired
%          kernel tensor K
%    bsen .. bound on sensitivity, computed in sensitCP(A,B,C) as
%          sum_r I_3 |a_r|^2|b_r|^2+I_2 |a_r|^2|c_r|^2+I_1 |c_r|^2|b_r|^2
%
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    K .... core tensor
%    A,B,C .... initial estimate of factor matrices (if available)
%
% Programmed by Petr Tichavsky, November 2020
% Modified for complex-valued data by Anh-Huy Phan, February 2021
%
[Ia,Ib,Ic]=size(T);
Iabc=Ia+Ib+Ic;
ia1=0; ib1=ia1+R*Ia; ic1=ib1+R*Ib; id1=ic1+R*Ic;
%K=zeros(R,R,R);
%K(1:1+R+R^2:R^3)=1;
if nargin<8
    i=sqrt(-1);
    A=randn(Ia,R);
    B=randn(Ib,R);
    C=randn(Ic,R);
end
sen=sensitCP(A,B,C);
fac=(bsen/sen)^(1/4);
A=A*fac; B=B*fac; C=C*fac;
iter=zeros(1,numit);
% initial parameter mu
nu=2;
err=chybaCP(T,A,B,C);              %%% computes the error of the current approximation
for it=1:numit
    E=T-reshape(A*krb(C,B).',Ia,Ib,Ic);
 %   gA=reshape(E,Ia,Ib*Ic)*reshape(multi2(K,[],B,C),R,Ib*Ic)';
 %   gB=reshape(permute(E,[2,3,1]),Ib,Ia*Ic)*reshape(permute(multi2(K,A,[],C),[2,3,1]),R,Ia*Ic)';
 %   gC=reshape(permute(E,[3,1,2]),Ic,Ia*Ib)*reshape(permute(multi2(K,A,B,[]),[3,1,2]),R,Ia*Ib)';
    gA=reshape(E,Ia,Ib*Ic)*conj(krb(C,B));
    gB=reshape(permute(E,[2,3,1]),Ib,Ia*Ic)*conj(krb(A,C));
    gC=reshape(permute(E,[3,1,2]),Ic,Ia*Ib)*conj(krb(B,A));
    g=[gA(:); gB(:); gC(:)];
    Y=zeros(id1,m+2); %  Krylov subspace
    ng=norm(g);
    Y(:,1)=g/ng;
    W=zeros(m+1,m+1);
    zz=g/ng;
    for i=1:m+1
        zz=Hx(A,B,C,zz);
        W(1:i,i)=Y(:,1:i)'*zz;
        Y(:,i+1)=zz-Y(:,1:i)*W(1:i,i);
        Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
        zz=Y(:,i+1);
    end
    W=W+W'-diag(diag(W));
    Z=zeros(id1,m+2); %  The second Krylov subspace
    u=GRADsenCP(A,B,C);
    nu=norm(u);  
    X=zeros(m+1,m+1);
    zz=u/nu;
    Z(:,1)=zz; 
    for i=1:m+1
        zz=Hx(A,B,C,zz);
        X(1:i,i)=Z(:,1:i)'*zz;
        Z(:,i+1)=zz-Z(:,1:i)*X(1:i,i);
        Z(:,i+1)=Z(:,i+1)/norm(Z(:,i+1));
        zz=Z(:,i+1);
    end
    X=X+X'-diag(diag(X));
    if it==1
        mu=max(real(diag(W)))/2;
    end
    d1=g/mu-Y(:,1:m+1)*(W*((eye(m+1)+W/mu)\(Y(:,1:m+1)'*g)))/mu^2;
    d2=u/mu-Z(:,1:m+1)*(X*((eye(m+1)+X/mu)\(Z(:,1:m+1)'*u)))/mu^2;
    lam=u'*d1/real(u'*d2);
    d=d1-lam*d2;
    A1=A; B1=B; C1=C; 
    A=A+reshape(d(ia1+1:ib1),Ia,R);
    B=B+reshape(d(1+ib1:ic1),Ib,R);
    C=C+reshape(d(1+ic1:id1),Ic,R);
    sen=sensitCP(A,B,C);
    fac=(bsen/sen)^(1/4);
    A=A*fac; B=B*fac; C=C*fac;
    err2=chybaCP(T,A,B,C);
    if err2>err                                     %%% step is not accepted
        mu=mu*nu; nu=2*nu; err2=err;
        A=A1; B=B1; C=C1;                  %%% return to the previous solution
    else
        rho=real((err-err2)/(d'*(g+mu*d)));
        err = err2;
        nu=2;
        mu=mu*max([1/3 1-(2*rho-1)^3]);
    end
    iter(it)=err2;
    if rem(it,10)==0
        [it err2]  % to monitor decrease of the cost function each 10 iterations
    end
end
semilogy(iter)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chybaCP(X,A,B,C)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
Y=A*krb(C,B).';
err=sum(abs(Y(:)-X(:)).^2); % corrected to work for complex-valued data
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
function y=Hx(A,B,C,x)
%
% to compute the matrix-vector product y=Hx for Tucker tensor decomposition
%
%[Ra,Rb,Rc]=size(K);
R=size(A,2);
%K=zeros(R,R,R);
%K(1:1+R+R^2:R^3)=1;
[Ia,Ra]=size(A);
[Ib,Rb]=size(B);
[Ic,Rc]=size(C);
%R=Ra;
ia1=0;
ib1=ia1+Ia*Ra;
ic1=ib1+Ib*Rb;
XA=reshape(x(ia1+1:ia1+Ia*Ra),Ia,Ra);
XB=reshape(x(ib1+1:ib1+Ib*Rb),Ib,Rb);
XC=reshape(x(ic1+1:ic1+Ic*Rc),Ic,Rc);
AA=A'*A; BB=B'*B; CC=C'*C;
YA=XA*conj(CC.*BB)+A*conj((XB'*B).*CC+(XC'*C).*BB);
YB=XB*conj(AA.*CC)+B*conj((XC'*C).*AA+(XA'*A).*CC);
YC=XC*conj(AA.*BB)+C*conj((XB'*B).*AA+(XA'*A).*BB);
y=[YA(:); YB(:); YC(:)];
end
%
function sen=sensitCP(A,B,C)
%
% Sensitivity of the structured Tucker tensor decomposition.
% Computed according to the expression in the paper
% "Krylov-Levenberg-Marquardt Algorithm for Structured Tucker Tensor Decompositions"
%  IEEE Journal of Selected Topics in Signal Processing, 2021
%
% Input:
%
%    K .... core tensor
%    A,B,C .... factor matrices
%
% Programmed by Petr Tichavsky, February 2021
%
Ia=size(A,1); Ib=size(B,1); Ic=size(C,1);
a1=sum(abs(A).^2); b1=sum(abs(B).^2); c1=sum(abs(C).^2);
sen=sum((Ib*a1+Ia*b1).*c1+Ic*a1.*b1);
end
%

function u=GRADsenCP(A,B,C)
%
% Gradient of sensitivity of the structured Tucker tensor decomposition.
% Input:
%
%    L  .. binary tensor indicationg zeros and nonzeros of the desired
%          kernel tensor K
%    K .... core tensor
%    A,B,C .... factor matrices
%
% Programmed by Petr Tichavsky, February 2021
%
Ia=size(A,1); Ib=size(B,1); Ic=size(C,1);
a1=sum(abs(A).^2); b1=sum(abs(B).^2); c1=sum(abs(C).^2);
u1=A.*repmat(Ic*b1+Ib*c1,Ia,1); u2=B.*repmat(Ic*a1+Ia*c1,Ib,1); u3=C.*repmat(Ia*b1+Ib*a1,Ic,1);
u=2*[u1(:); u2(:); u3(:)];
end

function T2=multi2(T,A,B,C)
T2 = T;
szT = size(T);
if ~isempty(A)
    szT(1)=size(A,1);
    T2=A*reshape(T2,[],szT(2)*szT(3));
end

if ~isempty(C)
    szT(3) = size(C,1);
    T2=reshape(T2,szT(1)*szT(2),[])*C.';% corrected to work for complex-valued data
end

if ~isempty(B)
    T2 = reshape(T2,szT);
    szT(2)= size(B,1);
    T2=permute(reshape(B*reshape(permute(T2,[2 1 3]),[],szT(1)*szT(3)),szT([2 1 3])),[2 1 3]);
end
T2 = reshape(T2,szT);
end