function [K,A,B,C,iter]=KLM3std(T,L,m,numit,K,A,B,C,tol)
%
% Krylov-Levenberg-Marquardt method for Tucker decomposition of order-3
% tensors. Input:
%
%    T... tensor to be decomposed
%    L .. binary tensor indicationg zeros and nonzeros of the desired
%          kernel tensor K
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    K .... core tensor
%    A,B,C .... initial estimate of factor matrices (if available)
%    tol ...... tolerance parameter to stop the algorithm if the cost
%    function does not decrease significantly
%
% Programmed by Petr Tichavsky, August 2020
%
[Ia,Ib,Ic]=size(T);
Iabc=Ia+Ib+Ic;
RR=size(L);
fx=find(L(:)==1);
lx=length(fx); rx=prod(RR);
R1=RR(1); R2=RR(2); R3=RR(3);
ia1=lx; ib1=ia1+R1*Ia; ic1=ib1+R2*Ib; id1=ic1+R3*Ic;
if nargin<8
    A=randn(Ia,R1);
    B=randn(Ib,R2);
    C=randn(Ic,R3);
    K=zeros(R1,R2,R3);
    K(fx)=randn(lx,1);
end
if nargin<9
    tol = 1e-5;
end
iter=zeros(1,numit);
na=sum(A.^2); nb=sum(B.^2); nc=sum(C.^2);
mu=max(K(:).^2)*max(max(nb)*[max(nc),max(na)+max(nc)*max(nb)]);
% initial parameter mu
nu=2;
err=chybaT(T,K,A,B,C)              %%% computes the error of the current approximation
for it=1:numit
    E=T-multi(K,A,B,C);
    gK=multi(E,A',B',C');
    gA=reshape(E,Ia,Ib*Ic)*reshape(multi(K,[],B,C),R1,Ib*Ic)';
    gB=reshape(permute(E,[2,3,1]),Ib,Ia*Ic)*reshape(permute(multi(K,A,[],C),[2,3,1]),R2,Ia*Ic)';
    gC=reshape(permute(E,[3,1,2]),Ic,Ia*Ib)*reshape(permute(multi(K,A,B,[]),[3,1,2]),R3,Ia*Ib)';
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
    d=g/mu-Y(:,1:m+1)/(inv(W)+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2;
    A1=A; B1=B; C1=C; K1=K;
    K(fx)=K(fx)+d(1:ia1);
    A=A+reshape(d(ia1+1:ib1),Ia,R1);
    B=B+reshape(d(1+ib1:ic1),Ib,R2);
    C=C+reshape(d(1+ic1:id1),Ic,R3);
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
        fprintf('%d  %d\n',[it err2])  % to monitor decrease of the cost function each 10 iterations
    end
    
     if (it)>5 && ((abs(iter(it)-iter(it-1))<tol*iter(it-1)) && std(abs(diff(iter(it-4:it))))<tol)
         break
     end
end
iter = iter(1:it);
%semilogy(iter)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Z=multi(XK,A,B,C)+multi(K,XA,B,C)+multi(K,A,XB,C)+multi(K,A,B,XC);
YK=multi(Z,A',B',C');
YA=reshape(Z,Ia,Ib*Ic)*reshape(multi(K,[],B,C),Ra,Ib*Ic)';
YB=reshape(permute(Z,[2,3,1]),Ib,Ia*Ic)*reshape(permute(multi(K,A,[],C),[2,3,1]),Rb,Ia*Ic)';
YC=reshape(permute(Z,[3,1,2]),Ic,Ia*Ib)*reshape(permute(multi(K,A,B,[]),[3,1,2]),Rc,Ia*Ib)';
y=[YK(:); YA(:); YB(:); YC(:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T2=multi(T,A,B,C)
T2 = T;
szT = size(T);
if ~isempty(A)
    szT(1)=size(A,1);
    T2=A*reshape(T2,[],szT(2)*szT(3));
end

if ~isempty(C)
    szT(3) = size(C,1);
    T2=reshape(T2,szT(1)*szT(2),[])*C.';
end

if ~isempty(B)
    T2 = reshape(T2,szT);    
    szT(2)= size(B,1);    
    T2=permute(reshape(B*reshape(permute(T2,[2 1 3]),[],szT(1)*szT(3)),szT([2 1 3])),[2 1 3]);    
end
T2 = reshape(T2,szT);    
end