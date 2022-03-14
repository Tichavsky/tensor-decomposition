function [A,iter]=KLMtc_sym(T,R,m,numit,A,tol)
%
% Krylov-Levenberg-Marquardt method for symmetric tensor chain decomposition
% of symmetric tensors
%
% Input:
%
%    T... tensor to be decomposed
%    R .... desired dimensions of core tensors: (R,sz,R)
%    m ... parameter controlling complexity of one iteration of KLM
%          i.e., dimension of the Krylov subspace
%    numit .... required number of iterations
%    A .... initial estimate of the core tensors (if available)
%    tol ...... tolerance parameter to stop the algorithm if the cost
%               function does not decrease significantly
%
% Programmed by Petr Tichavsky, January 2021
%
sz=size(T);
N=length(sz);
szn=prod(sz);                  %% number of the tensor elements
sz=sz(1);
if nargin<5
  A=randn(R,sz,R); %% random initialization
end
nel=sz*R^2;  %% total number of parameters
if nargin<6
    tol = 1e-8;
end
iter=zeros(1,numit);
nu=2;
err=chybaTCsym(N,T,A);       %%% computes the error of the current approximation
it=0;
while it<numit && err>1e-9
    it=it+1;
    E=T-fullTR_sym(N,A);     %%% error of the approximation (tensor)
    aux=(reshape(E,sz,sz^(N-1)))'; 
    g=reshape(aux,sz,sz^(N-1))*Unn(N,A)';
    g=g(:);
    Y=zeros(nel,m+2); %  Krylov subspace
    ng=norm(g);
    Y(:,1)=g/ng;
    W=zeros(m+1,m+1);
    zz=g/ng;
    for i=1:m+1
        zz=Hx_sym(N,A,zz,sz,R);
        W(1:i,i)=Y(:,1:i)'*zz;
        Y(:,i+1)=zz-Y(:,1:i)*W(1:i,i);
        Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
        zz=Y(:,i+1);
    end
    W=W+W'-diag(diag(W));
    if it==1
       mu=min(diag(W))/2;
    end   
%    d1=g/mu-Y(:,1:m+1)/(inv(W)+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2;
    d=g/mu-Y(:,1:m+1)*(W*((eye(m+1)+W/mu)\(Y(:,1:m+1)'*g)))/mu^2;
    A1=A; 
    aux=reshape(d,sz,R,R);
    A=A+permute(aux,[2,1,3]);
    err2=chybaTCsym(N,T,A);
    if err2>err                                     %%% step is not accepted
        mu=mu*nu; nu=2*nu; err2=err;
        A=A1;          %%% return to the previous solution
    else
        rho=(err-err2)/(d'*(g+mu*d));
        err = err2;
        nu=2;
        mu=mu*max([1/3 1-(2*rho-1)^3]);
    end
    iter(it)=err2;

  %  if rem(it,10)==0
        fprintf('%d  %d\n',[it err2])  % to monitor decrease of the cost function each 10 iterations
  %  end
     if err<1e-10 
         break
     end
 %    if (it)>5 && ((abs(iter(it)-iter(it-1))<tol*iter(it-1)) && std(abs(diff(iter(it-4:it))))<tol)
 %        break
 %    end
end
iter = iter(1:it);
%semilogy(iter)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chybaTCsym(N,X,A)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
Y=fullTR_sym(N,A);
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
function y=Hx_sym(N,A,x,sz,R)
%
% to compute the matrix-vector product y=Hx for Tucker tensor decomposition
%
%[Ra,Rb,Rc]=size(K);
ind=0;
Z=zeros(sz);
szn=prod(sz);
U=Unn(N,A);
Z0=reshape(x,sz,R^2)*U;
Z=Z0;
for n=1:N-1
    Z0=reshape(Z0',sz,sz^(N-1));
    Z=Z+Z0;
end   
y=Z*U';
y=y(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=Unn(d,tr)
%
% Computes the matrix Un
%
a = tr;    
[r0, n, r1] = size(a);
ns = zeros(1,d);
ns(1) = 1;

for k = 2:d-1
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr;
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
[r0, n, r1] = size(a);
a = reshape(permute(a,[3,1,2]), r0*r1, numel(a)/(r0*r1));
 
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [res]=fullTR_sym(d,tr)
%
%%% Converts a TR-representation into full format.
%%% Algorithm: multiply together the reshaped cores, then make this into
%%% an \alpha_0 , . , \alpha_0 tensor and sum over alpha_0 in a loop. Then
%%% reshape this into an n1 x ... x nd tensor.
%%% modified by Petr Tichavsky, October 2020

a = tr;    
[r0, n, r1] = size(a);
ns = zeros(1,d);
ns(2) = n;

for k = 3:d
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr;
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
b = tr;
[rold, n, rnew] = size(b);
ns(1)=n;
rr=rold*rnew;
b = reshape(permute(b,[2,1,3]),n,rr);
a = reshape(permute(a,[3,1,2]), rr, numel(a)/rr);
res=reshape(b*a,ns);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
