function [A,iter]=KLMtc(T,RR,m,numit,A,tol)
%
% Krylov-Levenberg-Marquardt method for tensor chain decomposition
% Input:
%
%    T... tensor to be decomposed
%    RR .... desired dimensions of core tensors
%    m ... parameter controlling complexity of one iteration of KLM
%          i.e., dimension of the Krylov subspace
%    numit .... required number of iterations
%    A .... initial estimate of the core tensors (if available)
%    tol ...... tolerance parameter to stop the algorithm if the cost
%               function does not decrease significantly
%
% Programmed by Petr Tichavsky, October 2020, revised February 2022.
%
sz=size(T);
N=length(sz);
szn=prod(sz);                  %% number of the tensor elements
if length(RR)==1
    RR=RR*ones(1,N);
end    
RR(N+1)=RR(1);
if nargin<5
  for n=1:N
      A{n}=randn(RR(n),sz(n),RR(n+1)); %% random initialization
  end 
  A=ALStc(T,RR,1,A);  %%% one iteration of ALS to initialize
  A=baltr(A);         %%% balance magnitudes of the core tensors
end
nel=sum(sz.*RR(1:N).*RR(2:N+1));  %% total number of parameters
if nargin<6
    tol = 1e-9;
end
iter=zeros(1,numit);
nu=2;
err=chybaTC(T,A);       %%% computes the error of the current approximation
for it=1:numit
    E=T-fullTR(A);     %%% error of the approximation (tensor)
    ind=0;
    g=zeros(nel,1);    %%% will be the error gradient
    for n=1:N
   %     aux=permute(E,[n:N,1:n-1]);
        aux=(reshape(E,prod(sz(1:n-1)),prod(sz(n:N))))'; 
        gA=reshape(aux,sz(n),szn/sz(n))*Unn(n,A)';
        nga=numel(gA);
        g(ind+1:ind+nga)=gA(:);
        ind=ind+nga;
    end
    Y=zeros(nel,m+2); %  Krylov subspace
    ng=norm(g);
    Y(:,1)=g/ng;
    W=zeros(m+1,m+1);
    zz=g/ng;
    for i=1:m+1
        zz=Hx2(A,zz);
    %    zz=Hx(A,zz,sz,RR);
        W(1:i,i)=Y(:,1:i)'*zz;
        Y(:,i+1)=zz-Y(:,1:i)*W(1:i,i);
        Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
        zz=Y(:,i+1);
    end
    W=W+W'-diag(diag(W));
    if it==1
       mu=min(diag(W))/10;
    end   
%    d1=g/mu-Y(:,1:m+1)/(inv(W)+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2;
    d=g/mu-Y(:,1:m+1)*(W*((eye(m+1)+W/mu)\(Y(:,1:m+1)'*g)))/mu^2;
    A1=A; ind=0;
    for n=1:N
        nga=numel(A{n});
        aux=reshape(d(ind+1:ind+nga,1),sz(n),RR(n),RR(n+1));
        A{n}=A{n}+permute(aux,[2,1,3]);
        ind=ind+nga;
    end
    err2=chybaTC(T,A);
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

    if rem(it,10)==0
        fprintf('%d  %d\n',[it err2])  % to monitor decrease of the cost function each 10 iterations
    end
     if err<1e-10 
         break
     end
  %   if (it)>5 && ((abs(iter(it)-iter(it-1))<tol*iter(it-1)) && std(abs(diff(iter(it-4:it))))<tol)
  %       break
  %   end
end
iter = iter(1:it);
%semilogy(iter)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chybaTC(X,A)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
Y=fullTR(A);
err=sum((Y(:)-X(:)).^2);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function y=Hx(A,x,sz,RR)
%
% to compute the matrix-vector product y=Hx for Tucker tensor decomposition
%
%[Ra,Rb,Rc]=size(K);
N=length(A);
ind=0;
Z=zeros(sz);
szn=prod(sz);
for n=1:N
  %  sza=size(A{n});
    nga=numel(A{n});
    A0n=A{n};
    aux=reshape(x(ind+1:ind+nga,1),sz(n),RR(n),RR(n+1));
    A{n}=permute(aux,[2,1,3]);
   % A{n}=reshape(x(ind+1:ind+nga),sza);
    Z=Z+fullTR(A);
    A{n}=A0n;
    ind=ind+nga;
end   
y=zeros(size(x));
ind=0;
for n=1:N
  %  aux=permute(Z,[n:N,1:n-1]);
    aux=(reshape(Z,prod(sz(1:n-1)),prod(sz(n:N))))';
    yA=reshape(aux,sz(n),szn/sz(n))*Unn(n,A)';
    nga=numel(yA);
    y(ind+1:ind+nga)=yA(:);
    ind=ind+nga;
end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=Unn(n,tr)
%
% Computes the matrix Un
%
d=length(tr);
i2=[1:d,1:d];
in=i2(n+1:n+d-1);
a = tr{in(1)};    
[~, n, ~] = size(a);
ns = zeros(1,d);
ns(1) = n;

for k = 2:d-1
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr{in(k)};
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
[r0, ~, r1] = size(a);
a = reshape(permute(a,[3,1,2]), r0*r1, numel(a)/(r0*r1));
 
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [res]=fullTR(tr)
%
%%% Converts a TR-representation into full format.
%%% Algorithm: multiply together the reshaped cores, then make this into
%%% an \alpha_0 , . , \alpha_0 tensor and sum over alpha_0 in a loop. Then
%%% reshape this into an n1 x ... x nd tensor.
%%% modified by Petr Tichavsky, October 2020

d=length(tr);
a = tr{2};    
[r0, n, r1] = size(a);
ns = zeros(1,d);
ns(2) = n;

for k = 3:d
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr{k};
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
b = tr{1};
[rold, n, rnew] = size(b);
ns(1)=n;
rr=rold*rnew;
b = reshape(permute(b,[2,1,3]),n,rr);
a = reshape(permute(a,[3,1,2]), rr, numel(a)/rr);
res=reshape(b*a,ns);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=baltr(A)
%
% balance core factors of the chain decomposition of a tensor
%
N=length(A);
val=zeros(N,1);
for n=1:N
    a=A{n};
    val(n,1)=sum(a(:).^2)/numel(a);
end
val=sqrt(val);
mval=(prod(val))^(1/N);
for n=1:N
    A{n}=(mval/val(n,1))*A{n};
end
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function y=Hx2(A,x)
%
% to compute the matrix-vector product y=Hx for TC tensor decomposition
%
N=length(A);
X=A; y=x;
ind=0;
for n=1:N
    nga=numel(A{n});
    sz=size(A{n});
    if length(sz)<3
        sz=[sz,1];
    end    
    X{n}=reshape(x(ind+1:ind+nga,1),sz([2,1,3]));
    ind=ind+nga;
end
i2=[1:N,1:N];
ind=0;
for n=1:N
    in=i2(n+1:n+N);
    a = A{in(1)}; 
    b = X{in(1)}; 
    [r0,ns,r1]=size(a);
    a=reshape(permute(a,[1,3,2]),r0*r1,ns);
    b=reshape(permute(b,[2,3,1]),r0*r1,ns);
    aa=reshape(permute(reshape(a*a',r0,r1,r0,r1),[1,3,2,4]),r0*r0,r1*r1);
    ax=reshape(permute(reshape(a*b',r0,r1,r0,r1),[1,3,2,4]),r0*r0,r1*r1);
    for k = 2:N-1
        c = A{in(k)};
        [r2, ns, r3] = size(c);
        c = reshape(permute(c,[1,3,2]),r2*r3, ns);
        cc=reshape(permute(reshape(c*c',r2,r3,r2,r3),[1,3,2,4]),r2*r2,r3*r3);
        d = X{in(k)};
        d = reshape(permute(d,[2,3,1]),r2*r3, ns);
        cx=reshape(permute(reshape(c*d',r2,r3,r2,r3),[1,3,2,4]),r2*r2,r3*r3);
        ax = ax*cc+aa*cx;
        aa = aa*cc;
    end 
    sz=size(A{in(N)});
    if length(sz)<3
        sz=[sz,1];
    end   
    aux=reshape(permute(reshape(aa,r0,r0,r3,r3),[3,1,4,2]),r3*r0,r3*r0)*reshape(permute(X{in(N)},[2,3,1]),sz(1)*sz(3),sz(2))...
    +reshape(permute(reshape(ax,r0,r0,r3,r3),[3,1,4,2]),r3*r0,r3*r0)*reshape(permute(A{in(N)},[1,3,2]),sz(1)*sz(3),sz(2));
    aux=aux';
    nga=prod(sz);
    y(ind+1:ind+nga,1)=aux(:);
    ind=ind+nga;
end
end