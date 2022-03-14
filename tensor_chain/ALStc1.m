function [A,iter]=ALStc1(T,RR,numit,A,tol)
%
% Alternating Least Squares method for tensor chain decomposition
% Version is specially suited for high tensor orders and ranks.
% Input:
%
%    T... tensor to be decomposed
%    RR .... desired dimensions of core tensors
%    numit .... required number of iterations
%    A .... initial estimate of the core tensors (if available)
%    tol ...... tolerance parameter to stop the algorithm if the cost
%               function does not decrease significantly
%
% Programmed by Petr Tichavsky, July 2021
%
sz=size(T);
N=length(sz);
if length(RR)==1
    RR=RR*ones(1,N);
end  
RR(N+1)=RR(1);
if nargin<4
  for n=1:N
      A{n}=randn(RR(n),sz(n),RR(n+1)); %% random initialization
  end    
end
% nel=sum(sz.*RR(1:N).*RR(2:N+1))  %% total number of parameters   
if nargin<5
    tol = 1e-10;
end
t0=sum(T(:).^2);
iter=zeros(1,numit);
for it=1:numit
    for n=1:N
        M=dejM(n,A);
        B=dejB(n,A,T);
   %    gA=(reshape(aux,sz(n),szn/sz(n))*aux2')/(aux2*aux2');
        gA=B/M;
        A{n}=permute(reshape(gA,sz(n),RR(n),RR(n+1)),[2,1,3]);
    end
    iter(it)=t0+sum(diag(gA*M*gA'))-2*sum(diag(gA'*B)); % =chybaTC(T,A);  
    it
    iter(it)
    if iter(it)<tol
         break
    end
    if (it)>5 && ((abs(iter(it)-iter(it-1))<tol*iter(it-1)) && std(abs(diff(iter(it-4:it))))<tol)
         break
    end
end
iter = iter(1:it);
%semilogy(iter)
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function a=dejM(n,tr)
%
% Computes the matrix product UnUn^T
%
d=length(tr);
i2=[1:d,1:d];
in=i2(n+1:n+d-1);
a = tr{in(1)};    
[r0, n, r1] = size(a);
a=reshape(permute(a,[1,3,2]),r0*r1,n);
a=reshape(permute(reshape(a*a',r0,r1,r0,r1),[1,3,2,4]),r0*r0,r1*r1);
for k = 2:d-1
    b = tr{in(k)};
    [r2, n, r3] = size(b);
    b = reshape(permute(b,[1,3,2]),r2*r3, n);
    b=reshape(permute(reshape(b*b',r2,r3,r2,r3),[1,3,2,4]),r2*r2,r3*r3);
    a = a*b;
end 
a=reshape(permute(reshape(a,r0,r0,r3,r3),[4,2,3,1]),r3*r0,r3*r0);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function a=dejB(n,tr,T)
%
% Computes the matrix product T_(n)Un^T
%
d=length(tr);
i2=[1:d,1:d];
in=i2(n+1:n+d);
a = tr{in(1)}; 
T=permute(T,[n+1:d,1:n]);
[rold, n, rnew] = size(a);
k=1;
while rold*rnew>n
    k=k+1;
    a = reshape(a, [rold*n, rnew]);
    b = tr{in(k)};
    [rnew, n, rnext] = size(b);
    b = reshape(b, [rnew, n*rnext]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
    [rold, n, rnew] = size(a);
end
[r0, n, r1] = size(a);
a = reshape(permute(a,[1,3,2]), r0*r1, n);
a=a*reshape(T,n,numel(T)/n);
mm=k;
for k = mm+1:d-1
    b = tr{in(k)};
    [rnew, n, r2] = size(b);
    b1=reshape(b,rnew*n,r2);
    sz3=numel(a)/(n*r1*r0);
    a=reshape(a,r0,r1,n,sz3);
    a=reshape(permute(a,[4,1,2,3]),r0*sz3,n*r1)*b1;
    a=reshape(permute(reshape(a,sz3,r0,r2),[2,3,1]),r0*r2,sz3); 
    r1=r2;
end
sz3=numel(a)/(r0*r1);
a=reshape(permute(reshape(a,[r0,r1,sz3]),[3,2,1]),sz3,r0*r1);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
