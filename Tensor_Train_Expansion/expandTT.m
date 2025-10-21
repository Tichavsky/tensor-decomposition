function [U,norms,err]=expandTT(x,bmax,M,eps,bmax2,n0)
%
% computes expansion of a tensor train as a sum of n TT's of bond dimension bmax
% 
% M ... maximum number of terms in the expansion
% eps ... minumum relative error in the expansion
% norms ... relative Frobenius norms of the terms
% err .... relative errors of the expansion
% bmax2 ... parameter used in ttsvdu
% n0 ...... max. number of items where ttsvdu0 is used instead of ttsvdu
%
if nargin<4
   eps=1e-5;
end
if nargin<3
   M=100;
end
if nargin<5
   bmax2=bmax;
end
if nargin<6
   n0=10;
end
sps=zeros(M,M);
xn=frobnorm(x);
U=cell(1,M); norms=zeros(1,M); err=norms;
U{1}=ttsvdtt3(x,eps,m);
V{1}=x;
V{2}=U{1};
V{2}{1}=-V{2}{1};
norms(1)=frobnorm(U{1})/xn;
sps(1,1:2)=[1 scalarprod(x,V{2})/xn^2];
sps(2,1:2)=[sps(1,2) norms(1)^2];
err(1)=sqrt(sps(1)+2*sps(1,2)+sps(2,2));
n=1;
while err(n)>eps && n<M
    n=n+1;
    if n<n0
       U{n}=ttsvdu0(V,eps,bmax);
    else
       U{n}=ttsvdu(V,eps,bmax,bmax2); 
    end
    V{n+1}=U{n};
    V{n+1}{1}=-V{n+1}{1};
    norms(n)=frobnorm(U{n})/xn;
    for j=1:n+1
        sps(j,n+1)=scalarprod2(V{n+1},V{j})/xn^2;
    end
    sps(n+1,:)=sps(:,n+1)';
    err(n)=sqrt(sum(reshape(sps(1:n+1,1:n+1),[],1)));
end    
end


