function [poradi,bd,npar]=ttopttt(F,tol,bmax)
%
% Optimization of the order of the modes in TT decomposition of a tensor.
%
% Input: F .... tensor in TT format (array of wagons)
%        tol .. toleration constant, typically small, such as 1e-5
% Output: poradi ... permutation of the wagons in which the bond dimensions
%                    should be the small as possible
%         bd ... vector of the bond dimensions in the new ordering
%         npar ... number of parameters (n. of elements of all wagons)
%
% The new wagons can be found using F, poradi, and function prepoctiTT.
%
N=length(F);
poradi=1:N;
volne=1:N;
bd=zeros(1,N-1);
[i0,bd0,konce]=prvniTT(F,tol,bmax);
poradi(1)=i0;
volne(i0)=[];
bd(1)=bd0;
npar=bd0*size(F{i0},2);
for n=2:N-1
    bdmin=1000; bvmin=1e10; pmin=0;
    for i=1:length(volne)
        ind=sort([poradi(1:n-1),volne(i)]);
        [bd0,bv0]=ttsvdaux(F,ind,tol,konce,bmax);
        if bd0<bdmin || (bd0==bdmin && bv0<bvmin)
           pmin=volne(i); 
           bdmin=bd0;
           bvmin=bv0;
        end
    end    
    poradi(n)=pmin;
    bd(n)=bdmin;
    npar=npar+bd(n)*bd(n-1)*size(F{pmin},2);
    volne(volne==pmin)=[];
end
poradi(N)=volne;
npar=npar+bd(N-1)*size(F{volne},2);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [s2,bv]=ttsvdaux(F,indx,tol,konce,bmax)
%
N=length(F);
M=length(indx);
G=cell(1,M);
xx=1;
for i=1:indx(1)-1
    [r1,sz,r2]=size(F{i});    
    xx=reshape(xx,[],r1)*reshape(F{i},r1,[]);
    xx=reshape(xx,[],r2)'*reshape(F{i},[],r2);
end
i0=indx(1);
[r1,sz,r2]=size(F{i0});  
xx3=permute(reshape(reshape(F{i0},[],r2)*konce{i0+1},[r1,sz,r2]),[2,3,1]);
xx4=permute(reshape(xx*reshape(F{i0},r1,[]),[r1,sz,r2]),[2,3,1]);
xx2=reshape(xx3,sz,[])*reshape(xx4,sz,[])';
[U,S,~]=svd(xx2);
sval=sqrt(diag(S));
s2=sum(sval>tol*sval(1));
s2=min([s2,bmax]);
G{1}=reshape(U(:,1:s2),1,[],s2);
F{i0}=reshape(reshape(permute(F{i0},[3,1,2]),[],sz)*U(:,1:s2),[r2,r1,s2]);
F{i0}=reshape(permute(F{i0},[2,3,1]),r1,s2,r2);
aux=reshape(F{i0},r1,[])'*xx*reshape(F{i0},r1,[]); 
for j=2:M
for i=indx(j-1)+1:indx(j)-1
    [r1,sz,r2]=size(F{i});
    aux=permute(reshape(aux,s2,r1,s2,r1),[1,3,2,4]);
    H=reshape(permute(F{i},[1,3,2]),[],sz);
    xx1=permute(reshape(H*H',r1,r2,r1,r2),[1,3,2,4]);
    aux=reshape(reshape(aux,[],r1^2)*reshape(xx1,r1^2,[]),[s2,s2,r2,r2]);
    aux=reshape(permute(aux,[1,3,2,4]),s2*r2,s2*r2);
end
i=indx(j); s1=s2;
[r1,sz,r2]=size(F{i});
    xx1=reshape(F{i},[],r2)*konce{i+1}*reshape(F{i},[],r2)';
    xx1=reshape(xx1,r1,sz,r1,sz);
    xx1=reshape(permute(xx1,[2,4,1,3]),[],r1^2);
    xx2=reshape(aux,s1,r1,s1,r1);  
    xx2=reshape(permute(xx2,[1,3,2,4]),[],r1^2);
    xx=reshape(xx1*xx2',[sz,sz,s1,s1]);
    xx=reshape(permute(xx,[3,1,4,2]),sz*s1,sz*s1);
    [U,S,~]=svd(xx);
    sval=sqrt(diag(S));
    s2=sum(sval>tol*sval(1));
    s2=min([s2,bmax]);
    G{j}=reshape(U(:,1:s2),[s1,sz,s2]);
    bux=reshape(permute(G{j},[1,3,2]),[],sz)*reshape(permute(F{i},[2,1,3]),sz,[]);
    bux=reshape(permute(reshape(bux,[s1,s2,r1,r2]),[1,3,2,4]),s1*r1,s2*r2);
    aux=bux'*aux*bux;
end
sval=sval(1:s2);
bv=sval(s2);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [i0,bd,konce]=prvniTT(F,tol,bmax)
%
N=length(F);
konce=cell(1,N); 
zacatky=cell(1,N);
xx0=1;
for j=N:-1:2
    [r1,~,r2]=size(F{j});
    xx0=reshape(F{j},[],r2)*reshape(xx0,r2,[]);
    xx0=reshape(xx0,r1,[])*reshape(F{j},r1,[])';
    konce{j}=xx0;
end
I2=1000; skr=1000;
xx=1;
zacatky{1}=1;
for i=1:N-1
    [r1,sz,r2]=size(F{i});
    xx3=permute(reshape(reshape(F{i},[],r2)*konce{i+1},[r1,sz,r2]),[2,3,1]);
    xx4=permute(reshape(xx*reshape(F{i},r1,[]),[r1,sz,r2]),[2,3,1]);
    xx2=reshape(xx3,sz,[])*reshape(xx4,sz,[])';
    [~,S,~]=svd(xx2);
    sval=sqrt(diag(S));
    I2a=sum(sval>tol*sval(1));
    I2a=min([I2a,bmax]);
    if I2a<I2 || (I2==I2a && S(I2,I2)<skr)
           i0=i; I2=I2a; skr=S(I2,I2);
    end
    xx=reshape(xx,[],r1)*reshape(F{i},r1,[]);
    xx=reshape(xx,[],r2)'*reshape(F{i},[],r2);
    zacatky{i+1}=xx;
end
bd=I2; 
konce{N+1}=1;
end