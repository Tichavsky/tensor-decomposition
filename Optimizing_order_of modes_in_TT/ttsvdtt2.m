function [G,bd,b]=ttsvdtt2(F,tol,bmax)
%
N=length(F);
if nargin<3
    bmax=10000;
end    
G=F; I1=1; b=zeros(1,N); bd=b;
mm=dejKonce(F);
for i=1:N-1
    xx=mm{i};
    [r1,sz,r2]=size(F{i});
    xx=reshape(F{i},[],r2)*xx*reshape(F{i},[],r2)';
    [U,S,V]=svd(xx);
    sval=sqrt(diag(S))
    I2=sum(sval>tol*sval(1));
    I2=min([I2,bmax]);
    if I2<1
       I2=1;
    end   
    b(i)=sval(I2); bd(i)=I2;
    G{i}=reshape(U(:,1:I2),[I1,sz,I2]);
    s=size(F{i+1});
    aux=U(:,1:I2)'*reshape(F{i},[],r2);
    F{i+1}=reshape(aux*reshape(F{i+1},s(1),[]),[I2,s(2:end)]);
    I1=I2;
end
G{N}=F{N};
if size(G{N},1)>size(G{N},2)
   [r1,sz,r2]=size(G{N-1});
   aux=reshape(G{N-1},[],r2)*reshape(G{N},r2,[]);
   [U,S,V]=svd(aux);
   I2=sum(diag(S)>tol*S(1,1));
   G{N-1}=reshape(U(:,1:I2),r1,sz,I2);
   G{N}=reshape(S(1:I2,1:I2)*V(:,1:I2)',I2,[]);
end    
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function mm=dejKonce(F)
%
N=length(F);
mm=cell(1,N);
xx0=1;
for j=N:-1:2
    [r1,~,r2]=size(F{j});
    xx0=reshape(F{j},[],r2)*reshape(xx0,r2,[]);
    xx0=reshape(xx0,r1,[])*reshape(F{j},r1,[])';
    mm{j-1}=xx0;
end
end
