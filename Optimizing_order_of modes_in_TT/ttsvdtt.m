function [G,b]=ttsvdtt(F,tol)
%
N=length(F);
G=F; I1=1; b=zeros(1,N);
for i=1:N-1
    xx=dejKonec(F,i);
    [r1,sz,r2]=size(F{i});
    xx=reshape(F{i},[],r2)*xx*reshape(F{i},[],r2)';
    [U,S,V]=svd(xx);
    sval=sqrt(diag(S));
    I2=sum(sval>tol*sval(1));
    b(i)=sval(I2);
    G{i}=reshape(U(:,1:I2),[I1,sz,I2]);
    s=size(F{i+1});
    aux=U(:,1:I2)'*reshape(F{i},[],r2);
    F{i+1}=reshape(aux*reshape(F{i+1},s(1),[]),[I2,s(2:end)]);
    I1=I2;
end
G{N}=F{N};
end
function xx0=dejKonec(F,i)
%
N=length(F);
xx0=1;
for j=N:-1:i+1
    [r1,~,r2]=size(F{j});
    xx0=reshape(F{j},[],r2)*reshape(xx0,r2,[]);
    xx0=reshape(xx0,r1,[])*reshape(F{j},r1,[])';
end
end