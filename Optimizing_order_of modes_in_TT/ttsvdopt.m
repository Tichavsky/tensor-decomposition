function [F,poradi]=ttsvdopt(T,tol)
if nargin<2
    tol=1e-6;
end
sz=size(T);
N=length(sz);
F=cell(1,N);
poradi=1:N;
I1=1;
X=T;
for i=1:N-1
    [U,S,V]=svd(reshape(X,sz(i)*I1,[]),'econ');
    I2=sum(diag(S)>tol*S(1,1));
    skr=S(I2,I2); pora=poradi;
    sza=sz;
    Xa=X;
    for j=i+1:N
        Xa=permute(reshape(Xa,I1,sza(j-1),[]),[1,3,2]);
        [Ua,Sa,Va]=svd(reshape(Xa,sz(i)*I1,[]),'econ');
        I2a=sum(diag(Sa)>tol*Sa(1,1));
        if I2a<I2 || (I2==I2a && Sa(I2,I2)<skr)
           U=Ua; S=Sa; V=Va; I2=I2a; skr=S(I2,I2);
           poradi=[pora(1:i-1),pora(j:N),pora(i:j-1)];
           sz=[sza(1:i-1),sza(j:N),sza(i:j-1)];
        end
    end    
    F{i}=reshape(U(:,1:I2),[I1,sz(i),I2]);
    X=S(1:I2,1:I2)*V(:,1:I2)';
    I1=I2;
end
F{N}=reshape(X,I2,sz(N),1);
end
