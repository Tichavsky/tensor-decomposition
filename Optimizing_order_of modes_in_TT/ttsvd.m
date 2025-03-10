function F=ttsvd(T,tol)
if nargin<2
    tol=1e-6;
end
sz=size(T);
N=length(sz);
F=cell(1,N);
X=T; I1=1;
for i=1:N-1
    [U,S,V]=svd(reshape(X,sz(i)*I1,[]),'econ');
    I2=sum(diag(S)>tol);
    F{i}=reshape(U(:,1:I2),[I1,sz(i),I2]);
    X=S(1:I2,1:I2)*V(:,1:I2)';
    I1=I2;
end
F{N}=reshape(X,I2,sz(N),1);
end
