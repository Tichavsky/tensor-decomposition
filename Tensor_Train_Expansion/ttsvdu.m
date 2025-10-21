function x=ttsvdu(U,tol,Bmax,Bmax2)
%
% SVD-based TT decomposition in SOTT format
% the algorithm has linear complexity wrt the number of the elements in SOTT
%
M=length(U); N=length(U{1});
ms=cell(1,N); b=cell(1,M); x=cell(1,N);
for m=1:M
    b{m}=1;
end
ms{N}=b;
for j=N:-1:2
    W=0;
    for m=1:M
        [r1,~,r2]=size(U{m}{j});
        b{m}=reshape(reshape(U{m}{j},[],r2)*reshape(b{m},r2,[]),r1,[]);
        W=W+b{m}'*b{m};
    end
    [V,S,~]=svd(W);
    B=min([Bmax2,sum(diag(S)>=tol*S(1,1))]);
    for m=1:M
        b{m}=b{m}*V(:,1:B);
    end
    ms{j-1}=b;
end
In=size(U{1}{1},2);
for m=1:M
    b{m}=reshape(U{m}{1},In,[]);
end
for i=1:N-1
    W=0; In=size(U{1}{i},2);
    for m=1:M
        W=W+b{m}*ms{i}{m};
    end
    [V,S,~]=svd(W,'econ');
    B=min([Bmax,sum(diag(S)>=tol*S(1,1))]);
    x{i}=reshape(V(:,1:B),[],In,B);
    for m=1:M
        aux=V(:,1:B)'*b{m};
        s0=size(aux,2);  s1=size(ms{i+1}{m},1);
        aux=aux*reshape(U{m}{i+1},s0,[]);
        b{m}=reshape(aux,[],s1);
    end
end
x{N}=b{1}; In=size(U{1}{N},2);
for m=2:M
    x{N}=x{N}+b{m};
end
x{N}=reshape(x{N},[],In);
end