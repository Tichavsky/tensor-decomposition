function F1=prepoctiTT(F,poradi,tol,bmax)
%
% computes wagons (cores) of a TT obtained by
% reordering wagons in F according to the order "poradi"
%
N=length(F);
F1=F; 
apor=1:N; 
for i=1:N-1
    ir=poradi(i);
    ip=find(apor==ir);
    for j=ip:-1:i+1
        F1=vymena(F1,j,tol,bmax);
    end
    apor(i:ip)=[ir,apor(i:ip-1)];
end
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function F=vymena(F,j,tol,bmax)
%
% switching wagons j and j-1
%
[r1,I,r2]=size(F{j});
[s1,J,~]=size(F{j-1});
aux=reshape(F{j-1},[],r1)*reshape(F{j},r1,[]);
aux=permute(reshape(aux,[s1,J,I,r2]),[1,3,2,4]);
[U,S,V]=svd(reshape(aux,s1*I,r2*J));
r=sum(diag(S)>S(1,1)*tol);
r=min([r,bmax]);
[j-1,j,r]
F{j-1}=reshape(U(:,1:r),[s1,I,r]);
F{j}=reshape(S(1:r,1:r)*V(:,1:r)',[r,J,r2]);
end