function sen=sensit_tc(A)
%
% Sensitivity of the tensor chain decomposition
% Input:

% Programmed by Petr Tichavsky, December 2020
%
N=length(A);
sen=0;
for n=1:N
    In=size(A{n},2); 
    aux=dejM(n,A);
    sen=sen+In*sum(diag(aux));
end 
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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