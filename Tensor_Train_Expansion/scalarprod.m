function s=scalarprod(x,y)
%
% Scalar product of two tensors in the tensor train format
%
   N=length(x); n = size(x{1},2);
   a=reshape(y{1},n,[])'*reshape(x{1},n,[]);
   for j=2:N-1
     r2 = size(x{j},1); [s2, n, ~] = size(y{j});
     a=reshape(a*reshape(x{j},r2,[]),s2*n,[]);
     a=reshape(y{j},s2*n,[])'*a;
   end
   s=a(:)'*reshape(y{N}*x{N}',[],1);
end