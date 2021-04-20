function T=MMtensor(n,m,k)
%
%creates a tensor that corresponds to multiplication of two matrices of the
%sizes n x m and m x k
%
T=zeros(n*m,m*k,n*k);
for j=1:m
    for im=1:n
        for l=1:k
            T((im-1)*m+j,(j-1)*k+l,(l-1)*n+im)=1;
        end
    end
end    
