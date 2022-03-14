function T=espig(N,D)
%
% The code generates an order-N tensor of
% dimension D x D x...x D
% by sampling the function
% f(x_1,...,x_N) = (1+sum_n x_n^2)^(-1/2)
%
osa=linspace(0,1,D);
X=reshape(repmat(osa.^2,D^(N-1),1),D*ones(1,N));
T=ones(size(X));
for n=1:N-1
    T=T+X;
    X=permute(X,[2:N,1]);
end    
T=T+X;
T=T.^(-1/2);
end