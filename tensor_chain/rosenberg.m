function T=rosenberg
%
% generates a sampled Rosenberg function of order 3
%
x=-2:0.05:2;
y=-1:0.05:3;
X1=reshape(repmat(x',1,81^2),81,81,81);
X2=reshape(repmat(x',1,81^2),81,81,81); X2=permute(X2,[2,3,1]);
X3=reshape(repmat(y',1,81^2),81,81,81); X3=permute(X3,[3,1,2]);
T=100*((X2-X1.^2).^2+(X2-X3).^2)+(X1-1).^2+(X2-1).^2;
end