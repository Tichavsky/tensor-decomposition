function F=plusTT(F1,F2,fac)
%
% computes sum of two tensor trains as a single TT
% The input tensors must have the same dimensions.
% The latter TT canbe multiplied by factor fac.
%
N=length(F1);
if nargin>2
    F2{1}=F2{1}*fac;
end
F=F1;
[r1,I,r2]=size(F1{1});
F{1}=reshape([reshape(F1{1},I,r2),reshape(F2{1},I,[])],1,I,[]);
F{N}=[F1{N};F2{N}];
for i=2:N-1
    [r1,I,r2]=size(F1{i});
    [r3,I,r4]=size(F2{i});
    F{i}=zeros(r1+r3,I,r2+r4);
    F{i}(1:r1,:,1:r2)=F1{i};
    F{i}(1+r1:r1+r3,:,1+r2:r2+r4)=F2{i};
end
end
