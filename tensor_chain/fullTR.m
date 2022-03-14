function [res]=fullTR(tr)
%%% Converts a TR-representation into full format.
%%% Algorithm: multiply together the reshaped cores, then make this into
%%% an \alpha_0 , . , \alpha_0 tensor and sum over alpha_0 in a loop. Then
%%% reshape this into an n1 x ... x nd tensor.
%%% modified by Petr Tichavsky, October 2020

d=length(tr);
a = tr{2};    
[r0, n, r1] = size(a);
ns = zeros(1,d);
ns(2) = n;

for k = 3:d
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr{k};
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
b = tr{1};
[rold, n, rnew] = size(b);
ns(1)=n;
rr=rold*rnew;
b = reshape(permute(b,[2,1,3]),n,rr);
a = reshape(permute(a,[3,1,2]), rr, numel(a)/rr);
res=reshape(b*a,ns);   
end