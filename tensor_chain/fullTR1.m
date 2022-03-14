function [res]=fullTR1(tr)
%
%fo Converts a TR-representation into full format.
% Coded by Petr Tichavsky, August 2021
%
d=length(tr);
d2=floor(d/2);
a = tr{1};    
[r0, n, r1] = size(a);
ns = zeros(1,d);
ns(1) = n;
for k = 2:d2
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr{k};
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
a1=reshape(permute(a,[2,1,3]),[],r0*rnext);
a = tr{d2+1};    
[r2, n, r1] = size(a);
ns(d2+1) = n;
for k = d2+2:d
    [rold, n, rnew] = size(a);
    a = reshape(a, [rold*n, rnew]);
    b = tr{k};
    [rnew, n, rnext] = size(b);
    ns(k) = n;
    b = reshape(b, [rnew, numel(b)/rnew]);
    tmp = a*b;
    a = reshape(tmp, [rold, numel(tmp)/(rold*rnext), rnext]);
end
a2=reshape(permute(a,[2,3,1]),[],r2*rnext);
res=reshape(a1*a2',ns);   
end