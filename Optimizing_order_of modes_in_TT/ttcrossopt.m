function [F,order]=ttcrossopt(getT,sz,tol)
%
% Optimizing order of modes by means of TT-cross
% Requires toolboox of Oseledets to be installed.
% Function getT can be defined through the function evalTT
% if the original tensor is given in TT form.
%
N=length(sz);
order=1:N; iord=1:N;
%z=greedy2_cross(sz,getT,tol);
z=dmrg_cross(N,sz(1),getT,1e-8);
F=core2cell(z)';
F=ttsvdtt(F,tol);
for i=1:N-1
    bmin=size(F{i},3);
    vmin=bondValue(F,i); 
    ord=order;
    for k=1:N-i
        ord=[ord(1:i-1),ord(i+1:N),ord(i)];
        iord(ord)=1:N;
        getTp=@(ix) getT(ix(iord));
     %   z=greedy2_cross(sz(ord),getTp,tol);
        z=dmrg_cross(N,sz(1),getTp,1e-4);
        Fa=core2cell(z)';
        Fa=ttsvdtt(Fa,tol);
        b=size(Fa{i},3);
        v=bondValue(Fa,i);
        if b<bmin || (b==bmin && v(b)<vmin(b))
           F=Fa; vmin=v; bmin=b; order=ord;
        end
    end
end    