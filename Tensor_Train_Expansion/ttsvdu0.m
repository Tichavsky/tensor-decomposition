function [G,bd,b]=ttsvdu0(U,tol,bmax)
%
M=length(U);
F=U{1};
for m=2:M
    F=plusTT(F,U{m});
end
[G,bd,b]=ttsvdtt3(F,tol,bmax);
end