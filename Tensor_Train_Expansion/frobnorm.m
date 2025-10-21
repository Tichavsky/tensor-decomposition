function val=frobnorm(tr)
%
% computes Frobenius norm of a tensor defined through wagons in 
% its tensor chain decomposition {{A{1},...A{N}}}
%
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
val=sqrt(scalarprod(tr,tr));
end
