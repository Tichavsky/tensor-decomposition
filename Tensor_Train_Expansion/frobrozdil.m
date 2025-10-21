function val=frobrozdil(tr,tr2)
%
val=sqrt(scalarprod(tr,tr)+scalarprod(tr2,tr2)-2*scalarprod(tr,tr2));
end
