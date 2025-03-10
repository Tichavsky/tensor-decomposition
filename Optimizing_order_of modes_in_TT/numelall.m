function num=numelall(F)
%
N=length(F);
num=0;
for i=1:N
    num=num+numel(F{i});
end
end