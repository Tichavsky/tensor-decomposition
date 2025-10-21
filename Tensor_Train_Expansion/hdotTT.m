function H=hdotTT(F,G)
%
% Hadamard product of two tensors in TT format
%
N=length(F);
H=cell(1,N);
for i=1:N
    sf=size(F{i});
    sg=size(G{i});
    if length(sf)<3
        sf(3)=1; 
    end 
    if length(sg)<3
        sg(3)=1;
    end 
    H{i}=zeros(sf(1)*sg(1),sf(2),sf(3)*sg(3));
    for k=1:sf(2)
        H{i}(:,k,:)=kron(reshape(F{i}(:,k,:),sf(1),[]),reshape(G{i}(:,k,:),sg(1),[]));
    end
end
end