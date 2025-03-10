function F=dejkvadtt(Q,sz,rastr)
%
N=size(Q,1);
F=cell(1,N);
M=floor(N/2);
F{1}=zeros(1,sz(1),3);
for i=1:sz(1)
    F{1}(1,i,:)=[1,rastr{1}(i),Q(1,1)*rastr{1}(i)^2];
end
for m=2:M
    F{m}=zeros(m+1,sz(m),m+2);
    for i=1:sz(m)
        xi=rastr{m}(i);
        F{m}(:,i,:)=[1,zeros(1,m-1),xi,Q(m,m)*xi^2; zeros(m-1,1),eye(m-1),zeros(m-1,1),2*Q(1:m-1,m)*xi; zeros(1,m+1),1];
    end
end
for m=M+1:N-1
    F{m}=zeros(N-m+3,sz(m),N-m+2);
    for i=1:sz(m)
        xi=rastr{m}(i);
        aux=[1,2*Q(m,m+1:N)*xi,Q(m,m)*xi^2; zeros(1,N-m+1),xi;zeros(N-m,1),eye(N-m),zeros(N-m,1); zeros(1,N-m+1),1];
        F{m}(:,i,:)=reshape(aux,[N-m+3,1,N-m+2]);
    end
end
F{N}=zeros(3,sz(N),1);
for i=1:sz(N)
    xi=rastr{N}(i);
    F{N}(:,i)=[Q(N,N)*xi^2; xi; 1];
end
if 2*M<N
   F{M+1}=zeros(M+2,sz(M+1),M+2); 
   for i=1:sz(M+1)
       xi=rastr{M+1}(i);
       aux=[1,2*Q(M+1,M+2:N)*xi,Q(M+1,M+1)*xi^2; zeros(M,1),2*Q(1:M,M+2:N),2*xi*Q(1:M,M+1);zeros(1,M+1),1];
       F{M+1}(:,i,:)=reshape(aux,[M+2,1,M+2]);
   end
else    
   H=2*Q(1:M,M+1:N);
   for i=1:sz(M+1)
      F{M+1}(2:M+1,i,2:M+1)=reshape(H*squeeze(F{M+1}(2:M+1,i,2:M+1)),M,1,M);
   end
end
end