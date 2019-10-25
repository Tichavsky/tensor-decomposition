function [A,B,C,D,iter]=KLM4cS(T,R,mez,m,numit,A,B,C,D)
%
% Krylov-Levenberg-Marquardt method for CP decomposition of order-4
% tensors with sensitivity limited by a constant (mez) 
%
%Input:
%
%    T... tensor to be decomposed
%    R ... tensor rank
%    mez ... constant to limit tr((A^TA)*(B^TB)*(C^TC)+(C^TC)*(B^TB)*(D^TD)
%                               +(A^TA)*(C^TC)*(D^TD)+(A^TA)*(B^TB)*(D^TD))
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    A,B,C,D .... initial estimate of factor matrices (if available)
%
% Programmed by Petr Tichavsky, October 2019
%
[Ia,Ib,Ic,Id]=size(T); 
Iab=Ia+Ib;
Iabc=Iab+Ic;
Iabcd=Iabc+Id;
if nargin<9
    A=randn(Ia,R);    %% random factor matrices ase default
    B=randn(Ib,R);
    C=randn(Ic,R);
    D=randn(Id,R);  
end 
if nargin<5
    numit=100;
end
if nargin<4
    m=20;
end        
a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2); d1=sum(D.^2);
p1=sum(a1.*b1.*(c1+d1)+c1.*d1.*(a1+b1));
nr=(mez/p1)^(1/6);
A=A*nr; B=B*nr; C=C*nr; D=D*nr;   
iter=zeros(1,numit); 
na=sum(A.^2); nb=sum(B.^2); nc=sum(C.^2);  nd=sum(D.^2);
mm=[nb.*nc.*nd na.*nc.*nd nd.*nb.*na na.*nb.*nc];   
it=1; 
nu=2; 
tau=1;
mu=tau*max(mm);
iter=chyba(T,A,B,C,D);  %%% computes the error of the current approximation
err=iter;
lam=0; %0.0947;
while it<numit 
   it=it+1;
   a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2); d1=sum(D.^2);
   f0=sum(a1.*b1.*(c1+d1)+c1.*d1.*(a1+b1));
   [dd,g,u]=krok4(T,A,B,C,D,m,lam,mu);    %% gradients and steps
   lam=(f0-mez-u'*dd(:,1))/(u'*dd(:,2));
   d=-dd(:,1)-lam*dd(:,2);            
   A1=A; B1=B; C1=C; D1=D;
   A=A+reshape(d(1:Ia*R),Ia,R); B=B+reshape(d(1+Ia*R:R*(Ia+Ib)),Ib,R); 
   C=C+reshape(d(1+(Ia+Ib)*R:Iabc*R),Ic,R); D=D+reshape(d(1+Iabc*R:Iabcd*R),Id,R); 
   a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2); d1=sum(D.^2);
   p1=sum(a1.*b1.*(c1+d1)+c1.*d1.*(a1+b1));
   nr=(mez/p1)^(1/6);
   A=A*nr; B=B*nr; C=C*nr; D=D*nr;  
   err2=chyba(T,A,B,C,D);
   rho=(err-err2)/(d'*(g+mu*d));
   if err2>err                               %%% step is not accepted
      mu=mu*nu; nu=2*nu; err2=err;    
      A=A1; B=B1; C=C1; D=D1;                %%% return to the previous solution
   else
      nu=2;                                        
      mu=mu*max([1/3 1-(2*rho-1)^3]);
      if rem(it,10)==0
       [log10(it) iter(end)]  % to monitor decrease of the cost function each 10 iterations
      end 
      err=err2;
   end   
   iter(it)=err2;
end   
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chyba(X,A,B,C,D)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B,C, and D
%
Ia=size(A,1);
Y=reshape(X,Ia,[])-A*krb(D,krb(C,B))';
err=sum(Y(:).^2);
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function AB = krb(A,B)

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in krb.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [dd,g,u,Y,W,Z,X]=krok4(T,A,B,C,D,m,lam,mu)
[Ia,R]=size(A);
[Ib R]=size(B);
[Ic R]=size(C);
[Id R]=size(D);
Iab=Ia+Ib; Iabc=Iab+Ic; Iabcd=Iabc+Id;
Y1=reshape(T,Ia,Ib*Ic*Id);
Y2=reshape(permute(T,[2,3,4,1]),Ib,Ia*Ic*Id);
Y3=reshape(permute(T,[3,4,1,2]),Ic,Ia*Ib*Id);
Y4=reshape(permute(T,[4,1,2,3]),Id,Ia*Ib*Ic);
AA=A'*A; BB=B'*B; CC=C'*C; DD=D'*D;
gA=Y1*krb(D,krb(C,B))-A*(BB.*CC.*DD);   % computes error gradient (mttkrp)
gB=Y2*krb(A,krb(D,C))-B*(AA.*CC.*DD);
gC=Y3*krb(B,krb(A,D))-C*(BB.*AA.*DD);
gD=Y4*krb(C,krb(B,A))-D*(BB.*AA.*CC);
g=[gA(:); gB(:); gC(:); gD(:)];
ng=norm(g);
a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2); d1=sum(D.^2);
a2=(b1+c1).*d1+b1.*c1; b2=(a1+c1).*d1+a1.*c1;
c2=(b1+a1).*d1+b1.*a1; d2=(b1+c1).*a1+b1.*c1;
u1=A.*repmat(b1.*c1+c1.*d1+b1.*d1,Ia,1); 
u2=B.*repmat(a1.*c1+c1.*d1+a1.*d1,Ib,1); 
u3=C.*repmat(b1.*a1+b1.*d1+a1.*d1,Ic,1); 
u4=D.*repmat(b1.*a1+b1.*c1+a1.*c1,Id,1);
u=[u1(:); u2(:); u3(:); u4(:)];
nu=norm(u);
Y=zeros(R*Iabcd,m+2); %  Krylov subspace 1 for vector g
W=zeros(m+1,m+1); 
Y(:,1)=g/ng;
gA=gA/ng; gB=gB/ng; gC=gC/ng; gD=gD/ng;
for i=1:m+1
    i1=i-1;
    if i==1
        i1=1;
    end
    %i1=max([1,i-1]);
    yA=gA*(CC.*BB.*DD)+A*((gB'*B).*CC.*DD+((gC'*C).*DD+(gD'*D).*CC).*BB);
    yB=gB*(AA.*CC.*DD)+B*(((gC'*C).*AA+(gA'*A).*CC).*DD+(gD'*D).*AA.*CC);
    yC=gC*(AA.*BB.*DD)+C*(((gB'*B).*AA+(gA'*A).*BB).*DD+(gD'*D).*BB.*AA);
    yD=gD*(AA.*BB.*CC)+D*(((gB'*B).*AA+(gA'*A).*BB).*CC+(gC'*C).*BB.*AA);
    %%
    ax1=sum(A.*gA); bx1=sum(B.*gB); cx1=sum(C.*gC); dx1=sum(D.*gD);
    h3a=repmat((c1+d1).*bx1+(b1+d1).*cx1+(b1+c1).*dx1,Ia,1).*A;
    h3b=repmat((c1+d1).*ax1+(a1+d1).*cx1+(a1+c1).*dx1,Ib,1).*B;
    h3c=repmat((a1+d1).*bx1+(b1+d1).*ax1+(b1+a1).*dx1,Ic,1).*C;
    h3d=repmat((c1+b1).*ax1+(a1+c1).*bx1+(a1+b1).*cx1,Id,1).*D;
    h3a=gA.*repmat(a2,Ia,1)+2*h3a; 
    h3b=gB.*repmat(b2,Ib,1)+2*h3b; 
    h3c=gC.*repmat(c2,Ic,1)+2*h3c; 
    h3d=gD.*repmat(d2,Id,1)+2*h3d;
    zz=[yA(:)+lam*h3a(:); yB(:)+lam*h3b(:); yC(:)+lam*h3c(:); yD(:)+lam*h3d(:)];
    W(i1:i,i)=Y(:,i1:i)'*zz;
    Y(:,i+1)=zz-Y(:,i1:i)*W(i1:i,i);
    Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
    gA=reshape(Y(1:R*Ia,i+1),Ia,R); 
    gB=reshape(Y(R*Ia+1:R*Iab,i+1),Ib,R); 
    gC=reshape(Y(R*Iab+1:R*Iabc,i+1),Ic,R);
    gD=reshape(Y(R*Iabc+1:R*Iabcd,i+1),Id,R);
end
W=W+W'-diag(diag(W)); 
iW=inv(W);
dd=[-g/mu+Y(:,1:m+1)/(iW+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2,g];
%
Z=zeros(R*Iabcd,m+2); %  Krylov subspace 2 for vector u (gradient of sensitivity)
X=zeros(m+1,m+1); 
uA=reshape(u(1:Ia*R)/nu,Ia,R); uB=reshape(u(1+Ia*R:Iab*R)/nu,Ib,R); 
uC=reshape(u(1+Iab*R:Iabc*R)/nu,Ic,R); uD=reshape(u(1+Iabc*R:Iabcd*R)/nu,Id,R);
Z(:,1)=u/nu;
for i=1:m+1
    yA=uA*(CC.*BB.*DD)+A*((uB'*B).*CC.*DD+((uC'*C).*DD+(uD'*D).*CC).*BB);
    yB=uB*(AA.*CC.*DD)+B*(((uC'*C).*AA+(uA'*A).*CC).*DD+(uD'*D).*AA.*CC);
    yC=uC*(AA.*BB.*DD)+C*(((uB'*B).*AA+(uA'*A).*BB).*DD+(uD'*D).*BB.*AA);
    yD=uD*(AA.*BB.*CC)+D*(((uB'*B).*AA+(uA'*A).*BB).*CC+(uC'*C).*BB.*AA);
    %
    ax1=sum(A.*uA); bx1=sum(B.*uB); cx1=sum(C.*uC); dx1=sum(D.*uD);
    h3a=repmat((c1+d1).*bx1+(b1+d1).*cx1+(b1+c1).*dx1,Ia,1).*A;
    h3b=repmat((c1+d1).*ax1+(a1+d1).*cx1+(a1+c1).*dx1,Ib,1).*B;
    h3c=repmat((a1+d1).*bx1+(b1+d1).*ax1+(b1+a1).*dx1,Ic,1).*C;
    h3d=repmat((c1+b1).*ax1+(a1+c1).*bx1+(a1+b1).*cx1,Id,1).*D;
    h3a=uA.*repmat(a2,Ia,1)+2*h3a; 
    h3b=uB.*repmat(b2,Ib,1)+2*h3b; 
    h3c=uC.*repmat(c2,Ic,1)+2*h3c; 
    h3d=uD.*repmat(d2,Id,1)+2*h3d;
    uu=[yA(:)+lam*h3a(:); yB(:)+lam*h3b(:); yC(:)+lam*h3c(:); yD(:)+lam*h3d(:)];
    X(1:i,i)=Z(:,1:i)'*uu;
    Z(:,i+1)=uu-Z(:,1:i)*X(1:i,i);
    Z(:,i+1)=Z(:,i+1)/norm(Z(:,i+1));
    uA=reshape(Z(1:R*Ia,i+1),Ia,R); 
    uB=reshape(Z(R*Ia+1:R*(Ia+Ib),i+1),Ib,R); 
    uC=reshape(Z(R*(Ia+Ib)+1:R*(Ia+Ib+Ic),i+1),Ic,R);
    uD=reshape(Z(1+Iabc*R:Iabcd*R,i+1),Id,R);
end
X=X+X'-diag(diag(X)); 
iX=inv(X);
dd(:,2)=u/mu-Z(:,1:m+1)/(iX+eye(m+1)/mu)*(Z(:,1:m+1)'*u)/mu^2;
end