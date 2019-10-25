function [A B C iter]=KLM3cS(X,r,mez,m,numit,A0,B0,C0)
%
%
% Krylov-Levenberg-Marquardt method for CP decomposition of order-3
% tensors with sensitivity bounded by constant "mez". Input:
%
%    T... tensor to be decomposed
%    r ... tensor rank
%    mez ... constant to limit tr((A^TA)*(B^TB)+(A^TA)*(C^TC)+(C^TC)*(B^TB))
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    A,B,C .... initial estimate of factor matrices (if available)
%
% Programmed by Petr Tichavsky, October 2019
%
%
if nargin<5
   numit=20;  %%% this is an ad hoc. It might be increased if "iter" does not converge to a constant.
end    
[Ia,Ib,Ic]=size(X);
X=reshape(X,Ia,Ib*Ic);
if nargin<7
   [A B]=zacni(X,r);   %%%  initial estimate of the factors, followed by one step of the alternating
   C=dejC(X,A,B);
   a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2);
   p1=(a1.*b1.*c1).^(1/3);
   p1=sqrt(mez/(3*sum(p1.^2)))*p1;
   A=A.*repmat(sqrt(p1./a1),Ia,1);
   B=B.*repmat(sqrt(p1./b1),Ib,1);
   C=C.*repmat(sqrt(p1./c1),Ic,1);
else
   A=A0; B=B0; C=C0;
   a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2);
   p1=sum(a1.*(b1+c1)+b1.*c1);
   nr=(mez/p1)^(1/4);
   A=A*nr; B=B*nr; C=C*nr;  
end    
it=1; incr=10; iter=0;
nu=2; Iabc=Ia+Ib+Ic; len=Iabc*r;
tau=1;
na=sum(A.^2); nb=sum(B.^2); nc=sum(C.^2);
mm=[nb.*nc na.*nc nb.*na];
mu=tau*max(mm);
iter(1)=chyba(X,A,B,C); 
lam=0.0947;
while it<numit && iter(end)>1e-12
   it=it+1;
   ar=sum(A.^2); br=sum(B.^2); cr=sum(C.^2); 
   f0=sum(ar.*(br+cr)+br.*cr);
   [dd,R,u]=krok3(reshape(X,Ia,Ib,Ic),A,B,C,m,lam,mu); %%% computes gradient of error, gradient of constraint, etc
   lam=(f0-mez-u'*dd(:,1))/(u'*dd(:,2));
   d=-dd(:,1)-lam*dd(:,2);
   err=chyba(X,A,B,C);              %%% computes the error of the current approximation
   A1=A; B1=B; C1=C;
   A=A+reshape(d(1:Ia*r),Ia,r); B=B+reshape(d(1+Ia*r:r*(Ia+Ib)),Ib,r); C=C+reshape(d(1+(Ia+Ib)*r:Iabc*r),Ic,r); 
   a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2);
   p1=sum(a1.*(b1+c1)+b1.*c1);
   nr=(mez/p1)^(1/4);
   A=A*nr; B=B*nr; C=C*nr;  
   err2=chyba(X,A,B,C);
   ar=sum(A.^2); br=sum(B.^2); cr=sum(C.^2); 
   f2=sum(ar.*(br+cr)+br.*cr);
   rho=(err-err2)/(d'*(R+mu*d));
   if err2>err                                     %%% step is not accepted
      mu=mu*nu; nu=2*nu; err2=err;    
      A=A1; B=B1; C=C1;                     %%% return to the previous solution
   else
      nu=2;                                        
      mu=mu*max([1/3 1-(2*rho-1)^3]);
  %    theta=reshape([A; B; C],r*Iabc,1);
      if rem(it,10)==0
       [log10(it) iter(end)]  % to monitor decrease of the cost function each 10 iterations
      end 
   end   
   iter(it)=err2;
end   
end
%semilogy(1:it,iter,'c:')

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [A B]=zacni(X,r)
%
% computes orthogonal bases of spaces of columns, rows
% and (??? how to call it)  in order to give an initial estimate of a
% tensor decomposition.
%
[Ia Ib Ic]=size(X);
A=randn(Ia,r);
B=randn(Ia,r);
end

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [C err]=dejC(X,A,B)
%
% computes vectors in A in a decomposition given B and C,
% respectively, in the mean square sense, in order to implement
% parallel factor analysis
%
[Ib r]=size(B);
[Ia r]=size(A);
[Ia Ibc]=size(X);
Ic=floor(Ibc/Ib);
M=inv((B'*B).*(A'*A));
H=A'*X;
S=[];
err=sum(X(:).^2);
for i=1:Ic
    ind=(i-1)*Ib;
    Aux=H(:,ind+1:ind+Ib);
    sloupec=diag(B'*Aux');
    S=[S sloupec];
    err=err-sloupec'*M*sloupec;
end
C=(M*S)';
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chyba(X,A,B,C)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
Y=X-A*krb(C,B)';
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

function [dd,g,u,Y,W,Z,X]=krok3(X,A,B,C,m,lam,mu)
[Ia,R]=size(A);
[Ib R]=size(B);
[Ic R]=size(C);
Iab=Ia+Ib; Iabc=Iab+Ic;
Y1=reshape(X,Ia,Ib*Ic);
Y2=reshape(permute(X,[2,1,3]),Ib,Ia*Ic);
Y3=reshape(permute(X,[3,1,2]),Ic,Ia*Ib);
AA=A'*A;
BB=B'*B;
CC=C'*C;
gA=Y1*krb(C,B)-A*(BB.*CC);   % computes error gradient (mttkrp)
gB=Y2*krb(C,A)-B*(AA.*CC);
gC=Y3*krb(B,A)-C*(BB.*AA);
g=[gA(:); gB(:); gC(:)];
ng=norm(g);
a1=sum(A.^2); b1=sum(B.^2); c1=sum(C.^2);
u1=A.*repmat(b1+c1,Ia,1); u2=B.*repmat(a1+c1,Ib,1); u3=C.*repmat(b1+a1,Ic,1);
u=[u1(:); u2(:); u3(:)];
nu=norm(u);
Y=zeros(R*(Ia+Ib+Ic),m+2); %  Krylov subspace 1
W=zeros(m+1,m+1); 
Y(:,1)=g/ng;
gA=gA/ng; gB=gB/ng; gC=gC/ng;
for i=1:m+1
    yA=gA*(CC.*BB)+A*((gB'*B).*CC+(gC'*C).*BB);
    yB=gB*(AA.*CC)+B*((gC'*C).*AA+(gA'*A).*CC);
    yC=gC*(AA.*BB)+C*((gB'*B).*AA+(gA'*A).*BB);
    ax1=sum(A.*gA); bx1=sum(B.*gB); cx1=sum(C.*gC);
    h3a=repmat(b1+c1,Ia,1).*gA+2*A.*repmat(bx1+cx1,Ia,1);
    h3b=repmat(a1+c1,Ib,1).*gB+2*B.*repmat(ax1+cx1,Ib,1);
    h3c=repmat(a1+b1,Ic,1).*gC+2*C.*repmat(ax1+bx1,Ic,1);
    zz=[yA(:)+lam*h3a(:); yB(:)+lam*h3b(:); yC(:)+lam*h3c(:)];
    W(1:i,i)=Y(:,1:i)'*zz;
    Y(:,i+1)=zz-Y(:,1:i)*W(1:i,i);
    Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
    gA=reshape(Y(1:R*Ia,i+1),Ia,R); 
    gB=reshape(Y(R*Ia+1:R*(Ia+Ib),i+1),Ib,R); 
    gC=reshape(Y(R*(Ia+Ib)+1:R*(Ia+Ib+Ic),i+1),Ic,R);
end
W=W+W'-diag(diag(W)); 
iW=inv(W);
dd=[-g/mu+Y(:,1:m+1)/(iW+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2,g];
%
Z=zeros(R*(Ia+Ib+Ic),m+2); %  Krylov subspace 2
X=zeros(m+1,m+1); 
uA=reshape(u(1:Ia*R)/nu,Ia,R); uB=reshape(u(1+Ia*R:Iab*R)/nu,Ib,R); uC=reshape(u(1+Iab*R:Iabc*R)/nu,Ic,R);
Z(:,1)=u/nu;
for i=1:m+1
    yA=uA*(CC.*BB)+A*((uB'*B).*CC+(uC'*C).*BB);
    yB=uB*(AA.*CC)+B*((uC'*C).*AA+(uA'*A).*CC);
    yC=uC*(AA.*BB)+C*((uB'*B).*AA+(uA'*A).*BB);
    ax1=sum(A.*uA); bx1=sum(B.*uB); cx1=sum(C.*uC);
    h3a=repmat(b1+c1,Ia,1).*uA+2*A.*repmat(bx1+cx1,Ia,1);
    h3b=repmat(a1+c1,Ib,1).*uB+2*B.*repmat(ax1+cx1,Ib,1);
    h3c=repmat(a1+b1,Ic,1).*uC+2*C.*repmat(ax1+bx1,Ic,1);
    uu=[yA(:)+lam*h3a(:); yB(:)+lam*h3b(:); yC(:)+lam*h3c(:)];
    X(1:i,i)=Z(:,1:i)'*uu;
    Z(:,i+1)=uu-Z(:,1:i)*X(1:i,i);
    Z(:,i+1)=Z(:,i+1)/norm(Z(:,i+1));
    uA=reshape(Z(1:R*Ia,i+1),Ia,R); 
    uB=reshape(Z(R*Ia+1:R*(Ia+Ib),i+1),Ib,R); 
    uC=reshape(Z(R*(Ia+Ib)+1:R*(Ia+Ib+Ic),i+1),Ic,R);
end
X=X+X'-diag(diag(X)); 
iX=inv(X);
dd(:,2)=u/mu-Z(:,1:m+1)/(iX+eye(m+1)/mu)*(Z(:,1:m+1)'*u)/mu^2;
end   