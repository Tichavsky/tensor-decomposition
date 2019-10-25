function [A,B,C,iter]=KLM3(T,R,m,numit,A,B,C)
%
% Krylov-Levenberg-Marquardt method for CP decomposition of order-3
% tensors. Input:
%
%    T... tensor to be decomposed
%    R ... tensor rank
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    A,B,C .... initial estimate of factor matrices (if available)
%
% Programmed by Petr Tichavsky, September 2018
%
[Ia,Ib,Ic]=size(T); 
Iabc=Ia+Ib+Ic;
if nargin<7
    A=randn(Ia,R);
    B=randn(Ib,R);
    C=randn(Ic,R);
end    
Y1=reshape(T,Ia,Ib*Ic);
Y2=reshape(permute(T,[2,1,3]),Ib,Ia*Ic);
Y3=reshape(permute(T,[3,1,2]),Ic,Ia*Ib);
iter=zeros(1,numit);
na=sum(A.^2); nb=sum(B.^2); nc=sum(C.^2);
mm=[nb.*nc na.*nc nb.*na];
mu=max(mm);  % initial parameter mu
nu=2;
err=chyba(Y1,A,B,C);              %%% computes the error of the current approximation
for it=1:numit
    AA=A'*A;
    BB=B'*B;
    CC=C'*C;
    gA=Y1*krb(C,B)-A*(BB.*CC);   % computes error gradient (mttkrp)
    gB=Y2*krb(C,A)-B*(AA.*CC);
    gC=Y3*krb(B,A)-C*(BB.*AA);
g=[gA(:); gB(:); gC(:)];
Y=zeros(R*(Ia+Ib+Ic),m+2); %  Krylov subspace
ng=norm(g);
Y(:,1)=g/ng;
gA=gA/ng; gB=gB/ng; gC=gC/ng;
W=zeros(m+1,m+1); 
for i=1:m+1
    yA=gA*(CC.*BB)+A*((gB'*B).*CC+(gC'*C).*BB);
    yB=gB*(AA.*CC)+B*((gC'*C).*AA+(gA'*A).*CC);
    yC=gC*(AA.*BB)+C*((gB'*B).*AA+(gA'*A).*BB);
    zz=[yA(:); yB(:); yC(:)];
    W(1:i,i)=Y(:,1:i)'*zz;
    Y(:,i+1)=zz-Y(:,1:i)*W(1:i,i);
    Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
    gA=reshape(Y(1:R*Ia,i+1),Ia,R); 
    gB=reshape(Y(R*Ia+1:R*(Ia+Ib),i+1),Ib,R); 
    gC=reshape(Y(R*(Ia+Ib)+1:R*(Ia+Ib+Ic),i+1),Ic,R);
end
W=W+W'-diag(diag(W));
d=g/mu-Y(:,1:m+1)/(inv(W)+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2;
A1=A; B1=B; C1=C;
A=A+reshape(d(1:Ia*R),Ia,R); 
B=B+reshape(d(1+Ia*R:(Ia+Ib)*R),Ib,R); 
C=C+reshape(d(1+(Ia+Ib)*R:end),Ic,R); 
err2=chyba(Y1,A,B,C);
if err2>err                                     %%% step is not accepted
      mu=mu*nu; nu=2*nu; err2=err;    
      A=A1; B=B1; C=C1;                     %%% return to the previous solution
   else
      rho=(err-err2)/(d'*(g+mu*d));
      err = err2;
      nu=2;
      mu=mu*max([1/3 1-(2*rho-1)^3]);
end   
iter(it)=err2;
if rem(it,10)==0
       [it err2]  % to monitor decrease of the cost function each 10 iterations
end    
end
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
%KRB Khatri-Rao product
%
% The columnwise Khatri-Rao-Bro product (Harshman, J.Chemom., 2002, 198-205)
% For two matrices with similar column dimension the khatri-Rao-Bro product
% is krb(A,B) = [kron(A(:,1),B(:,1)) .... kron(A(:,F),B(:,F))]
% 
% I/O AB = krb(A,B);
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.01 $ May 2001 $ Error in helpfile - A and B reversed $ RB $ Not compiled $

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
