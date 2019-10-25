function [A,B,C,D,iter]=KLM4(T,R,m,numit,A,B,C,D)
%
% Krylov-Levenberg-Marquardt method for CP decomposition of order-4
% tensors. Input:
%
%    T... tensor to be decomposed
%    R ... tensor rank
%    m ... parameter controlling complexity of one iteration of KLM
%    numit .... required number of iterations
%    A,B,C,D .... initial estimate of factor matrices (if available)
%
% Programmed by Petr Tichavsky, October 2018
%
[Ia,Ib,Ic,Id]=size(T); 
Iab=Ia+Ib;
Iabc=Iab+Ic;
Iabcd=Iabc+Id;
if nargin<8
    A=randn(Ia,R);
    B=randn(Ib,R);
    C=randn(Ic,R);
    D=randn(Id,R);
end    
% I2=reshape(1:Ia*Ib*Ic*Id,Ia,Ib,Ic,Id);  %%%% pridano
% I3=permute(I2,[3,4,1,2]);  I3=I3(:);
% I4=permute(I2,[4,1,2,3]);  I4=I4(:);
% I2=permute(I2,[2,3,4,1]);  I2=I2(:);
Y1=reshape(T,Ia,Ib*Ic*Id);
Y2=reshape(permute(T,[2,3,4,1]),Ib,Ia*Ic*Id);
Y3=reshape(permute(T,[3,4,1,2]),Ic,Ia*Ib*Id);
Y4=reshape(permute(T,[4,1,2,3]),Id,Ia*Ib*Ic);
iter=zeros(1,numit);
na=sum(A.^2); nb=sum(B.^2); nc=sum(C.^2);  nd=sum(D.^2);
mm=[nb.*nc.*nd na.*nc.*nd nd.*nb.*na na.*nb.*nc];
mu=max(mm);  % initial parameter mu
nu=2;
err=chyba(Y1,A,B,C,D);              %%% computes the error of the current approximation
for it=1:numit
    AA=A'*A;
    BB=B'*B;
    CC=C'*C;
    DD=D'*D;
    gA=reshape(T,Ia,Ib*Ic*Id)*krb(D,krb(C,B))-A*(BB.*CC.*DD);   % computes error gradient (mttkrp)
    gB=Y2*krb(A,krb(D,C))-B*(AA.*CC.*DD);
    gC=Y3*krb(B,krb(A,D))-C*(BB.*AA.*DD);
    gD=Y4*krb(C,krb(B,A))-D*(BB.*AA.*CC);
    g=[gA(:); gB(:); gC(:); gD(:)];
    Y=zeros(R*Iabcd,m+2); %  Krylov subspace
    Z=zeros(R*Iabcd,m+1);
    ng=norm(g);
    Y(:,1)=g/ng;
    gA=gA/ng; gB=gB/ng; gC=gC/ng; gD=gD/ng;
    W=zeros(m+1,m+1); 
    for i=1:m+1
        i1=max([1,i-1]);
        yA=gA*(CC.*BB.*DD)+A*((gB'*B).*CC.*DD+((gC'*C).*DD+(gD'*D).*CC).*BB);
        yB=gB*(AA.*CC.*DD)+B*(((gC'*C).*AA+(gA'*A).*CC).*DD+(gD'*D).*AA.*CC);
        yC=gC*(AA.*BB.*DD)+C*(((gB'*B).*AA+(gA'*A).*BB).*DD+(gD'*D).*BB.*AA);
        yD=gD*(AA.*BB.*CC)+D*(((gB'*B).*AA+(gA'*A).*BB).*CC+(gC'*C).*BB.*AA);
        Z(:,i)=[yA(:); yB(:); yC(:); yD(:)];
        W(i1:i,i)=Y(:,i1:i)'*Z(:,i);
        Y(:,i+1)=Z(:,i)-Y(:,i1:i)*W(i1:i,i);
        Y(:,i+1)=Y(:,i+1)/norm(Y(:,i+1));
        gA=reshape(Y(1:R*Ia,i+1),Ia,R); 
        gB=reshape(Y(R*Ia+1:R*Iab,i+1),Ib,R); 
        gC=reshape(Y(R*Iab+1:R*Iabc,i+1),Ic,R);
        gD=reshape(Y(R*Iabc+1:R*Iabcd,i+1),Id,R);
    end
    W(2:m+2:end)=W(m+2:m+2:end);
    %W=Y(:,1:m+1)'*Z;
    %b1=(H+mu*eye(size(H)))\g;
    d=g/mu-Y(:,1:m+1)/(inv(W)+eye(m+1)/mu)*(Y(:,1:m+1)'*g)/mu^2;
    A1=A; B1=B; C1=C;
    A=A+reshape(d(1:Ia*R),Ia,R); 
    B=B+reshape(d(1+Ia*R:Iab*R),Ib,R); 
    C=C+reshape(d(1+Iab*R:Iabc*R),Ic,R); 
    D=D+reshape(d(1+Iabc*R:Iabcd*R),Id,R); 
    err2=chyba(Y1,A,B,C,D);
    rho=(err-err2)/(d'*(g+mu*d));
    if err2>err                                     %%% step is not accepted
       mu=mu*nu; nu=2*nu; err2=err;    
       A=A1; B=B1; C=C1;                     %%% return to the previous solution
    else
      nu=2;                                        
      mu=mu*max([1/3 1-(2*rho-1)^3]);
      err=err2;
    end   
    iter(it)=err2;
    if rem(it,10)==0
        [it err2]
    end    
end
end
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
function [err,Y]=chyba(X,A,B,C,D)
%
% computes an error of approximation of X by sum
% of outer products of columns of A,B and C
%
Y=X-A*krb(D,krb(C,B))';
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
