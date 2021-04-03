function [ U,P,it,KKT ] = DecompPredPara(A,b,C,d,i0,u0,v0,p0,gamma,beta,rho_U,tol_U,itMax,tol )
% Decomposition par prediction version parallele for the quadratic case with affine inequality constraints
% 1/2<Au,u>-<b,u> such that Cu <= d
%   i0 : numero de la station avec allocation
%   u0 : initial point
%   v0 : initial allocation
%   p0 : initial price
%   0 < gamma <=1 pour allocation
%   0 < beta <=1 pour prix
%   rho_U : step of Uzawa
%   tol_U : tolerance for Uzawa

n=length(b); % number of variables
m=size(C,1); % number of constraints

P = zeros(m,itMax);
P(:,1)=p0;
U=zeros(n,itMax);
U(:,1)=u0;

it=1;
err=1;
C_temp=C;
C_temp(:,i0)=[];
KKT=zeros(itMax,1);
while err>tol && it<itMax && ~KKT(it)
    % i0
    [U(i0,it+1),lambda0, ~]=UzawaQuadraInequ(A(i0,i0),b(i0),C(:,i0),v0,u0(i0),P(:,it),rho_U,itMax,tol_U);
    % coordination for p
    P(:,it+1)=(1-beta)*P(:,it)+beta*lambda0;
    
    % else
    for i=1:n
        if i~=i0
            U(i,it+1)=A(i,i)\(b(i)-P(:,it)'*C(:,i));
        end
    end
    u1_temp=U(:,it+1);
    u1_temp(i0)=[];
    % coordination for v
    v0=(1-gamma)*v0+gamma*(d-C_temp*u1_temp);
    
    % verification KKT
    [ KKT1,KKT2,KKT3,KKT4 ] = testKKT( A,b,C,d,U(:,it+1),P(:,it+1));
    KKT(it+1) = all(abs(KKT1) < tol) & all(KKT2< tol) & all(KKT3 >=0) & all(abs(KKT4) <tol);
    err=norm(U(:,it+1)-U(:,it));
    it=it+1;
end
U = U(:,1:it);
P = P(:,1:it);
KKT = KKT(1:it);
fprintf('Decomposition par prediction parallel en %d iterations \n',it);
end

