function [ U,Lambda,it,KKT ] = DecompQuant( A,b,C,d,u0,v0,rho,rho_U,tol_U,itMax,tol )
% Algorithme decomposition par quantite for the quadratic case with affine inequality constraints
% 1/2<Au,u>-<b,u> such that Cu <= d
%   u0 : initial point
%   v0 : initial allocation
%   rho : step of projected gradient
%   rho_U : step of Uzawa
%   tol_U : tolerance for Uzawa

n=length(b); % number of variables
m=size(C,1); % number of constraints

Lambda = zeros(m,n,itMax);
U=zeros(n,itMax);
U(:,1)=u0;

it=1;
err=1;
KKT=zeros(itMax,1);
while err>tol && it<itMax && ~KKT(it)
    % subproblems
    for i=1:n
        [U(i,it+1),Lambda(:,i,it+1), ~]=UzawaQuadraInequ(A(i,i),b(i),C(:,i),v0(:,i),U(i,it),Lambda(:,i,it),rho_U,itMax,tol_U);
    end

    % coordination
    v1=v0+rho*(Lambda(:,:,it+1)-repmat(mean(Lambda(:,:,it+1),2),1,n));
    
    % verification KKT
    [ KKT1,KKT2,KKT3,KKT4 ] = testKKT( A,b,C,d,U(:,it+1),Lambda(:,n-1,it+1));
    KKT(it+1) = all(abs(KKT1) < tol) & all(KKT2< tol) & all(KKT3 >=0) & all(abs(KKT4) <tol);
    
    err= norm(v1-v0); 
    it=it+1;
    v0=v1;
end
U = U(:,1:it);
Lambda = Lambda(:,:,1:it);
KKT = KKT(1:it);
fprintf('Decomposition par quantite en %d iterations \n',it);
end

