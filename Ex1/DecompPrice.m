function [ U,P,it,KKT ] = DecompPrice( A,b,C,d,u0,p0,rho,itMax,tol )
% Algorithme decomposition par prix for the quadratic case with affine inequality constraints
% % 1/2<Au,u>-<b,u> such that Cu <= d
%   u0 : initial point
%   p0 : initial price
%   rho : pas de gradient

n=length(b); % number of variables
m=size(C,1); % number of constraints

U = zeros(n,itMax);
U(:,1) = u0;
P = zeros(m,itMax);
P(:,1) = p0;

it=1; 
err=1;
KKT=zeros(itMax,1);
while it<itMax && err>tol && ~KKT(it)
    % subproblems
    for i=1:n
        U(i,it+1)=A(i,i)\(b(i)-C(:,i)'*P(:,it)); 
    end
    
    % coordiantion
    P(:,it+1)= max(0,P(:,it)+rho*(C*U(:,it+1)-d));
    
    % verification KKT
    [ KKT1,KKT2,KKT3,KKT4 ] = testKKT( A,b,C,d,U(:,it+1),P(:,it+1));
    KKT(it+1) = all(abs(KKT1) < tol) & all(KKT2< tol) & all(KKT3 >=0) & all(abs(KKT4) <tol);
    
    err=norm(U(:,it+1)-U(:,it));
    it=it+1;
end

U = U(:,1:it);
P = P(:,1:it);
KKT = KKT(1:it);
fprintf('Decomposition par prix en %d iterations \n',it);
end